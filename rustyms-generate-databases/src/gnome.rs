use std::{collections::HashMap, io::Write, path::Path, sync::Arc};

use crate::{csv::parse_csv, glycan::*, SimpleModification};

use super::{
    obo::OboOntology, ontology_modification::OntologyModificationList, GnoComposition,
    GnoSubsumption, ModificationId, SimpleModificationInner,
};

use bincode::config::Configuration;
use thin_vec::ThinVec;

pub fn build_gnome_ontology(out_dir: &Path) {
    // Get all the basic info
    let mut mods = parse_gnome();
    let read_mods = mods.clone();
    let structures = parse_gnome_structures();

    // Fill all known info points
    for modification in mods.values_mut() {
        if modification.weight.is_none() {
            modification.weight = find_mass(&read_mods, modification.is_a.clone());
        }
        if let Some(structure) = structures.get(&modification.id.name) {
            modification.topology = Some(structure.structure.clone());
            if let Some(chebi) = structure.chebi {
                modification
                    .id
                    .cross_ids
                    .push(("ChEBI".to_string(), chebi.to_string()));
            }
            if let Some(pubchem) = structure.pubchem {
                modification
                    .id
                    .cross_ids
                    .push(("PubChemCID".to_string(), pubchem.to_string()));
            }
            modification.motif = structure.motif.clone();
            modification.taxonomy = structure.taxonomy.clone().into();
            modification.glycomeatlas = structure.glycomeatlas.clone().into();
        } else if let Some(id) = &modification.topology_id {
            modification.topology = structures.get(id).map(|s| s.structure.clone());
        }
        if let Some(composition_id) = &modification.composition_id {
            modification.composition = read_mods
                .get(composition_id)
                .and_then(|g| g.composition.clone());
        }
    }

    // Write out the data
    let dest_path = Path::new(&out_dir).join("gnome.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods
        .into_values()
        .filter(|m| m.weight.is_some())
        .sorted_unstable()
        .map(|m| (None, m.id.name.clone(), m.into_mod()))
        .collect::<OntologyModificationList>();
    println!("Found {} GNOme modifications", final_mods.len());
    file.write_all(
        &bincode::serde::encode_to_vec::<OntologyModificationList, Configuration>(
            final_mods,
            Configuration::default(),
        )
        .unwrap(),
    )
    .unwrap();
}

fn find_mass(mods: &HashMap<String, GNOmeModification>, mut name: String) -> Option<f64> {
    let mut mass = None;
    while mass.is_none() {
        mass = mods.get(&name)?.weight;
        name.clone_from(&mods[&name].is_a);
    }
    mass
}

#[expect(dead_code)]
mod gnome_terms {
    pub const SUBSUMPTION_MOLECULAR_WEIGHT: &str = "GNO:00000012";
    pub const SUBSUMPTION_BASECOMPOSITION: &str = "GNO:00000013";
    pub const SUBSUMPTION_COMPOSITION: &str = "GNO:00000014";
    pub const SUBSUMPTION_TOPOLOGY: &str = "GNO:00000015";
    pub const SUBSUMPTION_SACCHARIDE: &str = "GNO:00000016";

    pub const HAS_SUBSUMPTION_CATEGORY: &str = "GNO:00000021";
    pub const HAS_GLYTOUCAN_ID: &str = "GNO:00000022";
    pub const HAS_GLYTOUCAN_LINK: &str = "GNO:00000023";
    pub const IS_SUBSUMED_BY: &str = "GNO:00000024";
    pub const IS_RESTRICTION_MEMBER: &str = "GNO:00000025";
    /// Is the basic composition the same, meaning the same monosaccharides (without isomeric information)
    pub const HAS_BASECOMPOSITION: &str = "GNO:00000033";
    /// Is the composition the same, meaning the same monosaccharides
    pub const HAS_COMPOSITION: &str = "GNO:00000034";
    /// Is the basic structure the same (if anomeric and linkage information are thrown overboard)
    pub const HAS_TOPOLOGY: &str = "GNO:00000035";
    /// Is the linked structure the same (if anomeric and reducing end ring information are thrown overboard)
    pub const HAS_ARCHETYPE: &str = "GNO:00000036";
    pub const HAS_STRUCTURE_BROWSER_LINK: &str = "GNO:00000041";
    pub const HAS_COMPOSITION_BROWSER_LINK: &str = "GNO:00000042";
    pub const SHORTUCKB_COMPOSITION: &str = "GNO:00000101";
    /// Indicates the precision of the definition, lower is better, 0 is a fully defined glycan
    pub const HAS_STRUCTURE_CHARACTERISATION_SCORE: &str = "GNO:00000102";
    pub const HAS_BYONIC_NAME: &str = "GNO:00000202";
}

use gnome_terms::*;
use itertools::Itertools;

impl std::str::FromStr for GnoSubsumption {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            SUBSUMPTION_MOLECULAR_WEIGHT => Ok(Self::AverageWeight),
            SUBSUMPTION_BASECOMPOSITION => Ok(Self::BaseComposition),
            SUBSUMPTION_COMPOSITION => Ok(Self::Composition),
            SUBSUMPTION_TOPOLOGY => Ok(Self::Topology),
            SUBSUMPTION_SACCHARIDE => Ok(Self::Saccharide),
            _ => Err(()),
        }
    }
}

fn parse_gnome() -> HashMap<String, GNOmeModification> {
    let obo = OboOntology::from_file("rustyms-generate-databases/data/GNOme.obo")
        .expect("Not a valid obo file");
    let mut mods = HashMap::new();

    for obj in obo.objects {
        if obj.name != "Term" || !obj.lines.contains_key("is_a") {
            continue;
        }

        // name: glycan of molecular weight 40.03 Da
        let modification = GNOmeModification {
            id: ModificationId {
                ontology: super::Ontology::Gnome,
                name: obj.lines["id"][0][4..].to_ascii_lowercase(),
                id: None,
                description: obj.lines.get("def").map_or(String::new(), |d| {
                    d[0].trim_matches("\"[] ".chars().collect_vec().as_slice())
                        .to_string()
                }),
                synonyms: obj.lines.get("synonym").map_or(ThinVec::new(), |s| {
                    s.iter()
                        .filter_map(|s| s[1..].split_once('"').map(|(s, _)| s.to_string()))
                        .collect()
                }),
                cross_ids: obj
                    .property_values
                    .get(HAS_GLYTOUCAN_ID)
                    .map(|v| ("GlyTouCan".to_string(), v[0].to_string()))
                    .into_iter()
                    .chain(
                        obj.property_values
                            .get(HAS_GLYTOUCAN_LINK)
                            .map(|v| ("GlyTouCanURL".to_string(), v[0].to_string())),
                    )
                    .chain(
                        obj.property_values
                            .get(HAS_COMPOSITION_BROWSER_LINK)
                            .map(|v| ("CompositionBrowser".to_string(), v[0].to_string())),
                    )
                    .chain(
                        obj.property_values
                            .get(HAS_STRUCTURE_BROWSER_LINK)
                            .map(|v| ("StructureBrowser".to_string(), v[0].to_string())),
                    )
                    .collect(),
            },
            subsumption_level: obj
                .lines
                .get(HAS_SUBSUMPTION_CATEGORY)
                .map(|s| s[0].parse().unwrap())
                .unwrap_or_default(),
            structure_score: obj
                .lines
                .get(HAS_STRUCTURE_CHARACTERISATION_SCORE)
                .map(|s| s[0].parse().unwrap()),
            is_a: obj.lines["is_a"][0].trim()[4..]
                .split_once("!")
                .unwrap()
                .0
                .trim()
                .to_ascii_lowercase(),
            composition_id: obj
                .property_values
                .get(HAS_COMPOSITION)
                .map(|lines| lines[0].to_string()),
            topology_id: obj
                .property_values
                .get(HAS_TOPOLOGY)
                .map(|lines| lines[0].to_string()),
            weight: obj
                .lines
                .get("name")
                .map(|e| &e[0])
                .filter(|n| n.len() > 30)
                .and_then(|name| name[27..name.len() - 3].parse::<f64>().ok()),
            composition: obj
                .property_values
                .get(SHORTUCKB_COMPOSITION)
                .map(|lines| MonoSaccharide::from_composition(&lines[0].to_string()).unwrap()),
            topology: None, // Will be looked up later
            motif: None,
            taxonomy: ThinVec::new(),
            glycomeatlas: ThinVec::new(),
        };

        mods.insert(modification.id.name.clone(), modification);
    }

    mods
}

struct GlycosmosList {
    structure: GlycanStructure,
    motif: Option<(String, String)>,
    chebi: Option<usize>,
    pubchem: Option<usize>,
    taxonomy: Vec<(String, usize)>,
    glycomeatlas: Vec<(String, Vec<(String, String)>)>,
}

fn parse_gnome_structures() -> HashMap<String, GlycosmosList> {
    let mut glycans = HashMap::new();
    let mut errors = 0;
    for line in parse_csv(
        "rustyms-generate-databases/data/glycosmos_glycans_list.csv",
        b',',
        None,
    )
    .unwrap()
    .skip(1)
    {
        let line = line.unwrap();
        if !line.index_column("iupac condensed").unwrap().0.is_empty() {
            match GlycanStructure::from_short_iupac(
                line.line(),
                line.range(1).clone(),
                line.line_index() + 1,
            ) {
                Ok(glycan) => {
                    glycans.insert(
                        line.index_column("accession number")
                            .unwrap()
                            .0
                            .to_lowercase(),
                        GlycosmosList {
                            structure: glycan,
                            motif: line
                                .index_column("motif name(s)")
                                .ok()
                                .filter(|p| !p.0.is_empty())
                                .and_then(|p| p.0.split_once(':'))
                                .map(|(n, i)| (n.to_string(), i.to_ascii_lowercase())),
                            chebi: line
                                .index_column("chebi")
                                .ok()
                                .map(|p| p.0.to_string())
                                .filter(|p| !p.is_empty())
                                .and_then(|p| p.parse::<usize>().ok()),
                            pubchem: line
                                .index_column("pubchem cid")
                                .ok()
                                .map(|p| p.0.to_string())
                                .filter(|p| !p.is_empty())
                                .and_then(|p| p.parse::<usize>().ok()),
                            taxonomy: line
                                .index_column("taxonomy")
                                .ok()
                                .filter(|p| !p.0.is_empty())
                                .into_iter()
                                .flat_map(|p| {
                                    p.0.split(',').map(|s| {
                                        s.rsplit_once(':')
                                            .map(|(n, i)| {
                                                (n.to_string(), i.parse::<usize>().unwrap())
                                            })
                                            .unwrap()
                                    })
                                })
                                .collect(),
                            glycomeatlas: line
                                .index_column("glycomeatlas")
                                .ok()
                                .filter(|p| !p.0.is_empty())
                                .map(|p| p.0)
                                .into_iter()
                                .flat_map(|p| p.split(','))
                                .map(|p| p.split_once(':').unwrap())
                                .chunk_by(|(species, _)| species.to_string())
                                .into_iter()
                                .map(|(species, locations)| {
                                    (
                                        species.to_string(),
                                        locations
                                            .into_iter()
                                            .map(|location| {
                                                location
                                                    .1
                                                    .trim_end_matches(')')
                                                    .split_once('(')
                                                    .map(|(l, o)| (l.to_string(), o.to_string()))
                                                    .unwrap()
                                            })
                                            .collect(),
                                    )
                                })
                                .collect(),
                        },
                    );
                }
                Err(error) => {
                    if errors < 5 {
                        println!("{error}");
                    }
                    errors += 1;
                }
            }
        }
    }
    if errors > 0 {
        panic!(
            "Total glycan structure reading errors: {errors} total read {}",
            glycans.len()
        );
    }
    glycans
}

#[derive(Clone, Debug)]
struct GNOmeModification {
    id: ModificationId,
    /// subsumption level, indication to what extent this species is described
    subsumption_level: GnoSubsumption,
    /// id to prev in chain
    is_a: String,
    /// id of composition
    composition_id: Option<String>,
    /// id of topology
    topology_id: Option<String>,
    /// molecular weight if defined
    weight: Option<f64>,
    /// composition if defined
    composition: Option<Vec<(MonoSaccharide, isize)>>,
    /// structure if defined
    topology: Option<GlycanStructure>,
    /// the score for the structure (0 if fully defined)
    structure_score: Option<usize>,
    /// The underlying glycan motifs
    motif: Option<(String, String)>,
    /// Taxonomy of where the glycan exists
    taxonomy: ThinVec<(String, usize)>,
    /// Locations of where the glycan exists
    glycomeatlas: ThinVec<(String, Vec<(String, String)>)>,
}

impl GNOmeModification {
    fn into_mod(self) -> SimpleModification {
        Arc::new(SimpleModificationInner::Gno {
            composition: if let Some(structure) = self.topology {
                GnoComposition::Topology(structure)
            } else if let Some(composition) = self.composition {
                GnoComposition::Composition(composition)
            } else if let Some(mass) = self.weight {
                GnoComposition::Weight(crate::system::f64::da(mass).into())
            } else {
                panic!("unreachable")
            },
            id: self.id,
            structure_score: self.structure_score,
            subsumption_level: self.subsumption_level,
            motif: self.motif,
            taxonomy: self.taxonomy,
            glycomeatlas: self.glycomeatlas,
        })
    }
}

impl PartialEq for GNOmeModification {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
            && self.subsumption_level == other.subsumption_level
            && self.is_a == other.is_a
            && self.composition_id == other.composition_id
            && self.topology_id == other.topology_id
            && (self.weight.is_none() && other.weight.is_none()
                || self
                    .weight
                    .is_some_and(|sw| other.weight.is_some_and(|ow| sw.total_cmp(&ow).is_eq())))
            && self.composition == other.composition
            && self.topology == other.topology
            && self.structure_score == other.structure_score
    }
}
impl Eq for GNOmeModification {}

impl PartialOrd for GNOmeModification {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for GNOmeModification {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.id.name.cmp(&other.id.name)
    }
}
