//! Code to handle the GNOme ontology
use std::collections::HashMap;

use context_error::{
    BoxedError, CreateError, FullErrorContent, StaticErrorContent, combine_errors,
};
use thin_vec::ThinVec;

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType};

use crate::{
    glycan::{GlycanStructure, MonoSaccharide},
    ontology::Ontology,
    sequence::{GnoComposition, GnoSubsumption, ModificationId, SimpleModificationInner},
};

/// GNOme modifications
///
/// These integrate the GNOme ontology with the structures from glycosmos.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct Gnome {}

impl CVSource for Gnome {
    type Data = SimpleModificationInner;
    fn cv_name() -> &'static str {
        "GNOme"
    }

    fn files() -> &'static [CVFile] {
        &[
            CVFile {
                name: "GNOme",
                extension: "obo",
                url: Some("https://purl.obolibrary.org/obo/gno.obo"),
                compression: mzcv::CVCompression::None,
            },
            CVFile {
                name: "glycosmos_glycans_list",
                extension: "csv",
                url: Some("https://glycosmos.org/download/glycosmos_glycans_list.csv"),
                compression: mzcv::CVCompression::None,
            },
        ]
    }

    fn static_data() -> Option<(CVVersion, Vec<Self::Data>)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            use mzcv::CVCache;
            let cache: <SimpleModificationInner as mzcv::CVData>::Cache =
                bincode::decode_from_slice::<
                    <SimpleModificationInner as mzcv::CVData>::Cache,
                    Configuration,
                >(
                    include_bytes!("../databases/gnome.dat"),
                    Configuration::default(),
                )
                .unwrap()
                .0;
            Some(cache.deconstruct())
        }
        #[cfg(feature = "internal-no-data")]
        None
    }

    fn parse(
        mut readers: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
        let reader = readers.next().unwrap();
        let (version, mut mods) = OboOntology::from_raw(reader)
            .map_err(|e| {
                vec![
                    BoxedError::small(
                        CVError::FileCouldNotBeParsed,
                        e.get_short_description(),
                        e.get_long_description(),
                    )
                    .add_contexts(e.get_contexts().iter().cloned()),
                ]
            })
            .map(|obo| (obo.version(), parse_gnome(obo)))?;
        let read_mods = mods.clone();
        let structures = parse_gnome_structures(readers.next().unwrap());

        // Fill all known info points
        for modification in mods.values_mut() {
            if modification.weight.is_none() {
                modification.weight = find_mass(&read_mods, modification.is_a.clone());
            }
            if let Some(structure) = structures.get(&modification.id.name) {
                modification.topology = structure.structure.clone();
                if let Some(chebi) = structure.chebi {
                    modification
                        .id
                        .cross_ids
                        .push((Some("ChEBI".to_string()), chebi.to_string()));
                }
                if let Some(pubchem) = structure.pubchem {
                    modification
                        .id
                        .cross_ids
                        .push((Some("PubChemCID".to_string()), pubchem.to_string()));
                }
                modification.motif = structure.motif.clone();
                modification.taxonomy = structure.taxonomy.clone().into();
                modification.glycomeatlas = structure.glycomeatlas.clone().into();
            } else if let Some(id) = &modification.topology_id {
                modification.topology = structures.get(id).and_then(|s| s.structure.clone());
            }
            if modification.composition.is_none()
                && let Some(composition_id) = &modification.composition_id
            {
                modification.composition = read_mods
                    .get(composition_id)
                    .and_then(|g| g.composition.clone());
            }
        }

        Ok((
            version,
            mods.into_values().filter_map(|m| m.try_into().ok()),
        ))
    }
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
    pub(super) const SUBSUMPTION_MOLECULAR_WEIGHT: &str = "GNO:00000012";
    pub(super) const SUBSUMPTION_BASECOMPOSITION: &str = "GNO:00000013";
    pub(super) const SUBSUMPTION_COMPOSITION: &str = "GNO:00000014";
    pub(super) const SUBSUMPTION_TOPOLOGY: &str = "GNO:00000015";
    pub(super) const SUBSUMPTION_SACCHARIDE: &str = "GNO:00000016";

    pub(super) const HAS_SUBSUMPTION_CATEGORY: &str = "GNO:00000021";
    pub(super) const HAS_GLYTOUCAN_ID: &str = "GNO:00000022";
    pub(super) const HAS_GLYTOUCAN_LINK: &str = "GNO:00000023";
    pub(super) const IS_SUBSUMED_BY: &str = "GNO:00000024";
    pub(super) const IS_RESTRICTION_MEMBER: &str = "GNO:00000025";
    /// Is the basic composition the same, meaning the same monosaccharides (without isomeric information)
    pub(super) const HAS_BASECOMPOSITION: &str = "GNO:00000033";
    /// Is the composition the same, meaning the same monosaccharides
    pub(super) const HAS_COMPOSITION: &str = "GNO:00000034";
    /// Is the basic structure the same (if anomeric and linkage information are thrown overboard)
    pub(super) const HAS_TOPOLOGY: &str = "GNO:00000035";
    /// Is the linked structure the same (if anomeric and reducing end ring information are thrown overboard)
    pub(super) const HAS_ARCHETYPE: &str = "GNO:00000036";
    pub(super) const HAS_STRUCTURE_BROWSER_LINK: &str = "GNO:00000041";
    pub(super) const HAS_COMPOSITION_BROWSER_LINK: &str = "GNO:00000042";
    pub(super) const SHORTUCKB_COMPOSITION: &str = "GNO:00000101";
    /// Indicates the precision of the definition, lower is better, 0 is a fully defined glycan
    pub(super) const HAS_STRUCTURE_CHARACTERISATION_SCORE: &str = "GNO:00000102";
    pub(super) const HAS_BYONIC_NAME: &str = "GNO:00000202";
}

use gnome_terms::*;
use itertools::Itertools;

/// Get the GNO Subsumption level.
/// # Errors
/// If this is not a known level.
fn gno_subsumption_from_str(s: &str) -> Result<GnoSubsumption, ()> {
    match s {
        SUBSUMPTION_MOLECULAR_WEIGHT => Ok(GnoSubsumption::AverageWeight),
        SUBSUMPTION_BASECOMPOSITION => Ok(GnoSubsumption::BaseComposition),
        SUBSUMPTION_COMPOSITION => Ok(GnoSubsumption::Composition),
        SUBSUMPTION_TOPOLOGY => Ok(GnoSubsumption::Topology),
        SUBSUMPTION_SACCHARIDE => Ok(GnoSubsumption::Saccharide),
        _ => Err(()),
    }
}

/// Parse the GNOme ontology .obo file
/// # Errors
/// If the file is not valid.
fn parse_gnome(obo: OboOntology) -> HashMap<String, GNOmeModification> {
    let mut mods = HashMap::new();
    let mut errors = Vec::new();

    for obj in obo.objects {
        if obj.stanza_type != OboStanzaType::Term || !obj.lines.contains_key("is_a") {
            continue;
        }

        let modification = GNOmeModification {
            id: ModificationId {
                ontology: Ontology::Gnome,
                name: obj.id.1,
                id: None,
                description: obj.definition.map_or_else(String::new, |d| d.0),
                synonyms: obj.synonyms.iter().map(|s| s.synonym.clone()).collect(),
                cross_ids: obj
                    .property_values
                    .get(HAS_GLYTOUCAN_ID)
                    .map(|v| (Some("GlyTouCan".to_string()), v[0].0.to_string()))
                    .into_iter()
                    .chain(
                        obj.property_values
                            .get(HAS_GLYTOUCAN_LINK)
                            .map(|v| (Some("GlyTouCanURL".to_string()), v[0].0.to_string())),
                    )
                    .chain(
                        obj.property_values
                            .get(HAS_COMPOSITION_BROWSER_LINK)
                            .map(|v| (Some("CompositionBrowser".to_string()), v[0].0.to_string())),
                    )
                    .chain(
                        obj.property_values
                            .get(HAS_STRUCTURE_BROWSER_LINK)
                            .map(|v| (Some("StructureBrowser".to_string()), v[0].0.to_string())),
                    )
                    .collect(),
            },
            subsumption_level: obj
                .lines
                .get(HAS_SUBSUMPTION_CATEGORY)
                .map(|s| gno_subsumption_from_str(&s[0].0).unwrap())
                .unwrap_or_default(),
            structure_score: obj
                .lines
                .get(HAS_STRUCTURE_CHARACTERISATION_SCORE)
                .map(|s| s[0].0.parse().unwrap()),
            is_a: obj.lines["is_a"][0].0.trim()[4..]
                .trim()
                .to_ascii_lowercase(),
            composition_id: obj
                .property_values
                .get(HAS_COMPOSITION)
                .map(|lines| lines[0].0.to_string()),
            topology_id: obj
                .property_values
                .get(HAS_TOPOLOGY)
                .map(|lines| lines[0].0.to_string()),
            weight: obj
                .lines
                .get("name")
                .map(|e| &e[0])
                .filter(|(n, _, _)| n.len() > 30)
                .and_then(|(name, _, _)| name[27..name.len() - 3].parse::<f64>().ok()),
            composition: obj
                .property_values
                .get(SHORTUCKB_COMPOSITION)
                .and_then(|lines| {
                    match MonoSaccharide::pro_forma_composition::<false>(&lines[0].0.to_string()) {
                        Ok((v, _)) => Some(v),
                        Err(e) => {
                            combine_errors(
                                &mut errors,
                                e.into_iter().map(BoxedError::to_owned),
                                (),
                            );
                            None
                        }
                    }
                }),
            topology: None, // Will be looked up later
            motif: None,
            taxonomy: ThinVec::new(),
            glycomeatlas: ThinVec::new(),
        };

        mods.insert(modification.id.name.clone(), modification);
    }

    for error in errors {
        println!("{error}");
    }

    mods
}

#[derive(Debug)]
struct GlycosmosList {
    structure: Option<GlycanStructure>,
    motif: Option<(String, String)>,
    chebi: Option<usize>,
    pubchem: Option<usize>,
    taxonomy: Vec<(String, usize)>,
    glycomeatlas: Vec<(String, Vec<(String, String)>)>,
}

/// Parse the glycosmos glycan structures .csv file
/// # Panics
/// If the file is not valid.
fn parse_gnome_structures(
    file: HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>,
) -> HashMap<String, GlycosmosList> {
    let mut glycans = HashMap::new();
    let mut errors = 0;
    for line in crate::csv::parse_csv_raw(file, b',', None, None).unwrap() {
        let line = line.unwrap();

        glycans.insert(
            line.index_column("accession number").unwrap().0.to_string(),
            GlycosmosList {
                structure: line
                    .index_column("iupac condensed")
                    .ok()
                    .filter(|p| !p.0.is_empty())
                    .and_then(|(_, range)| {
                        match GlycanStructure::from_short_iupac(
                            line.line(),
                            range.clone(),
                            (line.line_index() + 1) as u32,
                        ) {
                            Ok(glycan) => Some(glycan),
                            Err(error) => {
                                if errors < 5 {
                                    println!("{error}");
                                }
                                errors += 1;
                                None
                            }
                        }
                    }),
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
                                .map(|(n, i)| (n.to_string(), i.parse::<usize>().unwrap()))
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
                    .chunk_by(|(species, _)| (*species).to_string())
                    .into_iter()
                    .map(|(species, locations)| {
                        (
                            (*species).to_string(),
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
    assert!(
        errors <= 0,
        "Total glycan structure reading errors: {errors} total read {}",
        glycans.len()
    );
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

impl TryFrom<GNOmeModification> for SimpleModificationInner {
    type Error = ();
    fn try_from(value: GNOmeModification) -> Result<Self, ()> {
        Ok(Self::Gno {
            composition: if let Some(structure) = value.topology {
                GnoComposition::Topology(structure)
            } else if let Some(composition) = value.composition {
                GnoComposition::Composition(composition)
            } else if let Some(mass) = value.weight {
                GnoComposition::Weight(crate::system::f64::da(mass).into())
            } else {
                return Err(());
            },
            id: value.id,
            structure_score: value.structure_score,
            subsumption_level: value.subsumption_level,
            motif: value.motif,
            taxonomy: value.taxonomy,
            glycomeatlas: value.glycomeatlas,
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
