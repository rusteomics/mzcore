//! Code to handle the PSI-MOD ontology
use std::{borrow::Cow, collections::HashMap};

use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};
use itertools::Itertools;

use mzcv::{
    CVData, CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType,
    SynonymScope, OboIdentifier
};

use crate::{
    chemistry::MolecularFormula,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        LinkerSpecificity, ModificationId, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};

/// PSI-MOD modifications
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct PsiMod {}

impl CVData for SimpleModificationInner {
    type Index = u32;
    fn index(&self) -> Option<u32> {
        self.description().and_then(ModificationId::id)
    }
    fn name(&self) -> Option<Cow<'_, str>> {
        self.description().map(|d| Cow::Borrowed(d.name.as_ref()))
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.description().into_iter().flat_map(|d| {
            d.synonyms
                .iter()
                .filter(|(s, _)| *s == SynonymScope::Exact)
                .map(|(_, s)| s.as_ref())
        })
    }
}

impl CVSource for PsiMod {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "PSI-MOD"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "PSI-MOD",
            extension: "obo",
            url: Some(
                "https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/refs/heads/master/PSI-MOD.obo",
            ),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/psimod.dat"),
                Configuration::default(),
            )
            .unwrap()
            .0;
            Some(cache)
        }
        #[cfg(feature = "internal-no-data")]
        None
    }

    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, Self::Structure), Vec<BoxedError<'static, CVError>>> {
        let reader = reader.next().unwrap();
        OboOntology::from_raw(reader)
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
            .map(|obo| {
                (obo.version(), {
                    let mut mods: Vec<SimpleModification> = Vec::new();
                    let mut ignored: HashMap<usize, Vec<OboIdentifier>> = HashMap::new();

                    'stanza: for obj in obo.objects {
                        if obj.stanza_type != OboStanzaType::Term
                            || obj.id == (Some("MOD".into()), "00000".into())
                            || obj.id == (Some("MOD".into()), "00004".into())
                            || obj.id == (Some("MOD".into()), "00008".into())
                        {
                            continue;
                        }
                        let mut modification = OntologyModification {
                            id: obj
                                .id
                                .1
                                .parse()
                                .expect("Incorrect psi mod id, should be numerical"),
                            name: obj.lines["name"][0].0.clone(),
                            ontology: Ontology::Psimod,
                            obsolete: obj.obsolete,
                            ..OntologyModification::default()
                        };
                        if let Some((description, cross_ids, _, _)) = obj.definition {
                            modification.description = description;
                            modification.cross_ids = cross_ids.into();
                        }
                        for synonym in obj.synonyms {
                            modification
                                .synonyms
                                .push((synonym.scope, synonym.synonym.clone())); // TODO: retain the scope
                        }

                        
                        let mut origins = Vec::new();
                        let mut term = None;
                        for (id, _values, _comment) in obj.xref {
                            match (id.0.as_deref(), id.1) {
                                (Some("DiffFormula"), s) if s.as_ref() != "\"none\"" => {
                                    modification.formula = MolecularFormula::psi_mod(
                                        s.trim_start_matches('\"').trim_end_matches('\"'),
                                    )
                                    .unwrap();
                                }
                                (Some("Origin"), value)  => {
                                    let commas = value.chars().filter(|c| *c == ',').count();
                                    let cross_link_numbers = obj.lines.get("comment").map_or(Vec::new(), |v| v[0].0.split(';').flat_map(|s| s.split('.')).filter_map(|tag| tag.trim().to_ascii_lowercase().strip_prefix("cross-link").map(|v| v.trim().trim_end_matches('.').parse::<usize>().map_err(|e| format!("Could not parse '{tag}': {e}")).unwrap())).collect());
                                    if cross_link_numbers.len() > 1 {
                                        println!("PSI-MOD: ERROR Found multiple cross-link numbers for MOD:{:05}: {}", modification.id, cross_link_numbers.iter().join(", "));
                                        continue 'stanza;
                                    }

                                    if let Some(amount) = cross_link_numbers.first() {
                                        if commas + 1 != *amount {
                                            println!("PSI-MOD: ERROR Found {value} at MOD:{:05}, cross-link: {amount} with {commas} commas", modification.id);
                                            continue 'stanza;
                                        }
                                    } else if value.contains(',') {
                                        println!("PSI-MOD: ERROR Found {value} at MOD:{:05} with {commas} commas but without a cross-link comment", modification.id);
                                        continue 'stanza;
                                    }
                                    if value.as_ref() != "\"none\"" {
                                        origins = value
                                            .trim_start_matches('\"')
                                            .trim_end_matches('\"')
                                            .split(',')
                                            .map(|s| s.trim().to_string())
                                            .collect();
                                    }
                                }
                                (Some("TermSpec"), value) => {
                                    let v = value.trim_start_matches('\"').trim_end_matches('\"');
                                    if v == "N-term" {
                                        term = Some(Position::AnyNTerm);
                                    } else if v == "C-term" {
                                        term = Some(Position::AnyCTerm);
                                    } else if v == "none" {
                                        term = Some(Position::Anywhere);
                                    } else {
                                        panic!("Invalid TermSpec: {v}");
                                    }
                                }
                                _ => (), // ignore
                            }
                        }
                        if origins.len() <= 1 {
                            let mut rules = Vec::new();
                            if let Some(origin) = origins.first() {
                                rules.push((
                                    vec![parse_rule(origin, term)],
                                    Vec::new(),
                                    Vec::new(),
                                ));
                            } else if let Some(term) = term {
                                rules.push((
                                    vec![PlacementRule::Terminal(term)],
                                    Vec::new(),
                                    Vec::new(),
                                ));
                            }
                            modification.data = ModData::Mod {
                                specificities: rules,
                            };
                        } else if origins.len() == 2 {
                            modification.data = ModData::Linker { length: crate::sequence::LinkerLength::Unknown, specificities: vec![LinkerSpecificity::Asymmetric { 
                                rules: (vec![parse_rule(&origins[0], term)], vec![parse_rule(&origins[1], term)]), 
                                stubs: Vec::new(), neutral_losses: Vec::new(), diagnostic: Vec::new() }] };
                        } else {
                            ignored.entry(origins.len()).or_default().push(obj.id.clone());
                        }
                        mods.push(modification.into());
                    }

                    for (sites, ids) in ignored.into_iter().sorted() {
                        println!(
                            "PSI-MOD: Ignored sites {sites}: {}",
                            ids.into_iter().sorted().join(","),
                        );
                    }

                    mods
                })
            })
    }
}

fn parse_rule(rule: &str, term: Option<Position>) -> PlacementRule {
    if rule.len() == 1 {
        if rule == "X" {
            term.map_or_else(|| PlacementRule::Anywhere, PlacementRule::Terminal)
        } else {
            PlacementRule::AminoAcid(
                vec![rule.try_into().unwrap()].into(),
                term.unwrap_or(Position::Anywhere),
            )
        }
    } else {
        let Some((_, id)) = rule.split_once(':') else {
            panic!("Invalid PSI-MOD in origin: '{rule}' does not contain a colon ':'")
        };
        let Ok(id) = id.parse() else {
            panic!("Invalid PSI-MOD in origin: '{rule}' the id is not numerical")
        };
        PlacementRule::PsiModification(id, term.unwrap_or(Position::Anywhere))
    }
}
