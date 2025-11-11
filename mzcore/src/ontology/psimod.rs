//! Code to handle the PSI-MOD ontology
use std::borrow::Cow;

use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVData, CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType,
    SynonymScope,
};

use crate::{
    chemistry::MolecularFormula,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        ModificationId, PlacementRule, Position, SimpleModification, SimpleModificationInner,
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

                    for obj in obo.objects {
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

                        let mut rules = Vec::new();
                        let mut origins = Vec::new();
                        let mut term = None;
                        for (id, values) in obj.property_values {
                            for (value, _, _) in values {
                                match (id.as_ref(), value) {
                                    ("DiffFormula", mzcv::OboValue::String(s)) => {
                                        modification.formula =
                                            MolecularFormula::psi_mod(&s).unwrap();
                                    }
                                    ("Origin", mzcv::OboValue::String(value)) => {
                                        origins = value
                                            .split(',')
                                            .map(|s| s.trim().to_string())
                                            .collect();
                                    }
                                    ("TermSpec", mzcv::OboValue::String(value)) => {
                                        if value == "N-term" {
                                            term = Some(Position::AnyNTerm);
                                        } else if value == "C-term" {
                                            term = Some(Position::AnyCTerm);
                                        } else {
                                            panic!("Invalid TermSpec: {value}");
                                        }
                                    }
                                    _ => (), // ignore
                                }
                            }
                        }
                        // If the list of possible origins contains "X" than the mod can be placed on any aminoacid
                        // But if there is a TermSpec definition that should still be accounted for
                        let all_aminoacids = origins.contains(&"X".to_string());
                        if !all_aminoacids {
                            for origin in &origins {
                                if origin.len() == 1 {
                                    rules.push((
                                        vec![PlacementRule::AminoAcid(
                                            vec![origin.try_into().unwrap()].into(),
                                            term.unwrap_or(Position::Anywhere),
                                        )],
                                        Vec::new(),
                                        Vec::new(),
                                    ));
                                } else {
                                    rules.push((
                                        vec![PlacementRule::PsiModification(
                                            origin
                                                .split_once(':')
                                                .expect(
                                                    "Incorrect psi mod id, should contain a colon",
                                                )
                                                .1
                                                .parse()
                                                .expect(
                                                    "Incorrect psi mod id, should be numerical",
                                                ),
                                            term.unwrap_or(Position::Anywhere),
                                        )],
                                        Vec::new(),
                                        Vec::new(),
                                    ));
                                }
                            }
                        }
                        if (origins.is_empty() || all_aminoacids)
                            && let Some(term) = term
                        {
                            rules.push((
                                vec![PlacementRule::Terminal(term)],
                                Vec::new(),
                                Vec::new(),
                            ));
                        }
                        modification.data = ModData::Mod {
                            specificities: rules,
                        };
                        mods.push(modification.into());
                    }

                    mods
                })
            })
    }
}
