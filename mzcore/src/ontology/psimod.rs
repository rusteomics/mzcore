//! Code to handle the PSI-MOD ontology
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVCacheSerde, CVData, CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology,
    OboStanzaType,
};

use crate::{
    chemistry::MolecularFormula,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{PlacementRule, Position, SimpleModificationInner},
};

pub struct PsiMod {}

impl CVData for SimpleModificationInner {
    type Index = usize;
    fn index(&self) -> Option<usize> {
        self.description().and_then(|d| d.id)
    }
    fn name(&self) -> Option<&str> {
        self.description().map(|d| d.name.as_str())
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.description()
            .into_iter()
            .flat_map(|d| d.synonyms.iter().map(String::as_str))
    }
    type Cache = CVCacheSerde<Self>;
}

impl CVSource for PsiMod {
    type Data = SimpleModificationInner;
    fn cv_name() -> &'static str {
        "MOD"
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
                    include_bytes!("../databases/psimod_new.dat"),
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
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
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
                    let mut mods: Vec<SimpleModificationInner> = Vec::new();

                    for obj in obo.objects {
                        if obj.stanza_type != OboStanzaType::Term
                            || obj.id == (Some("MOD".to_string()), "00000".to_string())
                            || obj.id == (Some("MOD".to_string()), "00004".to_string())
                            || obj.id == (Some("MOD".to_string()), "00008".to_string())
                        {
                            continue;
                        }
                        let mut modification = OntologyModification {
                            id: obj
                                .id
                                .1
                                .parse()
                                .expect("Incorrect psi mod id, should be numerical"),
                            name: obj.lines["name"][0].0.to_string(),
                            ontology: Ontology::Psimod,
                            ..OntologyModification::default()
                        };
                        if let Some((description, cross_ids, _, _)) = obj.definition {
                            modification.description = description;
                            modification.cross_ids = cross_ids.into();
                        }
                        for synonym in obj.synonyms {
                            modification.synonyms.push(synonym.synonym); // TODO: retain the scope
                        }

                        let mut rules = Vec::new();
                        let mut origins = Vec::new();
                        let mut term = None;
                        for (id, values) in obj.property_values {
                            for (value, _, _) in values {
                                match (id.as_str(), value) {
                                    ("DiffFormula", mzcv::OboValue::String(s)) => {
                                        modification.formula =
                                            MolecularFormula::from_psi_mod(&s, ..).unwrap();
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
                                            vec![origin.try_into().unwrap()],
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

                    mods.into_iter()
                })
            })
    }
}
