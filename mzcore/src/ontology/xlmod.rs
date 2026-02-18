//! Code to handle the XL-MOD ontology
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};
use std::sync::Arc;
use thin_vec::ThinVec;

use mzcv::{
    CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType, OboValue,
};

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        LinkerSpecificity, PlacementRule, Position, SimpleModification, SimpleModificationInner,
    },
};

/// XL-MOD modifications
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct XlMod {}

impl CVSource for XlMod {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "XLMOD"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "XLMOD",
            extension: "obo",
            url: Some(
                "https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/refs/heads/main/XLMOD.obo",
            ),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/xlmod.dat"),
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
                        if obj.stanza_type != OboStanzaType::Term {
                            continue;
                        }
                        let id = obj
                            .id
                            .1
                            .parse()
                            .expect("Incorrect XLMOD id, should be numerical");
                        let name = obj.lines["name"][0].0.clone();

                        let mut sites = None;
                        let mut length = None;
                        let mut mass = None;
                        let mut formula = None;
                        let mut origins = (Vec::new(), Vec::new());
                        let mut diagnostic_ions = Vec::new();
                        let mut neutral_losses = Vec::new();
                        let description = obj
                            .definition
                            .as_ref()
                            .map_or_else(Box::default, |d| d.0.clone());
                        let cross_ids = obj
                            .definition
                            .as_ref()
                            .map_or_else(ThinVec::new, |d| d.1.clone().into());
                        let synonyms = obj
                            .synonyms
                            .iter()
                            .map(|s| (s.scope, s.synonym.clone()))
                            .collect();

                        for (id, value) in &obj.property_values {
                            match id.as_ref() {
                                "reactionSites" => {
                                    sites = if let OboValue::Integer(n) = value[0].0 {
                                        Some(u8::try_from(n).unwrap())
                                    } else {
                                        dbg!(obj);
                                        unreachable!()
                                    }
                                }
                                "spacerLength" => {
                                    length = match &value[0].0 {
                                        OboValue::Float(n) => Some(
                                            ordered_float::OrderedFloat(*n)
                                                ..=ordered_float::OrderedFloat(*n),
                                        ),
                                        OboValue::FloatRange(n) => Some(
                                            ordered_float::OrderedFloat(*n.start())
                                                ..=ordered_float::OrderedFloat(*n.end()),
                                        ),
                                        _ => unreachable!(),
                                    }
                                }
                                "minSpacerLength" => {
                                    length = match &value[0].0 {
                                        OboValue::Float(n) => Some(
                                            ordered_float::OrderedFloat(*n)
                                                ..=length
                                                    .map_or(ordered_float::OrderedFloat(*n), |v| {
                                                        *v.end()
                                                    }),
                                        ),
                                        _ => unreachable!(),
                                    }
                                }
                                "maxSpacerLength" => {
                                    length = match &value[0].0 {
                                        OboValue::Float(n) => Some(
                                            length.map_or(ordered_float::OrderedFloat(*n), |v| {
                                                *v.start()
                                            })
                                                ..=ordered_float::OrderedFloat(*n),
                                        ),
                                        _ => unreachable!(),
                                    }
                                }
                                "monoIsotopicMass" => {
                                    mass = if let OboValue::Float(n) = value[0].0 {
                                        Some(ordered_float::OrderedFloat(n))
                                    } else {
                                        dbg!(obj);
                                        unreachable!()
                                    }
                                }
                                "deadEndFormula" => {
                                    sites = Some(1);
                                    formula = Some(
                                        MolecularFormula::xlmod(&value[0].0.to_string()).unwrap(),
                                    );
                                }
                                "neutralLossFormula" => {
                                    neutral_losses.push(
                                        MolecularFormula::xlmod(&value[0].0.to_string())
                                            .unwrap()
                                            .into(),
                                    );
                                }
                                "bridgeFormula" => {
                                    sites = Some(2);
                                    formula = Some(
                                        MolecularFormula::xlmod(&value[0].0.to_string()).unwrap(),
                                    );
                                }
                                "specificities" => {
                                    // specificities: "(C,U)" xsd:string
                                    // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                                    if let Some((l, r)) = value[0].0.to_string().split_once('&') {
                                        origins = (
                                            l.trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string())
                                                .collect(),
                                            r.trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string())
                                                .collect(),
                                        );
                                    } else {
                                        origins = (
                                            value[0]
                                                .0
                                                .to_string()
                                                .trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string())
                                                .collect(),
                                            Vec::new(),
                                        );
                                    }
                                }
                                "secondarySpecificities" => {
                                    // secondarySpecificities: "(S,T,Y)" xsd:string
                                    origins.0.extend(
                                        value[0]
                                            .0
                                            .to_string()
                                            .trim_matches(['(', ')'])
                                            .split(',')
                                            .map(|s| s.trim().to_string()),
                                    );
                                }
                                "reporterMass" | "CID_Fragment" => {
                                    // reporterMass: "555.2481" xsd:double
                                    // CID_Fragment: "828.5" xsd:double
                                    diagnostic_ions.push(DiagnosticIon(
                                        MolecularFormula::with_additional_mass(
                                            if let OboValue::Float(n) = value[0].0 {
                                                n
                                            } else {
                                                dbg!(obj);
                                                unreachable!()
                                            },
                                        ),
                                    ));
                                }
                                _ => {} // TODO: handle: stubDefinition, hydrophilicPEGchain?, maxAbsorption?, waveLengthRange?
                            }
                        }
                        let origins = (
                            read_placement_rules(&origins.0),
                            read_placement_rules(&origins.1),
                        );
                        if let Some(mass) = mass {
                            // Ignore the mass if a formula is set
                            if formula.is_none() {
                                formula = Some(MolecularFormula::with_additional_mass(mass.0));
                            }
                        }
                        if sites == Some(2) || !origins.1.is_empty() {
                            mods.push(
                                OntologyModification {
                                    formula: formula.unwrap_or_default(),
                                    name,
                                    description,
                                    cross_ids,
                                    synonyms,
                                    id,
                                    ontology: Ontology::Xlmod,
                                    obsolete: obj.obsolete,
                                    data: ModData::Linker {
                                        length,
                                        specificities: vec![if origins.1.is_empty() {
                                            LinkerSpecificity::Symmetric {
                                                rules: origins.0,
                                                stubs: Vec::new(),
                                                neutral_losses,
                                                diagnostic: diagnostic_ions,
                                            }
                                        } else {
                                            LinkerSpecificity::Asymmetric {
                                                rules: (origins.0, origins.1),
                                                stubs: Vec::new(),
                                                neutral_losses,
                                                diagnostic: diagnostic_ions,
                                            }
                                        }],
                                    },
                                }
                                .into(),
                            );
                        } else if sites == Some(3) {
                            // Ignore
                        } else {
                            mods.push(Arc::new(
                                OntologyModification {
                                    formula: formula.unwrap_or_default(),
                                    name,
                                    description,
                                    cross_ids,
                                    synonyms,
                                    ontology: Ontology::Xlmod,
                                    obsolete: obj.obsolete,
                                    id,
                                    data: ModData::Mod {
                                        specificities: vec![(
                                            origins.0,
                                            Vec::new(),
                                            diagnostic_ions,
                                        )],
                                    },
                                }
                                .into(),
                            ));
                        }
                    }

                    mods
                })
            })
    }
}

fn read_placement_rules(bricks: &[String]) -> Vec<PlacementRule> {
    if bricks.is_empty() {
        vec![PlacementRule::Anywhere]
    } else {
        bricks
            .iter()
            .filter_map(|brick| {
                if brick.len() == 1 {
                    Some(PlacementRule::AminoAcid(
                        vec![brick.try_into().unwrap()].into(),
                        Position::Anywhere,
                    ))
                } else if *brick == "Protein N-term" {
                    Some(PlacementRule::Terminal(Position::ProteinNTerm))
                } else if *brick == "Protein C-term" {
                    Some(PlacementRule::Terminal(Position::ProteinCTerm))
                } else if ["Thy"].contains(&brick.as_str()) {
                    None
                } else {
                    panic!("Invalid placement rule: '{brick}'")
                }
            })
            .collect()
    }
}
