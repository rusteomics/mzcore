//! Code to handle the XL-MOD ontology
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};
use itertools::Itertools;
use std::{collections::HashMap, sync::Arc};
use thin_vec::ThinVec;

use mzcv::{
    CVError, CVFile, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType, OboValue,
};

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        LinkerLength, LinkerSpecificity, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
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
                    let mut ignored: HashMap<Option<u8>, Vec<u32>> = HashMap::new();

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

                        let mut dna_linker = false;
                        let mut sites = None;
                        let mut length = LinkerLength::Unknown;
                        let mut mass = None;
                        let mut formula = None;
                        let mut origins = (Vec::new(), Vec::new());
                        let mut diagnostic_ions = Vec::new();
                        let mut neutral_losses = Vec::new();
                        let mut stubs = Vec::new();
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

                        for (property, value) in &obj.property_values {
                            match property.as_ref() {
                                "reactionSites" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 reactionSites definition for {id}")
                                    }
                                    sites = if let OboValue::Integer(n) = value[0].0 {
                                        Some(u8::try_from(n).unwrap())
                                    } else {
                                        dbg!(obj);
                                        unreachable!()
                                    }
                                }
                                "spacerLength" => {
                                    for (def, _, _) in value {
                                        match (&mut length, def) {
                                            (LinkerLength::Discreet(options), OboValue::Float(n)) => options.push((*n).into()),
                                            (l, OboValue::Float(n)) => *l = LinkerLength::Discreet(vec![(*n).into()]),
                                            _ => unreachable!(),
                                        }
                                    }
                                }
                                "minSpacerLength" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 minSpacerLength definition for {id}")
                                    }
                                    match (&mut length, &value[0].0) {
                                        (LinkerLength::InclusiveRange(start, _), OboValue::Float(n)) => *start = (*n).into(),
                                        (l, OboValue::Float(n)) => *l = LinkerLength::InclusiveRange((*n).into(),(*n).into(),),
                                        _ => unreachable!(),
                                    }
                                }
                                "maxSpacerLength" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 maxSpacerLength definition for {id}")
                                    }
                                    match (&mut length, &value[0].0) {
                                        (LinkerLength::InclusiveRange(_, end), OboValue::Float(n)) => *end = (*n).into(),
                                        (l, OboValue::Float(n)) => *l = LinkerLength::InclusiveRange((*n).into(),(*n).into(),),
                                        _ => unreachable!(),
                                    }
                                }
                                "monoIsotopicMass" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 monoIsotopicMass definition for {id}")
                                    }
                                    mass = if let OboValue::Float(n) = value[0].0 {
                                        Some(ordered_float::OrderedFloat(n))
                                    } else {
                                        dbg!(obj);
                                        unreachable!()
                                    }
                                }
                                "deadEndFormula" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 deadEndFormula definition for {id}")
                                    }
                                    sites = Some(1);
                                    formula = Some(
                                        MolecularFormula::xlmod(&value[0].0.to_string()).unwrap(),
                                    );
                                }
                                "neutralLossFormula" => {
                                    for (def, _, _) in value {
                                        neutral_losses.push(
                                            MolecularFormula::xlmod(&def.to_string())
                                                .unwrap()
                                                .into(),
                                        );
                                    }
                                }
                                "bridgeFormula" => {
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 bridgeFormula definition for {id}")
                                    }
                                    sites = sites.or(Some(2));
                                    formula = Some(
                                        MolecularFormula::xlmod(&value[0].0.to_string()).unwrap(),
                                    );
                                }
                                "specificities" => {
                                    // specificities: "(C,U)" xsd:string
                                    // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 specificities definition for {id}")
                                    }
                                    if let Some((l, r)) = value[0].0.to_string().split_once('&') {
                                        sites = sites.or(Some(2));
                                        origins.0.extend(
                                            l.trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string()),
                                        );
                                        origins.1.extend(
                                            r.trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string()),
                                        );
                                    } else {
                                        origins.0.extend(
                                            value[0]
                                                .0
                                                .to_string()
                                                .trim_matches(['(', ')'])
                                                .split(',')
                                                .map(|s| s.trim().to_string()),
                                        );
                                    }
                                }
                                "secondarySpecificities" => {
                                    // secondarySpecificities: "(S,T,Y)" xsd:string
                                    if value.len() > 1 {
                                        println!("XLMOD, more than 1 secondarySpecificities definition for {id}")
                                    }
                                    origins.0.extend(
                                        value[0]
                                            .0
                                            .to_string()
                                            .trim_matches(['(', ')'])
                                            .split(',')
                                            .map(|s| s.trim().to_string()),
                                    );
                                }
                                "baseSpecificities" | "secondarybaseSpecificities" => {
                                    dna_linker = true;
                                }
                                "reporterMass" | "CID_Fragment" => {
                                    // reporterMass: "555.2481" xsd:double
                                    // CID_Fragment: "828.5" xsd:double
                                    for (def, _, _) in value {
                                        diagnostic_ions.push(DiagnosticIon(
                                            MolecularFormula::with_additional_mass(
                                                if let OboValue::Float(n) = def {
                                                    *n
                                                } else {
                                                    dbg!(obj);
                                                    unreachable!()
                                                },
                                            ),
                                        ));
                                    }
                                }
                                "stubDefinition" | "stubFormula" => {
                                    // CID = H4 C3 O2 S1, -H2 -O1 : H2 C3 O1
                                    // ETD = -H1 :
                                    // TODO: Extend the stub logic to handle the techniques and neutral losses
                                    for (def, _, _) in value {
                                        if let OboValue::String(definition) = &def {
                                            let (techniques, definition) =
                                                definition.split_once('=').unwrap();
                                            let (first, second) =
                                                definition.split_once(':').unwrap();
                                            let mut split_first = first.split(',');
                                            let mut split_second = second.split(',');
                                            let formula_first = split_first
                                                .next()
                                                .map_or(MolecularFormula::default(), |v| {
                                                    MolecularFormula::xlmod(v).unwrap()
                                                });
                                            let losses_first: Vec<NeutralLoss> = split_first
                                                .map(|v| MolecularFormula::xlmod(v).unwrap().into())
                                                .collect();
                                            let formula_second = split_second
                                                .next()
                                                .map_or(MolecularFormula::default(), |v| {
                                                    MolecularFormula::xlmod(v).unwrap()
                                                });
                                            let losses_second: Vec<NeutralLoss> = split_second
                                                .map(|v| MolecularFormula::xlmod(v).unwrap().into())
                                                .collect();
                                            stubs.push((formula_first, formula_second));
                                        } else {
                                            dbg!(obj);
                                            unreachable!()
                                        }
                                    }
                                }
                                _ => {} // TODO: handle: hydrophilicPEGchain?, maxAbsorption?, waveLengthRange?
                            }
                        }
                        let mut origins = (
                            read_placement_rules(&origins.0),
                            read_placement_rules(&origins.1),
                        );
                        if origins.0.is_empty() {
                            origins.0 = vec![PlacementRule::Anywhere];
                        }
                        // Ignore the mass if a formula is set
                        formula =
                            formula.or(mass.map(|v| MolecularFormula::with_additional_mass(v.0)));

                        if sites == Some(2) {
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
                                                stubs,
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
                        } else if sites == Some(1) {
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
                                            neutral_losses,
                                            diagnostic_ions,
                                        )],
                                    },
                                }
                                .into(),
                            ));
                        } else if !dna_linker {
                            ignored.entry(sites).or_default().push(id);
                        }
                    }

                    for (sites, ids) in ignored.into_iter().sorted() {
                        println!(
                            "XLMOD: Ignored sites {}: {}",
                            sites.map_or("[not set]".to_string(), |v| v.to_string()),
                            ids.into_iter().sorted().join(","),
                        );
                    }

                    mods
                })
            })
    }
}

fn read_placement_rules(bricks: &[String]) -> Vec<PlacementRule> {
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
