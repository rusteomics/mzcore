use std::{io::Write, path::Path};

use bincode::config::Configuration;
use itertools::Itertools;
use rustyms::{
    chemistry::MolecularFormula,
    fragment::DiagnosticIon,
    ontology::{Ontology, OntologyModificationList},
    sequence::{LinkerSpecificity, PlacementRule, Position},
};
use thin_vec::ThinVec;

use super::{
    ModData,
    obo::{OboOntology, OboValue},
    ontology_modification::OntologyModification,
};

pub(crate) fn build_xlmod_ontology(out_dir: &Path) {
    let mods = parse_xlmod();

    let mut mods_file = std::fs::File::create(Path::new(&out_dir).join("xlmod.dat")).unwrap();
    let final_mods = mods
        .into_iter()
        .map(OntologyModification::into_mod)
        .sorted_unstable()
        .collect::<Vec<_>>();
    println!("Found {} XLMOD modifications", final_mods.len());
    mods_file
        .write_all(
            &bincode::serde::encode_to_vec::<OntologyModificationList, Configuration>(
                final_mods,
                Configuration::default(),
            )
            .unwrap(),
        )
        .unwrap();
}

fn parse_xlmod() -> Vec<OntologyModification> {
    let obo = OboOntology::from_file("rustyms-generate-databases/data/XLMOD.obo")
        .expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let id = obj.lines["id"][0]
            .split_once(':')
            .expect("Incorrect XLMOD id, should contain a colon")
            .1
            .parse()
            .expect("Incorrect XLMOD id, should be numerical");
        let name = obj.lines["name"][0].to_string();

        let mut sites = None;
        let mut length = None;
        let mut mass = None;
        let mut formula = None;
        let mut origins = (Vec::new(), Vec::new());
        let mut diagnostic_ions = Vec::new();
        let mut description = String::new();
        let mut cross_ids = ThinVec::new();
        let mut synonyms = ThinVec::new();
        if let Some(values) = obj.lines.get("def") {
            assert!(values.len() == 1);
            let line = values[0][1..].split_once('\"').unwrap();
            description = line.0.to_string();
            let ids = line.1.trim();
            cross_ids = ids[1..ids.len() - 1]
                .split(',')
                .map(|id| id.trim().split_once(':').unwrap())
                .map(|(r, i)| (r.to_string(), i.to_string()))
                .collect();
        }
        if let Some(values) = obj.lines.get("synonym") {
            for line in values {
                let line = line[1..].split_once('\"').unwrap();
                synonyms.push(line.0.to_string());
            }
        }
        for (id, value) in &obj.property_values {
            match id.as_str() {
                "reactionSites" => {
                    sites = if let OboValue::Integer(n) = value[0] {
                        Some(u8::try_from(n).unwrap())
                    } else {
                        dbg!(obj);
                        unreachable!()
                    }
                }
                "spacerLength" => {
                    length = if let OboValue::Float(n) = value[0] {
                        Some(ordered_float::OrderedFloat(n))
                    } else {
                        None // can contain ranges
                    }
                }
                "monoIsotopicMass" => {
                    mass = if let OboValue::Float(n) = value[0] {
                        Some(ordered_float::OrderedFloat(n))
                    } else {
                        dbg!(obj);
                        unreachable!()
                    }
                }
                "deadEndFormula" => {
                    sites = Some(1);
                    formula =
                        Some(MolecularFormula::from_xlmod(&value[0].to_string(), ..).unwrap());
                }
                "bridgeFormula" => {
                    sites = Some(2);
                    formula =
                        Some(MolecularFormula::from_xlmod(&value[0].to_string(), ..).unwrap());
                }
                "specificities" => {
                    // specificities: "(C,U)" xsd:string
                    // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                    if let Some((l, r)) = value[0].to_string().split_once('&') {
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
                            .to_string()
                            .trim_matches(['(', ')'])
                            .split(',')
                            .map(|s| s.trim().to_string()),
                    );
                }
                "reporterMass" | "CID_Fragment" => {
                    // reporterMass: "555.2481" xsd:double
                    // CID_Fragment: "828.5" xsd:double
                    diagnostic_ions.push(DiagnosticIon(MolecularFormula::with_additional_mass(
                        if let OboValue::Float(n) = value[0] {
                            n
                        } else {
                            dbg!(obj);
                            unreachable!()
                        },
                    )));
                }
                _ => {}
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
            mods.push(OntologyModification {
                formula: formula.unwrap_or_default(),
                name,
                description,
                cross_ids,
                synonyms,
                id,
                ontology: Ontology::Xlmod,
                data: ModData::Linker {
                    length,
                    specificities: vec![if origins.1.is_empty() {
                        LinkerSpecificity::Symmetric {
                            rules: origins.0,
                            stubs: Vec::new(),
                            neutral_losses: Vec::new(),
                            diagnostic: diagnostic_ions,
                        }
                    } else {
                        LinkerSpecificity::Asymmetric {
                            rules: (origins.0, origins.1),
                            stubs: Vec::new(),
                            neutral_losses: Vec::new(),
                            diagnostic: diagnostic_ions,
                        }
                    }],
                },
            });
        } else if sites == Some(3) {
            continue; // Ignore
        } else {
            mods.push(OntologyModification {
                formula: formula.unwrap_or_default(),
                name,
                description,
                cross_ids,
                synonyms,
                ontology: Ontology::Xlmod,
                id,
                data: ModData::Mod {
                    specificities: vec![(origins.0, Vec::new(), diagnostic_ions)],
                },
            });
        }
    }

    mods
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
                        vec![brick.try_into().unwrap()],
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
            .collect_vec()
    }
}
