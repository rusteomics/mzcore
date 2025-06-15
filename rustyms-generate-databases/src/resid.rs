use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::Path,
};

use super::{ModData, ontology_modification::OntologyModification};

use bincode::config::Configuration;
use roxmltree::*;
use rustyms::{
    chemistry::{MolecularFormula, MultiChemical},
    ontology::{Ontology, OntologyModificationList},
    sequence::{AminoAcid, LinkerSpecificity, PlacementRule, Position},
};

pub(crate) fn build_resid_ontology(out_dir: &Path) {
    let mods = parse_resid();

    let dest_path = Path::new(&out_dir).join("resid.dat");
    let mut file = File::create(dest_path).unwrap();
    let final_mods = mods
        .into_iter()
        .map(OntologyModification::into_mod)
        .collect::<Vec<_>>();
    println!("Found {} RESID modifications", final_mods.len());
    file.write_all(
        &bincode::serde::encode_to_vec::<OntologyModificationList, Configuration>(
            final_mods,
            Configuration::default(),
        )
        .unwrap(),
    )
    .unwrap();
}

fn parse_resid() -> Vec<OntologyModification> {
    let mut buf = String::new();
    let _ = BufReader::new(
        File::open("rustyms-generate-databases/data/RESID-RESIDUES.XML")
            .expect("Could not open RESID xml file 'RESID-RESIDUES.XML'"),
    )
    .read_to_string(&mut buf)
    .expect("Could not read RESID xml file");
    let document = Document::parse_with_options(
        &buf,
        ParsingOptions {
            allow_dtd: true,
            ..Default::default()
        },
    )
    .expect("Invalid xml in RESID xml");
    let mut modifications = Vec::new();

    'entry: for entry in document
        .root()
        .first_child()
        .expect("No Database node in RESID XML")
        .children()
    {
        if entry.has_tag_name("Entry") {
            let mut modification = OntologyModification::default();
            let mut corrections = Vec::new();
            let mut rules = Vec::new();

            modification.ontology = Ontology::Resid;
            modification.id = entry.attribute("id").unwrap()[2..].parse().unwrap();
            for data_block in entry.children() {
                match data_block.tag_name().name() {
                    "Names" => {
                        for name_node in data_block.children() {
                            match name_node.tag_name().name() {
                                "Name" => {
                                    modification.name =
                                        name_node.text().unwrap_or_default().to_string();
                                }
                                "AlternateName" | "SystematicName" => modification
                                    .synonyms
                                    .push(name_node.text().unwrap_or_default().to_string()),
                                "Xref" => {
                                    if let Some((a, b)) =
                                        name_node.text().unwrap_or_default().split_once(':')
                                    {
                                        modification.cross_ids.push((a.to_string(), b.to_string()));
                                    } else {
                                        panic!("Invalid Xref content")
                                    }
                                }
                                tag if tag.trim().is_empty() => (),
                                tag => panic!("RESID: Invalid Name tag: {tag} {name_node:?}"),
                            }
                        }
                    }
                    "FormulaBlock" => {
                        for formula_node in data_block.children() {
                            if formula_node.has_tag_name("Formula") {
                                modification.formula = MolecularFormula::from_resid(
                                    formula_node.text().unwrap_or_default(),
                                    ..,
                                )
                                .unwrap()
                                .to_vec()
                                .pop()
                                .unwrap(); // TODO: handle Multi cases, only used for B and Z (potentially just use those?)
                            }
                        }
                    }
                    "CorrectionBlock" => {
                        for formula_node in data_block.children() {
                            if formula_node.has_tag_name("Formula") {
                                corrections.push((
                                    formula_node.attribute("uids"),
                                    formula_node.attribute("link"),
                                    formula_node.attribute("label"),
                                    MolecularFormula::from_resid_single(
                                        formula_node.text().unwrap_or_default(),
                                        ..,
                                    )
                                    .unwrap(),
                                ));
                            }
                        }
                    } // Fixes on specific aas? or mods
                    "ReferenceBlock" => {
                        for ref_node in data_block.children() {
                            if ref_node.has_tag_name("Xref") {
                                if let Some((a, b)) =
                                    ref_node.text().unwrap_or_default().split_once(':')
                                {
                                    modification.cross_ids.push((a.to_string(), b.to_string()));
                                } else {
                                    panic!("Invalid Xref content")
                                }
                            }
                        }
                    }
                    "Comment" => modification.description += data_block.text().unwrap_or_default(),
                    "SequenceCode" => {
                        let mut rule =
                            (AminoAcid::Alanine, None, None, data_block.attribute("link"));
                        for rule_node in data_block.children() {
                            match rule_node.tag_name().name() {
                                "SequenceSpec" => {
                                    let txt = rule_node.text().unwrap_or_default();
                                    if let Some((a, b)) = txt.split_once(", ") {
                                        if b.contains(',') {
                                            println!(
                                                "RESID: Ignore SequenceSpec '{txt}' for {}",
                                                modification.id
                                            );
                                            continue 'entry; // Ignore any cross-link > 2
                                        }
                                        rule.0 = AminoAcid::try_from(a.trim())
                                            .unwrap_or_else(|()| panic!("Invalid AA: {a}"));
                                        rule.1 = Some(
                                            AminoAcid::try_from(b.trim())
                                                .unwrap_or_else(|()| panic!("Invalid AA: {b}")),
                                        );
                                    } else {
                                        rule.0 = AminoAcid::try_from(txt.trim())
                                            .unwrap_or_else(|()| panic!("Invalid AA: {txt}"));
                                    }
                                }
                                "Condition" => match rule_node.text().unwrap_or_default() {
                                    "amino-terminal" => rule.2 = Some(Position::AnyNTerm),
                                    "carboxyl-terminal" => rule.2 = Some(Position::AnyCTerm),
                                    "carboxamidine" => (),
                                    text if text.starts_with("cross-link")
                                        || text.starts_with("incidental")
                                        || text.starts_with("secondary") => {} // Ignore
                                    pos => panic!("Invalid condition position: {pos}"),
                                },
                                "Xref" => {
                                    if let Some((a, b)) =
                                        rule_node.text().unwrap_or_default().split_once(':')
                                    {
                                        modification.cross_ids.push((a.to_string(), b.to_string()));
                                    } else {
                                        panic!("Invalid Xref content")
                                    }
                                }
                                _ => (),
                            }
                        }
                        rules.push(rule);
                    } // Placement rules
                    _ => (),
                }
            }

            let mut shared_formula = None;
            let mut data = None;
            for rule in rules {
                if rule.0 == AminoAcid::AmbiguousAsparagine
                    || rule.0 == AminoAcid::AmbiguousGlutamine
                {
                    println!("RESID: B or Z used as target {}", modification.id);
                    continue 'entry;
                }
                let diff_formula = modification.formula.clone()
                    - rule.0.single_formula().expect("B or Z used as target")
                    - rule
                        .1
                        .map(|a| a.single_formula().expect("B or Z used as target"))
                        .unwrap_or_default();
                if shared_formula.is_some_and(|s| s != diff_formula) {
                    println!(
                        "RESID: Detected multiple diff formulas for {}",
                        modification.id
                    );
                    continue 'entry;
                }
                shared_formula = Some(diff_formula);

                if data.is_none() {
                    if rule.1.is_none() {
                        data = Some(ModData::Mod {
                            specificities: Vec::new(),
                        });
                    } else {
                        data = Some(ModData::Linker {
                            length: None,
                            specificities: Vec::new(),
                        });
                    }
                }
                if let (Some(ModData::Linker { specificities, .. }), Some(aa)) = (&mut data, rule.1)
                {
                    if rule.0 == aa {
                        specificities.push(LinkerSpecificity::Symmetric {
                            rules: vec![PlacementRule::AminoAcid(
                                vec![rule.0],
                                rule.2.unwrap_or(Position::Anywhere),
                            )],
                            stubs: Vec::new(),
                            neutral_losses: Vec::new(),
                            diagnostic: Vec::new(),
                        });
                    } else {
                        specificities.push(LinkerSpecificity::Asymmetric {
                            rules: (
                                vec![PlacementRule::AminoAcid(
                                    vec![rule.0],
                                    rule.2.unwrap_or(Position::Anywhere),
                                )],
                                vec![PlacementRule::AminoAcid(
                                    vec![aa],
                                    rule.2.unwrap_or(Position::Anywhere),
                                )],
                            ),
                            stubs: Vec::new(),
                            neutral_losses: Vec::new(),
                            diagnostic: Vec::new(),
                        });
                    }
                } else if let (Some(ModData::Mod { specificities }), None) = (&mut data, rule.1) {
                    specificities.push((
                        vec![PlacementRule::AminoAcid(
                            vec![rule.0],
                            rule.2.unwrap_or(Position::Anywhere),
                        )],
                        Vec::new(),
                        Vec::new(),
                    ));
                } else {
                    println!(
                        "RESID: Modification is both cross-linker and normal modification {}",
                        modification.id
                    );
                    continue 'entry;
                }
            }

            modification.data = data.unwrap_or_default();
            modifications.push(modification);
        }
    }

    modifications
}
