use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::Path,
};

use bincode::config::Configuration;
use roxmltree::*;
use rustyms::{
    chemistry::MolecularFormula,
    fragment::{DiagnosticIon, NeutralLoss},
    ontology::{Ontology, OntologyModificationList},
    sequence::PlacementRule,
};

use crate::ontology_modification::position_from_str;

use super::{ModData, ontology_modification::OntologyModification};

pub(crate) fn build_unimod_ontology(out_dir: &Path) {
    let mods = parse_unimod();

    let dest_path = Path::new(&out_dir).join("unimod.dat");
    let mut file = File::create(dest_path).unwrap();
    let final_mods = mods.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    println!("Found {} Unimod modifications", final_mods.len());
    file.write_all(
        &bincode::serde::encode_to_vec::<OntologyModificationList, Configuration>(
            final_mods,
            Configuration::default(),
        )
        .unwrap(),
    )
    .unwrap();
}

// TODO: use the other xml file because this one does not contain all mods
// TODO: see if cross-linkers can be recognised better
fn parse_unimod() -> Vec<OntologyModification> {
    let mut buf = String::new();
    BufReader::new(
        File::open("rustyms-generate-databases/data/unimod.xml")
            .expect("Could not open Unimod xml file"),
    )
    .read_to_string(&mut buf)
    .expect("Could not read Unimod xml file fully");
    let document = Document::parse_with_options(
        &buf,
        ParsingOptions {
            allow_dtd: true,
            ..Default::default()
        },
    )
    .expect("Invalid xml in Unimod xml");

    let mut errors = Vec::new();
    let mut modifications = Vec::new();

    for node in document.root().children() {
        if node.has_tag_name("unimod") {
            for node in node.children() {
                if node.has_tag_name("modifications") {
                    for node in node.children() {
                        if node.has_tag_name("mod") {
                            match parse_mod(&node) {
                                Ok(o) => modifications.push(o),
                                Err(e) => errors.push(e),
                            }
                        }
                    }
                }
            }
        }
    }
    modifications
}

fn parse_mod(node: &Node) -> Result<OntologyModification, String> {
    let mut formula = MolecularFormula::default();
    let full_name = node
        .attribute("full_name")
        .map(ToString::to_string)
        .ok_or("No defined description for modification")?;
    let mut diagnostics = Vec::new();
    let mut rules = Vec::new();
    let mut cross_ids = Vec::new();
    let mut synonyms = Vec::new();
    let mut description = full_name.clone();
    synonyms.push(full_name.clone());
    for child in node.children() {
        if child.has_tag_name("specificity") {
            let site = child.attribute("site").unwrap(); // Check if there is a way to recognise linkers
            let position = position_from_str(child.attribute("position").unwrap())?;
            let rule = match (site, position) {
                ("C-term" | "N-term", pos) => PlacementRule::Terminal(pos),
                (aa, pos) => PlacementRule::AminoAcid(
                    aa.chars()
                        .map(|c| c.try_into().unwrap_or_else(|_| panic!("Not an AA: {c}")))
                        .collect(),
                    pos,
                ),
            };
            let losses = child
                .children()
                .filter(|n| {
                    n.has_tag_name("NeutralLoss") && n.attribute("composition") != Some("0")
                })
                .map(|loss| {
                    NeutralLoss::Loss(
                        MolecularFormula::from_unimod(loss.attribute("composition").unwrap(), ..)
                            .map_err(|e| e.to_string())
                            .expect("Invalid composition in neutral loss"),
                    )
                })
                .collect();
            rules.push((rule, losses))
        }
        if child.has_tag_name("delta") {
            if let Some(composition) = child.attribute("composition") {
                formula =
                    MolecularFormula::from_unimod(composition, ..).map_err(|e| e.to_string())?;
            }
        }
        if child.has_tag_name("Ignore") {
            if let Some(composition) = child.attribute("composition") {
                diagnostics.push(DiagnosticIon(
                    MolecularFormula::from_unimod(composition, ..).map_err(|e| e.to_string())?,
                ));
            }
        }
        if child.has_tag_name("misc_notes") {
            description = child
                .text()
                .ok_or("Missing text for notes in modification")?
                .to_string();
        }
        if child.has_tag_name("alt_name") {
            synonyms.push(child.text().ok_or("Missing text for synonym")?.to_string())
        }
        if child.has_tag_name("xref") {
            let source = child
                .children()
                .find(|c| c.has_tag_name("source"))
                .and_then(|c| c.text())
                .unwrap()
                .to_string();
            let text = child
                .children()
                .find(|c| c.has_tag_name("text"))
                .and_then(|c| c.text())
                .unwrap()
                .to_string();
            let url = child
                .children()
                .find(|c| c.has_tag_name("url"))
                .and_then(|c| c.text())
                .map(|s| s.to_string())
                .unwrap_or_default();
            if url.is_empty() {
                cross_ids.push((source, text));
            } else {
                cross_ids.push((if source == "Misc. URL" { text } else { source }, url));
            }
        }
    }
    Ok(OntologyModification {
        name: node
            .attribute("title")
            .map(ToString::to_string)
            .ok_or("No defined name for modification")?,
        id: node
            .attribute("record_id")
            .ok_or("No defined id for modification")
            .and_then(|v| {
                v.parse::<usize>()
                    .map_err(|_| "Invalid id for modification")
            })?,
        ontology: Ontology::Unimod,
        description,
        synonyms: synonyms.into(),
        cross_ids: cross_ids.into(),
        formula,
        data: ModData::Mod {
            specificities: rules
                .into_iter()
                .map(|(r, l)| (vec![r], l, diagnostics.clone()))
                .collect(),
        },
    })
}
