//! Code to handle the Unimod ontology
use std::io::Read;

use context_error::{BoxedError, CreateError};
use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader, SynonymScope};
use roxmltree::*;

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{PlacementRule, SimpleModification, SimpleModificationInner},
};

/// Unimod modifications
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct Unimod {}

impl CVSource for Unimod {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "Unimod"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "Unimod",
            extension: "xml",
            url: Some("https://unimod.org/xml/unimod.xml"),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/unimod.dat"),
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
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, Self::Structure), Vec<BoxedError<'static, CVError>>> {
        let mut reader = reader.next().unwrap();
        let mut buf = String::new();
        reader.read_to_string(&mut buf).map_err(|e| {
            vec![BoxedError::small(
                CVError::FileCouldNotBeOpenend,
                "Could not read file",
                e.to_string(),
            )]
        })?;
        let document = Document::parse_with_options(
            &buf,
            ParsingOptions {
                allow_dtd: true,
                ..Default::default()
            },
        )
        .expect("Invalid xml in Unimod xml");

        let mut errors = Vec::new();
        let mut modifications: Vec<SimpleModification> = Vec::new();

        for node in document.root().children() {
            if node.has_tag_name("unimod") {
                for node in node.children() {
                    if node.has_tag_name("modifications") {
                        for node in node.children() {
                            if node.has_tag_name("mod") {
                                match parse_mod(&node) {
                                    Ok(o) => modifications.push(o.into()),
                                    Err(e) => errors.push(e),
                                }
                            }
                        }
                    }
                }
            }
        }
        if !errors.is_empty() {
            for err in &errors {
                println!("{err}");
            }
            panic!("{} Errors found while parsing Unimod", errors.len());
        }

        Ok((
            CVVersion {
                last_updated: None,
                version: None,
                hash: reader.hash().to_vec(),
            },
            modifications,
        ))
    }
}

fn parse_mod(node: &Node) -> Result<OntologyModification, String> {
    let mut formula = MolecularFormula::default();
    let full_name: Box<str> = node
        .attribute("full_name")
        .map(Into::into)
        .ok_or("No defined description for modification")?;
    let mut diagnostics = Vec::new();
    let mut rules = Vec::new();
    let mut cross_ids = Vec::new();
    let mut synonyms = Vec::new();
    let mut description = full_name.clone();
    synonyms.push((SynonymScope::Exact, full_name));
    for child in node.children() {
        if child.has_tag_name("specificity") {
            let site = child.attribute("site").unwrap(); // Check if there is a way to recognise linkers
            let position = child
                .attribute("position")
                .map_or(Ok(crate::sequence::Position::Anywhere), |p| {
                    p.parse().map_err(|()| format!("Invalid position '{p}'"))
                })?;
            let rule = match (site, position) {
                ("C-term" | "N-term", pos) => PlacementRule::Terminal(pos),
                (aa, pos) => PlacementRule::AminoAcid(
                    aa.chars()
                        .map(|c| c.try_into().unwrap_or_else(|()| panic!("Not an AA: {c}")))
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
                        1,
                        MolecularFormula::unimod(loss.attribute("composition").unwrap())
                            .map_err(|e| e.to_string())
                            .expect("Invalid composition in neutral loss"),
                    )
                })
                .collect();
            rules.push((rule, losses));
        }
        if child.has_tag_name("delta")
            && let Some(composition) = child.attribute("composition")
        {
            formula = MolecularFormula::unimod(composition).map_err(|e| e.to_string())?;
        }
        if child.has_tag_name("Ignore")
            && let Some(composition) = child.attribute("composition")
        {
            diagnostics.push(DiagnosticIon(
                MolecularFormula::unimod(composition).map_err(|e| e.to_string())?,
            ));
        }
        if child.has_tag_name("misc_notes") {
            description = child
                .text()
                .ok_or("Missing text for notes in modification")?
                .into();
        }
        if child.has_tag_name("alt_name") {
            synonyms.push((
                SynonymScope::Exact,
                child.text().ok_or("Missing text for synonym")?.into(),
            ));
        }
        if child.has_tag_name("xref") {
            let source: Box<str> = child
                .children()
                .find(|c| c.has_tag_name("source"))
                .and_then(|c| c.text())
                .unwrap()
                .into();
            let text: Box<str> = child
                .children()
                .find(|c| c.has_tag_name("text"))
                .and_then(|c| c.text())
                .unwrap()
                .into();
            let url: Box<str> = child
                .children()
                .find(|c| c.has_tag_name("url"))
                .and_then(|c| c.text())
                .map(Into::into)
                .unwrap_or_default();
            if url.is_empty() {
                cross_ids.push((Some(source), text));
            } else {
                cross_ids.push((
                    if *source == *"Misc. URL" {
                        Some(text)
                    } else {
                        Some(source)
                    },
                    url,
                ));
            }
        }
    }
    Ok(OntologyModification {
        name: node
            .attribute("title")
            .map(Into::into)
            .ok_or("No defined name for modification")?,
        id: node
            .attribute("record_id")
            .ok_or("No defined id for modification")
            .and_then(|v| v.parse::<u32>().map_err(|_| "Invalid id for modification"))?,
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
