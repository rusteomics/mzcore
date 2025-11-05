//! Code to handle the Unimod ontology
use std::io::Read;

use context_error::{BoxedError, CreateError};

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader};

use roxmltree::*;

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification, position_from_str},
    },
    sequence::{PlacementRule, SimpleModificationInner},
};

pub struct Unimod {}

impl CVSource for Unimod {
    type Data = SimpleModificationInner;
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
                    include_bytes!("../databases/unimod.dat"),
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
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
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
        let mut modifications = Vec::new();

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
            modifications.into_iter(),
        ))
    }
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
    synonyms.push(full_name);
    for child in node.children() {
        if child.has_tag_name("specificity") {
            let site = child.attribute("site").unwrap(); // Check if there is a way to recognise linkers
            let position = position_from_str(child.attribute("position").unwrap())?;
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
                        MolecularFormula::from_unimod(loss.attribute("composition").unwrap(), ..)
                            .map_err(|e| e.to_string())
                            .expect("Invalid composition in neutral loss"),
                    )
                })
                .collect();
            rules.push((rule, losses));
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
            synonyms.push(child.text().ok_or("Missing text for synonym")?.to_string());
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
                .map(ToString::to_string)
                .unwrap_or_default();
            if url.is_empty() {
                cross_ids.push((Some(source), text));
            } else {
                cross_ids.push((
                    if source == "Misc. URL" {
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
