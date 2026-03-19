//! Code to handle the Unimod ontology
use std::io::Read;

use context_error::{BoxedError, Context, CreateError, FullErrorContent};
use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader, SynonymScope};
use roxmltree::*;
use thin_vec::ThinVec;

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    helper_functions::explain_number_error,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    prelude::AminoAcid,
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
    ) -> Result<
        (
            CVVersion,
            Self::Structure,
            Vec<BoxedError<'static, CVError>>,
        ),
        Vec<BoxedError<'static, CVError>>,
    > {
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
        .map_err(|err| {
            vec![BoxedError::small(
                CVError::FileCouldNotBeParsed,
                "Invalid xml in Unimod xml",
                err.to_string(),
            )]
        })?;

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
        Ok((
            CVVersion {
                last_updated: None,
                version: None,
                hash: reader.hash(),
            },
            modifications,
            errors,
        ))
    }
}

/// Parse a XML node as a Unimod modification
/// # Errors
/// If it contains invalid data or is missing required data
fn parse_mod(node: &Node) -> Result<OntologyModification, BoxedError<'static, CVError>> {
    let mut formula = MolecularFormula::default();
    let mut diagnostics = Vec::new();
    let mut rules = Vec::new();
    let mut cross_ids = Vec::new();
    let mut synonyms = Vec::new();
    let mut description = node.attribute("full_name").map_or_else(
        || "".into(),
        |full_name| {
            synonyms.push((SynonymScope::Exact, full_name.into()));
            full_name.into()
        },
    );
    for child in node.children() {
        if child.has_tag_name("specificity") {
            let site = child.attribute("site").ok_or_else(|| {
                BoxedError::new(
                    CVError::ItemError,
                    "No defined site for modification",
                    "A Unimod modification must have a site set with 'site'",
                    Context::default().lines(0, format!("Byte range: {:?}", node.range())),
                )
            })?; // Check if there is a way to recognise linkers
            let position = child
                .attribute("position")
                .map_or(Ok(crate::sequence::Position::Anywhere), |p| {
                    p.parse().map_err(|()| BoxedError::new(CVError::ItemError, "Invalid position", "Position should be one of Anywhere, Any N-term, Protein N-term, Any C-term, or Protein C-term", Context::default().lines(0, p).to_owned()))
                })?;
            let rule = match (site, position) {
                ("C-term" | "N-term", pos) => PlacementRule::Terminal(pos),
                (aa, pos) => PlacementRule::AminoAcid(
                    aa.chars()
                        .map(|c| {
                            AminoAcid::try_from(c).map_err(|()| {
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Invalid amino acid",
                                    "Use any valid amino acid single character code",
                                    Context::default().lines(0, c.to_string()).to_owned(),
                                )
                            })
                        })
                        .collect::<Result<ThinVec<AminoAcid>, BoxedError<'static, CVError>>>()?,
                    pos,
                ),
            };
            let losses = child
                .children()
                .filter(|n| {
                    n.has_tag_name("NeutralLoss") && n.attribute("composition") != Some("0")
                })
                .map(|loss| {
                    Ok(NeutralLoss::Loss(
                        1,
                        MolecularFormula::unimod(loss.attribute("composition").ok_or_else(
                            || {
                                BoxedError::new(
                                    CVError::ItemError,
                                    "No defined composition for loss",
                                    "A Unimod loss must have a composition set with 'composition'",
                                    Context::default()
                                        .lines(0, format!("Byte range: {:?}", loss.range())),
                                )
                            },
                        )?)
                        .map_err(|e| {
                            e.to_owned()
                                .convert::<CVError, BoxedError<'static, CVError>>(|_| {
                                    CVError::ItemError
                                })
                        })?,
                    ))
                })
                .collect::<Result<Vec<NeutralLoss>, BoxedError<'static, CVError>>>()?;
            rules.push((rule, losses));
        }
        if child.has_tag_name("delta")
            && let Some(composition) = child.attribute("composition")
        {
            formula = MolecularFormula::unimod(composition).map_err(|e| {
                e.to_owned()
                    .convert::<CVError, BoxedError<'static, CVError>>(|_| CVError::ItemError)
            })?;
        }
        if child.has_tag_name("Ignore")
            && let Some(composition) = child.attribute("composition")
        {
            diagnostics.push(DiagnosticIon(
                MolecularFormula::unimod(composition).map_err(|e| {
                    e.to_owned()
                        .convert::<CVError, BoxedError<'static, CVError>>(|_| CVError::ItemError)
                })?,
            ));
        }
        if child.has_tag_name("misc_notes")
            && let Some(text) = child.text()
            && !text.is_empty()
        {
            description = text.into();
        }
        if child.has_tag_name("alt_name")
            && let Some(text) = child.text()
            && !text.is_empty()
        {
            synonyms.push((SynonymScope::Exact, text.into()));
        }
        if child.has_tag_name("xref") {
            let source: Option<Box<str>> = child
                .children()
                .find(|c| c.has_tag_name("source"))
                .and_then(|c| c.text())
                .filter(|t| !t.is_empty())
                .map(Into::into);
            let text: Option<Box<str>> = child
                .children()
                .find(|c| c.has_tag_name("text"))
                .and_then(|c| c.text())
                .filter(|t| !t.is_empty())
                .map(Into::into);
            let url: Option<Box<str>> = child
                .children()
                .find(|c| c.has_tag_name("url"))
                .and_then(|c| c.text())
                .filter(|t| !t.is_empty())
                .map(Into::into);
            if let Some(url) = url {
                cross_ids.push((
                    if source.as_deref() == Some("Misc. URL") {
                        text
                    } else {
                        source
                    },
                    url,
                ));
            } else if let Some(text) = text {
                cross_ids.push((source, text));
            }
        }
    }
    Ok(OntologyModification {
        name: node.attribute("title").map(Into::into).ok_or_else(|| {
            BoxedError::new(
                CVError::ItemError,
                "No defined name for modification",
                "A Unimod modification must have a name set with 'title'",
                Context::default().lines(0, format!("Byte range: {:?}", node.range())),
            )
        })?,
        id: node
            .attribute("record_id")
            .ok_or_else(|| {
                BoxedError::new(
                    CVError::ItemError,
                    "No defined ID for modification",
                    "A Unimod modification must have an ID set with 'record_id'",
                    Context::default().lines(0, format!("Byte range: {:?}", node.range())),
                )
            })
            .and_then(|v| {
                v.parse::<u32>().map_err(|err| {
                    BoxedError::new(
                        CVError::ItemError,
                        "Modification ID not numeric",
                        format!("The modification ID {}", explain_number_error(&err)),
                        Context::default().lines(0, format!("Byte range: {:?}", node.range())),
                    )
                })
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
        obsolete: false,
    })
}
