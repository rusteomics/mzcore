//! Code to handle the RESID ontology
use std::io::Read;

use context_error::{BoxedError, Context, CreateError, FullErrorContent, combine_error};

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader};

use roxmltree::*;

use crate::{
    chemistry::MolecularFormula,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        AminoAcid, LinkerSpecificity, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};

/// RESID modifications
///
/// Note that because RESID is stored on a ftp server the automatic downloading does not work.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct Resid {}

impl CVSource for Resid {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "RESID"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "RESID",
            extension: "xml",
            url: None, //FTP is not supported but for when it is: ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/RESIDUES.XML
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/resid.dat"),
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
            vec![BoxedError::new(
                CVError::FileCouldNotBeParsed,
                "Invalid XML",
                "The given RESID file does not contain valid XML",
                Context::default().lines(0, err.to_string()),
            )]
        })?;
        let mut errors = Vec::new();
        let mut modifications: Vec<SimpleModification> = Vec::new();
        let database = document.root().first_child().ok_or_else(|| {
            vec![BoxedError::new(
                CVError::FileCouldNotBeParsed,
                "No root node in RESID XML",
                "A Database node should be present as the root in a RESID XML file",
                Context::default(),
            )]
        })?;

        'entry: for entry in database.children() {
            if entry.has_tag_name("Entry") {
                let mut modification = OntologyModification::default();
                let mut rules = Vec::new();

                modification.ontology = Ontology::Resid;
                let Some(id) = entry.attribute("id").and_then(|id| id.parse().ok()) else {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid ID",
                            "A RESID ID should be present and of shape 'AA0000'",
                            Context::default().byte_range(entry.range()),
                        ),
                    );
                    continue;
                };
                modification.id = id;
                let mut formula_found = false;
                for data_block in entry.children() {
                    match data_block.tag_name().name() {
                        "Names" => {
                            for name_node in data_block.children() {
                                match name_node.tag_name().name() {
                                    "Name" => {
                                        modification.name =
                                            name_node.text().unwrap_or_default().into();
                                    }
                                    "AlternateName" | "SystematicName" => {
                                        modification.synonyms.push((
                                            mzcv::SynonymScope::Exact,
                                            name_node.text().unwrap_or_default().into(),
                                        ));
                                    }
                                    "Xref" => {
                                        if let Some((a, b)) =
                                            name_node.text().unwrap_or_default().split_once(':')
                                        {
                                            modification.cross_ids.push((Some(a.into()), b.into()));
                                        } else {
                                            combine_error(
                                                &mut errors,
                                                BoxedError::new(
                                                    CVError::ItemError,
                                                    "Invalid xref",
                                                    "A RESID xref should be of shape TYPE:VALUE",
                                                    Context::default()
                                                        .byte_range(name_node.range())
                                                        .lines(
                                                            0,
                                                            format!(
                                                                "RESID:AA{:04}",
                                                                modification.id
                                                            ),
                                                        ),
                                                ),
                                            );
                                        }
                                    }
                                    tag if tag.trim().is_empty() => (),
                                    _ => combine_error(
                                        &mut errors,
                                        BoxedError::new(
                                            CVError::ItemError,
                                            "Invalid name",
                                            "A RESID name should be Name, AlternateName, SystemicName, or Xref",
                                            Context::default().byte_range(name_node.range()).lines(
                                                0,
                                                format!("RESID:AA{:04}", modification.id),
                                            ),
                                        ),
                                    ),
                                }
                            }
                        }
                        "CorrectionBlock" => {
                            if formula_found {
                                combine_error(
                                    &mut errors,
                                    BoxedError::new(
                                        CVError::ItemWarning,
                                        "Multiple formulas",
                                        "Modifications with formulas (CorrectionBlock) that differ based on the used amino acid are ignored for now",
                                        Context::default()
                                            .byte_range(data_block.range())
                                            .lines(0, format!("RESID:AA{:04}", modification.id)),
                                    ),
                                );
                                continue 'entry;
                            }
                            for formula_node in data_block.children() {
                                if formula_node.has_tag_name("Formula") {
                                    match MolecularFormula::resid_single(
                                        formula_node.text().unwrap_or_default(),
                                    ) {
                                        Ok(f) => {
                                            formula_found = true;
                                            modification.formula = f;
                                        }
                                        Err(err) => combine_error(
                                            &mut errors,
                                            err.convert::<CVError, BoxedError<'_, CVError>>(|_| {
                                                CVError::ItemError
                                            })
                                            .to_owned(),
                                        ),
                                    }
                                }
                            }
                        }
                        "ReferenceBlock" => {
                            for ref_node in data_block.children() {
                                if ref_node.has_tag_name("Xref") {
                                    if let Some((a, b)) =
                                        ref_node.text().unwrap_or_default().split_once(':')
                                    {
                                        modification.cross_ids.push((Some(a.into()), b.into()));
                                    } else {
                                        combine_error(
                                            &mut errors,
                                            BoxedError::new(
                                                CVError::ItemError,
                                                "Invalid xref",
                                                "A RESID xref should be of shape TYPE:VALUE",
                                                Context::default()
                                                    .byte_range(ref_node.range())
                                                    .lines(
                                                        0,
                                                        format!("RESID:AA{:04}", modification.id),
                                                    ),
                                            ),
                                        );
                                    }
                                }
                            }
                        }
                        "Comment" => {
                            modification.description = format!(
                                "{}{}",
                                modification.description,
                                data_block.text().unwrap_or_default()
                            )
                            .into_boxed_str();
                        }
                        "SequenceCode" => {
                            let mut rule =
                                (AminoAcid::Alanine, None, None, data_block.attribute("link"));
                            for rule_node in data_block.children() {
                                match rule_node.tag_name().name() {
                                    "SequenceSpec" => {
                                        let txt = rule_node.text().unwrap_or_default();
                                        if let Some((a, b)) = txt.split_once(", ") {
                                            if b.contains(',') {
                                                combine_error(
                                                    &mut errors,
                                                    BoxedError::new(
                                                        CVError::ItemWarning,
                                                        "Higher order cross-linker",
                                                        "Only cross-linkers with two sites are supported for now",
                                                        Context::default()
                                                            .byte_range(rule_node.range())
                                                            .lines(
                                                                0,
                                                                format!(
                                                                    "RESID:AA{:04} {txt}",
                                                                    modification.id
                                                                ),
                                                            ),
                                                    ),
                                                );
                                                continue 'entry; // Ignore any cross-link > 2
                                            }
                                            if let Ok(a) = AminoAcid::try_from(a.trim()) {
                                                rule.0 = a;
                                            } else {
                                                combine_error(
                                                    &mut errors,
                                                    BoxedError::new(
                                                        CVError::ItemWarning,
                                                        "Invalid amino acid",
                                                        "Placement rules can only contain amino acids",
                                                        Context::default()
                                                            .byte_range(rule_node.range())
                                                            .lines(0, a.to_string()),
                                                    ),
                                                );
                                                continue 'entry;
                                            }
                                            if let Ok(b) = AminoAcid::try_from(b.trim()) {
                                                rule.1 = Some(b);
                                            } else {
                                                combine_error(
                                                    &mut errors,
                                                    BoxedError::new(
                                                        CVError::ItemWarning,
                                                        "Invalid amino acid",
                                                        "Placement rules can only contain amino acids",
                                                        Context::default()
                                                            .byte_range(rule_node.range())
                                                            .lines(0, b.to_string()),
                                                    ),
                                                );
                                                continue 'entry;
                                            }
                                        } else if let Ok(a) = AminoAcid::try_from(txt.trim()) {
                                            rule.0 = a;
                                        } else {
                                            combine_error(
                                                &mut errors,
                                                BoxedError::new(
                                                    CVError::ItemWarning,
                                                    "Invalid amino acid",
                                                    "Placement rules can only contain amino acids",
                                                    Context::default()
                                                        .byte_range(rule_node.range())
                                                        .lines(0, txt.to_string()),
                                                ),
                                            );
                                            continue 'entry;
                                        }
                                    }
                                    "Condition" => match rule_node.text().unwrap_or_default() {
                                        "amino-terminal" => rule.2 = Some(Position::AnyNTerm),
                                        "carboxyl-terminal" => rule.2 = Some(Position::AnyCTerm),
                                        "carboxamidine" => (),
                                        text if text.starts_with("cross-link")
                                            || text.starts_with("incidental")
                                            || text.starts_with("secondary") => {} // Ignore
                                        _ => combine_error(
                                            &mut errors,
                                            BoxedError::new(
                                                CVError::ItemError,
                                                "Invalid condition position",
                                                "A RESID condition position should be amino-terminal, carboxy-terminal, carboxamidine, cross-link, incidental, or secondary",
                                                Context::default()
                                                    .byte_range(rule_node.range())
                                                    .lines(
                                                        0,
                                                        format!("RESID:AA{:04}", modification.id),
                                                    ),
                                            ),
                                        ),
                                    },
                                    "Xref" => {
                                        if let Some((a, b)) =
                                            rule_node.text().unwrap_or_default().split_once(':')
                                        {
                                            modification.cross_ids.push((Some(a.into()), b.into()));
                                        } else {
                                            combine_error(
                                                &mut errors,
                                                BoxedError::new(
                                                    CVError::ItemError,
                                                    "Invalid xref",
                                                    "A RESID xref should be of shape TYPE:VALUE",
                                                    Context::default()
                                                        .byte_range(rule_node.range())
                                                        .lines(
                                                            0,
                                                            format!(
                                                                "RESID:AA{:04}",
                                                                modification.id
                                                            ),
                                                        ),
                                                ),
                                            );
                                        }
                                    }
                                    _ => (),
                                }
                            }
                            rules.push(rule);
                        } // Placement rules
                        "Features" => {
                            for feature in data_block.children() {
                                if feature.tag_name().name() == "Feature"
                                    && feature.attribute("type") == Some("UniProt")
                                    && let Some(text) = feature.text()
                                    && let Some(m) = text.strip_prefix("MOD_RES ")
                                {
                                    modification
                                        .cross_ids
                                        .push((Some("UniProt".into()), m.trim().into()));
                                }
                            }
                        }
                        _ => (),
                    }
                }

                let mut data = None;
                for rule in rules {
                    if data.is_none() {
                        if rule.1.is_none() {
                            data = Some(ModData::Mod {
                                specificities: Vec::new(),
                            });
                        } else {
                            data = Some(ModData::Linker {
                                length: crate::sequence::LinkerLength::Unknown,
                                specificities: Vec::new(),
                            });
                        }
                    }
                    if let (Some(ModData::Linker { specificities, .. }), Some(aa)) =
                        (&mut data, rule.1)
                    {
                        if rule.0 == aa {
                            specificities.push(LinkerSpecificity::Symmetric {
                                rules: vec![PlacementRule::AminoAcid(
                                    vec![rule.0].into(),
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
                                        vec![rule.0].into(),
                                        rule.2.unwrap_or(Position::Anywhere),
                                    )],
                                    vec![PlacementRule::AminoAcid(
                                        vec![aa].into(),
                                        rule.2.unwrap_or(Position::Anywhere),
                                    )],
                                ),
                                stubs: Vec::new(),
                                neutral_losses: Vec::new(),
                                diagnostic: Vec::new(),
                            });
                        }
                    } else if let (Some(ModData::Mod { specificities }), None) = (&mut data, rule.1)
                    {
                        specificities.push((
                            vec![PlacementRule::AminoAcid(
                                vec![rule.0].into(),
                                rule.2.unwrap_or(Position::Anywhere),
                            )],
                            Vec::new(),
                            Vec::new(),
                        ));
                    } else {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                CVError::ItemError,
                                "Both cross-linker and normal modification",
                                "This modification contains information for a normal modification and for a cross-linker",
                                Context::default()
                                    .byte_range(entry.range())
                                    .lines(0, format!("RESID:AA{:04}", modification.id)),
                            ),
                        );
                        continue 'entry;
                    }
                }

                modification.data = data.unwrap_or_default();
                modifications.push(modification.into());
            }
        }

        Ok((
            CVVersion {
                last_updated: database.attribute("date").map(|date| {
                    let mut res = date.split('-');
                    let day = res.next().and_then(|d| d.parse::<u8>().ok()).unwrap();
                    let month = res
                        .next()
                        .and_then(|m| {
                            let lower = m.to_ascii_lowercase();
                            [
                                "jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep",
                                "oct", "nov", "dec",
                            ]
                            .iter()
                            .position(|n| *n == lower)
                        })
                        .and_then(|p| u8::try_from(p).ok())
                        .unwrap();
                    let year = res.next().and_then(|d| d.parse::<u16>().ok()).unwrap();
                    (year, month, day, 0, 0)
                }),
                version: database.attribute("release").map(ToString::to_string),
                hash: reader.hash(),
            },
            modifications,
            errors,
        ))
    }
}
