//! Code to handle the RESID ontology
use std::io::Read;

use context_error::{BoxedError, CreateError};

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader};

use roxmltree::*;

use crate::{
    chemistry::{MolecularFormula, MultiChemical},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{AminoAcid, LinkerSpecificity, PlacementRule, Position, SimpleModificationInner},
};

/// RESID modifications
///
/// Note that because RESID is stored on a ftp server the automatic downloading does not work.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct Resid {}

impl CVSource for Resid {
    type Data = SimpleModificationInner;
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
                    include_bytes!("../databases/resid.dat"),
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
        .expect("Invalid xml in RESID xml");
        let mut modifications = Vec::new();
        let database = document
            .root()
            .first_child()
            .expect("No Database node in RESID XML");

        'entry: for entry in database.children() {
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
                                            modification
                                                .cross_ids
                                                .push((Some(a.to_string()), b.to_string()));
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
                                    modification.formula = MolecularFormula::resid(
                                        formula_node.text().unwrap_or_default(),
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
                                        MolecularFormula::resid_single(
                                            formula_node.text().unwrap_or_default(),
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
                                        modification
                                            .cross_ids
                                            .push((Some(a.to_string()), b.to_string()));
                                    } else {
                                        panic!("Invalid Xref content")
                                    }
                                }
                            }
                        }
                        "Comment" => {
                            modification.description += data_block.text().unwrap_or_default();
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
                                            modification
                                                .cross_ids
                                                .push((Some(a.to_string()), b.to_string()));
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
                    if let (Some(ModData::Linker { specificities, .. }), Some(aa)) =
                        (&mut data, rule.1)
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
                    } else if let (Some(ModData::Mod { specificities }), None) = (&mut data, rule.1)
                    {
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
                hash: reader.hash().to_vec(),
            },
            modifications.into_iter(),
        ))
    }
}
