use crate::modification::{
    AmbiguousLookup, CrossLinkLookup, CrossLinkName, Ontology, SimpleModification,
    SimpleModificationInner,
};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{
    num::NonZeroU16,
    ops::Range,
    sync::{Arc, OnceLock},
};

use regex::Regex;

use crate::{
    error::{Context, CustomError},
    glycan::{GlycanStructure, MonoSaccharide},
    helper_functions::*,
    ontologies::CustomDatabase,
    placement_rule::Position,
    system::{dalton, Mass, OrderedMass},
    AminoAcid, Element, MolecularFormula,
};

impl SimpleModificationInner {
    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `CustomError` explaining the error.
    pub fn try_from(
        line: &str,
        range: Range<usize>,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<ReturnModification, CustomError> {
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates over all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut last_result = Ok(None);
        let mut last_error = None;
        let mut offset = range.start;
        for part in line[range].split('|') {
            last_result = parse_single_modification(
                line,
                part,
                offset,
                ambiguous_lookup,
                cross_link_lookup,
                custom_database,
            );
            if let Ok(Some(m)) = last_result {
                return Ok(m);
            }
            if let Err(er) = &last_result {
                last_error = Some(er.clone());
            }
            offset += part.len() + 1;
        }
        last_error.map_or_else(
            || {
                last_result.map(|m| {
                    m.unwrap_or_else(|| {
                        ReturnModification::Defined(Arc::new(Self::Mass(OrderedMass::zero())))
                    })
                })
            },
            Err,
        )
    }
}

static MOD_REGEX: OnceLock<Regex> = OnceLock::new();

/// # Errors
/// It returns an error when the given line cannot be read as a single modification.
#[allow(clippy::missing_panics_doc)]
fn parse_single_modification(
    line: &str,
    full_modification: &str,
    offset: usize,
    ambiguous_lookup: &mut AmbiguousLookup,
    cross_link_lookup: &mut CrossLinkLookup,
    custom_database: Option<&CustomDatabase>,
) -> Result<Option<ReturnModification>, CustomError> {
    // Parse the whole intricate structure of the single modification (see here in action: https://regex101.com/r/pW5gsj/1)
    let regex = MOD_REGEX.get_or_init(|| {
        Regex::new(r"^(([^:#]*)(?::([^#]+))?)(?:#([0-9A-Za-z]+)(?:\((\d+\.\d+)\))?)?$").unwrap()
    });
    if let Some(groups) = regex.captures(full_modification) {
        // Capture the full mod name (head:tail), head, tail, ambiguous group, and localisation score
        let (full, head, tail, label_group, localisation_score) = (
            groups
                .get(1)
                .map(|m| (m.as_str(), m.start(), m.len()))
                .unwrap_or_default(),
            groups
                .get(2)
                .map(|m| (m.as_str().to_ascii_lowercase(), m.start(), m.len())),
            groups.get(3).map(|m| (m.as_str(), m.start(), m.len())),
            groups.get(4).map(|m| (m.as_str(), m.start(), m.len())),
            groups
                .get(5)
                .map(|m| {
                    m.as_str()
                        .parse::<f64>()
                        .map(OrderedFloat::from)
                        .map_err(|_| {
                            CustomError::error(
                        "Invalid modification localisation score",
                        "The ambiguous modification localisation score needs to be a valid number",
                        Context::line(None, line, offset + m.start(), m.len()),
                    )
                        })
                })
                .transpose()?,
        );

        let modification = if let (Some(head), Some(tail)) = (head.as_ref(), tail) {
            let basic_error = CustomError::error(
                "Invalid modification",
                "..",
                Context::line(None, line, offset + tail.1, tail.2),
            );
            match (head.0.as_str(), tail.0) {
                ("unimod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("Unimod accession number should be a number")
                    })?;
                    Ontology::Unimod
                        .find_id(id, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            basic_error.with_long_description(
                            "The supplied Unimod accession number is not an existing modification",
                        )
                        })
                }
                ("mod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("PSI-MOD accession number should be a number")
                    })?;
                    Ontology::Psimod
                        .find_id(id, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            basic_error.with_long_description(
                            "The supplied PSI-MOD accession number is not an existing modification",
                        )
                        })
                }
                ("resid", tail) => {
                    let id = tail[2..].parse::<usize>().map_err(|_| {
                        basic_error.with_long_description(
                            "RESID accession number should be a number prefixed with 'AA'",
                        )
                    })?;
                    Ontology::Resid
                        .find_id(id, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            basic_error.with_long_description(
                            "The supplied Resid accession number is not an existing modification",
                        )
                        })
                }
                ("xlmod", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("XLMOD accession number should be a number")
                    })?;
                    Ontology::Xlmod
                        .find_id(id, custom_database)
                        .map(Some)
                        .ok_or_else(|| {
                            basic_error.with_long_description(
                            "The supplied XLMOD accession number is not an existing modification",
                        )
                        })
                }
                ("custom", tail) => {
                    let id = tail.parse::<usize>().map_err(|_| {
                        basic_error
                            .with_long_description("Custom accession number should be a number")
                    })?;
                    Ontology::Custom
                    .find_id(id, custom_database)
                    .map(Some)
                    .ok_or_else(|| {
                        basic_error.with_long_description(
                                "The supplied Custom accession number is not an existing modification",
                            )
                    })
                }
                ("u", tail) => Ontology::Unimod
                    .find_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Unimod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("m", tail) => Ontology::Psimod
                    .find_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Psimod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("r", tail) => Ontology::Resid
                    .find_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Resid
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("x", tail) => Ontology::Xlmod
                    .find_name(tail, custom_database)
                    .ok_or_else(|| numerical_mod(tail))
                    .flat_err()
                    .map(Some)
                    .map_err(|_| {
                        Ontology::Xlmod
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("c", tail) => Ontology::Custom
                    .find_name(tail, custom_database)
                    .map(Some)
                    .ok_or_else(|| {
                        Ontology::Custom
                            .find_closest(tail, custom_database)
                            .with_context(basic_error.context().clone())
                    }),
                ("gno" | "g", tail) => Ontology::Gnome
                    .find_name(tail, custom_database)
                    .map(Some)
                    .ok_or_else(|| {
                        basic_error
                            .with_long_description("This modification cannot be read as a GNO name")
                    }),
                ("formula", tail) => Ok(Some(Arc::new(SimpleModificationInner::Formula(
                    MolecularFormula::from_pro_forma(tail, .., false, false, true).map_err(
                        |e| {
                            basic_error.with_long_description(format!(
                                "This modification cannot be read as a valid formula: {e}"
                            ))
                        },
                    )?,
                )))),
                ("glycan", tail) => Ok(Some(Arc::new(SimpleModificationInner::Glycan(
                    MonoSaccharide::from_composition(tail)
                        .map_err(|err| err.with_context(basic_error.context().clone()))?,
                )))),
                ("glycanstructure", _) => GlycanStructure::parse(
                    &line.to_ascii_lowercase(),
                    offset + tail.1..offset + tail.1 + tail.2,
                )
                .map(|g| Some(Arc::new(SimpleModificationInner::GlycanStructure(g)))),
                ("info", _) => Ok(None),
                ("obs", tail) => numerical_mod(tail).map(Some).map_err(|_| {
                    basic_error.with_long_description(
                        "This modification cannot be read as a numerical modification",
                    )
                }),
                (_, _) => Ontology::Unimod
                    .find_name(full.0, custom_database)
                    .or_else(|| Ontology::Psimod.find_name(full.0, custom_database))
                    .map(Some)
                    .ok_or_else(|| {
                        Ontology::find_closest_many(
                            &[
                                Ontology::Unimod,
                                Ontology::Psimod,
                                Ontology::Gnome,
                                Ontology::Xlmod,
                                Ontology::Resid,
                                Ontology::Custom,
                            ],
                            full.0,
                            custom_database,
                        )
                        .with_long_description(
                            "This modification cannot be read as a valid Unimod or PSI-MOD name.",
                        )
                        .with_context(Context::line(
                            None,
                            line,
                            offset + full.1,
                            full.2,
                        ))
                    }),
            }
        } else if full.0.is_empty() {
            Ok(None)
        } else {
            Ontology::Unimod.find_name(full.0 , custom_database)
                .or_else(|| Ontology::Psimod.find_name(full.0, custom_database))
                .ok_or_else(|| numerical_mod(full.0))
                .flat_err()
                .map(Some)
                .map_err(|_|
                    Ontology::find_closest_many(&[Ontology::Unimod, Ontology::Psimod, Ontology::Gnome, Ontology::Xlmod, Ontology::Resid, Ontology::Custom], full.0, custom_database)
                    .with_long_description("This modification cannot be read as a valid Unimod or PSI-MOD name, or as a numerical modification.")
                    .with_context(Context::line(None, line, offset+full.1, full.2))
                )
        };

        if let Some(group) = label_group {
            if group.0.to_ascii_lowercase() == "branch" {
                let index = cross_link_lookup
                    .iter()
                    .position(|c| c.0 == CrossLinkName::Branch)
                    .map_or_else(
                        || {
                            let index = cross_link_lookup.len();
                            cross_link_lookup.push((CrossLinkName::Branch, None));
                            index
                        },
                        |index| index,
                    );
                if let Some(linker) = modification? {
                    if cross_link_lookup[index]
                        .1
                        .as_ref()
                        .is_some_and(|l| *l != linker)
                    {
                        return Err(CustomError::error(
                                "Invalid branch definition", 
                                "A branch definition has to be identical at both sites, or only defined at one site.", 
                                Context::line(None, line, offset+full.1, full.2)
                            ));
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok(Some(ReturnModification::CrossLinkReferenced(index)))
            } else if let Some(name) = group.0.to_ascii_lowercase().strip_prefix("xl") {
                let name = CrossLinkName::Name(name.to_string());
                let index = cross_link_lookup
                    .iter()
                    .position(|c| c.0 == name)
                    .map_or_else(
                        || {
                            let index = cross_link_lookup.len();
                            cross_link_lookup.push((name, None));
                            index
                        },
                        |index| index,
                    );
                if let Some(linker) = modification? {
                    if cross_link_lookup[index]
                        .1
                        .as_ref()
                        .is_some_and(|l| *l != linker)
                    {
                        return Err(CustomError::error(
                            "Invalid cross-link definition", 
                            "A cross-link definition has to be identical at both sites, or only defined at one site.", 
                            Context::line(None, line, offset+full.1, full.2)
                        ));
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok(Some(ReturnModification::CrossLinkReferenced(index)))
            } else {
                handle_ambiguous_modification(
                    modification,
                    group,
                    localisation_score,
                    ambiguous_lookup,
                    Context::line(None, line, offset + full.1, full.2),
                )
            }
        } else {
            modification.map(|m| m.map(ReturnModification::Defined))
        }
    } else {
        Err(CustomError::error(
            "Invalid modification",
            "It does not match the ProForma definition for modifications",
            Context::line(None, line, offset, full_modification.len()),
        ))
    }
}

/// Handle the logic for an ambiguous modification
/// # Errors
/// If the content of the ambiguous modification was already defined
fn handle_ambiguous_modification(
    modification: Result<Option<SimpleModification>, CustomError>,
    group: (&str, usize, usize),
    localisation_score: Option<OrderedFloat<f64>>,
    ambiguous_lookup: &mut AmbiguousLookup,
    context: Context,
) -> Result<Option<ReturnModification>, CustomError> {
    let group_name = group.0.to_ascii_lowercase();
    // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
    let found_definition = ambiguous_lookup
        .iter()
        .enumerate()
        .find(|(_, (name, _))| name == &group_name)
        .map(|(index, (_, modification))| (index, modification.is_some()));
    // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
    match (modification, found_definition) {
        // Have a mod defined here and already in the lookup (error)
        (Ok(Some(_)), Some((_, true))) => Err(
            CustomError::error(
                "Invalid ambiguous modification",
                "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
                context,
            )),
        // Have a mod defined here, the name present in the lookup but not yet the mod
        (Ok(Some(m)), Some((index, false))) => {
            ambiguous_lookup[index].1 = Some(m);
            Ok(Some(ReturnModification::Ambiguous(index, localisation_score, true)))
        },
        // No mod defined, but the name is present in the lookup
        (Ok(None), Some((index, _))) => Ok(Some(ReturnModification::Ambiguous(index, localisation_score, false))),
        // The mod is not already in the lookup
        (Ok(m), None) => {
            let index = ambiguous_lookup.len();
            let preferred = m.is_some();
            ambiguous_lookup.push((group_name, m));
            Ok(Some(ReturnModification::Ambiguous(index, localisation_score, preferred)))
        },
        // Earlier error
        (Err(e), _) => Err(e),
    }
}

/// A modification as returned by the parser
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum ReturnModification {
    /// A fully self contained modification
    Defined(SimpleModification),
    /// A modification that references an ambiguous modification, (id, localisation score, preferred)
    Ambiguous(usize, Option<OrderedFloat<f64>>, bool),
    /// A modification that references a cross-link
    CrossLinkReferenced(usize),
}

impl ReturnModification {
    /// Force this modification to be defined
    #[must_use]
    pub fn defined(self) -> Option<SimpleModification> {
        match self {
            Self::Defined(modification) => Some(modification),
            _ => None,
        }
    }
}

/// Intermediate representation of a global modification
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum GlobalModification {
    /// A global isotope modification
    Isotope(Element, Option<NonZeroU16>),
    /// Can be placed on any place it fits, if that is the correct aminoacid and it fits according to the placement rules of the modification itself
    Fixed(Position, Option<AminoAcid>, SimpleModification),
}

/// # Errors
/// It returns an error when the text is not numerical
pub(super) fn numerical_mod(text: &str) -> Result<SimpleModification, String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| {
            Ok(Arc::new(SimpleModificationInner::Mass(
                Mass::new::<dalton>(n).into(),
            )))
        },
    )
}
