use std::{
    num::NonZeroU16,
    ops::Range,
    sync::{Arc, LazyLock},
};

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use regex::Regex;

use crate::{
    chemistry::{Element, MolecularFormula},
    error::{Context, CustomError},
    glycan::{GlycanStructure, MonoSaccharide},
    helper_functions::*,
    ontology::{CustomDatabase, Ontology},
    sequence::{
        AmbiguousLookup, AmbiguousLookupEntry, CrossLinkLookup, CrossLinkName, PlacementRule,
        SimpleModification, SimpleModificationInner,
    },
    system::{Mass, dalton},
};

impl SimpleModificationInner {
    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `CustomError` explaining the error.
    pub fn parse_pro_forma(
        line: &str,
        range: Range<usize>,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<(ReturnModification, MUPSettings), CustomError> {
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates over all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut modification = None;
        let mut settings = MUPSettings::default();
        let mut last_error = None;
        let mut offset = range.start;
        for part in line[range].split('|') {
            match parse_single_modification(
                line,
                part,
                offset,
                ambiguous_lookup,
                cross_link_lookup,
                custom_database,
            ) {
                Ok(SingleReturnModification::None) => (),
                Ok(SingleReturnModification::Modification(m)) => modification = Some(m),
                Ok(SingleReturnModification::Positions(p)) => settings.position = Some(p),
                Ok(SingleReturnModification::Limit(l)) => settings.limit = Some(l),
                Ok(SingleReturnModification::ColocalisePlacedModifications(s)) => {
                    settings.colocalise_placed_modifications = s;
                }
                Ok(SingleReturnModification::ColocaliseModificationsOfUnknownPosition(s)) => {
                    settings.colocalise_modifications_of_unknown_position = s;
                }
                Err(e) => last_error = Some(e),
            }
            offset += part.len() + 1;
        }
        if let Some(ReturnModification::Ambiguous(id, _, true)) = &modification {
            ambiguous_lookup[*id].copy_settings(&settings);
        }
        last_error.map_or_else(
            || {
                Ok((
                    modification.unwrap_or(ReturnModification::Defined(Arc::new(Self::Mass(
                        Mass::default().into(),
                    )))),
                    settings,
                ))
            },
            Err,
        )
    }
}

static MOD_REGEX: LazyLock<Regex> = LazyLock::new(|| {
    Regex::new(r"^(([^:#]*)(?::([^#]+))?)(?:#([0-9A-Za-z]+)(?:\((\d+\.\d+)\))?)?$").unwrap()
});

enum SingleReturnModification {
    None,
    Modification(ReturnModification),
    Positions(Vec<PlacementRule>),
    Limit(usize),
    ColocalisePlacedModifications(bool),
    ColocaliseModificationsOfUnknownPosition(bool),
}

/// Settings for a modification of unknown position
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MUPSettings {
    /// The additional placement rules
    pub(crate) position: Option<Vec<PlacementRule>>,
    /// The maximal number for a grouped mup (^x)
    pub(crate) limit: Option<usize>,
    /// Allow this mup to colocalise with placed modifications
    pub(crate) colocalise_placed_modifications: bool,
    /// Allow this mup to colocalise with other mups
    pub(crate) colocalise_modifications_of_unknown_position: bool,
}

impl Default for MUPSettings {
    fn default() -> Self {
        Self {
            position: None,
            limit: None,
            colocalise_placed_modifications: true,
            colocalise_modifications_of_unknown_position: true,
        }
    }
}

/// # Errors
/// It returns an error when the given line cannot be read as a single modification.
fn parse_single_modification(
    line: &str,
    full_modification: &str,
    offset: usize,
    ambiguous_lookup: &mut AmbiguousLookup,
    cross_link_lookup: &mut CrossLinkLookup,
    custom_database: Option<&CustomDatabase>,
) -> Result<SingleReturnModification, CustomError> {
    // Parse the whole intricate structure of the single modification (see here in action: https://regex101.com/r/pW5gsj/1)
    if let Some(groups) = MOD_REGEX.captures(full_modification) {
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
                    MolecularFormula::from_pro_forma(tail, .., true, false, true, true).map_err(|e| {
                        basic_error.with_long_description(format!(
                            "This modification cannot be read as a valid formula: {e}"
                        ))
                    })?,
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
                ("position", tail) => {
                    match super::parse::parse_placement_rules(tail, 0..tail.len()).map_err(
                        |error| basic_error.with_long_description(error.long_description()),
                    ) {
                        Ok(rules) => return Ok(SingleReturnModification::Positions(rules)),
                        Err(e) => Err(e),
                    }
                }
                ("limit", tail) => {
                    match tail.parse::<usize>().map_err(|error| {
                        basic_error.with_long_description(format!(
                            "Invalid limit for modification of unknown position, the number is {}",
                            explain_number_error(&error)
                        ))
                    }) {
                        Ok(l) => return Ok(SingleReturnModification::Limit(l)),
                        Err(e) => Err(e),
                    }
                }
                ("colocaliseplacedmodifications", tail) => {
                    match tail.parse::<bool>().map_err(|error| {
                        basic_error.with_long_description(format!(
                            "Invalid setting for colocalise placed modifications for modification of unknown position, the boolean is {error}",
                        ))
                    }) {
                        Ok(s) => return Ok(SingleReturnModification::ColocalisePlacedModifications(s)),
                        Err(e) => Err(e),
                    }
                }
                ("colocalisemodificationsofunknownposition", tail) => {
                    match tail.parse::<bool>().map_err(|error| {
                        basic_error.with_long_description(format!(
                            "Invalid setting for colocalise modifications of unknown position for modification of unknown position, the boolean is {error}",
                        ))
                    }) {
                        Ok(s) => return Ok(SingleReturnModification::ColocaliseModificationsOfUnknownPosition(s)),
                        Err(e) => Err(e),
                    }
                }
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
                            full.2
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
            if group.0.eq_ignore_ascii_case("branch") {
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
                            Context::line(None, line, offset + full.1, full.2),
                        ));
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok(SingleReturnModification::Modification(
                    ReturnModification::CrossLinkReferenced(index),
                ))
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
                            Context::line(None, line, offset + full.1, full.2),
                        ));
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok(SingleReturnModification::Modification(
                    ReturnModification::CrossLinkReferenced(index),
                ))
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
            modification.map(|m| {
                m.map_or(SingleReturnModification::None, |m| {
                    SingleReturnModification::Modification(ReturnModification::Defined(m))
                })
            })
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
) -> Result<SingleReturnModification, CustomError> {
    let group_name = group.0.to_ascii_lowercase();
    // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
    let found_definition = ambiguous_lookup
        .iter()
        .enumerate()
        .find(|(_, entry)| entry.name == group_name)
        .map(|(index, entry)| (index, entry.modification.is_some()));
    // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
    match (modification, found_definition) {
        // Have a mod defined here and already in the lookup (error)
        (Ok(Some(_)), Some((_, true))) => Err(CustomError::error(
            "Invalid ambiguous modification",
            "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
            context,
        )),
        // Have a mod defined here, the name present in the lookup but not yet the mod
        (Ok(Some(m)), Some((index, false))) => {
            ambiguous_lookup[index].modification = Some(m);
            Ok(SingleReturnModification::Modification(
                ReturnModification::Ambiguous(index, localisation_score, true),
            ))
        }
        // No mod defined, but the name is present in the lookup
        (Ok(None), Some((index, _))) => Ok(SingleReturnModification::Modification(
            ReturnModification::Ambiguous(index, localisation_score, false),
        )),
        // The mod is not already in the lookup
        (Ok(m), None) => {
            let index = ambiguous_lookup.len();
            let preferred = m.is_some();
            ambiguous_lookup.push(AmbiguousLookupEntry::new(group_name, m));
            Ok(SingleReturnModification::Modification(
                ReturnModification::Ambiguous(index, localisation_score, preferred),
            ))
        }
        // Earlier error
        (Err(e), _) => Err(e),
    }
}

/// A modification as returned by the parser
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
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
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum GlobalModification {
    /// A global isotope modification
    Isotope(Element, Option<NonZeroU16>),
    /// Can be placed on any place it fits, if that is the correct aminoacid and it fits according to the placement rules of the modification itself
    Fixed(PlacementRule, SimpleModification),
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
