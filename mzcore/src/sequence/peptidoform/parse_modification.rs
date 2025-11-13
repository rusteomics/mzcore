use std::{
    num::NonZeroU16,
    ops::Range,
    sync::{Arc, LazyLock},
};

use context_error::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    ParserResult,
    chemistry::{Element, MolecularFormula},
    glycan::{GlycanStructure, MonoSaccharide},
    helper_functions::*,
    ontology::{Ontologies, Ontology},
    sequence::{
        AmbiguousLookup, AmbiguousLookupEntry, CrossLinkLookup, CrossLinkName, MassTag,
        PlacementRule, SimpleModification, SimpleModificationInner,
    },
    system::{Mass, dalton},
};

impl SimpleModificationInner {
    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `BoxedError` explaining the error.
    pub fn pro_forma<'a>(
        line: &'a str,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, (ReturnModification, MUPSettings), BasicKind> {
        Self::pro_forma_inner(
            &Context::none().lines(0, line),
            line,
            0..line.len(),
            ambiguous_lookup,
            cross_link_lookup,
            ontologies,
        )
    }

    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `BoxedError` explaining the error.
    pub fn pro_forma_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, (ReturnModification, MUPSettings), BasicKind> {
        Self::pro_forma_main::<false>(
            base_context,
            line,
            range,
            ambiguous_lookup,
            cross_link_lookup,
            ontologies,
        )
    }

    /// Try to parse the modification. Any ambiguous modification will be numbered
    /// according to the lookup (which may be added to if necessary). The result
    /// is the modification, with, if applicable, its determined ambiguous group.
    /// # Errors
    /// If it is not a valid modification return a `BoxedError` explaining the error.
    pub(crate) fn pro_forma_main<'a, const STRICT: bool>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ambiguous_lookup: &mut AmbiguousLookup,
        cross_link_lookup: &mut CrossLinkLookup,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, (ReturnModification, MUPSettings), BasicKind> {
        let mut errors = Vec::new();
        // Because multiple modifications could be chained with the pipe operator
        // the parsing iterates through all links until it finds one it understands
        // it then returns that one. If no 'understandable' links are found it
        // returns the last link, if this is an info it returns a mass shift of 0,
        // but if any of the links returned an error it returns the last error.
        let mut modification = None;
        let mut settings = MUPSettings::default();
        let mut offset = range.start;
        for part in line[range].split('|') {
            match parse_single_modification::<STRICT>(
                base_context,
                line,
                part,
                offset,
                ambiguous_lookup,
                cross_link_lookup,
                ontologies,
            ) {
                Ok((result, w)) => {
                    combine_errors(&mut errors, w, ());
                    match result {
                        SingleReturnModification::None => (),
                        SingleReturnModification::Modification(m) => {
                            if modification.as_ref().is_none_or(|mo| m > *mo) {
                                modification = Some(m);
                            }
                        }
                        SingleReturnModification::Positions(p) => settings.position = Some(p),
                        SingleReturnModification::Limit(l) => settings.limit = Some(l),
                        SingleReturnModification::ColocalisePlacedModifications(s) => {
                            settings.colocalise_placed_modifications = s;
                        }
                        SingleReturnModification::ColocaliseModificationsOfUnknownPosition(s) => {
                            settings.colocalise_modifications_of_unknown_position = s;
                        }
                    }
                }
                Err(e) => combine_errors(&mut errors, e, ()),
            }
            offset += part.len() + 1;
        }
        if let Some(ReturnModification::Ambiguous(id, _, true)) = &modification {
            ambiguous_lookup[*id].copy_settings(&settings);
        }
        if errors.iter().any(|e| e.get_kind().is_error(())) {
            Err(errors)
        } else {
            Ok((
                (
                    modification.unwrap_or(ReturnModification::Defined(Arc::new(Self::Mass(
                        MassTag::None,
                        Mass::default().into(),
                        None,
                    )))),
                    settings,
                ),
                errors,
            ))
        }
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
    pub position: Option<Vec<PlacementRule>>,
    /// The maximal number for a grouped mup (^x)
    pub limit: Option<usize>,
    /// Allow this mup to colocalise with placed modifications
    pub colocalise_placed_modifications: bool,
    /// Allow this mup to colocalise with other mups
    pub colocalise_modifications_of_unknown_position: bool,
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
fn parse_single_modification<'error, const STRICT: bool>(
    base_context: &Context<'error>,
    line: &'error str,
    full_modification: &'error str,
    offset: usize,
    ambiguous_lookup: &mut AmbiguousLookup,
    cross_link_lookup: &mut CrossLinkLookup,
    ontologies: &Ontologies,
) -> ParserResult<'error, SingleReturnModification, BasicKind> {
    /// # Errors
    /// If the modification could not be found
    fn single_name_resolution<'a>(
        ontologies: &Ontologies,
        ontology: Ontology,
        name: &'a str,
        base_context: &Context<'a>,
        range: Range<usize>,
    ) -> ParserResult<'a, SimpleModification, BasicKind> {
        name_resolution(
            ontologies,
            &[ontology],
            &[ontology],
            name,
            base_context,
            range,
        )
    }

    /// # Errors
    /// If the modification could not be found
    fn name_resolution<'a>(
        ontologies: &Ontologies,
        match_selection: &[Ontology],
        search_selection: &[Ontology],
        name: &'a str,
        base_context: &Context<'a>,
        range: Range<usize>,
    ) -> ParserResult<'a, SimpleModification, BasicKind> {
        if let Some((by_name, modification)) =
            ontologies.get_by_name_or_synonym(match_selection, name)
        {
            if by_name {
                Ok((modification, Vec::new()))
            } else {
                Ok((
                    modification.clone(),
                    vec![BoxedError::new(
                        BasicKind::Warning,
                        "Used modification synonym",
                        "The name of the modification should be used to unambiguously identify a modification instead of a synonym",
                        base_context.clone().add_highlight((
                            0,
                            range,
                            format!("use {modification}"),
                        )),
                    )],
                ))
            }
        } else {
            Err(vec![
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid modification",
                    format!(
                        "The modification could not be found in {}",
                        if match_selection.len() == 1 {
                            match_selection[0].to_string()
                        } else if match_selection.len() == 2 {
                            format!("{} or {}", match_selection[0], match_selection[1])
                        } else if match_selection.is_empty() {
                            "Unimod, PSI-MOD, XL-MOD, GNOme, RESID, or Custom".to_string()
                        } else {
                            format!(
                                "{}, or {}",
                                match_selection[..match_selection.len() - 1]
                                    .iter()
                                    .map(ToString::to_string)
                                    .join(","),
                                match_selection[match_selection.len() - 1]
                            )
                        }
                    ),
                    base_context.clone().add_highlight((0, range)),
                )
                .suggestions(
                    ontologies
                        .search(search_selection, name)
                        .into_iter()
                        .map(|(m, s)| {
                            s.map_or_else(|| m.to_string(), |s| format!("`{s}` synonym for `{m}`"))
                        }),
                ),
            ])
        }
    }

    let mut errors = Vec::new();
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
                            vec![BoxedError::new(BasicKind::Error,
                        "Invalid modification localisation score",
                        "The ambiguous modification localisation score needs to be a valid number",
                        base_context.clone().add_highlight((0, offset + m.start(), m.len())),
                    )]
                        })
                })
                .transpose()?,
        );

        let modification = if let (Some(head), Some(tail)) = (head.as_ref(), tail) {
            let basic_error = BoxedError::new(
                BasicKind::Error,
                "Invalid modification",
                "..",
                base_context
                    .clone()
                    .add_highlight((0, offset + tail.1, tail.2)),
            );
            let sign_warning = BoxedError::new(
                BasicKind::Warning,
                "Improper modification",
                "A numerical modification should always be specified with a sign (+/-) to help it be recognised as a mass modification and not a modification index.",
                base_context
                    .clone()
                    .add_highlight((0, offset + tail.1, tail.2)),
            );
            match (head.0.as_str(), tail.0) {
                ("unimod", tail) => {
                    let id = tail.parse::<u32>().map_err(|_| {
                        vec![basic_error.clone()
                            .long_description("Unimod accession number should be a number")]
                    })?;
                    ontologies.unimod().get_by_index(&id)
                        .map(Some)
                        .ok_or_else(|| {
                            vec![basic_error.clone().long_description(
                            "The supplied Unimod accession number is not an existing modification",
                        )]
                        })
                }
                ("mod", tail) => {
                    let id = tail.parse::<u32>().map_err(|_| {
                        vec![basic_error.clone()
                            .long_description("PSI-MOD accession number should be a number")]
                    })?;
                    ontologies.psimod().get_by_index(&id)
                        .map(Some)
                        .ok_or_else(|| {
                            vec![basic_error.clone().long_description(
                            "The supplied PSI-MOD accession number is not an existing modification",
                        )]
                        })
                }
                ("resid", tail) => {
                    let id = tail[2..].parse::<u32>().map_err(|_| {
                        vec![basic_error.clone().long_description(
                            "RESID accession number should be a number prefixed with 'AA'",
                        )]
                    })?;
                    ontologies.resid().get_by_index(&id)
                        .map(Some)
                        .ok_or_else(|| {
                            vec![basic_error.clone().long_description(
                            "The supplied RESID accession number is not an existing modification",
                        )]
                        })
                }
                ("xlmod", tail) => {
                    let id = tail.parse::<u32>().map_err(|_| {
                        vec![basic_error.clone()
                            .long_description("XLMOD accession number should be a number")]
                    })?;
                    ontologies.xlmod().get_by_index(&id)
                        .map(Some)
                        .ok_or_else(|| {
                            vec![basic_error.clone().long_description(
                            "The supplied XLMOD accession number is not an existing modification",
                        )]
                        })
                }
                ("custom", tail) => {
                    let id = tail.parse::<u32>().map_err(|_| {
                        vec![basic_error.clone()
                            .long_description("Custom accession number should be a number")]
                    })?;
                    ontologies.custom().get_by_index(&id)
                    .map(Some)
                    .ok_or_else(|| {
                        vec![basic_error.clone().long_description(
                                "The supplied Custom accession number is not an existing modification",
                            )]
                    })
                }
                ("u", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Unimod), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Unimod, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("m", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Psimod), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Psimod, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("r", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Resid), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Resid, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("x", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Xlmod), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Xlmod, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("c", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Custom), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Custom, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("gno" | "g", name) => Ok(Some(handle!(errors, numerical_mod(MassTag::Ontology(Ontology::Gnome), name).map(|(sign, m)| (m, if STRICT && !sign {vec![sign_warning.clone()]} else {Vec::new()})).or_else(|_| 
                    single_name_resolution(ontologies, Ontology::Gnome, name, base_context, offset+tail.1..offset+tail.1+tail.2)
                )))),
                ("formula", _) => Ok(Some(Arc::new(SimpleModificationInner::Formula(
                    MolecularFormula::pro_forma_inner::<true, false>(base_context, line, offset + tail.1..offset + tail.1 + tail.2).map_err(|e| {
                       vec![e]
                    })?,
                )))),
                ("glycan", _) => Ok(Some(Arc::new(SimpleModificationInner::Glycan(
                    handle!(errors, MonoSaccharide::pro_forma_composition_inner::<STRICT>(base_context, line, offset + tail.1..offset + tail.1 + tail.2))
                )))),
                ("glycanstructure", _) => GlycanStructure::parse(
                    line,
                    offset + tail.1..offset + tail.1 + tail.2,
                )
                .map(|g| Some(Arc::new(SimpleModificationInner::GlycanStructure(g)))).map_err(|e| vec![e]),
                ("info", tail) => Ok(Some(Arc::new(SimpleModificationInner::Info(tail.to_string())))),
                ("obs", tail) => numerical_mod(MassTag::Observed,tail).map(|(sign, m)| {if STRICT && !sign {
                    combine_error(&mut errors, sign_warning.clone(), ());
                }
                Some(m)}).map_err(|_| {
                    vec![basic_error.long_description(
                        "This modification cannot be read as a numerical modification",
                    )]
                }),
                ("position", _) => {
                    match super::parse::parse_placement_rules(base_context, line, offset + tail.1..offset + tail.1 + tail.2) {
                        Ok(rules) => return Ok((SingleReturnModification::Positions(rules), errors)),
                        Err(e) => Err(vec![e]),
                    }
                }
                ("limit", tail) => {
                    match tail.parse::<usize>().map_err(|error| {
                        basic_error.long_description(format!(
                            "Invalid limit for modification of unknown position, the number is {}",
                            explain_number_error(&error)
                        ))
                    }) {
                        Ok(l) => return Ok((SingleReturnModification::Limit(l), errors)),
                        Err(e) => Err(vec![e]),
                    }
                }
                ("colocaliseplacedmodifications", tail) => {
                    match tail.parse::<bool>().map_err(|error| {
                        basic_error.long_description(format!(
                            "Invalid setting for colocalise placed modifications for modification of unknown position, the boolean is {error}",
                        ))
                    }) {
                        Ok(s) => return Ok((SingleReturnModification::ColocalisePlacedModifications(s), errors)),
                        Err(e) => Err(vec![e]),
                    }
                }
                ("colocalisemodificationsofunknownposition", tail) => {
                    match tail.parse::<bool>().map_err(|error| {
                        basic_error.long_description(format!(
                            "Invalid setting for colocalise modifications of unknown position for modification of unknown position, the boolean is {error}",
                        ))
                    }) {
                        Ok(s) => return Ok((SingleReturnModification::ColocaliseModificationsOfUnknownPosition(s), errors)),
                        Err(e) => Err(vec![e]),
                    }
                }
                (_, _) => Ok(Some(handle!(errors, name_resolution(
                        ontologies,
                        &[Ontology::Unimod, Ontology::Psimod],
                        &[],
                        full.0,
                        base_context,
                        offset + full.1..offset + full.1 + full.2
                    ))))
            }
        } else if full.0.is_empty() {
            Ok(None)
        } else {
            Ok(Some(handle!(
                errors,
                numerical_mod(MassTag::None, full.0)
                .map(|(sign, m)| (m, if STRICT && !sign {vec![BoxedError::new(
                    BasicKind::Warning,
                    "Improper modification",
                    "A numerical modification should always be specified with a sign (+/-) to help it be recognised as a mass modification and not a modification index.",
                    base_context
                        .clone()
                        .add_highlight((0, offset + full.1, full.2)),
                )]} else {Vec::new()}))
                    .or_else(|_| name_resolution(
                        ontologies,
                        &[Ontology::Unimod, Ontology::Psimod],
                        &[],
                        full.0,
                        base_context,
                        offset + full.1..offset + full.1 + full.2
                    ))
            )))
        };
        let modification = modification?;

        if let Some(modification) = &modification
            && let Some(description) = modification.description()
            && description.obsolete
        {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Warning,
                    "Obsolete modification",
                    "The used modification is marked obsolete",
                    base_context.clone().add_highlight((
                        0,
                        offset + full.1..offset + full.1 + full.2,
                        description.description.to_string(),
                    )),
                ),
                (),
            );
        }

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
                if let Some(linker) = modification {
                    if cross_link_lookup[index]
                        .1
                        .as_ref()
                        .is_some_and(|l| *l != linker)
                    {
                        return Err(vec![BoxedError::new(
                            BasicKind::Error,
                            "Invalid branch definition",
                            "A branch definition has to be identical at both sites, or only defined at one site.",
                            base_context
                                .clone()
                                .add_highlight((0, offset + full.1, full.2)),
                        )]);
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok((
                    SingleReturnModification::Modification(
                        ReturnModification::CrossLinkReferenced(index),
                    ),
                    errors,
                ))
            } else if let Some(name) = group.0.to_ascii_lowercase().strip_prefix("xl") {
                let name = CrossLinkName::Name(name.to_string().into_boxed_str());
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
                if let Some(linker) = modification {
                    if cross_link_lookup[index]
                        .1
                        .as_ref()
                        .is_some_and(|l| *l != linker)
                    {
                        return Err(vec![BoxedError::new(
                            BasicKind::Error,
                            "Invalid cross-link definition",
                            "A cross-link definition has to be identical at both sites, or only defined at one site.",
                            base_context
                                .clone()
                                .add_highlight((0, offset + full.1, full.2)),
                        )]);
                    }
                    cross_link_lookup[index].1 = Some(linker);
                }
                Ok((
                    SingleReturnModification::Modification(
                        ReturnModification::CrossLinkReferenced(index),
                    ),
                    errors,
                ))
            } else {
                handle_ambiguous_modification(
                    modification,
                    group,
                    localisation_score,
                    ambiguous_lookup,
                    base_context
                        .clone()
                        .add_highlight((0, offset + full.1 + 1, full.2)),
                )
                .map(|(m, w)| {
                    combine_errors(&mut errors, w, ());
                    (m, errors)
                })
            }
        } else {
            Ok((
                modification.map_or(SingleReturnModification::None, |m| {
                    SingleReturnModification::Modification(ReturnModification::Defined(m))
                }),
                errors,
            ))
        }
    } else {
        Err(vec![BoxedError::new(
            BasicKind::Error,
            "Invalid modification",
            "It does not match the ProForma definition for modifications",
            base_context
                .clone()
                .add_highlight((0, offset, full_modification.len())),
        )])
    }
}

/// Handle the logic for an ambiguous modification
/// # Errors
/// If the content of the ambiguous modification was already defined
fn handle_ambiguous_modification<'a>(
    modification: Option<SimpleModification>,
    group: (&str, usize, usize),
    localisation_score: Option<OrderedFloat<f64>>,
    ambiguous_lookup: &mut AmbiguousLookup,
    context: Context<'a>,
) -> ParserResult<'a, SingleReturnModification, BasicKind> {
    // Search for a previous definition of this name, store as Some((index, modification_definition_present)) or None if there is no definition in place
    let found_definition = ambiguous_lookup
        .iter()
        .enumerate()
        .find(|(_, entry)| entry.name.eq_ignore_ascii_case(group.0))
        .map(|(index, entry)| (index, entry.modification.as_ref()));
    // Handle all possible cases of having a modification found at this position and having a modification defined in the ambiguous lookup
    match (modification, found_definition) {
        // Have a mod defined here and already in the lookup (error)
        (Some(m), Some((index, Some(f)))) => {
            if *m == **f {
                Ok((
                    SingleReturnModification::Modification(ReturnModification::Ambiguous(
                        index,
                        localisation_score,
                        false,
                    )),
                    vec![BoxedError::new(
                        BasicKind::Warning,
                        "Invalid ambiguous modification",
                        "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
                        context,
                    )],
                ))
            } else {
                Err(vec![BoxedError::new(
                    BasicKind::Error,
                    "Invalid ambiguous modification",
                    "An ambiguous modification cannot be placed twice (for one of the modifications leave out the modification and only provide the group name)",
                    context,
                )])
            }
        }
        // Have a mod defined here, the name present in the lookup but not yet the mod
        (Some(m), Some((index, None))) => {
            ambiguous_lookup[index].modification = Some(m);
            Ok((
                SingleReturnModification::Modification(ReturnModification::Ambiguous(
                    index,
                    localisation_score,
                    true,
                )),
                Vec::new(),
            ))
        }
        // No mod defined, but the name is present in the lookup
        (None, Some((index, _))) => Ok((
            SingleReturnModification::Modification(ReturnModification::Ambiguous(
                index,
                localisation_score,
                false,
            )),
            Vec::new(),
        )),
        // The mod is not already in the lookup
        (m, None) => {
            let index = ambiguous_lookup.len();
            let preferred = m.is_some();
            ambiguous_lookup.push(AmbiguousLookupEntry::new(group.0.to_string(), m));
            Ok((
                SingleReturnModification::Modification(ReturnModification::Ambiguous(
                    index,
                    localisation_score,
                    preferred,
                )),
                Vec::new(),
            ))
        }
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

/// Returns a mass modification is this text could be parsed as a number.
/// It also returns a boolean indicating if a sign was present (`true`) or not present to be able to create a warning for strict ProForma handling.
/// # Errors
/// It returns an error when the text is not numerical
pub(super) fn numerical_mod(
    tag: MassTag,
    text: &str,
) -> Result<(bool, SimpleModification), String> {
    text.parse().map_or_else(
        |_| Err("Invalid number".to_string()),
        |n| {
            Ok((
                text.starts_with('-') || text.starts_with('+'),
                Arc::new(SimpleModificationInner::Mass(
                    tag,
                    Mass::new::<dalton>(n).into(),
                    float_digits(text),
                )),
            ))
        },
    )
}
