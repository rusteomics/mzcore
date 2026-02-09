use std::{collections::BTreeMap, num::NonZeroU16, ops::Range};

use context_error::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::{
    ParserResult,
    chemistry::{Element, MolecularCharge, MolecularFormula},
    helper_functions::*,
    ontology::Ontologies,
    sequence::{
        AmbiguousLookup, AmbiguousLookupEntry, AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon,
        CrossLinkLookup, Linked, MUPSettings, Peptidoform, PeptidoformIon, PlacementRule, Position,
        SequenceElement, SequencePosition, SimpleModification, SimpleModificationInner,
    },
};

use super::{GlobalModification, Linear, ReturnModification, SemiAmbiguous};

#[derive(Debug, Eq, PartialEq)]
enum End {
    Empty,
    CrossLink,
    Chimeric,
}

struct LinearPeptideResult {
    peptide: Peptidoform<Linear>,
    index: usize,
    ending: End,
    cross_links: Vec<(usize, SequencePosition)>,
}

impl Peptidoform<Linked> {
    /// Convenience wrapper to parse a linear peptide in ProForma notation, to handle all possible ProForma sequences look at [`CompoundPeptidoformIon::pro_forma`].
    /// # Errors
    /// It gives an error when the peptide is not correctly formatted. (Also see the `CompoundPeptidoformIon` main function for this.)
    /// It additionally gives an error if the peptide specified was chimeric (see [`CompoundPeptidoformIon::singular`] and [`PeptidoformIon::singular`]).
    pub fn pro_forma<'a>(
        value: &'a str,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_inner(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
            ontologies,
        )
    }

    /// This parses a substring of the given string as a [ProForma](https://psidev.info/proforma) definition. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string. Or when the string has multiple peptidoforms.
    pub fn pro_forma_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        let (pep, mut warnings) =
            PeptidoformIon::pro_forma_inner(base_context, line, range.clone(), ontologies)?;
        if let Some(pep) = pep.singular() {
            if warnings.iter().any(|e| e.get_kind().is_error(())) {
                Err(warnings)
            } else {
                Ok((pep, warnings))
            }
        } else {
            combine_error(
                &mut warnings,
                BoxedError::new(
                    BasicKind::Error,
                    "Peptidoform ion found",
                    "A linear peptidoform was expected but a cross linked peptidoform ion was found.",
                    base_context.clone().add_highlight((0, range)),
                ),
            );
            Err(warnings)
        }
    }
}

impl PeptidoformIon {
    /// Parse a peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string. Or when the string has multiple peptidoforms.
    pub fn pro_forma<'a>(
        value: &'a str,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_inner(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
            ontologies,
        )
    }

    /// This parses a substring of the given string as a [ProForma](https://psidev.info/proforma) definition. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string. Or when the string has multiple peptidoforms.
    pub fn pro_forma_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        let (pep, mut warnings) =
            CompoundPeptidoformIon::pro_forma_inner(base_context, line, range.clone(), ontologies)?;
        if let Some(pep) = pep.singular() {
            Ok((pep, warnings))
        } else {
            combine_error(
                &mut warnings,
                BoxedError::new(
                    BasicKind::Error,
                    "Compound peptidoform ion found",
                    "A peptide ion was expected but chimeric peptidoforms were found.",
                    base_context.clone().add_highlight((0, range)),
                ),
            );
            Err(warnings)
        }
    }
}

impl CompoundPeptidoformIon {
    /// Parse a compound peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    pub fn pro_forma<'a>(
        value: &'a str,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_inner(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
            ontologies,
        )
    }

    /// This parses a substring of the given string as a [ProForma](https://psidev.info/proforma) definition. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    pub fn pro_forma_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_main::<false>(base_context, line, range, ontologies)
    }

    /// Parse a compound peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    /// Generate warnings for any and all things that are not exactly to specification (not that
    /// not all warnings are in yet at the time of writing.)
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    pub fn pro_forma_strict<'a>(
        value: &'a str,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_inner_strict(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
            ontologies,
        )
    }

    /// This parses a substring of the given string as a [ProForma](https://psidev.info/proforma) definition. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    /// Generate warnings for any and all things that are not exactly to specification (not that
    /// not all warnings are in yet at the time of writing.)
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    pub fn pro_forma_inner_strict<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        Self::pro_forma_main::<true>(base_context, line, range, ontologies)
    }

    fn pro_forma_main<'a, const STRICT: bool>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
    ) -> ParserResult<'a, Self, BasicKind> {
        let mut peptidoforms = Vec::new();
        let mut errors = Vec::new();
        // Global modifications
        let (mut start, global_modifications) = handle!(
            errors,
            global_modifications::<STRICT>(base_context, line, range.clone(), ontologies)
        );
        let (peptidoform, tail) = handle!(
            errors,
            Self::parse_peptidoform::<STRICT>(
                base_context,
                line,
                start..range.end,
                &global_modifications,
                ontologies,
            )
        );
        start = tail;
        peptidoforms.push(peptidoform);

        // Parse any following chimeric species
        while start < line.len() && start < range.end {
            let (peptidoform, tail) = handle!(
                errors,
                Self::parse_peptidoform::<STRICT>(
                    base_context,
                    line,
                    start..range.end,
                    &global_modifications,
                    ontologies,
                )
            );
            peptidoforms.push(peptidoform);
            start = tail;
        }

        if peptidoforms.is_empty() {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Error,
                    "No peptide found",
                    "The peptide definition is empty",
                    base_context.clone().add_highlight((0, range)),
                ),
            );
            Err(errors)
        } else {
            Ok((Self(peptidoforms), errors))
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    fn parse_peptidoform<'a, const STRICT: bool>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        global_modifications: &[GlobalModification],
        ontologies: &Ontologies,
    ) -> ParserResult<'a, (PeptidoformIon, usize), BasicKind> {
        let mut index: usize = range.start;
        let mut peptides = Vec::new();
        let mut ending = End::CrossLink;
        let mut cross_link_lookup = Vec::new();
        // Grouped on cross link id and stores peptide id, sequence index
        let mut cross_links_found = BTreeMap::new();
        let mut errors = Vec::new();

        // Parse any following cross-linked species
        while index < line.len() && ending == End::CrossLink && index < range.end {
            let mut result = handle!(
                errors,
                Self::parse_linear_peptide::<STRICT>(
                    base_context,
                    line,
                    index..range.end,
                    ontologies,
                    &mut cross_link_lookup,
                )
            );
            if !result
                .peptide
                .apply_global_modifications(global_modifications)
            {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global isotope modification",
                        "There is an invalid global isotope modification",
                        base_context.clone().add_highlight((0, range.clone())),
                    ),
                );
            } else if result.peptide.is_empty() {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "No peptide found",
                        "The peptidoform definition is empty",
                        base_context.clone().add_highlight((0, range)),
                    ),
                );
                return Err(errors);
            }
            peptides.push(result.peptide);
            index = result.index;
            ending = result.ending;
            for cross_link in result.cross_links {
                cross_links_found
                    .entry(cross_link.0)
                    .or_insert(Vec::new())
                    .push((peptides.len() - 1, cross_link.1));
            }
        }

        if peptides.is_empty() {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Error,
                    "No peptide found",
                    "The peptidoform definition is empty",
                    base_context.clone().add_highlight((0, range)),
                ),
            );
            Err(errors)
        } else {
            let peptidoform = handle!(
                errors,
                super::validate::cross_links(peptides, cross_links_found, &cross_link_lookup, line)
            );

            if errors.iter().any(|e| e.get_kind().is_error(())) {
                Err(errors)
            } else {
                Ok(((peptidoform, index), errors))
            }
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    #[expect(clippy::missing_panics_doc)] // Can not panic
    fn parse_linear_peptide<'a, const STRICT: bool>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        ontologies: &Ontologies,
        cross_link_lookup: &mut CrossLinkLookup,
    ) -> ParserResult<'a, LinearPeptideResult, BasicKind> {
        let mut errors = Vec::new();
        let index: usize = range.start;
        if line.trim().is_empty() {
            return Err(vec![BoxedError::new(
                BasicKind::Error,
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                base_context.clone().add_highlight((0, index, 1)),
            )]);
        }
        let mut peptide = Peptidoform::default();
        let chars: &[u8] = line.as_bytes();
        let mut c_term = false;
        let mut ambiguous_aa_counter = std::num::NonZeroU32::MIN;
        let mut ambiguous_aa = None;
        let mut ambiguous_lookup = Vec::new();
        let mut cross_link_found_positions: Vec<(usize, SequencePosition)> = Vec::new();
        let mut ambiguous_found_positions: Vec<(
            SequencePosition,
            bool,
            usize,
            Option<OrderedFloat<f64>>,
        )> = Vec::new();
        let mut ranged_unknown_position_modifications = Vec::new();
        let mut ending = End::Empty;

        // Unknown position mods
        let (index, unknown_position_modifications) = handle!(
            errors,
            global_unknown_position_mods::<STRICT>(
                base_context,
                index..range.end,
                line,
                ontologies,
                &mut ambiguous_lookup
            )
        );

        // Labile modification(s)
        let (mut index, labile) = handle!(
            errors,
            labile_modifications::<STRICT>(base_context, line, index, ontologies)
        );
        peptide = peptide.labile(labile);

        // N term modifications
        let (end_index, n_term_mods) = handle!(
            errors,
            multiple_mods::<STRICT>(
                base_context,
                line,
                index..range.end,
                ontologies,
                &mut ambiguous_lookup,
                cross_link_lookup,
            )
        );
        index = end_index;

        if chars.get(index) == Some(&b'-') {
            index += 1;
            for (modification, settings, _r) in n_term_mods {
                place_modification(
                    modification,
                    settings,
                    SequencePosition::NTerm,
                    &mut peptide,
                    &mut ambiguous_found_positions,
                    &mut cross_link_found_positions,
                );
            }
        }

        // Rest of the sequence
        let mut braces_start = None; // Sequence index where the last unopened braces started
        while index < chars.len() && index < range.end {
            match (c_term, chars[index]) {
                (false, b'(') if chars.get(index + 1) == Some(&b'?') => {
                    if braces_start.is_some() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ambiguous amino acid set",
                                "Ambiguous amino acid sets cannot be nested within ranged ambiguous modifications",
                                base_context.clone().add_highlight((0, index, 1)),
                            ),
                        );
                        return Err(errors);
                    }
                    if ambiguous_aa.is_some() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ambiguous amino acid set",
                                "Ambiguous amino acid sets cannot be nested within ambiguous amino acid sets",
                                base_context.clone().add_highlight((0, index, 1)),
                            ),
                        );
                        return Err(errors);
                    }
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter = handle!(single errors, ambiguous_aa_counter.checked_add(1).ok_or_else(|| {BoxedError::new(BasicKind::Error,
                        "Invalid ambiguous amino acid set",
                        format!("There are too many ambiguous amino acid sets, there can only be {} in one linear peptide", std::num::NonZeroU32::MAX),
                        base_context.clone().add_highlight((0, index, 1)))}));
                    index += 2;
                }
                (false, b')') if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
                    index += 1;
                }
                (false, b'(') => {
                    if braces_start.is_some() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ranged ambiguous modification",
                                "Ranged ambiguous modifications cannot be nested within ranged ambiguous modifications",
                                base_context.clone().add_highlight((0, index, 1)),
                            ),
                        );
                        return Err(errors);
                    }
                    if ambiguous_aa.is_some() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ranged ambiguous modification",
                                "Ranged ambiguous modifications cannot be nested within ambiguous amino acid sets",
                                base_context.clone().add_highlight((0, index, 1)),
                            ),
                        );
                        return Err(errors);
                    }
                    braces_start = Some(peptide.len());
                    index += 1;
                }
                (false, b')') if braces_start.is_some() => {
                    let start = braces_start.unwrap();
                    if start == peptide.len() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ranged ambiguous modification",
                                "The ranged ambiguous modification is placed on an empty range",
                                base_context
                                    .clone()
                                    .add_highlight((0, index - 1..index + 2)),
                            ),
                        );
                        return Err(errors);
                    }
                    braces_start = None;
                    index += 1;
                    let (end_index, mods) = handle!(
                        errors,
                        multiple_mods::<STRICT>(
                            base_context,
                            line,
                            index..range.end,
                            ontologies,
                            &mut ambiguous_lookup,
                            cross_link_lookup
                        )
                    );
                    index = end_index;
                    for (m, _settings, r) in mods {
                        match m.defined() {
                            Some(m) => ranged_unknown_position_modifications.push((
                            start,
                            peptide.len().saturating_sub(1),
                            m,
                        )),
                        None => errors.push(BoxedError::new(BasicKind::Error,
                            "Invalid ranged ambiguous modification",
                            "A ranged ambiguous modification has to be fully defined, so no ambiguous modification is allowed",
                            base_context.clone().add_highlight((0, r)))),
                        }
                    }
                }
                (false, b'/') => {
                    // Chimeric peptide
                    if chars.get(index + 1) == Some(&b'/') {
                        index += 2; // Potentially this can be followed by another peptide
                        ending = End::CrossLink;
                    } else {
                        let (buf, charge_carriers) =
                            handle!(errors, parse_charge_state_2_0(base_context, line, index));
                        index = buf;
                        peptide = peptide.charge_carriers(Some(charge_carriers));
                        if index < chars.len() && chars[index] == b'+' {
                            index += 1; // Potentially this can be followed by another peptide
                            ending = End::Chimeric;
                        }
                    }
                    break;
                }
                (is_c_term, b'[') => {
                    let end_index = handle!(single errors, end_of_enclosure(line, index + 1, b'[', b']')
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid modification",
                            "No valid closing delimiter",
                            base_context.clone().add_highlight((0, index, 1)),
                        )
                    }));
                    let (modification, settings) = handle!(
                        errors,
                        SimpleModificationInner::pro_forma_main::<STRICT>(
                            base_context,
                            line,
                            index + 1..end_index,
                            &mut ambiguous_lookup,
                            cross_link_lookup,
                            ontologies,
                        )
                    );
                    let start_index = index + 1;
                    index = end_index + 1;
                    if is_c_term {
                        place_modification(
                            modification,
                            settings,
                            SequencePosition::CTerm,
                            &mut peptide,
                            &mut ambiguous_found_positions,
                            &mut cross_link_found_positions,
                        );

                        let (end_index, mods) = handle!(
                            errors,
                            multiple_mods::<STRICT>(
                                base_context,
                                line,
                                index..range.end,
                                ontologies,
                                &mut ambiguous_lookup,
                                cross_link_lookup
                            )
                        );
                        index = end_index;

                        for (modification, settings, _r) in mods {
                            place_modification(
                                modification,
                                settings,
                                SequencePosition::CTerm,
                                &mut peptide,
                                &mut ambiguous_found_positions,
                                &mut cross_link_found_positions,
                            );
                        }

                        if index + 1 < chars.len()
                            && chars[index] == b'/'
                            && chars[index + 1] != b'/'
                        {
                            let (buf, charge_carriers) =
                                handle!(errors, parse_charge_state_2_0(base_context, line, index));
                            index = buf;
                            peptide = peptide.charge_carriers(Some(charge_carriers));
                        }
                        if index < chars.len() && chars[index] == b'+' {
                            index += 1; // If a peptide in a chimeric definition contains a C terminal modification
                            ending = End::Chimeric;
                        } else if index + 1 < chars.len() && chars[index..=index + 1] == *b"//" {
                            index += 2; // If a peptide in a cross-linked definition contains a C terminal modification
                            ending = End::CrossLink;
                        }
                        c_term = false; // Fix false negative errors on single ending hyphen
                        break;
                    }

                    if let Some((sequence_index, _aa)) =
                        peptide.sequence_mut().iter_mut().enumerate().next_back()
                    {
                        place_modification(
                            modification,
                            settings,
                            SequencePosition::Index(sequence_index),
                            &mut peptide,
                            &mut ambiguous_found_positions,
                            &mut cross_link_found_positions,
                        );
                    } else {
                        handle!(
                        fail errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid modification",
                                "A modification cannot be placed before any amino acid, did you want to use an N terminal modification ('[mod]-AA..')? or did you want a modification of unknown position ('[mod]?AA..')?",
                                base_context.clone().add_highlight((0, start_index, index - start_index - 1)),
                            )
                        );
                    }
                }
                (false, b'-') => {
                    c_term = true;
                    index += 1;
                }
                (false, b'+') => {
                    // Chimeric spectrum stop for now, remove the plus
                    index += 1;
                    ending = End::Chimeric;
                    break;
                }
                (false, ch) => {
                    peptide.sequence_mut().push(SequenceElement::new(
                        handle!(single errors, CheckedAminoAcid::<SemiAmbiguous>::try_from(ch)
                        .map_err(|()| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid amino acid",
                                "This character is not a valid amino acid",
                                base_context.clone().add_highlight((0, index, 1)),
                            )
                        }))
                        .into(),
                        ambiguous_aa,
                    ));
                    if STRICT && ch.is_ascii_lowercase() {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Warning,
                                "Lowercase amino acid",
                                "Use uppercase for amino acids",
                                base_context.clone().add_highlight((0, index, 1)),
                            ),
                        );
                    }
                    index += 1;
                }
                (true, _) => {
                    handle!(
                        fail
                        errors,
                        BoxedError::new(
                            BasicKind::Error,
                            "Parsing error",
                            "A singular hyphen cannot exist ('-'), if a C terminal modification is intended use 'SEQ-[MOD]'",
                            base_context.clone().add_highlight((0, index, 1)),
                        )
                    );
                }
            }
        }
        if c_term {
            handle!(
                fail
                errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid peptide",
                    "A single hyphen cannot end the definition, if a C terminal modification is intended use 'SEQ-[MOD]'",
                    base_context.clone().add_highlight((0, line.len().saturating_sub(2), 1)),
                )
            );
        }
        if let Some(pos) = braces_start {
            handle!(
                fail
                errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid peptide",
                    format!("Unclosed brace at amino acid position {pos}"),
                    base_context.clone().add_highlight((0, range)),
                )
            );
        }
        if ambiguous_aa.is_some() {
            handle!(
                fail
                errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid peptide",
                    "Unclosed ambiguous amino acid group",
                    base_context.clone().add_highlight((0, range)),
                )
            );
        }
        if peptide.is_empty() {
            handle!(
                fail
                errors,
                BoxedError::new(
                    BasicKind::Error,
                    "No amino acids found",
                    "The peptide definition is empty",
                    base_context.clone().add_highlight((0, range)),
                )
            );
        }

        // Fill in ambiguous positions, ambiguous contains (index, preferred, id, localisation_score)
        for (id, ambiguous) in ambiguous_found_positions
            .into_iter()
            .into_group_map_by(|aa| aa.2)
        {
            let positions = ambiguous
                .iter()
                .map(|(index, _, _, score)| (*index, *score))
                .collect_vec();
            let preferred = ambiguous.iter().find_map(|p| p.1.then_some(p.0));
            if let Some(modification) = ambiguous_lookup[id].modification.clone() {
                if !peptide.add_ambiguous_modification(
                    modification,
                    Some(ambiguous_lookup[id].name.clone()),
                    &positions,
                    preferred,
                    None,
                    true,
                ) {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            BasicKind::Error,
                            "Modification of unknown position cannot be placed",
                            format!(
                                "There is no position where this ambiguous modification {} can be placed based on the placement rules in the database.",
                                ambiguous_lookup[id].name
                            ),
                            base_context.clone().add_highlight((0, range.clone())),
                        ),
                    );
                }
            } else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ambiguous modification",
                        format!(
                            "Ambiguous modification {} did not have a definition for the actual modification",
                            ambiguous_lookup[id].name
                        ),
                        base_context.clone().add_highlight((0, range.clone())),
                    ),
                );
            }
        }

        if let Err(errs) = peptide
            .apply_unknown_position_modification(&unknown_position_modifications, &ambiguous_lookup)
        {
            combine_errors(&mut errors, errs);
        }
        if let Err(errs) = peptide
            .apply_ranged_unknown_position_modification(&ranged_unknown_position_modifications)
        {
            combine_errors(&mut errors, errs);
        }
        combine_errors(
            &mut errors,
            peptide.enforce_modification_rules_with_context(base_context),
        );

        if errors.iter().any(|e| e.get_kind().is_error(())) {
            Err(errors)
        } else {
            Ok((
                LinearPeptideResult {
                    peptide,
                    index,
                    ending,
                    cross_links: cross_link_found_positions,
                },
                errors,
            ))
        }
    }
}

/// Parse multiple modifications at the given position
/// # Errors
/// When no end of a mod enclosure could be found.
fn multiple_mods<'a, const STRICT: bool>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    ontologies: &Ontologies,
    ambiguous_lookup: &mut AmbiguousLookup,
    cross_link_lookup: &mut CrossLinkLookup,
) -> ParserResult<'a, (usize, Vec<(ReturnModification, MUPSettings, Range<usize>)>), BasicKind> {
    let mut errors = Vec::new();
    let mut mods = Vec::new();
    let mut temp_index = range.start;
    while line.as_bytes().get(temp_index) == Some(&b'[') {
        let end_index = handle!(single errors, end_of_enclosure(line, temp_index + 1, b'[', b']').ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid modification",
                    "No valid closing delimiter",
                    base_context.clone().add_highlight((0, temp_index, 1))
                )
        }));
        let (m, mup) = handle!(
            errors,
            SimpleModificationInner::pro_forma_main::<STRICT>(
                base_context,
                line,
                temp_index + 1..end_index,
                ambiguous_lookup,
                cross_link_lookup,
                ontologies,
            )
        );
        mods.push((m, mup, temp_index..end_index + 1));

        temp_index = end_index + 1;
    }
    Ok(((temp_index, mods), errors))
}

/// Place modification on the right location and handle all bookkeeping
fn place_modification(
    modification: ReturnModification,
    _settings: MUPSettings,
    position: SequencePosition,
    peptidoform: &mut Peptidoform<Linear>,
    ambiguous_found_positions: &mut Vec<(SequencePosition, bool, usize, Option<OrderedFloat<f64>>)>,
    cross_link_found_positions: &mut Vec<(usize, SequencePosition)>,
) {
    match modification {
        ReturnModification::Defined(simple) => {
            peptidoform.add_simple_modification(position, simple);
        }
        ReturnModification::CrossLinkReferenced(id) => {
            cross_link_found_positions.push((id, position));
        }
        ReturnModification::Ambiguous(id, localisation_score, preferred) => {
            ambiguous_found_positions.push((position, preferred, id, localisation_score));
        }
    }
}

/// Parse global modifications
/// # Errors
/// If the global modifications are not defined to the specification
pub(super) fn global_modifications<'a, const STRICT: bool>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    ontologies: &Ontologies,
) -> ParserResult<'a, (usize, Vec<GlobalModification>), BasicKind> {
    let mut index = range.start;
    let mut errors = Vec::new();
    let chars = line.as_bytes();
    let mut global_modifications = Vec::new();
    while index < chars.len() && chars[index] == b'<' && index < range.end {
        let Some(end_index) = end_of_enclosure_with_brackets(line, index + 1, b'<', b'>') else {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Global modification not closed",
                    "A global modification should be closed with a closing angle bracket '>'",
                    base_context.clone().add_highlight((0, index, 1)),
                ),
            );
            return Err(errors);
        };
        let start_index = index;
        index = end_index + 1; // Already set the index for if the parsing fails
        if let Some(at_index) = next_char(chars, start_index, b'@') {
            let at_index = at_index + 1;
            if at_index > end_index {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        "A global modification should have an at '@' sign inside the enclosing angle brackets '<>'",
                        base_context.clone().add_highlight((
                            0,
                            start_index + 1,
                            at_index - start_index - 1,
                        )),
                    ),
                );
                continue;
            }
            if chars[start_index + 1] != b'[' || chars[at_index - 2] != b']' {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        "A global modification should always be enclosed in square brackets '[]'",
                        base_context.clone().add_highlight((
                            0,
                            start_index + 1,
                            at_index - start_index - 2,
                        )),
                    ),
                );
                continue;
            }
            let modification = handle!(
                errors,
                SimpleModificationInner::pro_forma_main::<STRICT>(
                    base_context,
                    line,
                    start_index + 2..at_index - 2,
                    &mut Vec::new(),
                    &mut Vec::new(),
                    ontologies,
                )
                .map(|((m, _mup_settings), mut warnings)| {
                    let w = warnings.clone();
                    m.defined().map(|m| (m, w)).ok_or_else(|| {
                        warnings.push(BoxedError::new(
                            BasicKind::Error,
                            "Invalid global modification",
                            "A global modification cannot be ambiguous or a cross-linker",
                            base_context.clone().add_highlight((
                                0,
                                start_index + 2,
                                at_index - start_index - 4,
                            )),
                        ));
                        warnings
                    })
                })
                .flat_err()
            );
            let rules = match parse_placement_rules(base_context, line, at_index..end_index) {
                Ok(rules) => rules,
                Err(err) => {
                    combine_error(&mut errors, err);
                    continue;
                }
            };
            global_modifications.extend(
                rules
                    .into_iter()
                    .map(|r| GlobalModification::Fixed(r, modification.clone())),
            );
        } else if &line[start_index + 1..end_index].to_ascii_lowercase() == "d" {
            global_modifications.push(GlobalModification::Isotope(Element::H, NonZeroU16::new(2)));
        } else {
            let num = &line[start_index + 1..end_index]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            let el = &line[start_index + 1 + num.len()..end_index];
            let el: Element = if let Ok(el) = el.try_into() {
                el
            } else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        "Could not determine the element",
                        base_context.clone().add_highlight((
                            0,
                            start_index + num.len(),
                            end_index - (start_index + 1 + num.len()),
                        )),
                    ),
                );
                continue;
            };
            let num = Some(match num.parse::<NonZeroU16>() {
                Ok(num) => num,
                Err(err) => {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid global modification",
                            format!("The isotope number is {}", explain_number_error(&err)),
                            base_context.clone().add_highlight((
                                0,
                                start_index + 1,
                                end_index - start_index,
                            )),
                        ),
                    );
                    continue;
                }
            });
            if el.is_valid(num) {
                global_modifications.push(GlobalModification::Isotope(el, num));
            } else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        format!(
                            "This element {el} does not have a defined weight {}",
                            num.map_or_else(String::new, |num| format!("for isotope {num}"))
                        ),
                        base_context.clone().add_highlight((
                            0,
                            start_index + 1,
                            end_index - start_index,
                        )),
                    ),
                );
            }
        }
    }
    if errors.iter().any(|e| e.get_kind().is_error(())) {
        Err(errors)
    } else {
        Ok(((index, global_modifications), errors))
    }
}

/// Parse a set of placement rules.
/// # Errors
/// When any rule is invalid.
pub(super) fn parse_placement_rules<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<Vec<PlacementRule>, BoxedError<'a, BasicKind>> {
    let mut result = Vec::new();
    for aa in line[range.clone()].split(',') {
        if aa.to_ascii_lowercase().starts_with("n-term") {
            if let Some((_, aa)) = aa.split_once(':') {
                result.push(PlacementRule::AminoAcid(
                    vec![TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid global modification",
                            "The location could not be read as an amino acid",
                            base_context.clone().add_highlight((0, range.clone())),
                        )
                    })?]
                    .into(),
                    Position::AnyNTerm,
                ));
            } else {
                result.push(PlacementRule::Terminal(Position::AnyNTerm));
            }
        } else if aa.to_ascii_lowercase().starts_with("c-term") {
            if let Some((_, aa)) = aa.split_once(':') {
                result.push(PlacementRule::AminoAcid(
                    vec![TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid global modification",
                            "The location could not be read as an amino acid",
                            base_context.clone().add_highlight((0, range.clone())),
                        )
                    })?]
                    .into(),
                    Position::AnyCTerm,
                ));
            } else {
                result.push(PlacementRule::Terminal(Position::AnyCTerm));
            }
        } else {
            result.push(PlacementRule::AminoAcid(
                vec![TryInto::<AminoAcid>::try_into(aa).map_err(|()| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        "The location could not be read as an amino acid",
                        base_context.clone().add_highlight((0, range.clone())),
                    )
                })?]
                .into(),
                Position::Anywhere,
            ));
        }
    }
    Ok(result)
}

/// MUP: `[Mod][Mod]?AAA`
/// If the text is recognised as an unknown mods list it is Some(...), if it has errors during parsing Some(Err(...))
/// The returned happy path contains the mods and the index from where to continue parsing.
/// # Errors
/// Give all errors when the text cannot be read as mods of unknown position.
pub(super) fn global_unknown_position_mods<'a, const STRICT: bool>(
    base_context: &Context<'a>,
    range: Range<usize>,
    line: &'a str,
    ontologies: &Ontologies,
    ambiguous_lookup: &mut AmbiguousLookup,
) -> ParserResult<'a, (usize, Vec<usize>), BasicKind> {
    let mut index = range.start;
    let mut modifications = Vec::new();
    let mut errs = Vec::new();
    let mut cross_link_lookup = Vec::new();

    // Parse until no new modifications are found
    while line.as_bytes().get(index) == Some(&b'[') {
        let start_index = index;
        index = end_of_enclosure(line, index + 1, b'[', b']').ok_or_else(|| {
            vec![BoxedError::new(
                BasicKind::Error,
                "Global unknown position modification not closed",
                "All global unknown position modifications should be closed with a closing square bracket ']'",
                base_context.clone().add_highlight((0, index, 1)),
            )]
        })? + 1;
        let id = match SimpleModificationInner::pro_forma_main::<STRICT>(
            base_context,
            line,
            start_index + 1..index - 1,
            ambiguous_lookup,
            &mut cross_link_lookup,
            ontologies,
        ) {
            Ok(((ReturnModification::Defined(m), settings), warnings)) => {
                combine_errors(&mut errs, warnings);
                let id = ambiguous_lookup.len();
                ambiguous_lookup.push(AmbiguousLookupEntry::new(format!("u{id}"), Some(m)));
                ambiguous_lookup[id].copy_settings(&settings);
                id
            }
            Ok(((ReturnModification::Ambiguous(id, _, _), settings), warnings)) => {
                combine_errors(&mut errs, warnings);
                ambiguous_lookup[id].copy_settings(&settings);
                id
            }
            Ok(((ReturnModification::CrossLinkReferenced(_), _), warnings)) => {
                combine_errors(&mut errs, warnings);
                errs.push(BoxedError::new(
                    BasicKind::Error,
                    "Invalid unknown position modification",
                    "A modification of unknown position cannot be a cross-link",
                    base_context
                        .clone()
                        .add_highlight((0, (start_index + 1)..index)),
                ));
                continue;
            }
            Err(e) => {
                combine_errors(&mut errs, e);
                continue;
            }
        };
        let number = if line.as_bytes().get(index) == Some(&b'^') {
            if let Some((len, num)) = next_num(line.as_bytes(), index + 1, false) {
                index += len + 1;
                if num < 0 {
                    errs.push(
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid unknown position modification",
                            "A modification of unknown position with multiple copies cannot have more a negative number of copies",
                            base_context.clone().add_highlight((0, index, 1))));
                    0
                } else if num > i16::MAX as isize {
                    errs.push(
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid unknown position modification",
                            format!("A modification of unknown position with multiple copies cannot have more then {} copies", i16::MAX),
                             base_context.clone().add_highlight((0, index, 1))));
                    0
                } else {
                    num as usize
                }
            } else {
                errs.push(
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid unknown position modification",
                        "A modification of unknown position with multiple copies needs the copy number after the caret ('^') symbol",
                         base_context.clone().add_highlight((0, index, 1))));
                0
            }
        } else {
            1
        };
        if number > 1 {
            let group_name = ambiguous_lookup.len();
            ambiguous_lookup[id].group = Some(group_name);
            modifications.push(id);
            for _ in 0..number {
                modifications.push(ambiguous_lookup.len());
                ambiguous_lookup.push(ambiguous_lookup[id].clone());
            }
        } else {
            modifications.push(id);
        }
    }
    if line.as_bytes().get(index) == Some(&b'?') {
        if errs.iter().any(|e| e.get_kind().is_error(())) {
            Err(errs)
        } else {
            Ok(((index + 1, modifications), errs))
        }
    } else {
        ambiguous_lookup.clear(); // Any ambiguous N terminal modification was incorrectly already added to the lookup
        Ok(((range.start, Vec::new()), errs))
    }
}

/// Parse labile modifications `{mod}{mod2}`. These are assumed to fall off from the peptide in the MS.
/// # Errors
/// If the mods are not followed by a closing brace. Or if the mods are ambiguous.
fn labile_modifications<'a, const STRICT: bool>(
    base_context: &Context<'a>,
    line: &'a str,
    mut index: usize,
    ontologies: &Ontologies,
) -> ParserResult<'a, (usize, Vec<SimpleModification>), BasicKind> {
    let mut errors = Vec::new();
    let chars = line.as_bytes();
    let mut labile = Vec::new();
    while chars.get(index) == Some(&b'{') {
        let Some(end_index) = end_of_enclosure(line, index + 1, b'{', b'}') else {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid labile modification",
                    "No valid closing delimiter, a labile modification should be closed by '}'",
                    base_context.clone().add_highlight((0, index, 1)),
                ),
            );
            return Err(errors);
        };

        match SimpleModificationInner::pro_forma_main::<STRICT>(
            base_context,
            line,
            index + 1..end_index,
            &mut Vec::new(),
            &mut Vec::new(),
            ontologies,
        ) {
            Ok(((ReturnModification::Defined(m), _), warnings)) => {
                combine_errors(&mut errors, warnings);
                labile.push(m);
            }
            Ok((_, warnings)) => {
                combine_errors(&mut errors, warnings);
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid labile modification",
                        "A labile modification cannot be ambiguous or a cross-linker",
                        base_context
                            .clone()
                            .add_highlight((0, index + 1, end_index - 1 - index)),
                    ),
                );
            }
            Err(e) => combine_errors(&mut errors, e),
        }
        index = end_index + 1;
    }
    if errors.iter().any(|e| e.get_kind().is_error(())) {
        Err(errors)
    } else {
        Ok(((index, labile), errors))
    }
}

/// Parse a charge state `/2` or more complex ones like `/2[+2Na+]`.
/// Assumes the text starts with `/`.
/// # Errors
/// If the charge state is not following the specification.
/// # Panics
/// Panics if the text is not UTF-8.
pub(super) fn parse_charge_state_2_0<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    index: usize,
) -> ParserResult<'a, (usize, MolecularCharge), BasicKind> {
    let mut errors = Vec::new();
    let chars = line.as_bytes();
    let (charge_len, total_charge) = handle!(single errors, next_num(chars, index + 1, false)
    .ok_or_else(|| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid peptide charge state",
            "There should be a number dictating the total charge of the peptide",
            base_context.clone().add_highlight((0, index + 1, 1))
        )
    }));
    if chars.get(index + 1 + charge_len) == Some(&b'[') {
        let end_index = handle!(single errors, end_of_enclosure(line, index + 2 + charge_len, b'[', b']')
        .ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid adduct ion",
                "No valid closing delimiter",
                base_context.clone().add_highlight((0, index + 2 + charge_len, 1))
            )
        }));
        let mut offset = index + 2 + charge_len;
        let mut charge_carriers: Vec<(isize, MolecularFormula)> = Vec::new();
        let mut found_charge: isize = 0;

        for set in chars[index + 2 + charge_len..end_index].split(|c| *c == b',') {
            // num
            let (count_len, count) = handle!(single errors, next_num(chars, offset, true)
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid adduct ion",
                    "Invalid adduct ion count",
                    base_context.clone().add_highlight((0, index, 1))
                )
            }));

            // charge
            let charge_len = set.iter().rev().take_while(|c| c.is_ascii_digit()).count();
            let charge = if charge_len == 0 {
                1
            } else {
                handle!(single errors, line[offset + set.len() - charge_len..offset + set.len()]
                .parse::<i32>()
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid adduct ion",
                        format!("The adduct ion number {err}"),
                        base_context.clone().add_highlight((0, offset + set.len() - charge_len, charge_len))
                    )
                }))
            };
            let (charge_len, charge) = match (set.len() - charge_len)
                .checked_sub(1)
                .and_then(|i| set.get(i))
            {
                Some(b'+') => (charge_len + 1, charge),
                Some(b'-') => (charge_len + 1, -charge),
                _ => {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid adduct ion",
                            "The adduct ion number should be preceded by a sign",
                            base_context.clone().add_highlight((
                                0,
                                offset + set.len() - charge_len - 1,
                                1,
                            )),
                        ),
                    );
                    return Err(errors);
                }
            };

            // Check for empty formula
            if count_len + charge_len == set.len() {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid adduct ion",
                        "The adduct ion should have a formula defined",
                        base_context.clone().add_highlight((0, offset, set.len())),
                    ),
                );
                return Err(errors);
            }

            // formula
            let formula_range = offset + count_len..offset + set.len() - charge_len;
            let formula = if !formula_range.is_empty()
                && line[formula_range].eq_ignore_ascii_case("e")
            {
                MolecularFormula::new(&[(Element::Electron, None, -charge)], &[]).unwrap()
            } else {
                let mut formula = handle!(single errors, MolecularFormula::pro_forma_inner::<true, false>(
                    base_context,
                    line,
                    offset + count_len..offset + set.len() - charge_len,
                ));
                let _ = formula.add((
                    Element::Electron,
                    None,
                    formula.charge().value as i32 - charge,
                ));
                formula
            };

            // Deduplicate
            if let Some((amount, _)) = charge_carriers.iter_mut().find(|(_, f)| *f == formula) {
                *amount = handle!(single errors, amount
                .checked_add(count)
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid peptide charge amount",
                        "The peptide charge amount is too big to store inside an isize",
                        base_context.clone().add_highlight((0, index, offset))
                    )
                }));
            } else {
                charge_carriers.push((count, formula));
            }

            offset += set.len() + 1;
            found_charge = handle!(single errors, found_charge
            .checked_add(
                handle!(single errors, count
                    .checked_mul(charge as isize)
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid peptide charge state",
                            "The peptide charge state is too big to store inside an isize",
                            base_context.clone().add_highlight((0, index, offset))
                        )
                    })),
            )
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid peptide charge state",
                    "The peptide charge state is too big to store inside an isize",
                    base_context.clone().add_highlight((0, index, offset))
                )
            }));
        }
        if total_charge == found_charge {
            Ok((
                (end_index + 1, MolecularCharge::new(&charge_carriers)),
                errors,
            ))
        } else {
            combine_error(
                &mut errors,
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid peptide charge state",
                    "The peptide charge state number has to be equal to the sum of all separate adduct ions",
                    base_context.clone().add_highlight((0, index, offset)),
                ),
            );
            Err(errors)
        }
    } else {
        // If no adduct ions are provided assume it is just protons
        Ok((
            (
                index + charge_len + 1,
                MolecularCharge::proton(crate::system::isize::Charge::new::<crate::system::e>(
                    total_charge,
                )),
            ),
            errors,
        ))
    }
}
