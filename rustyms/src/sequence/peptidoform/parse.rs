use std::{collections::BTreeMap, num::NonZeroU16};

use context_error::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;

use crate::{
    chemistry::{Element, MolecularCharge, MolecularFormula},
    helper_functions::*,
    ontology::CustomDatabase,
    sequence::{
        AmbiguousLookup, AmbiguousLookupEntry, AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon,
        CrossLinkLookup, Linked, Modification, Peptidoform, PeptidoformIon, PlacementRule,
        Position, SequenceElement, SequencePosition, SimpleModification, SimpleModificationInner,
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
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        CompoundPeptidoformIon::pro_forma(value, custom_database)?
            .singular()
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Complex peptide found",
                    "A linear peptide was expected but a chimeric peptide was found.",
                    Context::show(value),
                )
            })
            .and_then(|p| {
                p.singular().ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Complex peptide found",
                        "A linear peptide was expected but a cross linked peptidoform was found.",
                        Context::show(value),
                    )
                })
            })
    }
}

impl PeptidoformIon {
    /// Parse a peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string. Or when the string has multiple peptidoforms.
    pub fn pro_forma<'a>(
        value: &'a str,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        CompoundPeptidoformIon::pro_forma(value, custom_database)?
            .singular()
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Complex peptide found",
                    "A linear peptide was expected but a chimeric peptide was found.",
                    Context::show(value),
                )
            })
    }
}

impl CompoundPeptidoformIon {
    /// Parse a compound peptidoform in the [ProForma specification](https://github.com/HUPO-PSI/ProForma).
    ///
    /// # Errors
    /// It fails when the string is not a valid ProForma string.
    pub fn pro_forma<'a>(
        value: &'a str,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let mut peptidoforms = Vec::new();
        // Global modification(s)
        let (mut start, global_modifications) = global_modifications(value, 0, custom_database)?;
        let (peptidoform, tail) =
            Self::parse_peptidoform(value, start, &global_modifications, custom_database)?;
        start = tail;
        peptidoforms.push(peptidoform);

        // Parse any following chimeric species
        while start < value.len() {
            let (peptidoform, tail) =
                Self::parse_peptidoform(value, start, &global_modifications, custom_database)?;
            peptidoforms.push(peptidoform);
            start = tail;
        }

        if peptidoforms.is_empty() {
            Err(BoxedError::new(
                BasicKind::Error,
                "No peptide found",
                "The peptide definition is empty",
                Context::full_line(0, value),
            ))
        } else {
            Ok(Self(peptidoforms))
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    fn parse_peptidoform<'a>(
        line: &'a str,
        mut index: usize,
        global_modifications: &[GlobalModification],
        custom_database: Option<&CustomDatabase>,
    ) -> Result<(PeptidoformIon, usize), BoxedError<'a, BasicKind>> {
        let mut peptides = Vec::new();
        let mut ending = End::CrossLink;
        let mut cross_link_lookup = Vec::new();
        // Grouped on cross link id and stores peptide id, sequence index
        let mut cross_links_found = BTreeMap::new();

        // Parse any following cross-linked species
        while index < line.len() && ending == End::CrossLink {
            let mut result =
                Self::parse_linear_peptide(line, index, custom_database, &mut cross_link_lookup)?;
            if !result
                .peptide
                .apply_global_modifications(global_modifications)
            {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid global isotope modification",
                    "There is an invalid global isotope modification",
                    Context::full_line(0, line),
                ));
            } else if result.peptide.is_empty() {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "No peptide found",
                    "The peptide definition is empty",
                    Context::full_line(0, line),
                ));
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
            Err(BoxedError::new(
                BasicKind::Error,
                "No peptide found",
                "The peptidoform definition is empty",
                Context::full_line(0, line),
            ))
        } else {
            let peptidoform = super::validate::cross_links(
                peptides,
                cross_links_found,
                &cross_link_lookup,
                line,
            )?;
            Ok((peptidoform, index))
        }
    }

    /// # Errors
    /// It returns an error if the line is not a supported ProForma line.
    #[expect(clippy::missing_panics_doc)] // Can not panic
    fn parse_linear_peptide<'a>(
        line: &'a str,
        mut index: usize,
        custom_database: Option<&CustomDatabase>,
        cross_link_lookup: &mut CrossLinkLookup,
    ) -> Result<LinearPeptideResult, BoxedError<'a, BasicKind>> {
        if line.trim().is_empty() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, index, 1),
            ));
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
        let mut unknown_position_modifications = Vec::new();
        let mut ranged_unknown_position_modifications = Vec::new();
        let mut ending = End::Empty;

        // Unknown position mods
        if let Some(result) =
            global_unknown_position_mods(chars, index, line, custom_database, &mut ambiguous_lookup)
        {
            let (buf, mods) = result.map_err(|errors| {
                BoxedError::new(
                    BasicKind::Error,
                    "Some unknown position modifications are invalid",
                    "See the underlying errors for more details.",
                    Context::show(line),
                )
                .add_underlying_errors(errors)
            })?;
            index = buf;

            unknown_position_modifications = mods;
        }

        // Labile modification(s)
        let (mut index, labile) = labile_modifications(line, index, custom_database)?;
        peptide = peptide.labile(labile);

        // N term modifications
        let mut n_term_mods = Vec::new();
        let mut n_end_seen = false;
        let mut temp_index = index;
        while chars.get(temp_index) == Some(&b'[') {
            let end_index =
                end_of_enclosure(line, temp_index + 1, b'[', b']').ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid N term modification",
                        "No valid closing delimiter",
                        Context::line(None, line, temp_index, 1),
                    )
                })?;
            n_term_mods.push(SimpleModificationInner::parse_pro_forma(
                line,
                temp_index + 1..end_index,
                &mut ambiguous_lookup,
                cross_link_lookup,
                custom_database,
            )?);

            temp_index = end_index + 1;

            if chars.get(temp_index) == Some(&b'-') {
                index = temp_index + 1;
                n_end_seen = true;
                break;
            }
        }

        if n_end_seen {
            for (m, _mup) in n_term_mods {
                match m {
                    ReturnModification::Defined(simple) => peptide.add_simple_n_term(simple),
                    ReturnModification::CrossLinkReferenced(id) => {
                        cross_link_found_positions.push((id, SequencePosition::NTerm));
                    }
                    ReturnModification::Ambiguous(id, localisation_score, preferred) => {
                        ambiguous_found_positions.push((
                            SequencePosition::NTerm,
                            preferred,
                            id,
                            localisation_score,
                        ));
                    }
                }
            }
        }

        // Rest of the sequence
        let mut braces_start = None; // Sequence index where the last unopened braces started
        while index < chars.len() {
            match (c_term, chars[index]) {
                (false, b'(') if chars.get(index + 1) == Some(&b'?') => {
                    if braces_start.is_some() {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ambiguous amino acid set",
                            "Ambiguous amino acid sets cannot be nested within ranged ambiguous modifications",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    if ambiguous_aa.is_some() {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ambiguous amino acid set",
                            "Ambiguous amino acid sets cannot be nested within ambiguous amino acid sets",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    ambiguous_aa = Some(ambiguous_aa_counter);
                    ambiguous_aa_counter = ambiguous_aa_counter.checked_add(1).ok_or_else(|| BoxedError::new(BasicKind::Error,
                        "Invalid ambiguous amino acid set",
                        format!("There are too many ambiguous amino acid sets, there can only be {} in one linear peptide", std::num::NonZeroU32::MAX),
                        Context::line(None, line, index, 1),
                    ))?;
                    index += 2;
                }
                (false, b')') if ambiguous_aa.is_some() => {
                    ambiguous_aa = None;
                    index += 1;
                }
                (false, b'(') => {
                    if braces_start.is_some() {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ranged ambiguous modification",
                            "Ranged ambiguous modifications cannot be nested within ranged ambiguous modifications",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    if ambiguous_aa.is_some() {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ranged ambiguous modification",
                            "Ranged ambiguous modifications cannot be nested within ambiguous amino acid sets",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    braces_start = Some(peptide.len());
                    index += 1;
                }
                (false, b')') if braces_start.is_some() => {
                    let start = braces_start.unwrap();
                    braces_start = None;
                    index += 1;
                    while chars.get(index) == Some(&b'[') {
                        let end_index =
                            end_of_enclosure(line, index + 1, b'[', b']').ok_or_else(|| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid ranged ambiguous modification",
                                    "No valid closing delimiter",
                                    Context::line(None, line, index, 1),
                                )
                            })?;
                        let modification = SimpleModificationInner::parse_pro_forma(
                            line, index + 1..end_index,
                            &mut ambiguous_lookup, cross_link_lookup, custom_database,
                        )?.0.defined().ok_or_else(|| BoxedError::new(BasicKind::Error,
                            "Invalid ranged ambiguous modification",
                            "A ranged ambiguous modification has to be fully defined, so no ambiguous modification is allowed",
                            Context::line(None, line, index, 1),
                        ))?;
                        index = end_index + 1;
                        ranged_unknown_position_modifications.push((
                            start,
                            peptide.len().saturating_sub(1),
                            modification,
                        ));
                    }
                }
                (false, b'/') => {
                    // Chimeric peptide
                    if chars.get(index + 1) == Some(&b'/') {
                        index += 2; // Potentially this can be followed by another peptide
                        ending = End::CrossLink;
                    } else {
                        let (buf, charge_carriers) = parse_charge_state(line, index)?;
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
                    let end_index =
                        end_of_enclosure(line, index + 1, b'[', b']').ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid modification",
                                "No valid closing delimiter",
                                Context::line(None, line, index, 1),
                            )
                        })?;
                    let (modification, _) = SimpleModificationInner::parse_pro_forma(
                        line,
                        index + 1..end_index,
                        &mut ambiguous_lookup,
                        cross_link_lookup,
                        custom_database,
                    )?;
                    let start_index = index + 1;
                    index = end_index + 1;
                    if is_c_term {
                        match modification {
                            ReturnModification::Defined(simple) => {
                                peptide.add_simple_c_term(simple);
                            }
                            ReturnModification::CrossLinkReferenced(id) => {
                                cross_link_found_positions.push((id, SequencePosition::CTerm));
                            }
                            ReturnModification::Ambiguous(id, localisation_score, preferred) => {
                                ambiguous_found_positions.push((
                                    SequencePosition::CTerm,
                                    preferred,
                                    id,
                                    localisation_score,
                                ));
                            }
                        }

                        while chars.get(index) == Some(&b'[') {
                            let end_index = end_of_enclosure(line, index + 1, b'[', b']')
                                .ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid C term modification",
                                        "No valid closing delimiter",
                                        Context::line(None, line, index, 1),
                                    )
                                })?;
                            let modification = SimpleModificationInner::parse_pro_forma(
                                line,
                                index + 1..end_index,
                                &mut ambiguous_lookup,
                                cross_link_lookup,
                                custom_database,
                            )?;

                            match modification.0 {
                                ReturnModification::Defined(simple) => {
                                    peptide.add_simple_c_term(simple);
                                }
                                ReturnModification::CrossLinkReferenced(id) => {
                                    cross_link_found_positions.push((id, SequencePosition::CTerm));
                                }
                                ReturnModification::Ambiguous(
                                    id,
                                    localisation_score,
                                    preferred,
                                ) => {
                                    ambiguous_found_positions.push((
                                        SequencePosition::CTerm,
                                        preferred,
                                        id,
                                        localisation_score,
                                    ));
                                }
                            }

                            index = end_index + 1;
                        }

                        if index + 1 < chars.len()
                            && chars[index] == b'/'
                            && chars[index + 1] != b'/'
                        {
                            let (buf, charge_carriers) = parse_charge_state(line, index)?;
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

                    if let Some((sequence_index, aa)) =
                        peptide.sequence_mut().iter_mut().enumerate().next_back()
                    {
                        match modification {
                            ReturnModification::Defined(m) => {
                                aa.modifications.push(Modification::Simple(m));
                            }
                            ReturnModification::Ambiguous(id, localisation_score, preferred) => {
                                ambiguous_found_positions.push((
                                    SequencePosition::Index(sequence_index),
                                    preferred,
                                    id,
                                    localisation_score,
                                ));
                            }
                            ReturnModification::CrossLinkReferenced(id) => {
                                cross_link_found_positions
                                    .push((id, SequencePosition::Index(sequence_index)));
                            }
                        }
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid modification",
                            "A modification cannot be placed before any amino acid, did you want to use an N terminal modification ('[mod]-AA..')? or did you want a modification of unknown position ('[mod]?AA..')?",
                            Context::line(None, line, start_index, index - start_index - 1),
                        ));
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
                        CheckedAminoAcid::<SemiAmbiguous>::try_from(ch)
                            .map_err(|()| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid amino acid",
                                    "This character is not a valid amino acid",
                                    Context::line(None, line, index, 1),
                                )
                            })?
                            .into(),
                        ambiguous_aa,
                    ));
                    index += 1;
                }
                (true, _) => {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Parsing error",
                        "A singular hyphen cannot exist ('-'), if this is part of a c-terminus follow the format 'AA-[modification]'",
                        Context::line(None, line, index, 1),
                    ));
                }
            }
        }
        if c_term {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid peptide",
                "A single hyphen cannot end the definition, if a C terminal modification is intended use 'SEQ-[MOD]'",
                Context::line(None, line, line.len().saturating_sub(2), 1),
            ));
        }
        if let Some(pos) = braces_start {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid peptide",
                format!("Unclosed brace at amino acid position {pos}"),
                Context::full_line(0, line),
            ));
        }
        if ambiguous_aa.is_some() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid peptide",
                "Unclosed ambiguous amino acid group",
                Context::full_line(0, line),
            ));
        }
        if peptide.is_empty() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "No amino acids found",
                "The peptide definition is empty",
                Context::full_line(0, line),
            ));
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
            if !peptide.add_ambiguous_modification(ambiguous_lookup[id].modification.clone().ok_or_else(||
                BoxedError::new(BasicKind::Error,
                    "Invalid ambiguous modification",
                    format!("Ambiguous modification {} did not have a definition for the actual modification", ambiguous_lookup[id].name),
                    Context::full_line(0, line),
                )
                )?, Some(ambiguous_lookup[id].name.clone()), &positions, preferred, None,  true) {
                return Err(BoxedError::new(BasicKind::Error,
                    "Modification of unknown position cannot be placed",
                    format!("There is no position where this ambiguous modification {} can be placed based on the placement rules in the database.", ambiguous_lookup[id].name),
                    Context::full_line(0, line),
                    ));
            }
        }

        peptide.apply_unknown_position_modification(
            &unknown_position_modifications,
            &ambiguous_lookup,
        )?;
        peptide
            .apply_ranged_unknown_position_modification(&ranged_unknown_position_modifications)?;
        peptide.enforce_modification_rules()?;

        Ok(LinearPeptideResult {
            peptide,
            index,
            ending,
            cross_links: cross_link_found_positions,
        })
    }
}

/// Parse global modifications
/// # Errors
/// If the global modifications are not defined to the specification
pub(super) fn global_modifications<'a>(
    line: &'a str,
    mut index: usize,
    custom_database: Option<&CustomDatabase>,
) -> Result<(usize, Vec<GlobalModification>), BoxedError<'a, BasicKind>> {
    let chars = line.as_bytes();
    let mut global_modifications = Vec::new();
    while index < chars.len() && chars[index] == b'<' {
        let end_index =
            end_of_enclosure_with_brackets(line, index + 1, b'<', b'>').ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Global modification not closed",
                    "A global modification should be closed with a closing angle bracket '>'",
                    Context::line(None, line, index, 1),
                )
            })?;
        if let Some(offset) = next_char(chars, index, b'@') {
            let at_index = index + 1 + offset;
            if at_index > end_index {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid global modification",
                    "A global modification should have an at '@' sign inside the enclosing angle brackets '<>'",
                    Context::line(None, line, index + 1, at_index - index - 2),
                ));
            }
            if chars[index + 1] != b'[' || chars[at_index - 2] != b']' {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid global modification",
                    "A global modification should always be enclosed in square brackets '[]'",
                    Context::line(None, line, index + 1, at_index - index - 2),
                ));
            }
            let modification = SimpleModificationInner::parse_pro_forma(
                line,
                index + 2..at_index - 2,
                &mut Vec::new(),
                &mut Vec::new(),
                custom_database,
            )
            .map(|m| {
                m.0.defined().ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid global modification",
                        "A global modification cannot be ambiguous or a cross-linker",
                        Context::line(None, line, index + 2, at_index - index - 4),
                    )
                })
            })
            .flat_err()?;
            let rules = parse_placement_rules(line, at_index..end_index)?;
            global_modifications.extend(
                rules
                    .into_iter()
                    .map(|r| GlobalModification::Fixed(r, modification.clone())),
            );
        } else if &line[index + 1..end_index].to_ascii_lowercase() == "d" {
            global_modifications.push(GlobalModification::Isotope(Element::H, NonZeroU16::new(2)));
        } else {
            let num = &line[index + 1..end_index]
                .chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            let el = &line[index + 1 + num.len()..end_index];
            let el: Element = el.try_into().map_err(|()| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid global modification",
                    "Could not determine the element",
                    Context::line(
                        None,
                        line,
                        index + num.len(),
                        end_index - (index + 1 + num.len()),
                    ),
                )
            })?;
            let num = Some(num.parse::<NonZeroU16>().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid global modification",
                    format!("The isotope number is {}", explain_number_error(&err)),
                    Context::line(None, line, index + 1, end_index - index),
                )
            })?);
            if !el.is_valid(num) {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid global modification",
                    format!(
                        "This element {el} does not have a defined weight {}",
                        num.map_or_else(String::new, |num| format!("for isotope {num}"))
                    ),
                    Context::line(None, line, index + 1, end_index - index),
                ));
            }
            global_modifications.push(GlobalModification::Isotope(el, num));
        }

        index = end_index + 1;
    }
    Ok((index, global_modifications))
}

/// Parse a set of placement rules.
/// # Errors
/// When any rule is invalid.
pub(super) fn parse_placement_rules(
    line: &str,
    range: std::ops::Range<usize>,
) -> Result<Vec<PlacementRule>, BoxedError<'_, BasicKind>> {
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
                            Context::line(None, line, range.start, range.len()),
                        )
                    })?],
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
                            Context::line(None, line, range.start, range.len()),
                        )
                    })?],
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
                        Context::line(None, line, range.start, range.len()),
                    )
                })?],
                Position::Anywhere,
            ));
        }
    }
    Ok(result)
}

/// If the text is recognised as an unknown mods list it is Some(...), if it has errors during parsing Some(Err(...))
/// The returned happy path contains the mods and the index from where to continue parsing.
/// # Errors
/// Give all errors when the text cannot be read as mods of unknown position.
pub(super) fn global_unknown_position_mods<'a>(
    bytes: &'a [u8],
    start: usize,
    line: &'a str,
    custom_database: Option<&CustomDatabase>,
    ambiguous_lookup: &mut AmbiguousLookup,
) -> Option<Result<(usize, Vec<usize>), Vec<BoxedError<'a, BasicKind>>>> {
    let mut index = start;
    let mut modifications = Vec::new();
    let mut errs = Vec::new();
    let mut cross_link_lookup = Vec::new();

    // Parse until no new modifications are found
    while bytes.get(index) == Some(&b'[') {
        let start_index = index;
        index = next_char(bytes, index + 1, b']')? + 1;
        let id = match SimpleModificationInner::parse_pro_forma(
            line,
            start_index + 1..index - 1,
            ambiguous_lookup,
            &mut cross_link_lookup,
            custom_database,
        ) {
            Ok((ReturnModification::Defined(m), settings)) => {
                let id = ambiguous_lookup.len();
                ambiguous_lookup.push(AmbiguousLookupEntry::new(format!("u{id}"), Some(m)));
                ambiguous_lookup[id].copy_settings(&settings);
                id
            }
            Ok((ReturnModification::Ambiguous(id, _, _), settings)) => {
                ambiguous_lookup[id].copy_settings(&settings);
                id
            }
            Ok((ReturnModification::CrossLinkReferenced(_), _)) => {
                errs.push(BoxedError::new(
                    BasicKind::Error,
                    "Invalid unknown position modification",
                    "A modification of unknown position cannot be a cross-link",
                    Context::line_range(None, line, (start_index + 1)..index),
                ));
                continue;
            }
            Err(e) => {
                errs.push(e);
                continue;
            }
        };
        let number = if bytes.get(index) == Some(&b'^') {
            if let Some((len, num)) = next_num(bytes, index + 1, false) {
                index += len + 1;
                if num < 0 {
                    errs.push(
                        BoxedError::new(BasicKind::Error,"Invalid unknown position modification", "A modification of unknown position with multiple copies cannot have more a negative number of copies", Context::line(None, line, index, 1)));
                    0
                } else if num > i16::MAX as isize {
                    errs.push(
                        BoxedError::new(BasicKind::Error,"Invalid unknown position modification", format!("A modification of unknown position with multiple copies cannot have more then {} copies", i16::MAX), Context::line(None, line, index, 1)));
                    0
                } else {
                    num as usize
                }
            } else {
                errs.push(
                    BoxedError::new(BasicKind::Error,"Invalid unknown position modification", "A modification of unknown position with multiple copies needs the copy number after the caret ('^') symbol", Context::line(None, line, index, 1)));
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
    if bytes.get(index) == Some(&b'?') {
        Some(if errs.is_empty() {
            Ok((index + 1, modifications))
        } else {
            Err(errs)
        })
    } else {
        ambiguous_lookup.clear(); // Any ambiguous N terminal modification was incorrectly already added to the lookup
        None
    }
}

/// Parse labile modifications `{mod}{mod2}`. These are assumed to fall off from the peptide in the MS.
/// # Errors
/// If the mods are not followed by a closing brace. Or if the mods are ambiguous.
fn labile_modifications<'a>(
    line: &'a str,
    mut index: usize,
    custom_database: Option<&CustomDatabase>,
) -> Result<(usize, Vec<SimpleModification>), BoxedError<'a, BasicKind>> {
    let chars = line.as_bytes();
    let mut labile = Vec::new();
    while chars.get(index) == Some(&b'{') {
        let end_index = end_of_enclosure(line, index + 1, b'{', b'}').ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid labile modification",
                "No valid closing delimiter, a labile modification should be closed by '}'",
                Context::line(None, line, index, 1),
            )
        })?;

        labile.push(
            SimpleModificationInner::parse_pro_forma(
                line,
                index + 1..end_index,
                &mut Vec::new(),
                &mut Vec::new(),
                custom_database,
            )
            .and_then(|m| {
                m.0.defined().ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid labile modification",
                        "A labile modification cannot be ambiguous or a cross-linker",
                        Context::line(None, line, index + 1, end_index - 1 - index),
                    )
                })
            })?,
        );
        index = end_index + 1;
    }
    Ok((index, labile))
}

/// Parse a charge state `/2` or more complex ones like `/2[+2Na+]`.
/// Assumes the text starts with `/`.
/// # Errors
/// If the charge state is not following the specification.
/// # Panics
/// Panics if the text is not UTF-8.
pub(super) fn parse_charge_state(
    line: &str,
    index: usize,
) -> Result<(usize, MolecularCharge), BoxedError<'_, BasicKind>> {
    let chars = line.as_bytes();
    let (charge_len, total_charge) = next_num(chars, index + 1, false).ok_or_else(|| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid peptide charge state",
            "There should be a number dictating the total charge of the peptide",
            Context::line(None, line, index + 1, 1),
        )
    })?;
    if chars.get(index + 1 + charge_len) == Some(&b'[') {
        let end_index =
            end_of_enclosure(line, index + 2 + charge_len, b'[', b']').ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid adduct ion",
                    "No valid closing delimiter",
                    Context::line(None, line, index + 2 + charge_len, 1),
                )
            })?;
        let mut offset = index + 2 + charge_len;
        let mut charge_carriers = Vec::new();
        let mut found_charge: isize = 0;

        for set in chars[index + 2 + charge_len..end_index].split(|c| *c == b',') {
            // num
            let (count_len, count) = next_num(chars, offset, true).ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid adduct ion",
                    "Invalid adduct ion count",
                    Context::line(None, line, offset, 1),
                )
            })?;

            // charge
            let charge_len = set.iter().rev().take_while(|c| c.is_ascii_digit()).count();
            let charge = if charge_len == 0 {
                1
            } else {
                line[offset + set.len() - charge_len..offset + set.len()]
                    .parse::<i32>()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid adduct ion",
                            format!("The adduct ion number {err}"),
                            Context::line(None, line, offset + set.len() - charge_len, charge_len),
                        )
                    })?
            };
            let (charge_len, charge) = match (set.len() - charge_len)
                .checked_sub(1)
                .and_then(|i| set.get(i))
            {
                Some(b'+') => (charge_len + 1, charge),
                Some(b'-') => (charge_len + 1, -charge),
                _ => {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid adduct ion",
                        "The adduct ion number should be preceded by a sign",
                        Context::line(None, line, offset + set.len() - charge_len - 1, 1),
                    ));
                }
            };

            // Check for empty formula
            if count_len + charge_len == set.len() {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid adduct ion",
                    "The adduct ion should have a formula defined",
                    Context::line(None, line, offset, set.len()),
                ));
            }

            // formula
            let mut formula = MolecularFormula::from_pro_forma(
                line,
                offset + count_len..offset + set.len() - charge_len,
                true,
                false,
                true,
                true,
            )?;
            let _ = formula.add((
                Element::Electron,
                None,
                formula.charge().value as i32 - charge,
            ));

            // Deduplicate
            if let Some((amount, _)) = charge_carriers.iter_mut().find(|(_, f)| *f == formula) {
                *amount += count;
            } else {
                charge_carriers.push((count, formula));
            }

            offset += set.len() + 1;
            found_charge = found_charge
                .checked_add(count.checked_mul(charge as isize).ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid peptide charge state",
                        "The peptide charge state is too big to store inside an isize",
                        Context::line(None, line, index, offset),
                    )
                })?)
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid peptide charge state",
                        "The peptide charge state is too big to store inside an isize",
                        Context::line(None, line, index, offset),
                    )
                })?;
        }
        if total_charge == found_charge {
            Ok((end_index + 1, MolecularCharge::new(&charge_carriers)))
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid peptide charge state",
                "The peptide charge state number has to be equal to the sum of all separate adduct ions",
                Context::line(None, line, index, offset),
            ))
        }
    } else {
        // If no adduct ions are provided assume it is just protons
        Ok((
            index + charge_len + 1,
            MolecularCharge::proton(total_charge),
        ))
    }
}
