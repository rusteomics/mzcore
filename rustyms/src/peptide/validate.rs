use itertools::Itertools;
use std::collections::BTreeMap;

use crate::{
    error::{Context, CustomError},
    modification::{CrossLinkName, SimpleModification},
    LinearPeptide, Modification, Peptidoform, SequencePosition,
};

use super::{GlobalModification, Linear};

/// Validate all cross links
/// # Errors
/// If there is a cross link with more then 2 locations. Or if there never is a definition for this cross link.
/// Or if there are peptides that cannot be reached from the first peptide.
pub fn cross_links(
    peptides: Vec<LinearPeptide<Linear>>,
    cross_links_found: BTreeMap<usize, Vec<(usize, SequencePosition)>>,
    cross_link_lookup: &[(CrossLinkName, Option<SimpleModification>)],
    line: &str,
) -> Result<Peptidoform, CustomError> {
    let mut peptidoform = Peptidoform(peptides.into_iter().map(Into::into).collect());
    for (id, locations) in cross_links_found {
        let definition = &cross_link_lookup[id];
        if let Some(linker) = &definition.1 {
            match locations.len() {
                0 => {return Err(CustomError::error(
                    "Invalid cross-link",
                    format!("The cross-link named '{}' has no listed locations, this is an internal error please report this", definition.0),
                    Context::full_line(0, line),
                ))},
                1 => {
                    let (index, position) = locations[0];
                    if linker.is_possible(&peptidoform.0[index][position], position).any_possible() {
                        peptidoform.0[index].add_simple_modification(position, linker.clone());
                    } else {
                        let rules = linker.placement_rules();
                        return Err(CustomError::error(
                            "Modification incorrectly placed",
                            format!(
                                "Modification {linker} is not allowed on {}{}",
                                match position {
                                    SequencePosition::NTerm => "the N-terminus".to_string(),
                                    SequencePosition::CTerm => "the C-terminus".to_string(),
                                    SequencePosition::Index(index) =>
                                        format!("the side chain of {} at index {index}", peptidoform.0[index][position].aminoacid),
                                },
                                if rules.is_empty() {
                                    String::new()
                                } else {
                                    format!(", this modification is only allowed at the following locations: {}", rules.join(", "))
                                }
                            ),
                            Context::full_line(0, line),
                        ));
                    }
                },
                2 => {
                    if !peptidoform.add_cross_link(locations[0], locations[1], linker.clone(), definition.0.clone()) {
                        return Err(CustomError::error(
                            "Invalid cross-link",
                            format!("The cross-link named '{}' cannot be placed according to its location specificities", definition.0),
                            Context::full_line(0, line),
                        ))
                    }
                },
                _ => {return Err(CustomError::error(
                    "Invalid cross-link",
                    format!("The cross-link named '{}' has more than 2 attachment locations, only cross-links spanning two locations are allowed", definition.0),
                    Context::full_line(0, line),
                ))}
            }
        } else {
            let (c, name, description) = if definition.0 == CrossLinkName::Branch {
                ("MOD", "00134", " N6-glycyl-L-lysine")
            } else {
                ("X", "DSS", "")
            };
            return Err(CustomError::error(
                "Invalid cross-link",
                format!("The cross-link named '{0}' is never defined, for example for {name}{description} define it like: '[{c}:{name}{0}]'", definition.0),
                Context::full_line(0, line),
            ));
        }
    }

    // Check if all peptides can be reached from the first one
    let mut found_peptides = Vec::new();
    let mut stack = vec![0];

    while let Some(index) = stack.pop() {
        found_peptides.push(index);
        if let Some(Modification::CrossLink { peptide, .. }) = &peptidoform.0[index].get_n_term() {
            if !found_peptides.contains(peptide) && !stack.contains(peptide) {
                stack.push(*peptide);
            }
        }
        if let Some(Modification::CrossLink { peptide, .. }) = &peptidoform.0[index].get_c_term() {
            if !found_peptides.contains(peptide) && !stack.contains(peptide) {
                stack.push(*peptide);
            }
        }
        for seq in peptidoform.0[index].sequence() {
            for modification in &seq.modifications {
                if let Modification::CrossLink { peptide, .. } = modification {
                    if !found_peptides.contains(peptide) && !stack.contains(peptide) {
                        stack.push(*peptide);
                    }
                }
            }
        }
    }

    if found_peptides.len() != peptidoform.peptides().len() {
        return Err(CustomError::error(
            "Unconnected peptidoform",
            "Not all peptides in this peptidoform are connected with cross-links or branches, if separate peptides were intended use the chimeric notation `+` instead of the peptidoform notation `//`.",
             Context::full_line(0, line),
        ));
    }

    Ok(peptidoform)
}

impl LinearPeptide<Linear> {
    /// Apply a global modification if this is a global isotope modification with invalid isotopes it returns false
    #[must_use]
    pub(super) fn apply_global_modifications(
        &mut self,
        global_modifications: &[GlobalModification],
    ) -> bool {
        for modification in global_modifications {
            match modification {
                GlobalModification::Fixed(pos, aa, modification) => {
                    let positions = self
                        .iter(..)
                        .filter(|(position, seq)| {
                            pos.is_possible(position.sequence_index)
                                && aa.map_or(true, |aa| aa == seq.aminoacid.aminoacid())
                                && modification
                                    .is_possible(seq, position.sequence_index)
                                    .any_possible()
                        })
                        .map(|(position, _)| position)
                        .collect_vec();
                    for position in positions {
                        self.add_simple_modification(position.sequence_index, modification.clone());
                    }
                }
                GlobalModification::Isotope(el, isotope) if el.is_valid(*isotope) => {
                    let _ = self.add_global((*el, *isotope)); // Already validated
                }
                GlobalModification::Isotope(..) => return false,
            }
        }
        true
    }

    /// Place all global unknown positions at all possible locations as ambiguous modifications
    /// # Errors
    /// When a mod cannot be placed anywhere
    pub(super) fn apply_unknown_position_modification(
        &mut self,
        unknown_position_modifications: &[SimpleModification],
    ) -> Result<(), CustomError> {
        for modification in unknown_position_modifications {
            if !self.add_unknown_position_modification(modification.clone(), None, ..) {
                return Err(CustomError::error(
                    "Modification of unknown position cannot be placed", 
                    "There is no position where this modification can be placed based on the placement rules in the database.",
                     Context::show(modification)
                    ));
            }
        }
        Ok(())
    }

    /// Place all ranged unknown positions at all possible locations as ambiguous modifications
    /// # Errors
    /// When a mod cannot be placed anywhere
    /// # Panics
    /// It panics when information for an ambiguous modification is missing (name/mod).
    pub(super) fn apply_ranged_unknown_position_modification(
        &mut self,
        ranged_unknown_position_modifications: &[(usize, usize, SimpleModification)],
    ) -> Result<(), CustomError> {
        for (start, end, modification) in ranged_unknown_position_modifications {
            if !self.add_unknown_position_modification(modification.clone(), None, start..=end) {
                return Err(CustomError::error(
                    "Modification of unknown position on a range cannot be placed", 
                    "There is no position where this modification can be placed based on the placement rules in the database.", 
                    Context::show(modification)
                ));
            }
        }
        Ok(())
    }
}

impl<T> LinearPeptide<T> {
    /// # Errors
    /// If a modification rule is broken it returns an error.
    pub(crate) fn enforce_modification_rules(&self) -> Result<(), CustomError> {
        for (position, seq) in self.iter(..) {
            seq.enforce_modification_rules(position.sequence_index)?;
        }
        Ok(())
    }
}
