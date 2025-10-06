#![warn(dead_code)]

use std::{
    collections::{HashMap, HashSet},
    fmt::Write,
    marker::PhantomData,
    num::NonZeroU32,
};

use context_error::*;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, MultiChemical},
    glycan::{BackboneFragmentKind, GlycanAttachement},
    quantities::Multi,
    sequence::{
        AtLeast, CheckedAminoAcid, CrossLinkName, Linked, LinkerSpecificity, Modification,
        Peptidoform, PlacementRule, RulePossible, SequencePosition, SimpleModification,
        SimpleModificationInner,
    },
};

/// One block in a sequence meaning an aminoacid and it's accompanying modifications
#[derive(Debug, Default, Deserialize, Ord, PartialOrd, Serialize)]
pub struct SequenceElement<T> {
    /// The aminoacid
    pub aminoacid: CheckedAminoAcid<T>,
    /// All present modifications
    pub modifications: ThinVec<Modification>,
    /// If this aminoacid is part of an ambiguous sequence group `(QA)?` in ProForma
    pub ambiguous: Option<NonZeroU32>,
    /// The marker indicating which level of complexity this sequence element uses as higher bound
    marker: PhantomData<T>,
}

impl<T> Clone for SequenceElement<T> {
    fn clone(&self) -> Self {
        Self {
            aminoacid: self.aminoacid,
            modifications: self.modifications.clone(),
            ambiguous: self.ambiguous,
            marker: PhantomData,
        }
    }
}

impl<A, B> PartialEq<SequenceElement<B>> for SequenceElement<A> {
    fn eq(&self, other: &SequenceElement<B>) -> bool {
        self.aminoacid == other.aminoacid
            && self.modifications == other.modifications
            && self.ambiguous == other.ambiguous
    }
}

impl<T> std::hash::Hash for SequenceElement<T> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.aminoacid.hash(state);
        self.modifications.hash(state);
        self.ambiguous.hash(state);
    }
}

impl<T> Eq for SequenceElement<T> {}

impl<T> SequenceElement<T> {
    /// Mark this sequence element as the following complexity level, the level is not validated
    pub(super) fn mark<M>(self) -> SequenceElement<M> {
        SequenceElement {
            aminoacid: self.aminoacid.mark::<M>(),
            modifications: self.modifications,
            ambiguous: self.ambiguous,
            marker: PhantomData,
        }
    }

    /// Cast a sequence element into a more complex sequence element. This does not change the
    /// content of the sequence element. It only allows to pass this as higher complexity if needed.
    pub fn cast<OtherComplexity: AtLeast<T>>(self) -> SequenceElement<OtherComplexity> {
        self.mark()
    }

    /// Create a new aminoacid without any modifications
    pub fn new(aminoacid: CheckedAminoAcid<T>, ambiguous: Option<NonZeroU32>) -> Self {
        Self {
            aminoacid,
            modifications: ThinVec::new(),
            ambiguous,
            marker: PhantomData,
        }
    }

    /// Add a modification to this sequence element
    #[must_use]
    pub fn with_simple_modification(mut self, modification: SimpleModification) -> Self {
        self.modifications.push(Modification::Simple(modification));
        self
    }

    /// Add a modification to this sequence element
    pub fn add_simple_modification(&mut self, modification: SimpleModification) {
        self.modifications.push(Modification::Simple(modification));
    }
}

impl<T> SequenceElement<T> {
    /// # Errors
    /// If the underlying formatter errors.
    pub(crate) fn display(
        &self,
        f: &mut impl Write,
        placed_ambiguous: &[usize],
        preferred_ambiguous_location: &[Option<SequencePosition>],
        index: usize,
        last_ambiguous: Option<NonZeroU32>,
        specification_compliant: bool,
    ) -> Result<Vec<usize>, std::fmt::Error> {
        let mut extra_placed = Vec::new();
        if last_ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, ")")?;
        }
        if self.ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, "(?")?;
        }
        write!(f, "{}", self.aminoacid)?;
        for m in &self.modifications {
            let mut display_ambiguous = false;
            if let Modification::Ambiguous { id, .. } = m
                && (!placed_ambiguous.contains(id) && preferred_ambiguous_location[*id].is_none()
                    || preferred_ambiguous_location[*id]
                        .is_some_and(|p| p == SequencePosition::Index(index)))
            {
                display_ambiguous = true;
                extra_placed.push(*id);
            }
            write!(f, "[")?;
            m.display(f, specification_compliant, display_ambiguous)?;
            write!(f, "]")?;
        }
        Ok(extra_placed)
    }

    /// Get the molecular formulas for this position without any the ambiguous modifications
    pub(crate) fn formulas_base(
        &self,
        all_peptidoforms: &[Peptidoform<Linked>],
        visited_peptidoforms: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
        glycan_model: &impl GlycanAttachement,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        self.formulas_generic(
            &mut |_| false,
            all_peptidoforms,
            visited_peptidoforms,
            applied_cross_links,
            allow_ms_cleavable,
            sequence_index,
            peptidoform_index,
            peptidoform_ion_index,
            glycan_model,
        )
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotope modifications
    #[expect(clippy::too_many_arguments)]
    pub(crate) fn formulas_greedy(
        &self,
        placed: &mut [bool],
        all_peptidoforms: &[Peptidoform<Linked>],
        visited_peptidoforms: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
        glycan_model: &impl GlycanAttachement,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        self.formulas_generic(
            &mut |id| (!placed[id]).then(|| placed[id] = true).is_some(),
            all_peptidoforms,
            visited_peptidoforms,
            applied_cross_links,
            allow_ms_cleavable,
            sequence_index,
            peptidoform_index,
            peptidoform_ion_index,
            glycan_model,
        )
    }

    /// Get the molecular formulas for this position with all ambiguous modifications, without any global isotope modifications
    pub(crate) fn formulas_all(
        &self,
        all_peptidoforms: &[Peptidoform<Linked>],
        visited_peptidoforms: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
        glycan_model: &impl GlycanAttachement,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        self.formulas_generic(
            &mut |_| true,
            all_peptidoforms,
            visited_peptidoforms,
            applied_cross_links,
            allow_ms_cleavable,
            sequence_index,
            peptidoform_index,
            peptidoform_ion_index,
            glycan_model,
        )
    }

    /// Get the molecular formulas for this position with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotope modifications
    #[expect(clippy::too_many_arguments)]
    pub(crate) fn formulas_generic(
        &self,
        place_ambiguous: &mut impl FnMut(usize) -> bool,
        all_peptidoforms: &[Peptidoform<Linked>],
        visited_peptidoforms: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
        glycan_model: &impl GlycanAttachement,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        let (formula, specific, seen) = self
            .modifications
            .iter()
            .filter_map(|m| {
                if let Modification::Ambiguous {
                    id, modification, ..
                } = m
                {
                    place_ambiguous(*id).then(|| {
                        let default = glycan_model.get_default_fragments();
                        let specific =
                            glycan_model.get_specific_fragments(Some(self.aminoacid.aminoacid()));
                        (
                            modification.formula_inner(
                                sequence_index,
                                peptidoform_index,
                                default,
                                Some(self.aminoacid.aminoacid()),
                            ),
                            specific
                                .into_iter()
                                .map(|(k, setting)| {
                                    (
                                        k,
                                        modification.formula_inner(
                                            sequence_index,
                                            peptidoform_index,
                                            setting,
                                            Some(self.aminoacid.aminoacid()),
                                        ),
                                    )
                                })
                                .collect(),
                            HashSet::default(),
                        )
                    })
                } else {
                    let (formula, specific, seen) = m.formula_inner(
                        all_peptidoforms,
                        visited_peptidoforms,
                        applied_cross_links,
                        allow_ms_cleavable,
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                        glycan_model,
                        Some(self.aminoacid.aminoacid()),
                    );
                    Some((formula, specific, seen))
                }
            })
            .fold(
                (Multi::default(), HashMap::new(), HashSet::new()),
                |(am, asp, av), (m, sp, v)| {
                    (
                        am * m,
                        crate::helper_functions::merge_hashmap(asp, sp),
                        av.union(&v).cloned().collect(),
                    )
                },
            );
        let own =
            self.aminoacid
                .formulas_inner(sequence_index, peptidoform_index, peptidoform_ion_index);
        (
            formula * &own,
            specific.into_iter().map(|(k, v)| (k, v * &own)).collect(),
            seen,
        )
    }

    /// Enforce the placement rules of predefined modifications.
    /// # Errors
    /// If a rule has been broken.
    /// # Panics
    /// If any placement rule is placement on a PSI modification that does not exist.
    pub(crate) fn enforce_modification_rules(
        &self,
        position: SequencePosition,
    ) -> Result<(), BoxedError<'static, BasicKind>> {
        for modification in &self.modifications {
            if modification.is_possible(self, position) == RulePossible::No {
                let rules = modification
                    .simple()
                    .map(|s| s.placement_rules())
                    .unwrap_or_default();
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Modification incorrectly placed",
                    format!(
                        "Modification {modification} is not allowed on {}{}",
                        match position {
                            SequencePosition::NTerm => "the N-terminus".to_string(),
                            SequencePosition::CTerm => "the C-terminus".to_string(),
                            SequencePosition::Index(index) =>
                                format!("the side chain of {} at index {index}", self.aminoacid),
                        },
                        if rules.is_empty() {
                            String::new()
                        } else {
                            format!(
                                ", this modification is only allowed at the following locations: {}",
                                rules.join(", ")
                            )
                        }
                    ),
                    Context::none(),
                ));
            }
        }
        Ok(())
    }

    /// Get all possible diagnostic ions
    pub(crate) fn diagnostic_ions(
        &self,
        position: SequencePosition,
        n_term: &[Modification],
        c_term: &[Modification],
    ) -> Vec<DiagnosticIon> {
        let mut diagnostic_ions = Vec::new();
        let modifications = match position {
            SequencePosition::NTerm => n_term,
            SequencePosition::Index(_) => &self.modifications,
            SequencePosition::CTerm => c_term,
        };
        for modification in modifications {
            match modification {
                Modification::CrossLink { linker, side, .. } => {
                    diagnostic_ions.extend_from_slice(&side.allowed_rules(linker).2);
                }
                Modification::Simple(modification)
                | Modification::Ambiguous { modification, .. } => match &**modification {
                    SimpleModificationInner::Database { specificities, .. } => {
                        for (rules, _, ions) in specificities {
                            if PlacementRule::any_possible(rules, self, position) {
                                diagnostic_ions.extend_from_slice(ions);
                            }
                        }
                    }
                    SimpleModificationInner::Linker { specificities, .. } => {
                        for rule in specificities {
                            match rule {
                                LinkerSpecificity::Symmetric {
                                    rules, diagnostic, ..
                                } => {
                                    if PlacementRule::any_possible(rules, self, position) {
                                        diagnostic_ions.extend_from_slice(diagnostic);
                                    }
                                }
                                LinkerSpecificity::Asymmetric {
                                    rules: (rules_left, rules_right),

                                    diagnostic,
                                    ..
                                } => {
                                    if PlacementRule::any_possible(rules_left, self, position)
                                        || PlacementRule::any_possible(rules_right, self, position)
                                    {
                                        diagnostic_ions.extend_from_slice(diagnostic);
                                    }
                                }
                            }
                        }
                    }
                    _ => (),
                },
            }
        }
        diagnostic_ions
    }
}

impl<C, T> From<T> for SequenceElement<C>
where
    T: Into<CheckedAminoAcid<C>>,
{
    fn from(value: T) -> Self {
        Self::new(value.into(), None)
    }
}
