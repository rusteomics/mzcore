use std::{collections::BTreeSet, fmt::Write};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{MolecularCharge, MolecularFormula},
    glycan::FullGlycan,
    quantities::Multi,
    sequence::{
        CrossLinkName, CrossLinkSide, Linked, Modification, Peptidoform, RulePossible,
        SequencePosition, SimpleModification, SimpleModificationInner,
    },
};

/// A single peptidoform ion, can contain multiple peptidoforms
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct PeptidoformIon {
    name: String,
    pub(crate) peptidoforms: Vec<Peptidoform<Linked>>,
}

impl PeptidoformIon {
    /// Create a new [`PeptidoformIon`] from many [`Peptidoform`]s. This returns None if the
    /// global isotope modifications or the charge carriers of all peptides are not identical.
    pub fn new<Complexity>(
        name: String,
        iter: impl IntoIterator<Item = Peptidoform<Complexity>>,
    ) -> Option<Self> {
        let result = Self {
            name,
            peptidoforms: iter.into_iter().map(Peptidoform::mark).collect(),
        };
        let global_and_charge_equal = result.peptidoforms().iter().tuple_windows().all(|(a, b)| {
            a.get_global() == b.get_global() && a.get_charge_carriers() == b.get_charge_carriers()
        });
        global_and_charge_equal.then_some(result)
    }

    /// Create a new [`PeptidoformIon`] from many [`Peptidoform`]s. This returns None if the
    /// global isotope modifications or the charge carriers of all peptides are not identical.
    pub fn from_vec(mut iter: Vec<Peptidoform<Linked>>) -> Option<Self> {
        iter.shrink_to_fit();
        let result = Self {
            name: String::new(),
            peptidoforms: iter,
        };
        let global_and_charge_equal = result.peptidoforms().iter().tuple_windows().all(|(a, b)| {
            a.get_global() == b.get_global() && a.get_charge_carriers() == b.get_charge_carriers()
        });
        global_and_charge_equal.then_some(result)
    }

    /// Shrink to fit on all peptidoforms
    pub fn shrink_to_fit(&mut self) {
        self.name.shrink_to_fit();
        self.peptidoforms.shrink_to_fit();
        for p in &mut self.peptidoforms {
            p.shrink_to_fit();
        }
    }

    /// Gives all possible formulas for this peptidoform (including breakage of cross-links that can break).
    /// Includes the full glycan, if there are any glycans.
    /// Assumes all peptides in this peptidoform are connected.
    /// If there are no peptides in this peptidoform it returns [`Multi::default`].
    pub fn formulas(&self) -> Multi<MolecularFormula> {
        self.formulas_inner(0)
    }

    /// Gives all possible formulas for this peptidoform (including breakage of cross-links that can break).
    /// Includes the full glycan, if there are any glycans.
    /// Assumes all peptides in this peptidoform are connected.
    /// If there are no peptides in this peptidoform it returns [`Multi::default`].
    pub(crate) fn formulas_inner(&self, peptidoform_ion_index: usize) -> Multi<MolecularFormula> {
        self.peptidoforms
            .first()
            .map(|p| {
                p.formulas_inner(
                    0,
                    peptidoform_ion_index,
                    &self.peptidoforms,
                    &[],
                    &mut Vec::new(),
                    true,
                    &FullGlycan {},
                )
                .0
            })
            .unwrap_or_default()
    }

    /// Assume there is exactly one peptide in this collection.
    pub fn singular(mut self) -> Option<Peptidoform<Linked>> {
        (self.peptidoforms.len() == 1)
            .then(|| self.peptidoforms.pop())
            .flatten()
    }

    /// Assume there is exactly one peptide in this collection.
    pub fn singular_ref(&self) -> Option<&Peptidoform<Linked>> {
        (self.peptidoforms.len() == 1).then(|| &self.peptidoforms[0])
    }

    /// Get all peptides making up this peptidoform
    pub fn peptidoforms(&self) -> &[Peptidoform<Linked>] {
        &self.peptidoforms
    }

    /// Get all peptides making up this peptidoform
    pub fn peptidoforms_mut(&mut self) -> &mut [Peptidoform<Linked>] {
        &mut self.peptidoforms
    }

    /// Set the charge carriers
    #[expect(clippy::needless_pass_by_value)]
    pub fn set_charge_carriers(&mut self, charge_carriers: Option<MolecularCharge>) {
        for peptide in &mut self.peptidoforms {
            peptide.set_charge_carriers(charge_carriers.clone());
        }
    }

    /// Get the charge carriers
    pub fn get_charge_carriers(&self) -> Option<&MolecularCharge> {
        // Take the charge of the first peptidoform, as all are required to have the same charge
        self.peptidoforms
            .first()
            .and_then(|p| p.get_charge_carriers())
    }

    /// Get the name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the name
    pub const fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }

    /// Add a cross-link to this peptidoform and check if it is placed according to its placement rules.
    /// The positions are first the peptide index and second the sequence index.
    pub fn add_cross_link(
        &mut self,
        position_1: (usize, SequencePosition),
        position_2: (usize, SequencePosition),
        linker: SimpleModification,
        name: CrossLinkName,
    ) -> bool {
        let pos_1 = self
            .peptidoforms
            .get(position_1.0)
            .map(|seq| &seq[position_1.1]);
        let pos_2 = self
            .peptidoforms
            .get(position_2.0)
            .map(|seq| &seq[position_2.1]);
        if let (Some(pos_1), Some(pos_2)) = (pos_1, pos_2) {
            let left = linker.is_possible(pos_1, position_1.1);
            let right = linker.is_possible(pos_2, position_1.1);
            let (left, right, according_to_rules) = if matches!(
                &*linker,
                SimpleModificationInner::Formula(_)
                    | SimpleModificationInner::Glycan(_)
                    | SimpleModificationInner::GlycanStructure(_)
                    | SimpleModificationInner::Gno { .. }
                    | SimpleModificationInner::Mass(_, _, _)
            ) {
                (
                    CrossLinkSide::Symmetric(BTreeSet::default()),
                    CrossLinkSide::Symmetric(BTreeSet::default()),
                    true,
                )
            } else {
                match (left, right) {
                    (RulePossible::Symmetric(a), RulePossible::Symmetric(b)) => {
                        let intersection: BTreeSet<usize> = a.intersection(&b).copied().collect();
                        if intersection.is_empty() {
                            (
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                false,
                            )
                        } else {
                            (
                                CrossLinkSide::Symmetric(intersection.clone()),
                                CrossLinkSide::Symmetric(intersection),
                                true,
                            )
                        }
                    }
                    (
                        RulePossible::AsymmetricLeft(a),
                        RulePossible::AsymmetricRight(b) | RulePossible::Symmetric(b),
                    ) => {
                        let intersection: BTreeSet<usize> = a.intersection(&b).copied().collect();
                        if intersection.is_empty() {
                            (
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                false,
                            )
                        } else {
                            (
                                CrossLinkSide::Left(intersection.clone()),
                                CrossLinkSide::Right(intersection),
                                true,
                            )
                        }
                    }
                    (
                        RulePossible::AsymmetricRight(a),
                        RulePossible::AsymmetricLeft(b) | RulePossible::Symmetric(b),
                    ) => {
                        let intersection: BTreeSet<usize> = a.intersection(&b).copied().collect();
                        if intersection.is_empty() {
                            (
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                CrossLinkSide::Symmetric(BTreeSet::default()),
                                false,
                            )
                        } else {
                            (
                                CrossLinkSide::Right(intersection.clone()),
                                CrossLinkSide::Left(intersection),
                                true,
                            )
                        }
                    }
                    _ => (
                        CrossLinkSide::Symmetric(BTreeSet::default()),
                        CrossLinkSide::Symmetric(BTreeSet::default()),
                        false,
                    ),
                }
            };
            self.peptidoforms[position_1.0].add_modification(
                position_1.1,
                Modification::CrossLink {
                    peptide: position_2.0,
                    sequence_index: position_2.1,
                    linker: linker.clone(),
                    name: name.clone(),
                    side: left,
                },
            );
            self.peptidoforms[position_2.0].add_modification(
                position_2.1,
                Modification::CrossLink {
                    peptide: position_1.0,
                    sequence_index: position_1.1,
                    linker,
                    name,
                    side: right,
                },
            );
            according_to_rules
        } else {
            false // TODO: maybe generate better error on invalid positions
        }
    }

    /// Display this peptidoform.
    /// `specification_compliant` Displays this peptidoform either normalised to the internal representation or as fully spec compliant ProForma
    /// (no glycan structure or custom modifications).
    /// # Panics
    /// When some peptides do not have the same global isotope modifications.
    /// # Errors
    /// If the underlying formatter errors.
    pub fn display(
        &self,
        f: &mut impl Write,
        show_global_mods: bool,
        specification_compliant: bool,
    ) -> std::fmt::Result {
        if show_global_mods {
            let global_equal = self
                .peptidoforms()
                .iter()
                .map(Peptidoform::get_global)
                .tuple_windows()
                .all(|(a, b)| a == b);
            assert!(
                global_equal,
                "Not all global isotope modifications on all peptides on this peptidoform are identical"
            );
            let empty = Vec::new();
            let global = self
                .peptidoforms()
                .first()
                .map_or(empty.as_slice(), |p| p.get_global());
            for (element, isotope) in global {
                write!(
                    f,
                    "<{}{}>",
                    isotope.map(|i| i.to_string()).unwrap_or_default(),
                    element
                )?;
            }
        }
        if !self.name.is_empty() {
            write!(f, "(>>{})", self.name)?;
        }

        let mut first = true;
        for (index, p) in self.peptidoforms().iter().enumerate() {
            if !first {
                write!(f, "//")?;
            }
            p.display(
                f,
                false,
                index == self.peptidoforms().len() - 1,
                specification_compliant,
            )?;
            first = false;
        }
        Ok(())
    }
}

impl std::fmt::Display for PeptidoformIon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true, true)
    }
}

impl<Complexity> From<Peptidoform<Complexity>> for PeptidoformIon {
    fn from(value: Peptidoform<Complexity>) -> Self {
        Self {
            name: String::new(),
            peptidoforms: vec![value.mark()],
        }
    }
}

impl crate::space::Space for PeptidoformIon {
    fn space(&self) -> crate::space::UsedSpace {
        self.name.space() + self.peptidoforms.space()
    }
}
