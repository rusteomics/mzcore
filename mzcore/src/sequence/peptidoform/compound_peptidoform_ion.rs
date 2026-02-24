use std::fmt::{Display, Write};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    chemistry::MolecularFormula,
    prelude::Chemical,
    quantities::Multi,
    sequence::{Linked, Peptidoform, PeptidoformIon},
};

/// A single full ProForma entry. This entry can contain multiple sets of cross-linked peptides.
/// A single set of cross-linked peptides is a [`PeptidoformIon`]. A ProForma entry with two chimeric
/// peptides will be saved as one [`CompoundPeptidoformIon`] with two [`PeptidoformIon`]s that each
/// contain one of the [`Peptidoform`]s.
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct CompoundPeptidoformIon {
    name: String,
    pub(super) peptidoforms: Vec<PeptidoformIon>,
}

impl CompoundPeptidoformIon {
    /// Create a new [`CompoundPeptidoformIon`] from many [`Peptidoform`]s. This returns None if the
    /// global isotope modifications of all peptidoforms are not identical.
    pub fn new(name: String, iter: impl IntoIterator<Item = PeptidoformIon>) -> Option<Self> {
        let result = Self {
            name,
            peptidoforms: iter.into_iter().collect(),
        };
        let global_equal = result
            .peptidoform_ions()
            .iter()
            .flat_map(PeptidoformIon::peptidoforms)
            .tuple_windows()
            .all(|(a, b)| a.get_global() == b.get_global());
        global_equal.then_some(result)
    }

    /// Get all possible formulas for this compound peptidoform
    pub fn formulas(&self) -> Multi<MolecularFormula> {
        self.peptidoforms
            .iter()
            .enumerate()
            .flat_map(|(i, p)| {
                (p.formulas_inner(i)
                    + p.get_charge_carriers()
                        .map(Chemical::formula)
                        .unwrap_or_default())
                .to_vec()
            })
            .collect()
    }

    /// Get all possible neutral formulas for this compound peptidoform
    pub fn neutral_formulas(&self) -> Multi<MolecularFormula> {
        self.peptidoforms
            .iter()
            .enumerate()
            .flat_map(|(i, p)| p.formulas_inner(i).to_vec())
            .collect()
    }

    /// Assume there is exactly one peptidoform in this compound peptidoform.
    #[doc(alias = "assume_linear")]
    pub fn singular(mut self) -> Option<PeptidoformIon> {
        (self.peptidoforms.len() == 1)
            .then(|| self.peptidoforms.pop())
            .flatten()
    }

    /// Assume there is exactly one peptidoform in this compound peptidoform.
    pub fn singular_ref(&self) -> Option<&PeptidoformIon> {
        (self.peptidoforms.len() == 1).then(|| &self.peptidoforms[0])
    }

    /// Assume there is exactly one peptide in this compound peptidoform.
    pub fn singular_peptidoform(self) -> Option<Peptidoform<Linked>> {
        self.singular().and_then(PeptidoformIon::singular)
    }

    /// Assume there is exactly one peptide in this compound peptidoform.
    pub fn singular_peptidoform_ref(&self) -> Option<&Peptidoform<Linked>> {
        self.singular_ref().and_then(PeptidoformIon::singular_ref)
    }

    /// Get all peptidoform ions making up this compound peptidoform.
    pub fn peptidoform_ions(&self) -> &[PeptidoformIon] {
        &self.peptidoforms
    }

    /// Get all peptidoform ions making up this compound peptidoform.
    pub fn into_peptidoform_ions(self) -> Vec<PeptidoformIon> {
        self.peptidoforms
    }

    /// Get all peptidoforms making up this compound peptidoform.
    pub fn peptidoforms(&self) -> impl Iterator<Item = &Peptidoform<Linked>> {
        self.peptidoforms
            .iter()
            .flat_map(PeptidoformIon::peptidoforms)
    }

    /// Display this compound peptidoform.
    /// `specification_compliant` Displays this compound peptidoform either normalised to the
    /// internal representation (with false) or as fully spec compliant ProForma (no glycan
    /// structure or custom modifications) (with true).
    /// # Errors
    /// Only if the underlying formatter (`f`) errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        // The global isotope modifications are guaranteed to be identical, so take the first
        let empty = Vec::new();
        if !self.name.is_empty() {
            write!(f, "(>>>{})", self.name)?;
        }
        let global = self
            .peptidoform_ions()
            .iter()
            .flat_map(PeptidoformIon::peptidoforms)
            .next()
            .map_or(empty.as_slice(), |p| p.get_global());
        for (element, isotope) in global {
            write!(
                f,
                "<{}{}>",
                isotope.map(|i| i.to_string()).unwrap_or_default(),
                element
            )?;
        }

        let mut first = true;
        for p in self.peptidoform_ions() {
            if !first {
                write!(f, "+")?;
            }
            p.display(f, false, specification_compliant)?;
            first = false;
        }
        Ok(())
    }

    /// Get the name
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Get the name
    pub fn name_mut(&mut self) -> &mut String {
        &mut self.name
    }
}

impl Display for CompoundPeptidoformIon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl<Complexity> From<Peptidoform<Complexity>> for CompoundPeptidoformIon {
    fn from(value: Peptidoform<Complexity>) -> Self {
        Self {
            name: String::new(),
            peptidoforms: vec![value.into()],
        }
    }
}

impl From<PeptidoformIon> for CompoundPeptidoformIon {
    fn from(value: PeptidoformIon) -> Self {
        Self {
            name: String::new(),
            peptidoforms: vec![value],
        }
    }
}

impl<Complexity> From<Vec<Peptidoform<Complexity>>> for CompoundPeptidoformIon {
    fn from(value: Vec<Peptidoform<Complexity>>) -> Self {
        Self {
            name: String::new(),
            peptidoforms: value.into_iter().map(Into::into).collect(),
        }
    }
}

impl<Complexity> From<ThinVec<Peptidoform<Complexity>>> for CompoundPeptidoformIon {
    fn from(value: ThinVec<Peptidoform<Complexity>>) -> Self {
        Self {
            name: String::new(),
            peptidoforms: value.into_iter().map(Into::into).collect(),
        }
    }
}

impl FromIterator<PeptidoformIon> for CompoundPeptidoformIon {
    fn from_iter<T: IntoIterator<Item = PeptidoformIon>>(iter: T) -> Self {
        Self {
            name: String::new(),
            peptidoforms: iter.into_iter().collect(),
        }
    }
}

impl<Complexity> FromIterator<Peptidoform<Complexity>> for CompoundPeptidoformIon {
    fn from_iter<T: IntoIterator<Item = Peptidoform<Complexity>>>(iter: T) -> Self {
        Self {
            name: String::new(),
            peptidoforms: iter.into_iter().map(|p| p.into()).collect(),
        }
    }
}

impl crate::space::Space for CompoundPeptidoformIon {
    fn space(&self) -> crate::space::UsedSpace {
        self.name.space() + self.peptidoforms.space()
    }
}
