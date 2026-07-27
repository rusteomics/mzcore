use crate::{chemistry::MolecularFormula, quantities::Multi, sequence::SequencePosition};

/// Any item that has a clearly defined single molecular formula
pub trait Chemical {
    /// Get the molecular formula
    fn formula(&self) -> MolecularFormula {
        self.formula_inner(SequencePosition::default(), 0)
    }

    /// Get the molecular formula while keeping track of all ambiguous labels
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula;
}

impl<T: Chemical> Chemical for &[T] {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

impl<T: Chemical> Chemical for &Vec<T> {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

/// Any item that has a number of potential chemical formulas
pub trait MultiChemical {
    /// Get all possible molecular formulas
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.formulas_inner(SequencePosition::default(), 0, 0)
    }

    /// Get all possible molecular formulas while keeping track of all ambiguous labels
    fn formulas_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Multi<MolecularFormula>;

    /// Get the charge of this chemical, it returns None if no charge is defined.
    fn charge(&self) -> Option<crate::system::isize::Charge> {
        self.formulas()
            .first()
            .map(MolecularFormula::charge)
            .filter(|c| c.value != 0)
    }

    /// Return a single formula if this `MultiChemical` has only one possible formula
    fn single_formula(&self) -> Option<MolecularFormula> {
        let formulas = self.formulas();
        (formulas.len() == 1).then_some(formulas.to_vec().pop().unwrap())
    }

    /// Return a single formula if this `MultiChemical` has only one possible formula while keeping
    /// track of all ambiguous labels
    fn single_formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Option<MolecularFormula> {
        let formulas =
            self.formulas_inner(sequence_index, peptidoform_index, peptidoform_ion_index);
        (formulas.len() == 1).then_some(formulas.to_vec().pop().unwrap())
    }
}
