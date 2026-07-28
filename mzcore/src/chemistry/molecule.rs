use crate::{
    chemistry::MolecularFormula, quantities::Multi, sequence::SequencePosition, system::Mass,
};

/// A trait to generalise over different ways of calculating total masses.
pub trait MassOutputMode {
    /// The output type
    type Output: std::iter::Sum + Clone;

    fn from_mass(mass: Mass) -> Self::Output;
    fn from_formula(formula: MolecularFormula) -> Self::Output;
}

/// A trait to signify that only the final monoisotopic mass is interesting, note that this will be
/// less precise than going through a [`MolecularFormula`] first.
pub struct OutputMonoIsotopic;

impl MassOutputMode for OutputMonoIsotopic {
    type Output = Mass;

    fn from_mass(mass: Mass) -> Self::Output {
        mass
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula.monoisotopic_mass()
    }
}

/// A trait to signify that only the final average weight is interesting, note that this will be
/// less precise than going through a [`MolecularFormula`] first.
pub struct OutputAverageWeight;

impl MassOutputMode for OutputAverageWeight {
    type Output = Mass;

    fn from_mass(mass: Mass) -> Self::Output {
        mass
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula.average_weight()
    }
}

pub struct OutputMolecularFormula;

impl MassOutputMode for OutputMolecularFormula {
    type Output = MolecularFormula;

    fn from_mass(mass: Mass) -> Self::Output {
        MolecularFormula::with_additional_mass(mass.value)
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula
    }
}

/// Any molecule that has a clearly defined single molecular formula
pub trait Molecule<Mode: MassOutputMode> {
    /// Get the molecular formula
    fn formula(&self) -> Mode::Output {
        self.formula_inner(SequencePosition::default(), 0)
    }

    /// Get the molecular formula while keeping track of all ambiguous labels
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output;
}

impl<T: Molecule<Mode>, Mode: MassOutputMode> Molecule<Mode> for &[T] {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

impl<T: Molecule<Mode>, Mode: MassOutputMode> Molecule<Mode> for &Vec<T> {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

/// Any molecule that has a number of potential chemical formulas
pub trait AmbiguousMolecule<Mode: MassOutputMode> {
    /// Get all possible molecular formulas
    fn formulas(&self) -> Multi<Mode::Output> {
        self.formulas_inner(SequencePosition::default(), 0, 0)
    }

    /// Get all possible molecular formulas while keeping track of all ambiguous labels
    fn formulas_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Multi<Mode::Output>;

    /// Get the charge of this chemical, it returns None if no charge is defined.
    fn charge(&self) -> Option<crate::system::isize::Charge>;

    /// Return a single formula if this `MultiChemical` has only one possible formula
    fn single_formula(&self) -> Option<Mode::Output> {
        let formulas = self.formulas();
        (formulas.len() == 1).then_some(formulas[0].clone())
    }

    /// Return a single formula if this `MultiChemical` has only one possible formula while keeping
    /// track of all ambiguous labels
    fn single_formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Option<Mode::Output> {
        let formulas =
            self.formulas_inner(sequence_index, peptidoform_index, peptidoform_ion_index);
        (formulas.len() == 1).then_some(formulas[0].clone())
    }
}
