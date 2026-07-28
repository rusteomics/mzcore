use std::num::NonZeroU16;

use crate::{
    chemistry::{AmbiguousLabel, Element, MassMode, MolecularFormula},
    quantities::Multi,
    sequence::SequencePosition,
    system::{Mass, OrderedMass},
};

/// A trait to generalise over different ways of calculating total masses.
pub trait MassOutputMode: Sized {
    /// The output type
    type Output: MassOutputType;

    /// Create this output from a simple mass
    fn from_mass(mass: Mass) -> Self::Output;

    /// Create this output from a molecular formula, for example from the values as defined in the
    /// modification databases.
    fn from_ref_formula(formula: &MolecularFormula) -> Self::Output;

    /// Create this output from a molecular formula, for example from the values as defined in the
    /// source code.
    fn from_formula(formula: MolecularFormula) -> Self::Output;

    /// Pick one of the three supplied values, to prevent having to allocate for the formula if that
    /// is not necessary.
    fn pick(
        monoisotopic: MassOutput,
        average_weight: MassOutput,
        formula: impl FnOnce() -> MolecularFormula,
    ) -> Self::Output;
}

/// A type that can be used as mass output in calculations.
pub trait MassOutputType
where
    Self: std::iter::Sum<Self>
        + Clone
        + Default
        + std::ops::Add<Self, Output = Self>
        + std::ops::Sub<Self, Output = Self>
        + std::ops::Mul<i32, Output = Self>
        + Eq
        + std::hash::Hash,
{
    /// Add the following labels to this formula
    #[must_use]
    fn with_labels(self, labels: &[AmbiguousLabel]) -> Self;

    /// Add the following label to this formula
    #[must_use]
    fn with_label(self, label: AmbiguousLabel) -> Self;

    /// The labels of sources of ambiguity/multiplicity
    fn labels(&self) -> &[AmbiguousLabel];

    /// Apply the given global isotope modifications. It can return None if any of the substitutents
    /// are invalid.
    #[must_use]
    fn with_global_isotope_modifications(
        self,
        isotopes: &[(Element, Option<NonZeroU16>)],
    ) -> Option<Self>;
}

/// A trait to signify that only the final monoisotopic mass is interesting, note that this will be
/// less precise than going through a [`MolecularFormula`] first. And this will miss any global
/// isotope modifications.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct OutputMonoIsotopic;

impl MassOutputMode for OutputMonoIsotopic {
    type Output = MassOutput;

    fn from_mass(mass: Mass) -> Self::Output {
        mass.into()
    }

    fn from_ref_formula(formula: &MolecularFormula) -> Self::Output {
        formula.monoisotopic_mass().into()
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula.monoisotopic_mass().into()
    }

    fn pick(
        monoisotopic: MassOutput,
        _average_weight: MassOutput,
        _formula: impl FnOnce() -> MolecularFormula,
    ) -> Self::Output {
        monoisotopic
    }
}

/// A trait to signify that only the final average weight is interesting, note that this will be
/// less precise than going through a [`MolecularFormula`] first. And this will miss any global
/// isotope modifications.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct OutputAverageWeight;

impl MassOutputMode for OutputAverageWeight {
    type Output = MassOutput;

    fn from_mass(mass: Mass) -> Self::Output {
        mass.into()
    }

    fn from_ref_formula(formula: &MolecularFormula) -> Self::Output {
        formula.average_weight().into()
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula.average_weight().into()
    }

    fn pick(
        _monoisotopic: MassOutput,
        average_weight: MassOutput,
        _formula: impl FnOnce() -> MolecularFormula,
    ) -> Self::Output {
        average_weight
    }
}

/// A mass output, contains the mass, charge, and any ambiguous labels.
#[derive(Debug, Clone, Default)]
pub struct MassOutput {
    pub(crate) mass: Mass,
    pub(crate) charge: crate::system::isize::Charge,
    pub(crate) labels: thin_vec::ThinVec<AmbiguousLabel>,
}

impl MassOutput {
    /// Get the mass
    pub fn mass(&self) -> Mass {
        self.mass
    }

    /// Get the charge
    pub fn charge(&self) -> crate::system::isize::Charge {
        self.charge
    }

    /// Get the mass over charge
    pub fn mz(&self) -> crate::system::MassOverCharge {
        crate::system::MassOverCharge::new::<crate::system::thomson>(
            self.mass.value / self.charge.value as f64,
        )
    }
}

impl From<Mass> for MassOutput {
    fn from(value: Mass) -> Self {
        Self {
            mass: value,
            ..Default::default()
        }
    }
}

impl From<MassOutput> for Mass {
    fn from(value: MassOutput) -> Self {
        value.mass
    }
}

impl From<MassOutput> for OrderedMass {
    fn from(value: MassOutput) -> Self {
        value.mass.into()
    }
}

impl From<&MassOutput> for Mass {
    fn from(value: &MassOutput) -> Self {
        value.mass
    }
}

impl From<&MassOutput> for OrderedMass {
    fn from(value: &MassOutput) -> Self {
        value.mass.into()
    }
}

impl PartialEq for MassOutput {
    fn eq(&self, other: &Self) -> bool {
        self.mass.value.total_cmp(&other.mass.value).is_eq()
            && self.charge.value == other.charge.value
    }
}

impl Eq for MassOutput {}

impl Ord for MassOutput {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.mass
            .value
            .total_cmp(&other.mass.value)
            .then(self.charge.value.cmp(&other.charge.value))
    }
}

impl PartialOrd for MassOutput {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl std::hash::Hash for MassOutput {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        // Implementation from ordered float crate
        if self.mass.value.is_nan() {
            0x7ff8000000000000_u64.hash(state);
        } else {
            (self.mass.value + 0.0).to_bits().hash(state);
        }
        self.charge.value.hash(state);
    }
}

impl std::ops::Mul<i32> for MassOutput {
    type Output = Self;

    fn mul(self, rhs: i32) -> Self::Output {
        Self {
            mass: self.mass * rhs as f64,
            charge: self.charge * rhs as isize,
            labels: self.labels,
        }
    }
}

impl std::ops::Add<MassOutput> for MassOutput {
    type Output = Self;

    fn add(self, rhs: MassOutput) -> Self::Output {
        let mut labels = self.labels;
        labels.extend_from_slice(&rhs.labels);
        Self {
            mass: self.mass + rhs.mass,
            charge: self.charge + rhs.charge,
            labels,
        }
    }
}

impl std::ops::Sub<MassOutput> for MassOutput {
    type Output = Self;

    fn sub(self, rhs: MassOutput) -> Self::Output {
        Self {
            mass: self.mass - rhs.mass,
            charge: self.charge - rhs.charge,
            labels: self.labels,
        }
    }
}

impl std::iter::Sum<MassOutput> for MassOutput {
    fn sum<I: Iterator<Item = MassOutput>>(iter: I) -> Self {
        iter.fold(Self::default(), |acc, i| acc + i)
    }
}

impl MassOutputType for MassOutput {
    fn labels(&self) -> &[AmbiguousLabel] {
        &self.labels
    }

    fn with_label(mut self, label: AmbiguousLabel) -> Self {
        self.labels.push(label);
        self
    }

    fn with_labels(mut self, labels: &[AmbiguousLabel]) -> Self {
        self.labels.extend_from_slice(labels);
        self
    }

    fn with_global_isotope_modifications(
        self,
        substitutions: &[(Element, Option<NonZeroU16>)],
    ) -> Option<Self> {
        if substitutions.is_empty() || substitutions.iter().all(|(_, i)| i.is_none()) {
            Some(self)
        } else {
            None // No elements are stored so this cannot be applied
        }
    }
}

/// Calculate masses as molecular formula. This allows for maximally precise values, most abundant
/// isotope calculation, isotopic pattern generation and more.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct OutputMolecularFormula;

impl MassOutputMode for OutputMolecularFormula {
    type Output = MolecularFormula;

    fn from_mass(mass: Mass) -> Self::Output {
        MolecularFormula::with_additional_mass(mass.value)
    }

    fn from_ref_formula(formula: &MolecularFormula) -> Self::Output {
        formula.clone()
    }

    fn from_formula(formula: MolecularFormula) -> Self::Output {
        formula
    }

    fn pick(
        _monoisotopic: MassOutput,
        _average_weight: MassOutput,
        formula: impl FnOnce() -> MolecularFormula,
    ) -> Self::Output {
        formula()
    }
}

/// Any molecule that has a clearly defined single molecular formula
pub trait Molecule {
    /// Get the molecular formula
    fn calculate_mass<Mode: MassOutputMode>(&self) -> Mode::Output {
        self.calculate_mass_inner::<Mode>(SequencePosition::default(), 0)
    }

    /// Get the molecular formula while keeping track of all ambiguous labels
    fn calculate_mass_inner<Mode: MassOutputMode>(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output;

    /// Get the mass of this molecule, this takes a shortcut when calculating monoisotopic mass or
    /// average weight to be faster but less precise.
    fn mass(&self, mass_mode: MassMode) -> MassOutput {
        match mass_mode {
            MassMode::Monoisotopic => self.calculate_mass::<OutputMonoIsotopic>(),
            MassMode::Average => self.calculate_mass::<OutputAverageWeight>(),
            #[cfg(feature = "isotopes")]
            MassMode::MostAbundant => self.monoisotopic_mass(),
        }
    }

    /// Get the molecular formula for this molecule
    fn formula(&self) -> MolecularFormula {
        self.calculate_mass::<OutputMolecularFormula>()
    }

    /// Get the monoisotopic mass for this molecule, less precise then going through a formula first
    /// but faster
    fn monoisotopic_mass(&self) -> MassOutput {
        self.calculate_mass::<OutputMonoIsotopic>()
    }

    /// Get the average weight for this molecule, less precise then going through a formula first
    /// but faster
    fn average_weight(&self) -> MassOutput {
        self.calculate_mass::<OutputAverageWeight>()
    }

    /// Get the most abundant mass for this molecule, less precise then going through a formula
    /// first but faster
    #[cfg(feature = "isotopes")]
    fn most_abundant_mass(&self) -> MassOutput {
        let f = self.calculate_mass::<OutputMolecularFormula>();

        MassOutput {
            mass: f.mass(MassMode::MostAbundant),
            charge: f.charge(),
            labels: f.labels,
        }
    }
}

impl<T: Molecule> Molecule for &[T] {
    fn calculate_mass_inner<Mode: MassOutputMode>(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output {
        self.iter()
            .map(|f| f.calculate_mass_inner::<Mode>(sequence_index, peptidoform_index))
            .sum()
    }
}

impl<T: Molecule> Molecule for &Vec<T> {
    fn calculate_mass_inner<Mode: MassOutputMode>(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output {
        self.iter()
            .map(|f| f.calculate_mass_inner::<Mode>(sequence_index, peptidoform_index))
            .sum()
    }
}

/// Any molecule that has a number of potential chemical formulas
pub trait AmbiguousMolecule {
    /// Get all possible molecular formulas
    fn calculate_masses<Mode: MassOutputMode>(&self) -> Multi<Mode::Output> {
        self.calculate_masses_inner::<Mode>(SequencePosition::default(), 0, 0)
    }

    /// Get all possible molecular formulas while keeping track of all ambiguous labels
    fn calculate_masses_inner<Mode: MassOutputMode>(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Multi<Mode::Output>;

    /// Get the mass of this molecule, this takes a shortcut when calculating monoisotopic mass or
    /// average weight to be faster but less precise.
    fn masses(&self, mass_mode: MassMode) -> Multi<MassOutput> {
        match mass_mode {
            MassMode::Monoisotopic => self.calculate_masses::<OutputMonoIsotopic>(),
            MassMode::Average => self.calculate_masses::<OutputAverageWeight>(),
            #[cfg(feature = "isotopes")]
            MassMode::MostAbundant => self.monoisotopic_masses(),
        }
    }

    /// Get the molecular formula for this molecule
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.calculate_masses::<OutputMolecularFormula>()
    }

    /// Get the monoisotopic mass for this molecule, less precise then going through a formula first
    /// but faster
    fn monoisotopic_masses(&self) -> Multi<MassOutput> {
        self.calculate_masses::<OutputMonoIsotopic>()
    }

    /// Get the average weight for this molecule, less precise then going through a formula first
    /// but faster
    fn average_weights(&self) -> Multi<MassOutput> {
        self.calculate_masses::<OutputAverageWeight>()
    }

    /// Get the most abundant mass for this molecule, less precise then going through a formula
    /// first but faster
    #[cfg(feature = "isotopes")]
    fn most_abundant_masses(&self) -> Multi<MassOutput> {
        self.calculate_masses::<OutputMolecularFormula>()
            .to_vec()
            .into_iter()
            .map(|f| MassOutput {
                mass: f.mass(MassMode::MostAbundant),
                charge: f.charge(),
                labels: f.labels,
            })
            .collect()
    }
}
