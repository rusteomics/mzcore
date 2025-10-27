//! Handle fragment related issues, access provided if you want to dive deeply into fragments in your own code.

use std::fmt::{Debug, Display};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use mzcore::{
    chemistry::{AmbiguousLabel, CachedCharge, ChargeRange, NeutralLoss},
    molecular_formula,
    prelude::*,
    quantities::{Multi, Tolerance},
    system::{self, MassOverCharge, OrderedMassOverCharge, Ratio, isize::Charge},
};
use thin_vec::ThinVec;

use crate::{annotation::model::PossiblePrimaryIons, fragment::FragmentType};

/// A theoretical fragment
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Fragment {
    /// The theoretical composition
    pub formula: Option<MolecularFormula>,
    /// The charge
    pub charge: Charge,
    /// The annotation for this fragment
    pub ion: FragmentType,
    /// The peptidoform this fragment comes from, saved as the index into the list of peptidoform in the overarching [`mzcore::sequence::CompoundPeptidoformIon`] struct
    pub peptidoform_ion_index: Option<usize>,
    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching [`mzcore::sequence::PeptidoformIon`] struct
    pub peptidoform_index: Option<usize>,
    /// Any neutral losses applied
    pub neutral_loss: ThinVec<NeutralLoss>,
    /// m/z deviation, if known (from mzPAF)
    pub deviation: Option<Tolerance<OrderedMassOverCharge>>,
    /// Confidence in this annotation (from mzPAF)
    pub confidence: Option<OrderedFloat<f64>>,
    /// If this is an auxiliary fragment (from mzPAF)
    pub auxiliary: bool,
}

impl Fragment {
    /// Get the mz
    pub fn mz(&self, mode: MassMode) -> Option<MassOverCharge> {
        self.formula.as_ref().map(|f| {
            f.mass(mode) / system::f64::Charge::new::<system::charge::e>(self.charge.value as f64)
        })
    }

    /// Get the ppm difference between two fragments
    pub fn ppm(&self, other: &Self, mode: MassMode) -> Option<Ratio> {
        self.mz(mode)
            .and_then(|mz| other.mz(mode).map(|omz| (mz, omz)))
            .map(|(mz, omz)| mz.ppm(omz))
    }

    /// Create a new fragment
    #[must_use]
    pub fn new(
        theoretical_mass: MolecularFormula,
        charge: Charge,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        ion: FragmentType,
    ) -> Self {
        Self {
            formula: Some(theoretical_mass),
            charge,
            ion,
            peptidoform_ion_index: Some(peptidoform_ion_index),
            peptidoform_index: Some(peptidoform_index),
            neutral_loss: ThinVec::new(),
            deviation: None,
            confidence: None,
            auxiliary: false,
        }
    }

    /// Generate a list of possible fragments from the list of possible preceding termini and neutral losses
    /// # Panics
    /// When the charge range results in a negative charge
    #[expect(clippy::too_many_arguments)]
    #[must_use]
    pub fn generate_all(
        theoretical_mass: &Multi<MolecularFormula>,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        annotation: &FragmentType,
        termini: &Multi<MolecularFormula>,
        neutral_losses: &[Vec<NeutralLoss>],
        charge_carriers: &mut CachedCharge,
        charge_range: ChargeRange,
    ) -> Vec<Self> {
        termini
            .iter()
            .cartesian_product(theoretical_mass.iter())
            .cartesian_product(charge_carriers.range(charge_range))
            .cartesian_product(std::iter::once(None).chain(neutral_losses.iter().map(Some)))
            .map(|(((term, mass), charge), losses)| Self {
                formula: Some(
                    term + mass
                        + charge.formula_inner(SequencePosition::default(), peptidoform_index)
                        + losses
                            .iter()
                            .flat_map(|l| l.iter())
                            .sum::<MolecularFormula>(),
                ),
                charge: Charge::new::<system::e>(charge.charge().value),
                ion: annotation.clone(),
                peptidoform_ion_index: Some(peptidoform_ion_index),
                peptidoform_index: Some(peptidoform_index),
                neutral_loss: losses.cloned().unwrap_or_default().into(),
                deviation: None,
                confidence: None,
                auxiliary: false,
            })
            .collect()
    }

    /// Generate a list of possible fragments from the list of possible preceding termini and neutral losses
    /// # Panics
    /// When the charge range results in a negative charge
    #[must_use]
    #[expect(clippy::too_many_arguments)] // Needs many different pieces of information
    pub fn generate_series(
        theoretical_mass: &Multi<MolecularFormula>,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        annotation: &FragmentType,
        termini: &Multi<MolecularFormula>,
        neutral_losses: &[Vec<NeutralLoss>],
        charge_carriers: &mut CachedCharge,
        settings: &PossiblePrimaryIons,
    ) -> Vec<Self> {
        termini
            .iter()
            .cartesian_product(theoretical_mass.iter())
            .cartesian_product(charge_carriers.range(settings.1))
            .cartesian_product(
                std::iter::once(None)
                    .chain(settings.0.iter().map(Some))
                    .chain(neutral_losses.iter().map(Some)),
            )
            .cartesian_product(settings.2.iter())
            .map(|((((term, mass), charge), losses), variant)| Self {
                formula: Some(
                    term + mass
                        + charge.formula_inner(SequencePosition::default(), peptidoform_index)
                        + losses
                            .iter()
                            .flat_map(|l| l.iter())
                            .sum::<MolecularFormula>()
                        + molecular_formula!(H 1) * variant,
                ),
                charge: Charge::new::<system::e>(charge.charge().value),
                ion: annotation.with_variant(*variant),
                peptidoform_ion_index: Some(peptidoform_ion_index),
                peptidoform_index: Some(peptidoform_index),
                neutral_loss: losses.cloned().unwrap_or_default().into(),
                deviation: None,
                confidence: None,
                auxiliary: false,
            })
            .collect()
    }

    /// Create a copy of this fragment with the given charge
    #[must_use]
    fn with_charge(&self, charge: &MolecularCharge) -> Self {
        let formula = charge
            .formula()
            .with_labels(&[AmbiguousLabel::ChargeCarrier(charge.formula())]);
        let c = Charge::new::<system::charge::e>(formula.charge().value);
        Self {
            formula: Some(self.formula.clone().unwrap_or_default() + &formula),
            charge: c,
            ..self.clone()
        }
    }

    /// Create a copy of this fragment with the given charges
    pub fn with_charge_range(
        self,
        charge_carriers: &mut CachedCharge,
        charge_range: ChargeRange,
    ) -> impl Iterator<Item = Self> {
        charge_carriers
            .range(charge_range)
            .into_iter()
            .map(move |c| self.with_charge(&c))
    }

    /// Create a copy of this fragment with the given charges
    pub fn with_charge_range_slice(
        self,
        charges: &[MolecularCharge],
    ) -> impl Iterator<Item = Self> {
        charges.iter().map(move |c| self.with_charge(c))
    }

    /// Create a copy of this fragment with the given neutral loss
    #[must_use]
    pub fn with_neutral_loss(&self, neutral_loss: &NeutralLoss) -> Self {
        let mut new_neutral_loss = self.neutral_loss.clone();
        new_neutral_loss.push(neutral_loss.clone());
        Self {
            formula: Some(self.formula.clone().unwrap_or_default() + neutral_loss),
            neutral_loss: new_neutral_loss,
            ..self.clone()
        }
    }

    /// Create copies of this fragment with the given neutral losses (and a copy of this fragment itself)
    #[must_use]
    pub fn with_neutral_losses(&self, neutral_losses: &[NeutralLoss]) -> Vec<Self> {
        let mut output = Vec::with_capacity(neutral_losses.len() + 1);
        output.push(self.clone());
        output.extend(
            neutral_losses
                .iter()
                .map(|loss| self.with_neutral_loss(loss)),
        );
        output
    }
}

impl Display for Fragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}@{}{:+}{}",
            self.ion,
            self.mz(MassMode::Monoisotopic)
                .map_or(String::new(), |mz| mz.value.to_string()),
            self.charge.value,
            self.neutral_loss.iter().map(ToString::to_string).join("")
        )
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use mzcore::sequence::PeptidePosition;

    use super::*;

    #[test]
    fn neutral_loss() {
        let a = Fragment::new(
            AminoAcid::AsparticAcid.formulas()[0].clone(),
            Charge::new::<system::charge::e>(1),
            0,
            0,
            FragmentType::Precursor,
        );
        let loss = a.with_neutral_losses(&[NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]);
        dbg!(&a, &loss);
        assert_eq!(a.formula, loss[0].formula);
        assert_eq!(
            a.formula.unwrap(),
            &loss[1].formula.clone().unwrap() + &molecular_formula!(H 2 O 1)
        );
    }

    #[test]
    fn flip_terminal() {
        let n0 = PeptidePosition::n(SequencePosition::Index(0), 2);
        let n1 = PeptidePosition::n(SequencePosition::Index(1), 2);
        let n2 = PeptidePosition::n(SequencePosition::Index(2), 2);
        let c0 = PeptidePosition::c(SequencePosition::Index(0), 2);
        let c1 = PeptidePosition::c(SequencePosition::Index(1), 2);
        let c2 = PeptidePosition::c(SequencePosition::Index(2), 2);
        assert_eq!(n0.flip_terminal(), c0);
        assert_eq!(n1.flip_terminal(), c1);
        assert_eq!(n2.flip_terminal(), c2);
    }
}
