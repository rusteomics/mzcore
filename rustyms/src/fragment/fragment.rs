//! Handle fragment related issues, access provided if you want to dive deeply into fragments in your own code.

use std::{
    fmt::{Debug, Display, Write},
    sync::LazyLock,
};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::{ChargeRange, PossiblePrimaryIons},
    chemistry::{
        AmbiguousLabel, CachedCharge, Chemical, MassMode, MolecularCharge, MolecularFormula,
        MultiChemical,
    },
    fragment::{BackboneCFragment, BackboneNFragment, FragmentType, NeutralLoss},
    molecular_formula,
    quantities::{Multi, Tolerance},
    sequence::{BACKBONE, IsAminoAcid, SequencePosition},
    system::{
        OrderedMassOverCharge,
        f64::{MassOverCharge, Ratio},
        usize::Charge,
    },
};

/// A theoretical fragment
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Fragment {
    /// The theoretical composition
    pub formula: Option<MolecularFormula>,
    /// The charge
    pub charge: Charge,
    /// The annotation for this fragment
    pub ion: FragmentType,
    /// The peptidoform this fragment comes from, saved as the index into the list of peptidoform in the overarching [`crate::sequence::CompoundPeptidoformIon`] struct
    pub peptidoform_ion_index: Option<usize>,
    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching [`crate::sequence::PeptidoformIon`] struct
    pub peptidoform_index: Option<usize>,
    /// Any neutral losses applied
    pub neutral_loss: Vec<NeutralLoss>,
    /// m/z deviation, if known (from mzPAF)
    pub deviation: Option<Tolerance<OrderedMassOverCharge>>,
    /// Confidence in this annotation (from mzPAF)
    pub confidence: Option<OrderedFloat<f64>>,
    /// If this is an auxiliary fragment (from mzPAF)
    pub auxiliary: bool,
}

impl Fragment {
    /// Write the fragment as an mzPAF string
    #[allow(non_snake_case)]
    pub fn to_mzPAF(&self) -> String {
        let mut output = String::new();
        if self.auxiliary {
            output.push('&');
        }
        // Push the ion type info (plus maybe some neutral losses if needed)
        match &self.ion {
            FragmentType::a(pos, variant)
            | FragmentType::b(pos, variant)
            | FragmentType::c(pos, variant)
            | FragmentType::x(pos, variant)
            | FragmentType::y(pos, variant) => write!(
                &mut output,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 0 {
                    String::new()
                } else {
                    format!("{variant:+}H")
                }
            )
            .unwrap(),
            FragmentType::z(pos, variant) => write!(
                &mut output,
                "{}{}{}",
                self.ion.kind(),
                pos.series_number,
                if *variant == 1 {
                    String::new()
                } else {
                    format!("{:+}H", variant - 1)
                }
            )
            .unwrap(),
            FragmentType::d(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "d{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else if let Some(loss) = aa
                    .satellite_ion_fragments(
                        pos.sequence_index,
                        self.peptidoform_index.unwrap_or_default(),
                        self.peptidoform_ion_index.unwrap_or_default(),
                    )
                    .and_then(|fragments| {
                        fragments
                            .iter()
                            .find(|f| f.0 == *label)
                            .map(|(_, loss)| loss.clone())
                    })
                {
                    write!(
                        &mut output,
                        "a{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::v(pos, aa, distance, variant) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "v{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(
                        &mut output,
                        "y{}-{}{}",
                        pos.series_number,
                        aa.formulas()
                            .first()
                            .map(|f| f - LazyLock::force(&BACKBONE))
                            .unwrap_or_default(),
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                }
            }
            FragmentType::w(pos, aa, distance, variant, label) => {
                if *distance == 0 {
                    write!(
                        &mut output,
                        "w{label}{}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else if let Some(loss) = aa
                    .satellite_ion_fragments(
                        pos.sequence_index,
                        self.peptidoform_index.unwrap_or_default(),
                        self.peptidoform_ion_index.unwrap_or_default(),
                    )
                    .and_then(|fragments| {
                        fragments
                            .iter()
                            .find(|f| f.0 == *label)
                            .map(|(_, loss)| loss.clone())
                    })
                {
                    write!(
                        &mut output,
                        "z{}-{loss}{}",
                        pos.series_number,
                        if *variant == 0 {
                            String::new()
                        } else {
                            format!("{variant:+}H")
                        }
                    )
                    .unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::Precursor => write!(&mut output, "p").unwrap(),
            FragmentType::PrecursorSideChainLoss(_, aa) => {
                write!(&mut output, "p-r[sidechain_{aa}]").unwrap();
            }
            FragmentType::Immonium(_, seq) => write!(
                &mut output,
                "I{}{}",
                seq.aminoacid,
                seq.modifications.iter().map(|m| format!("[{m}]")).join("") // TODO: how to handle ambiguous mods? maybe store somewhere which where applied for this fragment
            )
            .unwrap(),
            FragmentType::Unknown(num) => write!(
                &mut output,
                "?{}",
                num.map_or(String::new(), |u| u.to_string())
            )
            .unwrap(),
            FragmentType::Diagnostic(_)
            | FragmentType::B { .. }
            | FragmentType::BComposition(_, _)
            | FragmentType::Y(_)
            | FragmentType::YComposition(_, _) => {
                if let Some(formula) = &self.formula {
                    // TODO: better way of storing?
                    write!(&mut output, "f{{{formula}}}",).unwrap();
                } else {
                    write!(&mut output, "?",).unwrap();
                }
            }
            FragmentType::Internal(Some(name), a, b) => write!(
                &mut output,
                "m{}:{}{}",
                a.sequence_index + 1,
                b.sequence_index + 1,
                match name {
                    (BackboneNFragment::a, BackboneCFragment::x)
                    | (BackboneNFragment::b, BackboneCFragment::y)
                    | (BackboneNFragment::c, BackboneCFragment::z) => "",
                    (BackboneNFragment::a, BackboneCFragment::y) => "-CO",
                    (BackboneNFragment::a, BackboneCFragment::z) => "-CHNO",
                    (BackboneNFragment::b, BackboneCFragment::x) => "+CO",
                    (BackboneNFragment::b, BackboneCFragment::z) => "-NH",
                    (BackboneNFragment::c, BackboneCFragment::x) => "+CHNO",
                    (BackboneNFragment::c, BackboneCFragment::y) => "+NH",
                }
            )
            .unwrap(),
            FragmentType::Internal(None, a, b) => write!(
                &mut output,
                "m{}:{}",
                a.sequence_index + 1,
                b.sequence_index + 1
            )
            .unwrap(),
        }
        // More losses
        for loss in &self.neutral_loss {
            match loss {
                NeutralLoss::SideChainLoss(_, aa) => {
                    write!(&mut output, "-r[sidechain_{aa}]").unwrap();
                }
                l => write!(&mut output, "{l}").unwrap(),
            }
        }
        // Isotopes: not handled
        // Charge state
        if self.charge.value != 1 {
            write!(&mut output, "^{}", self.charge.value).unwrap();
        }
        // Deviation
        match self.deviation {
            Some(Tolerance::Absolute(abs)) => write!(&mut output, "/{}", abs.value).unwrap(),
            Some(Tolerance::Relative(ppm)) => write!(&mut output, "/{}ppm", ppm.value).unwrap(),
            None => (),
        }
        // Confidence
        if let Some(confidence) = self.confidence {
            write!(&mut output, "*{confidence}").unwrap();
        }
        output
    }

    /// Get the mz
    pub fn mz(&self, mode: MassMode) -> Option<MassOverCharge> {
        self.formula.as_ref().map(|f| {
            f.mass(mode)
                / crate::system::f64::Charge::new::<crate::system::charge::e>(
                    self.charge.value as f64,
                )
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
            neutral_loss: Vec::new(),
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
                charge: Charge::new::<crate::system::e>(charge.charge().value.try_into().unwrap()),
                ion: annotation.clone(),
                peptidoform_ion_index: Some(peptidoform_ion_index),
                peptidoform_index: Some(peptidoform_index),
                neutral_loss: losses.cloned().unwrap_or_default(),
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
    pub fn generate_series(
        theoretical_mass: &Multi<MolecularFormula>,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        annotation: &FragmentType,
        termini: &Multi<MolecularFormula>,
        charge_carriers: &mut CachedCharge,
        settings: &PossiblePrimaryIons,
    ) -> Vec<Self> {
        termini
            .iter()
            .cartesian_product(theoretical_mass.iter())
            .cartesian_product(charge_carriers.range(settings.1))
            .cartesian_product(std::iter::once(None).chain(settings.0.iter().map(Some)))
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
                charge: Charge::new::<crate::system::e>(charge.charge().value.try_into().unwrap()),
                ion: annotation.with_variant(*variant),
                peptidoform_ion_index: Some(peptidoform_ion_index),
                peptidoform_index: Some(peptidoform_index),
                neutral_loss: losses.cloned().unwrap_or_default(),
                deviation: None,
                confidence: None,
                auxiliary: false,
            })
            .collect()
    }

    /// Create a copy of this fragment with the given charge
    /// # Panics
    /// If the charge is negative.
    #[must_use]
    fn with_charge(&self, charge: &MolecularCharge) -> Self {
        let formula = charge
            .formula()
            .with_labels(&[AmbiguousLabel::ChargeCarrier(charge.formula())]);
        let c = Charge::new::<crate::system::charge::e>(
            usize::try_from(formula.charge().value).unwrap(),
        );
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
    use crate::{fragment::PeptidePosition, prelude::AminoAcid};

    use super::*;

    #[test]
    fn neutral_loss() {
        let a = Fragment::new(
            AminoAcid::AsparticAcid.formulas()[0].clone(),
            Charge::new::<crate::system::charge::e>(1),
            0,
            0,
            FragmentType::Precursor,
        );
        let loss = a.with_neutral_losses(&[NeutralLoss::Loss(molecular_formula!(H 2 O 1))]);
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
