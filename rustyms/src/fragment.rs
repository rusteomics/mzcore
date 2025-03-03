//! Handle fragment related issues, access provided if you want to dive deeply into fragments in your own code.

use std::{
    borrow::Cow,
    fmt::{Debug, Display},
};

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    glycan::MonoSaccharide,
    model::ChargeRange,
    molecular_charge::{CachedCharge, MolecularCharge},
    system::{
        f64::{MassOverCharge, Ratio},
        usize::Charge,
        OrderedMassOverCharge,
    },
    AmbiguousLabel, AminoAcid, Chemical, MassMode, Modification, MolecularFormula, Multi,
    NeutralLoss, SemiAmbiguous, SequenceElement, SequencePosition, Tolerance,
};

/// A theoretical fragment of a peptide
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize, Default)]
pub struct Fragment {
    /// The theoretical composition
    pub formula: Option<MolecularFormula>,
    /// The charge
    pub charge: Charge,
    /// All possible annotations for this fragment saved as a tuple of peptide index and its type
    pub ion: FragmentType,
    /// The peptidoform this fragment comes from, saved as the index into the list of peptidoform in the overarching [`crate::CompoundPeptidoform`] struct
    pub peptidoform_ion_index: Option<usize>,
    /// The peptide this fragment comes from, saved as the index into the list of peptides in the overarching [`crate::Peptidoform`] struct
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
        neutral_losses: &[NeutralLoss],
        charge_carriers: &mut CachedCharge,
        charge_range: ChargeRange,
    ) -> Vec<Self> {
        termini
            .iter()
            .cartesian_product(theoretical_mass.iter())
            .cartesian_product(charge_carriers.range(charge_range))
            .cartesian_product(std::iter::once(None).chain(neutral_losses.iter().map(Some)))
            .map(|(((term, mass), charge), loss)| Self {
                formula: Some(
                    term + mass
                        + charge.formula_inner(SequencePosition::default(), peptidoform_index)
                        + loss.unwrap_or(&NeutralLoss::Gain(MolecularFormula::default())),
                ),
                charge: Charge::new::<crate::system::e>(charge.charge().value.try_into().unwrap()),
                ion: annotation.clone(),
                peptidoform_ion_index: Some(peptidoform_ion_index),
                peptidoform_index: Some(peptidoform_index),
                neutral_loss: loss.map(|l| vec![l.clone()]).unwrap_or_default(),
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
            self.neutral_loss
                .iter()
                .map(std::string::ToString::to_string)
                .join("")
        )
    }
}

// /// An isotope annotation.
// #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
// pub struct MatchedIsotopeDistribution {
//     /// The index of the matched peak in the spectrum, if found
//     pub peak_index: Option<usize>,
//     /// The isotope offset in whole daltons from the monoisotopic peak
//     pub isotope_offset: usize,
//     /// The theoretical abundance of this isotope (normalised to 1 for the whole distribution)
//     pub theoretical_isotope_abundance: OrderedFloat<f64>,
// }

/// The definition of the position of an ion
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
#[non_exhaustive]
pub struct PeptidePosition {
    /// The sequence index (0 based into the peptide sequence)
    pub sequence_index: SequencePosition,
    /// The series number (1 based from the ion series terminal)
    pub series_number: usize,
    /// The length of the whole sequence
    pub sequence_length: usize,
}

impl PeptidePosition {
    /// Generate a position for N terminal ion series
    pub const fn n(sequence_index: SequencePosition, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: match sequence_index {
                SequencePosition::NTerm => 0,
                SequencePosition::Index(i) => i + 1,
                SequencePosition::CTerm => length,
            },
            sequence_length: length,
        }
    }
    /// Generate a position for C terminal ion series
    pub const fn c(sequence_index: SequencePosition, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: match sequence_index {
                SequencePosition::NTerm => length,
                SequencePosition::Index(i) => length - i,
                SequencePosition::CTerm => 0,
            },
            sequence_length: length,
        }
    }
    /// Check if this position is on the N terminus
    pub fn is_n_terminal(&self) -> bool {
        self.sequence_index == SequencePosition::NTerm
    }
    /// Check if this position is on the C terminus
    pub fn is_c_terminal(&self) -> bool {
        self.sequence_index == SequencePosition::CTerm
    }
    /// Flip to the other series (N->C and C->N)
    #[must_use]
    pub const fn flip_terminal(self) -> Self {
        Self {
            sequence_index: self.sequence_index,
            series_number: self.sequence_length + 1 - self.series_number,
            sequence_length: self.sequence_length,
        }
    }
}

/// The definition of the position of an ion inside a glycan
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct GlycanPosition {
    /// The depth starting at the amino acid
    pub inner_depth: usize,
    /// The series number (from the ion series terminal)
    pub series_number: usize,
    /// The branch naming
    pub branch: Vec<usize>,
    /// The aminoacid index where this glycan is attached
    pub attachment: Option<(AminoAcid, usize)>,
}

impl GlycanPosition {
    /// Get the branch names
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn branch_names(&self) -> String {
        self.branch
            .iter()
            .enumerate()
            .map(|(i, b)| {
                if i == 0 {
                    char::from_u32(
                        (0x03B1..=0x03C9)
                            .chain(0x0391..=0x03A9)
                            .nth(*b)
                            .expect("Too many branches in glycan, out of greek letters"),
                    )
                    .unwrap()
                    .to_string()
                } else if i == 1 {
                    "\'".repeat(*b)
                } else {
                    format!(",{b}")
                }
            })
            .collect::<String>()
    }
    /// Generate the label for this glycan position, example: `1α'`
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn label(&self) -> String {
        format!("{}{}", self.series_number, self.branch_names())
    }
    /// Generate the label for this glycan attachment eg N1 (1 based numbering) or an empty string if the attachment is unknown
    pub fn attachment(&self) -> String {
        self.attachment
            .map(|(aa, pos)| format!("{aa}{pos}"))
            .unwrap_or_default()
    }
}

/// Any position on a glycan or a peptide
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum DiagnosticPosition {
    /// A position on a glycan
    Glycan(GlycanPosition, MonoSaccharide),
    /// A position on a compositional glycan (attachment AA + sequence index + the sugar)
    GlycanCompositional(MonoSaccharide, Option<(AminoAcid, usize)>),
    /// A position on a peptide
    Peptide(PeptidePosition, AminoAcid),
    /// Labile modification
    Labile(Modification),
    /// Reporter ion
    Reporter,
}

/// The possible types of fragments
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize, Default)]
#[expect(non_camel_case_types)]
pub enum FragmentType {
    /// a
    a(PeptidePosition),
    /// b
    b(PeptidePosition),
    /// c
    c(PeptidePosition),
    /// d
    d(PeptidePosition),
    /// v
    v(PeptidePosition),
    /// w
    w(PeptidePosition),
    /// x
    x(PeptidePosition),
    /// y
    y(PeptidePosition),
    /// z
    z(PeptidePosition),
    /// z·
    z·(PeptidePosition),
    // glycan A fragment (Never generated)
    //A(GlycanPosition),
    /// glycan B fragment
    B(GlycanPosition),
    // glycan C fragment (Never generated)
    //C(GlycanPosition),
    // glycan X fragment (Never generated)
    //X(GlycanPosition),
    /// glycan Y fragment, generated by one or more branches broken
    Y(Vec<GlycanPosition>),
    // glycan Z fragment (Never generated)
    // Z(GlycanPosition),
    /// Internal glycan fragment, meaning both a B and Y breakages (and potentially multiple of both), resulting in a set of monosaccharides
    Oxonium(Vec<GlycanBreakPos>),
    /// A B or internal glycan fragment for a glycan where only the composition is known, also saves the attachment (AA + sequence index)
    OxoniumComposition(Vec<(MonoSaccharide, isize)>, Option<(AminoAcid, usize)>),
    /// A B or internal glycan fragment for a glycan where only the composition is known, also saves the attachment (AA + sequence index)
    YComposition(Vec<(MonoSaccharide, isize)>, Option<(AminoAcid, usize)>),
    /// Immonium ion
    Immonium(PeptidePosition, SequenceElement<SemiAmbiguous>),
    /// Precursor with amino acid side chain loss
    PrecursorSideChainLoss(PeptidePosition, AminoAcid),
    /// Diagnostic ion for a given position
    Diagnostic(DiagnosticPosition),
    /// An internal fragment, potentially with the named bonds that resulted in this fragment
    Internal(
        Option<(BackboneNFragment, BackboneCFragment)>,
        PeptidePosition,
        PeptidePosition,
    ),
    /// An unknown series, with potentially the series number
    Unknown(Option<usize>),
    /// precursor
    #[default]
    Precursor,
}

impl FragmentType {
    /// Get the position of this ion (or None if it is a precursor ion)
    pub const fn position(&self) -> Option<&PeptidePosition> {
        match self {
            Self::a(n)
            | Self::b(n)
            | Self::c(n)
            | Self::d(n)
            | Self::v(n)
            | Self::w(n)
            | Self::x(n)
            | Self::y(n)
            | Self::z(n)
            | Self::z·(n)
            | Self::Diagnostic(DiagnosticPosition::Peptide(n, _))
            | Self::Immonium(n, _)
            | Self::PrecursorSideChainLoss(n, _) => Some(n),
            _ => None,
        }
    }

    /// Get the glycan position of this ion (or None not applicable)
    pub const fn glycan_position(&self) -> Option<&GlycanPosition> {
        match self {
            Self::B(n) | Self::Diagnostic(DiagnosticPosition::Glycan(n, _)) => Some(n),
            _ => None,
        }
    }

    /// Get the position label, unless it is a precursor ion
    pub fn position_label(&self) -> Option<String> {
        match self {
            Self::a(n)
            | Self::b(n)
            | Self::c(n)
            | Self::d(n)
            | Self::v(n)
            | Self::w(n)
            | Self::x(n)
            | Self::y(n)
            | Self::z(n)
            | Self::z·(n)
            | Self::Diagnostic(DiagnosticPosition::Peptide(n, _))
            | Self::Immonium(n, _)
            | Self::PrecursorSideChainLoss(n, _) => Some(n.series_number.to_string()),
            Self::B(n) | Self::Diagnostic(DiagnosticPosition::Glycan(n, _)) => Some(n.label()),
            Self::Y(bonds) => Some(bonds.iter().map(GlycanPosition::label).join("")),
            Self::Oxonium(breakages) => Some(
                breakages
                    .iter()
                    .map(std::string::ToString::to_string)
                    .join(""),
            ),
            Self::YComposition(sugars, _) | Self::OxoniumComposition(sugars, _) => Some(
                sugars
                    .iter()
                    .map(|(sugar, amount)| format!("{sugar}{amount}"))
                    .join(""),
            ),
            Self::Internal(_, pos1, pos2) => {
                Some(format!("{}:{}", pos1.sequence_index, pos2.sequence_index,))
            }
            Self::Precursor
            | Self::Unknown(_)
            | Self::Diagnostic(
                DiagnosticPosition::Labile(_)
                | DiagnosticPosition::GlycanCompositional(_, _)
                | DiagnosticPosition::Reporter,
            ) => None,
        }
    }

    /// Get the label for this fragment type
    pub fn label(&self) -> Cow<str> {
        match self {
            Self::a(_) => Cow::Borrowed("a"),
            Self::b(_) => Cow::Borrowed("b"),
            Self::c(_) => Cow::Borrowed("c"),
            Self::d(_) => Cow::Borrowed("d"),
            Self::v(_) => Cow::Borrowed("v"),
            Self::w(_) => Cow::Borrowed("w"),
            Self::x(_) => Cow::Borrowed("x"),
            Self::y(_) => Cow::Borrowed("y"),
            Self::z(_) => Cow::Borrowed("z"),
            Self::z·(_) => Cow::Borrowed("z·"),
            Self::B(_) => Cow::Borrowed("B"),
            Self::Y(_) | Self::YComposition(_, _) => Cow::Borrowed("Y"),
            Self::Diagnostic(DiagnosticPosition::Peptide(_, aa)) => Cow::Owned(format!("d{aa}")),
            Self::Diagnostic(DiagnosticPosition::Reporter) => Cow::Borrowed("r"),
            Self::Diagnostic(DiagnosticPosition::Labile(m)) => Cow::Owned(format!("d{m}")),
            Self::Diagnostic(
                DiagnosticPosition::Glycan(_, sug)
                | DiagnosticPosition::GlycanCompositional(sug, _),
            ) => Cow::Owned(format!("d{sug}")),
            Self::Oxonium(_) | Self::OxoniumComposition(_, _) => Cow::Borrowed("oxonium"),
            Self::Immonium(_, aa) => Cow::Owned(format!("i{}", aa.aminoacid)),
            Self::PrecursorSideChainLoss(_, aa) => Cow::Owned(format!("p-s{aa}")),
            Self::Precursor => Cow::Borrowed("p"),
            Self::Internal(fragmentation, _, _) => Cow::Owned(format!(
                "m{}",
                fragmentation.map_or(String::new(), |(n, c)| format!("{n}:{c}")),
            )),
            Self::Unknown(series) => Cow::Owned(format!(
                "?{}",
                series.map_or(String::new(), |s| s.to_string()),
            )),
        }
    }

    /// Get the kind of fragment, easier to match against
    pub const fn kind(&self) -> FragmentKind {
        match self {
            Self::a(_) => FragmentKind::a,
            Self::b(_) => FragmentKind::b,
            Self::c(_) => FragmentKind::c,
            Self::d(_) => FragmentKind::d,
            Self::v(_) => FragmentKind::v,
            Self::w(_) => FragmentKind::w,
            Self::x(_) => FragmentKind::x,
            Self::y(_) => FragmentKind::y,
            Self::z(_) | Self::z·(_) => FragmentKind::z,
            Self::Y(_) | Self::YComposition(_, _) => FragmentKind::Y,
            Self::Diagnostic(
                DiagnosticPosition::Glycan(_, _) | DiagnosticPosition::GlycanCompositional(_, _),
            )
            | Self::B(_)
            | Self::Oxonium(_)
            | Self::OxoniumComposition(_, _) => FragmentKind::Oxonium,
            Self::Diagnostic(_) => FragmentKind::diagnostic,
            Self::Immonium(_, _) => FragmentKind::immonium,
            Self::PrecursorSideChainLoss(_, _) => FragmentKind::precursor_side_chain_loss,
            Self::Precursor => FragmentKind::precursor,
            Self::Internal(_, _, _) => FragmentKind::internal,
            Self::Unknown(_) => FragmentKind::unknown,
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self.label(),
            self.position_label().unwrap_or_default()
        )
    }
}

/// The possible kinds of N terminal backbone fragments.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
#[expect(non_camel_case_types)]
pub enum BackboneNFragment {
    /// a
    a,
    /// b
    b,
    /// c
    c,
}

impl Display for BackboneNFragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a => "a",
                Self::b => "b",
                Self::c => "c",
            }
        )
    }
}

/// The possible kinds of C terminal backbone fragments.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
#[expect(non_camel_case_types)]
pub enum BackboneCFragment {
    /// x
    x,
    /// y
    y,
    /// z and z·
    z,
}

impl Display for BackboneCFragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::x => "x",
                Self::y => "y",
                Self::z => "z",
            }
        )
    }
}

/// The possible kinds of fragments, same options as [`FragmentType`] but without any additional data
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
#[expect(non_camel_case_types)]
pub enum FragmentKind {
    /// a
    a,
    /// b
    b,
    /// c
    c,
    /// d
    d,
    /// v
    v,
    /// w
    w,
    /// x
    x,
    /// y
    y,
    /// z and z·
    z,
    /// glycan Y fragment, generated by one or more branches broken
    Y,
    /// B or glycan diagnostic ion or Internal glycan fragment, meaning both a B and Y breakages (and potentially multiple of both), resulting in a set of monosaccharides
    Oxonium,
    /// Immonium ion
    immonium,
    /// Precursor with amino acid side chain loss
    precursor_side_chain_loss,
    /// Diagnostic ion for a given position
    diagnostic,
    /// Internal ion
    internal,
    /// precursor
    precursor,
    /// unknown fragment
    unknown,
}

impl Display for FragmentKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a => "a",
                Self::b => "b",
                Self::c => "c",
                Self::d => "d",
                Self::x => "x",
                Self::y => "y",
                Self::v => "v",
                Self::w => "w",
                Self::z => "z",
                Self::Y => "Y",
                Self::Oxonium => "oxonium",
                Self::immonium => "immonium",
                Self::precursor_side_chain_loss => "precursor side chain loss",
                Self::diagnostic => "diagnostic",
                Self::internal => "m",
                Self::precursor => "precursor",
                Self::unknown => "unknown",
            }
        )
    }
}

/// All positions where a glycan can break
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum GlycanBreakPos {
    /// No breaks just until the end of a chain
    End(GlycanPosition),
    /// Break at a Y position
    Y(GlycanPosition),
    /// Break at a B position
    B(GlycanPosition),
}

impl GlycanBreakPos {
    /// Get the position of this breaking position
    pub const fn position(&self) -> &GlycanPosition {
        match self {
            Self::B(p) | Self::End(p) | Self::Y(p) => p,
        }
    }

    /// Get the label for this breaking position
    pub const fn label(&self) -> &str {
        match self {
            Self::End(_) => "End",
            Self::Y(_) => "Y",
            Self::B(_) => "B",
        }
    }
}

impl std::fmt::Display for GlycanBreakPos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.label(), self.position().label())
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {

    use crate::{AminoAcid, MultiChemical};

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
