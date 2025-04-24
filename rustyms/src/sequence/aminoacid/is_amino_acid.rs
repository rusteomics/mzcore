//! Module used to create the [`IsAminoAcid`] trait

use crate::{
    chemistry::{MassMode, MolecularFormula, MultiChemical},
    fragment::SatelliteLabel,
    quantities::Multi,
    sequence::SequencePosition,
    system::Mass,
};

use std::borrow::Cow;

/// A general trait to define amino acids.
pub trait IsAminoAcid: MultiChemical {
    /// The full name for this amino acid.
    fn name(&self) -> Cow<'_, str>;
    /// The three letter code for this amino acid. Or None if there is no common three letter
    /// definition for this amino acid.
    fn three_letter_code(&self) -> Option<Cow<'_, str>>;
    /// The one letter code for this amino acid. Or None if there is no common single character
    /// definition for this amino acid.
    #[doc(alias = "code")]
    fn one_letter_code(&self) -> Option<char>;
    /// The ProForma definition for this amino acid. If this is not a simple amino acid it can be
    /// defined as an amino acid with an additional modification. For example `X[H9C2N2]` could be
    /// used if Arginine was not defined as `R` in ProForma.
    fn pro_forma_definition(&self) -> Cow<'_, str>;
    /// The monoisotopic mass of this amino acid. Should be redefined for better performance.
    fn monoisotopic_mass(&self) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(
            self.formulas()
                .iter()
                .map(MolecularFormula::monoisotopic_mass)
                .collect(),
        )
    }
    /// The average weight of this amino acid. Should be redefined for better performance.
    fn average_weight(&self) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(
            self.formulas()
                .iter()
                .map(MolecularFormula::average_weight)
                .collect(),
        )
    }
    /// The mass with a given mass mode for this amino acid. Should be redefined for better performance.
    fn mass(&self, mode: MassMode) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(self.formulas().iter().map(|f| f.mass(mode)).collect())
    }
    /// The molecular formula of the side chain of the amino acid. The `sequence_index` and
    /// `peptidoform_index` are used to keep track of ambiguous amino acids.
    fn side_chain(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Cow<'_, Multi<MolecularFormula>>;
    /// The molecular formulas that can fragment for satellite ions (d and w). Commonly the fragment
    /// after the second carbon into the side chain. `MolecularFormula::default()` can be returned
    /// if no satellite ions are possible. The `sequence_index` and `peptidoform_index` are used to
    /// keep track of ambiguous amino acids.
    fn satellite_ion_fragments(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Option<Cow<'_, Vec<(SatelliteLabel, MolecularFormula)>>>;
}
