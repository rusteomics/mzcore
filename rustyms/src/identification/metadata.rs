use std::{borrow::Cow, ops::Range};

use crate::{
    identification::{FastaIdentifier, KnownFileFormat, SpectrumIds},
    prelude::CompoundPeptidoformIon,
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

/// Generalised access to meta data of identified peptidoforms
pub trait MetaData {
    /// Get the compound peptidoform ion, if present
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>>;

    /// Get the format and version for this peptidoform
    fn format(&self) -> KnownFileFormat;

    /// Get the PSM identifier
    fn id(&self) -> String;

    /// Get the confidence, a score between -1 and 1 describing the confidence in the entire PSM
    fn confidence(&self) -> Option<f64>;

    /// Get the local confidence, a score between -1 and 1 for each amino acid in the peptide
    fn local_confidence(&self) -> Option<&[f64]>;

    /// Get the original confidence
    fn original_confidence(&self) -> Option<f64>;

    /// Get the original local confidence, a score for each amino acid in the peptide
    fn original_local_confidence(&self) -> Option<&[f64]>;

    /// The charge of the precursor/PSM, if known
    fn charge(&self) -> Option<Charge>;

    /// Which fragmentation mode was used, if known
    fn mode(&self) -> Option<&str>; // TODO: should create an enum or use mzdata formats at some point

    /// The retention time, if known
    fn retention_time(&self) -> Option<Time>;

    /// The scans per rawfile that are at the basis for this identified peptide
    fn scans(&self) -> SpectrumIds;

    /// Get the mz as experimentally determined
    fn experimental_mz(&self) -> Option<MassOverCharge>;

    /// Get the mass as experimentally determined
    fn experimental_mass(&self) -> Option<Mass>;

    /// Get the absolute ppm error between the experimental and theoretical precursor mass, if there are multiple masses possible returns the smallest ppm
    fn ppm_error(&self) -> Option<Ratio> {
        let exp_mass = self.experimental_mass()?;
        self.compound_peptidoform_ion().and_then(|f| {
            f.formulas()
                .iter()
                .map(|theo_mass| theo_mass.monoisotopic_mass().ppm(exp_mass))
                .min_by(|a, b| a.value.total_cmp(&b.value))
        })
    }

    /// Get the absolute mass error between the experimental and theoretical precursor mass, if there are multiple masses possible returns the smallest difference
    fn mass_error(&self) -> Option<Mass> {
        let exp_mass = self.experimental_mass()?;
        self.compound_peptidoform_ion().and_then(|f| {
            f.formulas()
                .iter()
                .map(|theo_mass| (exp_mass - theo_mass.monoisotopic_mass()).abs())
                .min_by(|a, b| a.value.total_cmp(&b.value))
        })
    }

    /// Get the protein name if this was database matched data
    fn protein_name(&self) -> Option<FastaIdentifier<String>>;

    /// Get the protein id if this was database matched data
    fn protein_id(&self) -> Option<usize>;

    /// Get the protein location if this was database matched data
    fn protein_location(&self) -> Option<Range<usize>>;

    // Get the matched fragments, potentially with m/z and intensity
    // #[doc(hidden)]
    // pub fn matched_fragments(
    //     &self,
    // ) -> Option<Vec<(Option<MassOverCharge>, Option<f64>, Fragment)>> {
    //     // OPair, MaxQuant, PLGS
    //     None
    // }
}
