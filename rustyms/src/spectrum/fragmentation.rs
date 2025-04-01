use crate::{
    model::MatchingParameters, system::MassOverCharge, CompoundPeptidoformIon, Fragment, MassMode,
};

use super::AnnotatedSpectrum;

/// A spectrum that can be annotated. Within rustyms this is implemented for the build in
/// [mgf reader](crate::rawfile::mgf) and for mzdata [`SpectrumLike`](mzdata::prelude::SpectrumLike).
/// For up to date information see that crate, but at the moment of writing this supports mgf, mzML,
/// indexed mzML, mzMLb, and Thermo RAW. Note that any 'Missing' and
/// [`RawData`](mzdata::spectrum::RawSpectrum) from mzdata result in an empty annotated spectrum.
/// Also note that the feature `mzdata` is required for the mzdata spectra to work.
pub trait AnnotatableSpectrum {
    /// The tolerance type that should be used for searching peaks in the spectrum.
    type Tolerance: From<crate::Tolerance<MassOverCharge>> + Copy;

    /// Create an empty annotated spectrum, which is required to fill the spectrum vector with
    /// [`blank`](crate::spectrum::AnnotatedPeak::background) annotated peaks.
    fn empty_annotated(&self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum;

    /// Search for a specific mz within the tolerance. Has to return the index in the annotated
    /// spectrum vector for closest peak (if there is any).
    fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize>;

    /// Annotate this spectrum with the given peptidoform and given fragments see
    /// [`crate::CompoundPeptidoform::generate_theoretical_fragments`].
    fn annotate(
        &self,
        peptide: CompoundPeptidoformIon,
        theoretical_fragments: &[Fragment],
        parameters: &MatchingParameters,
        mode: MassMode,
    ) -> AnnotatedSpectrum {
        let tolerance = parameters.tolerance.into();
        let mut annotated = Self::empty_annotated(self, peptide);

        for fragment in theoretical_fragments {
            // Determine fragment mz and see if it is within the model range.
            if let Some(mz) = fragment.mz(mode) {
                if !parameters.mz_range.contains(&mz) {
                    continue;
                }

                // Get the index of the element closest to this value
                if let Some(index) = Self::search(self, mz, tolerance) {
                    // Keep the theoretical fragments sorted to have the highest theoretical likelihood on top
                    match annotated.spectrum[index].annotation.binary_search(fragment) {
                        Ok(ai) | Err(ai) => annotated.spectrum[index]
                            .annotation
                            .insert(ai, fragment.clone()),
                    }
                }
            }
        }

        annotated
    }
}
