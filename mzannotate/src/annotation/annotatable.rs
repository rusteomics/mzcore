use mzcore::prelude::{CompoundPeptidoformIon, MassMode};
use mzdata::mzpeaks::PeakCollection;

use crate::prelude::{AnnotatedSpectrum, Fragment, MatchingParameters};

/// A spectrum that can be annotated. Within rustyms this is implemented for the build in
/// [mgf reader](crate::spectrum::mgf) and for mzdata [`SpectrumLike`](mzdata::prelude::SpectrumLike).
/// For up to date information see that crate, but at the moment of writing this supports mgf, mzML,
/// indexed mzML, mzMLb, Thermo RAW, and Bruker TDF. Note that any 'Missing' and
/// [`RawData`](mzdata::spectrum::RawSpectrum) from mzdata result in an empty annotated spectrum.
/// Also note that the feature `mzdata` is required for the mzdata spectra to work.
pub trait AnnotatableSpectrum: Sized {
    /// Create an empty annotated spectrum, which is required to fill the spectrum vector with
    /// [`blank`](crate::annotation::AnnotatedPeak::background) annotated peaks.
    fn empty_annotated(self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum;

    /// Annotate this spectrum with the given peptidoform and given fragments see
    /// [`crate::sequence::CompoundPeptidoformIon::generate_theoretical_fragments`]
    /// to generate the fragments.
    fn annotate(
        self,
        peptide: CompoundPeptidoformIon,
        theoretical_fragments: &[Fragment],
        parameters: &MatchingParameters,
        mode: MassMode,
    ) -> AnnotatedSpectrum {
        let tolerance = match parameters.tolerance {
            mzcore::quantities::Tolerance::Absolute(mz) => mzdata::prelude::Tolerance::Da(mz.value),
            mzcore::quantities::Tolerance::Relative(ratio) => {
                mzdata::prelude::Tolerance::PPM(ratio.get::<mzcore::system::ratio::ppm>())
            }
        };
        let mut annotated = Self::empty_annotated(self, peptide);

        for fragment in theoretical_fragments {
            // Determine fragment mz and see if it is within the model range.
            if let Some(mz) = fragment.mz(mode) {
                if !parameters.mz_range.contains(&mz) {
                    continue;
                }

                // Get the index of the element closest to this value
                if let Some(index) = annotated
                    .peaks
                    .search(mz.get::<mzcore::system::thomson>(), tolerance)
                {
                    // Keep the theoretical fragments sorted to have the highest theoretical likelihood on top
                    match annotated.peaks[index].annotations.binary_search(fragment) {
                        Ok(ai) | Err(ai) => annotated.peaks[index]
                            .annotations
                            .insert(ai, fragment.clone()),
                    }
                }
            }
        }

        annotated
    }
}

impl<T: Into<AnnotatedSpectrum>> AnnotatableSpectrum for T {
    fn empty_annotated(self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum {
        let mut spectrum: AnnotatedSpectrum = self.into();
        for (index, peptidoform_ion) in peptide.into_peptidoform_ions().into_iter().enumerate() {
            spectrum.add_analyte(crate::mzspeclib::Analyte::new(
                index as u32,
                crate::mzspeclib::AnalyteTarget::PeptidoformIon(peptidoform_ion),
                Vec::new(),
            ));
        }

        spectrum
    }
}
