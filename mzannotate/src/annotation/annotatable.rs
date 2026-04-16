use std::num::NonZeroU32;

use mzcore::{
    molecular_formula,
    prelude::{MassMode, PeptidoformIonSet},
    system::{MassOverCharge, thomson},
};
use mzdata::mzpeaks::PeakCollection;

use crate::prelude::{AnnotatedSpectrum, Fragment, MatchingParameters};

/// A spectrum that can be annotated. The best way to use this is with mzdata
/// [`SpectrumLike`](mzdata::prelude::SpectrumLike). For up to date information see that crate, but
/// at the moment of writing this supports mgf, mzML, indexed mzML, mzMLb, Thermo RAW, and Bruker
/// TDF. Note this only takes the centroided data and that any 'Missing' and
/// [`RawData`](mzdata::spectrum::RawSpectrum) from mzdata result in an empty annotated spectrum.
/// Also note that the feature `mzdata` is required for the mzdata spectra to work.
pub trait AnnotatableSpectrum: Sized {
    /// Create an empty annotated spectrum. This spectrum is assumed to contain all peaks but
    /// without any annotations.
    fn empty_annotated(self, peptide: PeptidoformIonSet) -> AnnotatedSpectrum;

    /// Annotate this spectrum with the given peptidoform and given fragments see
    /// [`crate::prelude::PeptidoformFragmentation::generate_theoretical_fragments`]
    /// to generate the fragments.
    fn annotate(
        self,
        peptide: PeptidoformIonSet,
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
        let isotope_tolerance = match parameters.isotope_tolerance {
            mzcore::quantities::Tolerance::Absolute(mz) => mzdata::prelude::Tolerance::Da(mz.value),
            mzcore::quantities::Tolerance::Relative(ratio) => {
                mzdata::prelude::Tolerance::PPM(ratio.get::<mzcore::system::ratio::ppm>())
            }
        };
        let mut annotated = Self::empty_annotated(self, peptide);
        let isotope_shift = molecular_formula!([13 C 1] [12 C -1]).mass(MassMode::Monoisotopic);

        for fragment in theoretical_fragments {
            // Determine fragment mz and see if it is within the model range.
            if let Some(mz) = fragment.mz(mode) {
                if !parameters.mz_range.contains(&mz) {
                    continue;
                }

                // Get the index of the element closest to this value
                if let Some(index) = annotated.peaks.search(mz.get::<thomson>(), tolerance) {
                    // #[cfg(feature = "mzdata/isotopes")]
                    if parameters.match_isotopes
                        && let Some(formula) = &fragment.formula
                    {
                        let offset = annotated.peaks[index].mz - mz;
                        let envelope = formula.isotopic_distribution(0.01);
                        let base = if mode == MassMode::Monoisotopic {
                            mz
                        } else {
                            formula.mass(MassMode::Monoisotopic) / fragment.charge.to_float()
                        } + offset;
                        let mut matches = Vec::with_capacity(envelope.len());
                        for (index, intensity) in envelope.into_iter().enumerate() {
                            let isotope_mz = base
                                + MassOverCharge::new::<thomson>(
                                    (index as f64 * isotope_shift).value,
                                );
                            matches.push((
                                intensity,
                                annotated
                                    .peaks
                                    .search(isotope_mz.get::<thomson>(), isotope_tolerance),
                            ));
                        }
                        let similarity = todo!();
                        if similarity >= parameters.isotope_filter {
                            for (index, (_, peak)) in matches.into_iter().enumerate() {
                                if let Some(peak) = peak {
                                    let frag = fragment.clone();
                                    // Store the isotope info somewhere and add to the annotated peaks
                                    annotated.peaks[peak].annotations.insert(index, element);
                                }
                            }
                        }
                    } else {
                        // Keep the theoretical fragments sorted to have the highest theoretical likelihood on top
                        match annotated.peaks[index].annotations.binary_search(fragment) {
                            Ok(ai) | Err(ai) => annotated.peaks[index]
                                .annotations
                                .insert(ai, fragment.clone()),
                        }
                    }
                }
            }
        }

        annotated
    }
}

impl<T: Into<AnnotatedSpectrum>> AnnotatableSpectrum for T {
    fn empty_annotated(self, peptide: PeptidoformIonSet) -> AnnotatedSpectrum {
        let mut spectrum: AnnotatedSpectrum = self.into();
        for (index, peptidoform_ion) in peptide.into_peptidoform_ions().into_iter().enumerate() {
            spectrum.analytes.push(crate::mzspeclib::Analyte::new(
                NonZeroU32::new(index as u32 + 1).unwrap(),
                crate::mzspeclib::AnalyteTarget::PeptidoformIon(peptidoform_ion),
            ));
        }

        spectrum
    }
}
