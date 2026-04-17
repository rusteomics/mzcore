use std::num::NonZeroU32;

use mzcore::{
    prelude::{MassMode, PeptidoformIonSet},
    system::thomson,
};
use mzdata::mzpeaks::PeakCollection;

use crate::prelude::{AnnotatedSpectrum, Fragment, MatchingParameters};

#[cfg(feature = "isotopes")]
use crate::fragment::Isotope;
#[cfg(feature = "isotopes")]
use mzcore::{molecular_formula, system::MassOverCharge};

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
        #[cfg(feature = "isotopes")]
        let isotope_tolerance = match parameters.isotope_tolerance {
            mzcore::quantities::Tolerance::Absolute(mz) => mzdata::prelude::Tolerance::Da(mz.value),
            mzcore::quantities::Tolerance::Relative(ratio) => {
                mzdata::prelude::Tolerance::PPM(ratio.get::<mzcore::system::ratio::ppm>())
            }
        };
        #[cfg(feature = "isotopes")]
        let isotope_shift = molecular_formula!([13 C 1] [12 C -1]).mass(MassMode::Monoisotopic);
        let mut annotated = Self::empty_annotated(self, peptide);

        for fragment in theoretical_fragments {
            // Determine fragment mz and see if it is within the model range.
            if let Some(mz) = fragment.mz(mode) {
                if !parameters.mz_range.contains(&mz) {
                    continue;
                }

                // Get the index of the element closest to this value
                if let Some(index) = annotated.peaks.search(mz.get::<thomson>(), tolerance) {
                    #[cfg(feature = "isotopes")]
                    if parameters.match_isotopes
                        && let Some(formula) = &fragment.formula
                    {
                        let envelope =
                            formula.isotopic_distribution(parameters.isotope_minimum_probability);
                        let base = if mode == MassMode::Monoisotopic {
                            mz
                        } else {
                            formula.mass(MassMode::Monoisotopic) / fragment.charge.to_float()
                        };
                        let mut matched_envelope = Vec::with_capacity(envelope.len());
                        let mut matches = Vec::with_capacity(envelope.len());
                        for (index, intensity) in envelope.into_iter().enumerate() {
                            let isotope_mz =
                                base + (index as f64 * isotope_shift) / fragment.charge.to_float();
                            let peak = annotated
                                .peaks
                                .search(isotope_mz.get::<thomson>(), isotope_tolerance);
                            matched_envelope.push((
                                isotope_mz.value,
                                peak.map_or(0.0, |i| annotated.peaks[i].mz.value),
                                intensity,
                                peak.map_or(0.0, |i| f64::from(annotated.peaks[i].intensity)),
                            ));
                            matches.push(peak);
                        }
                        let similarity = cosine_similarity(&matched_envelope);
                        if similarity >= parameters.isotope_filter {
                            for (index, peak) in matches.into_iter().enumerate() {
                                if let Some(peak) = peak {
                                    let frag = fragment
                                        .clone()
                                        .with_isotope(&[(index as i32, Isotope::General)]);
                                    let frag_index = annotated.peaks[peak]
                                        .annotations
                                        .binary_search(&frag)
                                        .unwrap_or_else(|i| i);
                                    annotated.peaks[peak].annotations.insert(frag_index, frag);
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

                    #[cfg(not(feature = "isotopes"))]
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

// Weighted cosine spectral similarity = https://mzmine.github.io/mzmine_documentation/module_docs/id_spectral_library_search/spectral-similarity-measures.html
/// (mz1, mz2, i1, i2)
#[cfg(feature = "isotopes")]
fn cosine_similarity(peaks: &[(f64, f64, f64, f64)]) -> f64 {
    let sumx = peaks
        .iter()
        .map(|(m1, m2, i1, i2)| m1 * m2 * i1 * i2)
        .sum::<f64>();
    let sumx2 = peaks
        .iter()
        .map(|(m1, _, i1, _)| (m1 * i1).powi(2))
        .sum::<f64>()
        .sqrt()
        * peaks
            .iter()
            .map(|(_, m2, _, i2)| (m2 * i2).powi(2))
            .sum::<f64>()
            .sqrt();
    sumx / sumx2
}
