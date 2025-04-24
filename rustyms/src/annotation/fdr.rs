use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::{model::MatchingParameters, AnnotatedSpectrum},
    chemistry::MassMode,
    fragment::Fragment,
    quantities::WithinTolerance,
    system::{MassOverCharge, Ratio},
};

impl AnnotatedSpectrum {
    /// Get a false discovery rate estimation for this annotation. See the [`Fdr`] struct for all
    /// statistics that can be retrieved. The returned tuple has the FDR for all peptides combined
    /// as first item and as second item a vector with for each peptide its individual scores.
    ///
    /// # Estimation
    /// It estimated the FDR by permutation. It takes all MZs for all fragments (or only the subset
    /// for each separate peptide) and counts how many match with the spectrum with offset
    /// `-25..=25 + π`. The returned result give insight in the average number of matches and
    /// standard deviation. The `π` offset is needed to guarantee non integer offsets preventing
    /// spurious matches from 1 Da isotopes.
    pub fn fdr(
        &self,
        fragments: &[Fragment],
        parameters: &MatchingParameters,
        mass_mode: MassMode,
    ) -> (Fdr, Vec<Vec<Fdr>>) {
        let mzs = fragments
            .iter()
            .filter_map(|f| {
                f.mz(mass_mode)
                    .map(|mz| (mz, f.peptidoform_ion_index, f.peptidoform_index))
            })
            .filter(|(mz, _, _)| parameters.mz_range.contains(mz))
            .collect_vec();

        let individual_peptides = self
            .peptide
            .peptidoform_ions()
            .iter()
            .enumerate()
            .map(|(peptidoform_ion_index, peptidoform)| {
                peptidoform
                    .peptidoforms()
                    .iter()
                    .enumerate()
                    .map(|(peptidoform_index, _)| {
                        self.internal_fdr(
                            mzs.iter()
                                .filter_map(|(mz, pi, ppi)| {
                                    (*pi == Some(peptidoform_ion_index)
                                        && *ppi == Some(peptidoform_index))
                                    .then_some(*mz)
                                })
                                .collect_vec()
                                .as_slice(),
                            parameters,
                        )
                    })
                    .collect()
            })
            .collect();
        (
            self.internal_fdr(
                mzs.iter().map(|(mz, _, _)| *mz).collect_vec().as_slice(),
                parameters,
            ),
            individual_peptides,
        )
    }

    fn internal_fdr(&self, mzs: &[MassOverCharge], parameters: &MatchingParameters) -> Fdr {
        let mut results = Vec::with_capacity(51);
        let total_intensity = self.spectrum.iter().map(|s| s.intensity.0).sum::<f64>();

        for offset in -25..=25 {
            let peaks = self
                .spectrum
                .iter()
                .map(|p| {
                    p.experimental_mz
                        + MassOverCharge::new::<crate::system::mz>(
                            std::f64::consts::PI + f64::from(offset),
                        )
                })
                .collect_vec();
            let mut peak_annotated = vec![false; peaks.len()];
            let mut number_peaks_annotated = 0;
            let mut intensity_annotated = 0.0;
            for mass in mzs {
                // Get the index of the element closest to this value (spectrum is defined to always be sorted)
                let index = peaks
                    .binary_search_by(|p| p.value.total_cmp(&mass.value))
                    .unwrap_or_else(|i| i);

                // Check index-1, index and index+1 (if existing) to find the one with the lowest ppm
                let mut closest = (0, Ratio::new::<crate::system::ratio::ppm>(f64::INFINITY));
                #[expect(clippy::needless_range_loop)] // I like this better
                for i in if index == 0 { 0 } else { index - 1 }
                    ..=(index + 1).min(self.spectrum.len().saturating_sub(1))
                {
                    let ppm = peaks[i].ppm(*mass);
                    if ppm < closest.1 {
                        closest = (i, ppm);
                    }
                }

                if parameters
                    .tolerance
                    .within(&self.spectrum[closest.0].experimental_mz, mass)
                    && !peak_annotated[closest.0]
                {
                    number_peaks_annotated += 1;
                    intensity_annotated += self.spectrum[closest.0].intensity.0;
                    peak_annotated[closest.0] = true;
                }
            }
            results.push((
                f64::from(number_peaks_annotated) / self.spectrum.len() as f64,
                intensity_annotated / total_intensity,
            ));
        }
        let peaks_average = results.iter().map(|r| r.0).sum::<f64>() / results.len() as f64;
        let peaks_st_dev = (results
            .iter()
            .map(|x| (x.0 - peaks_average).powi(2))
            .sum::<f64>()
            / results.len() as f64)
            .sqrt();
        let intensity_average = results.iter().map(|r| r.1).sum::<f64>() / results.len() as f64;
        let intensity_st_dev = (results
            .iter()
            .map(|x| (x.1 - intensity_average).powi(2))
            .sum::<f64>()
            / results.len() as f64)
            .sqrt();
        let actual: (u32, f64) = self
            .spectrum
            .iter()
            .filter(|p| !p.annotation.is_empty())
            .fold((0, 0.0), |acc, p| (acc.0 + 1, acc.1 + p.intensity.0));

        Fdr {
            peaks_actual: f64::from(actual.0) / self.spectrum.len() as f64,
            peaks_average_false: peaks_average,
            peaks_standard_deviation_false: peaks_st_dev,
            intensity_actual: actual.1 / total_intensity,
            intensity_average_false: intensity_average,
            intensity_standard_deviation_false: intensity_st_dev,
        }
    }
}

/// A false discovery rate for an annotation to a spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct Fdr {
    /// The fraction of the total (assumed to be true) peaks that could be annotated
    pub peaks_actual: f64,
    /// The average fraction of the false peaks that could be annotated
    pub peaks_average_false: f64,
    /// The standard deviation of the false peaks that could be annotated
    pub peaks_standard_deviation_false: f64,
    /// The fraction of the total (assumed to be true) intensity that could be annotated
    pub intensity_actual: f64,
    /// The average fraction of the false intensity that could be annotated
    pub intensity_average_false: f64,
    /// The standard deviation of the false intensity that could be annotated
    pub intensity_standard_deviation_false: f64,
}

impl Fdr {
    /// Get the false discovery rate (as a fraction).
    /// The average number of false peaks annotated divided by the average number of annotated peaks.
    pub fn peaks_fdr(&self) -> f64 {
        self.peaks_average_false / self.peaks_actual
    }

    /// Get the number of standard deviations the number of annotated peaks is from the average number of false annotations.
    pub fn peaks_sigma(&self) -> f64 {
        (self.peaks_actual - self.peaks_average_false) / self.peaks_standard_deviation_false
    }

    /// Get the peaks score of this annotation. Defined as the log2 of the sigma.
    pub fn peaks_score(&self) -> f64 {
        self.peaks_sigma().log2()
    }

    /// Get the false discovery rate (as a fraction).
    /// The average number of false intensity annotated divided by the average number of annotated intensity.
    pub fn intensity_fdr(&self) -> f64 {
        self.intensity_average_false / self.intensity_actual
    }

    /// Get the number of standard deviations the annotated intensity is from the average false annotations.
    pub fn intensity_sigma(&self) -> f64 {
        (self.intensity_actual - self.intensity_average_false)
            / self.intensity_standard_deviation_false
    }

    /// Get the intensity score of this annotation. Defined as the log2 of the sigma.
    pub fn intensity_score(&self) -> f64 {
        self.intensity_sigma().log2()
    }
}
