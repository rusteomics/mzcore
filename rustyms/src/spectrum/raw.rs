//! Raw spectra (not annotated)

use std::cmp::Ordering;

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    spectrum::PeakSpectrum,
    system::{
        f64::{Mass, MassOverCharge, Ratio, Time},
        isize::Charge,
    },
};

/// A raw spectrum (meaning not annotated yet)
#[derive(Clone, Debug, Default, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct RawSpectrum {
    /// The title (as used in MGF)
    pub title: String,
    /// The number of scans
    pub num_scans: u64,
    /// The retention time
    pub rt: Option<Time>,
    /// The found precursor charge
    pub charge: Option<Charge>,
    /// The found precursor mass
    pub mass: Option<Mass>,
    /// The found precursor intensity
    pub intensity: Option<f64>,
    /// The peaks of which this spectrum consists
    pub(crate) spectrum: Vec<RawPeak>,
    /// MGF: if present the SEQUENCE line
    pub sequence: Option<String>,
    /// MGF TITLE: if present the raw file where this mgf was made from
    pub raw_file: Option<String>,
    /// MGF TITLE: if present the raw file scan number
    pub raw_scan_number: Option<usize>,
    /// MGF TITLE: index number
    pub raw_index: Option<usize>,
    /// MGF TITLE: sample number
    pub sample: Option<usize>,
    /// MGF TITLE: period number
    pub period: Option<usize>,
    /// MGF TITLE: cycle number
    pub cycle: Option<usize>,
    /// MGF TITLE: experiment number
    pub experiment: Option<usize>,
    /// MGF TITLE: controllerType number
    pub controller_type: Option<usize>,
    /// MGF TITLE: controllerNumber number
    pub controller_number: Option<usize>,
}

impl RawSpectrum {
    /// Filter the spectrum to retain all with an intensity above `filter_threshold` times the maximal intensity.
    ///
    /// # Panics
    /// It panics if any peaks has an intensity that is NaN.
    pub fn relative_noise_filter(&mut self, filter_threshold: f64) {
        let max = self
            .spectrum
            .iter()
            .map(|p| *p.intensity)
            .reduce(f64::max)
            .unwrap_or(f64::INFINITY);
        self.spectrum
            .retain(|p| *p.intensity >= max * filter_threshold);
        self.spectrum.shrink_to_fit();
    }

    /// Filter the spectrum to retain all with an intensity above `filter_threshold`.
    pub fn absolute_noise_filter(&mut self, filter_threshold: f64) {
        self.spectrum.retain(|p| *p.intensity >= filter_threshold);
        self.spectrum.shrink_to_fit();
    }

    /// Filter a spectrum by dividing it in windows and within each window only retain the `top` number of peaks.
    #[expect(clippy::missing_panics_doc)] // Cannot panic as it checks with peek first
    pub fn top_x_filter(&mut self, window_size: f64, top: usize) {
        let mut new_spectrum = Vec::with_capacity(
            self.spectrum
                .last()
                .and_then(|l| self.spectrum.first().map(|f| (f, l)))
                .map(|(f, l)| ((l.mz.value - f.mz.value) / window_size).round() as usize * top)
                .unwrap_or_default(),
        );
        let mut spectrum = self.spectrum.iter().copied().peekable();
        let mut window = 1;
        let mut peaks = Vec::new();

        while spectrum.peek().is_some() {
            while let Some(peek) = spectrum.peek() {
                if peek.mz.value <= f64::from(window) * window_size {
                    peaks.push(spectrum.next().unwrap());
                } else {
                    break;
                }
            }
            new_spectrum.extend(
                peaks
                    .iter()
                    .copied()
                    .k_largest_by(top, |a, b| a.mz.value.total_cmp(&b.mz.value)),
            );
            peaks.clear();
            window += 1;
        }

        self.spectrum = new_spectrum;
    }
}

impl Extend<RawPeak> for RawSpectrum {
    fn extend<T: IntoIterator<Item = RawPeak>>(&mut self, iter: T) {
        self.spectrum.extend(iter);
        self.spectrum.sort_unstable();
    }
}

impl IntoIterator for RawSpectrum {
    type Item = RawPeak;
    type IntoIter = std::vec::IntoIter<RawPeak>;
    fn into_iter(self) -> Self::IntoIter {
        self.spectrum.into_iter()
    }
}

impl std::ops::Index<usize> for RawSpectrum {
    type Output = RawPeak;
    fn index(&self, index: usize) -> &Self::Output {
        &self.spectrum[index]
    }
}

impl PeakSpectrum for RawSpectrum {
    type PeakType = RawPeak;
    type Iter<'a> = std::slice::Iter<'a, Self::PeakType>;

    /// Return the slice of peaks that is within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[RawPeak] {
        let left_idx = match self
            .spectrum
            .binary_search_by(|a| a.mz.value.total_cmp(&low.value))
        {
            Ok(idx) | Err(idx) => {
                let mut idx = idx.saturating_sub(1);
                while idx > 0 && self.spectrum[idx].mz.value.total_cmp(&low.value) != Ordering::Less
                {
                    idx -= 1;
                }
                idx
            }
        };

        let right_idx = match self.spectrum[left_idx..]
            .binary_search_by(|a| a.mz.value.total_cmp(&high.value))
        {
            Ok(idx) | Err(idx) => {
                let mut idx = idx + left_idx;
                while idx < self.spectrum.len()
                    && self.spectrum[idx].mz.value.total_cmp(&high.value) != Ordering::Greater
                {
                    idx = idx.saturating_add(1);
                }
                idx.min(self.spectrum.len())
            }
        };
        &self.spectrum[left_idx..right_idx]
    }

    fn spectrum(&self) -> Self::Iter<'_> {
        self.spectrum.iter()
    }

    fn add_peak(&mut self, item: Self::PeakType) {
        let index = self.spectrum.binary_search(&item).unwrap_or_else(|i| i);
        self.spectrum.insert(index, item);
    }
}

/// A raw peak
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct RawPeak {
    /// The mz value of this peak
    pub mz: MassOverCharge,
    /// The intensity of this peak
    pub intensity: OrderedFloat<f64>,
}

impl PartialOrd for RawPeak {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for RawPeak {
    /// Use `f64::total_cmp` on `self.mz`
    fn cmp(&self, other: &Self) -> Ordering {
        self.mz.value.total_cmp(&other.mz.value)
    }
}

impl PartialEq for RawPeak {
    /// Use `f64::total_cmp` on all fields to detect total equality
    fn eq(&self, other: &Self) -> bool {
        self.mz.value.total_cmp(&other.mz.value) == Ordering::Equal
            && self.intensity.total_cmp(&other.intensity) == Ordering::Equal
    }
}

impl Eq for RawPeak {}

impl RawPeak {
    /// Determine the ppm error for the given fragment
    pub fn ppm(&self, mz: MassOverCharge) -> Ratio {
        self.mz.ppm(mz)
    }
}
