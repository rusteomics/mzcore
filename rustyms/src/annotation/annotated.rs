//! Annotated spectra

use std::cmp::Ordering;

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    fragment::Fragment,
    sequence::CompoundPeptidoformIon,
    spectrum::{PeakSpectrum, RawPeak},
    system::{
        f64::{Mass, MassOverCharge, Time},
        usize::Charge,
    },
};

/// An annotated spectrum
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct AnnotatedSpectrum {
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
    /// The peptide with which this spectrum was annotated
    pub peptide: CompoundPeptidoformIon,
    /// The spectrum
    pub(super) spectrum: Vec<AnnotatedPeak>,
}

impl Extend<AnnotatedPeak> for AnnotatedSpectrum {
    fn extend<T: IntoIterator<Item = AnnotatedPeak>>(&mut self, iter: T) {
        self.spectrum.extend(iter);
        self.spectrum.sort_unstable();
    }
}

impl IntoIterator for AnnotatedSpectrum {
    type Item = AnnotatedPeak;
    type IntoIter = std::vec::IntoIter<AnnotatedPeak>;
    fn into_iter(self) -> Self::IntoIter {
        self.spectrum.into_iter()
    }
}

impl std::ops::Index<usize> for AnnotatedSpectrum {
    type Output = AnnotatedPeak;
    fn index(&self, index: usize) -> &Self::Output {
        &self.spectrum[index]
    }
}

impl PeakSpectrum for AnnotatedSpectrum {
    type PeakType = AnnotatedPeak;
    type Iter<'a> = std::slice::Iter<'a, Self::PeakType>;

    /// Return the slice of peaks that have experimental mz values within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[AnnotatedPeak] {
        let left_idx = match self
            .spectrum
            .binary_search_by(|a| a.experimental_mz.value.total_cmp(&low.value))
        {
            Ok(idx) | Err(idx) => {
                let mut idx = idx.saturating_sub(1);
                while idx > 0
                    && self.spectrum[idx]
                        .experimental_mz
                        .value
                        .total_cmp(&low.value)
                        != Ordering::Less
                {
                    idx -= 1;
                }
                idx
            }
        };

        let right_idx = match self.spectrum[left_idx..]
            .binary_search_by(|a| a.experimental_mz.value.total_cmp(&high.value))
        {
            Ok(idx) | Err(idx) => {
                let mut idx = idx + left_idx;
                while idx < self.spectrum.len()
                    && self.spectrum[idx]
                        .experimental_mz
                        .value
                        .total_cmp(&high.value)
                        != Ordering::Greater
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

/// An annotated peak
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct AnnotatedPeak {
    /// The experimental mz
    pub experimental_mz: MassOverCharge,
    /// The experimental intensity
    pub intensity: OrderedFloat<f64>,
    /// The annotation, if present
    pub annotation: Vec<Fragment>, // Could become Vec<(Fragment, Vec<MatchedIsotopeDistribution>)> when isotope matching is finally in place
    /// Any annotation as isotope from a given fragment
    pub isotope_annotation: Vec<(usize, usize)>,
}

impl AnnotatedPeak {
    /// Make a new annotated peak with the given annotation
    pub fn new(peak: &RawPeak, annotation: Fragment) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            annotation: vec![annotation],
            isotope_annotation: Vec::new(),
        }
    }

    /// Make a new annotated peak if no annotation is possible
    pub fn background(peak: &RawPeak) -> Self {
        Self {
            experimental_mz: peak.mz,
            intensity: peak.intensity,
            annotation: Vec::new(),
            isotope_annotation: Vec::new(),
        }
    }
}

impl PartialOrd for AnnotatedPeak {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for AnnotatedPeak {
    /// Use `f64::total_cmp` on `self.mz`
    fn cmp(&self, other: &Self) -> Ordering {
        self.experimental_mz
            .value
            .total_cmp(&other.experimental_mz.value)
    }
}

impl PartialEq for AnnotatedPeak {
    /// Use `f64::total_cmp` on all fields to detect total equality
    fn eq(&self, other: &Self) -> bool {
        self.experimental_mz
            .value
            .total_cmp(&other.experimental_mz.value)
            == Ordering::Equal
            && self.intensity.total_cmp(&other.intensity) == Ordering::Equal
            && self.annotation == other.annotation
    }
}

impl Eq for AnnotatedPeak {}
