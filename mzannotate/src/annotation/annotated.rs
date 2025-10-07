//! Annotated spectra

use std::cmp::Ordering;

use mzcore::{
    prelude::CompoundPeptidoformIon,
    system::{Mass, MassOverCharge, Time, isize::Charge},
};
use serde::{Deserialize, Serialize};

use crate::{
    fragment::Fragment,
    spectrum::{AnnotatedPeak, PeakSpectrum},
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
    pub(super) spectrum: Vec<AnnotatedPeak<Fragment>>,
}

impl Extend<AnnotatedPeak<Fragment>> for AnnotatedSpectrum {
    fn extend<T: IntoIterator<Item = AnnotatedPeak<Fragment>>>(&mut self, iter: T) {
        self.spectrum.extend(iter);
        self.spectrum.sort_unstable();
    }
}

impl IntoIterator for AnnotatedSpectrum {
    type Item = AnnotatedPeak<Fragment>;
    type IntoIter = std::vec::IntoIter<AnnotatedPeak<Fragment>>;
    fn into_iter(self) -> Self::IntoIter {
        self.spectrum.into_iter()
    }
}

impl std::ops::Index<usize> for AnnotatedSpectrum {
    type Output = AnnotatedPeak<Fragment>;
    fn index(&self, index: usize) -> &Self::Output {
        &self.spectrum[index]
    }
}

impl PeakSpectrum for AnnotatedSpectrum {
    type PeakType = AnnotatedPeak<Fragment>;
    type Iter<'a> = std::slice::Iter<'a, Self::PeakType>;

    /// Return the slice of peaks that have experimental mz values within the given tolerance bounds.
    fn binary_search(
        &self,
        low: MassOverCharge,
        high: MassOverCharge,
    ) -> &[AnnotatedPeak<Fragment>] {
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
