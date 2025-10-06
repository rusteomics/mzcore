//! General trait for a spectrum

use std::iter::FusedIterator;

use mzcore::system::f64::MassOverCharge;

/// The trait for all spectra that contain peaks.
pub trait PeakSpectrum:
    Extend<Self::PeakType>
    + IntoIterator<Item = Self::PeakType>
    + std::ops::Index<usize, Output = Self::PeakType>
{
    /// The type of peaks this spectrum contains
    type PeakType;
    /// The type of spectrum iterator this spectrum generates
    type Iter<'a>: DoubleEndedIterator + ExactSizeIterator + FusedIterator
    where
        Self: 'a;
    /// Return the slice of peaks that is within the given tolerance bounds.
    fn binary_search(&self, low: MassOverCharge, high: MassOverCharge) -> &[Self::PeakType];
    /// Get the full spectrum
    fn spectrum(&self) -> Self::Iter<'_>;
    /// Add a single peak
    fn add_peak(&mut self, item: Self::PeakType);
}
