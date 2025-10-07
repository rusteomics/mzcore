use std::fmt::Display;

use mzcore::system::MassOverCharge;
use mzdata::mzpeaks::prelude::*;
use serde::{Deserialize, Serialize};

/// An annotated peak. So a peak that contains some interpretation(s).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotatedPeak<Annotation> {
    /// The mass over charge of the peaks
    pub mz: MassOverCharge,
    /// The intensity of a peak
    pub intensity: f32,
    /// The index of a peak
    pub index: mzdata::mzpeaks::IndexType,
    /// The annotations or interpretations
    pub annotations: Vec<Annotation>,
    /// The aggregations for this peak
    pub aggregations: Vec<String>,
}

impl<A> Default for AnnotatedPeak<A> {
    fn default() -> Self {
        Self {
            mz: MassOverCharge::default(),
            intensity: f32::default(),
            index: mzdata::mzpeaks::IndexType::default(),
            annotations: Vec::new(),
            aggregations: Vec::new(),
        }
    }
}

impl<Annotation> AnnotatedPeak<Annotation> {
    /// Create a new peak
    pub const fn new(
        mz: MassOverCharge,
        intensity: f32,
        index: mzdata::mzpeaks::IndexType,
        annotations: Vec<Annotation>,
        aggregations: Vec<String>,
    ) -> Self {
        Self {
            mz,
            intensity,
            index,
            annotations,
            aggregations,
        }
    }
}

impl<B> AnnotatedPeak<B> {
    /// Convert an annotated peak into a different annotation type using `From`.
    pub fn from<A>(value: AnnotatedPeak<A>) -> Self
    where
        B: From<A>,
    {
        Self {
            mz: value.mz,
            intensity: value.intensity,
            index: value.index,
            annotations: value.annotations.into_iter().map(From::from).collect(),
            aggregations: value.aggregations,
        }
    }
}

impl<B> AnnotatedPeak<B> {
    /// Convert an annotated peak into a different annotation type using `TryFrom`.
    /// # Errors
    /// If the `TryFrom` failed.
    pub fn try_from<A>(value: AnnotatedPeak<A>) -> Result<Self, B::Error>
    where
        B: TryFrom<A>,
    {
        Ok(Self {
            mz: value.mz,
            intensity: value.intensity,
            index: value.index,
            annotations: value
                .annotations
                .into_iter()
                .map(B::try_from)
                .collect::<Result<Vec<_>, _>>()?,
            aggregations: value.aggregations,
        })
    }
}

impl<A, T: CentroidLike> PartialEq<T> for AnnotatedPeak<A> {
    #[inline]
    fn eq(&self, other: &T) -> bool {
        (self.mz.value - other.coordinate()).abs() < 1e-3
            && (self.intensity - other.intensity()).abs() < 1e-3
    }
}

impl<A> Eq for AnnotatedPeak<A> {}

impl<A, T: CentroidLike> PartialOrd<T> for AnnotatedPeak<A> {
    #[inline]
    fn partial_cmp(&self, other: &T) -> Option<std::cmp::Ordering> {
        Some(
            self.mz
                .value
                .total_cmp(&other.coordinate())
                .then_with(|| self.intensity.total_cmp(&other.intensity())),
        )
    }
}

impl<A> Ord for AnnotatedPeak<A> {
    #[inline]
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.mz
            .value
            .total_cmp(&other.mz.value)
            .then_with(|| self.intensity.total_cmp(&other.intensity))
    }
}

impl<A> CoordinateLike<mzdata::mzpeaks::MZ> for AnnotatedPeak<A> {
    #[inline]
    fn coordinate(&self) -> f64 {
        self.mz.value
    }
}

impl<A> CoordinateLikeMut<mzdata::mzpeaks::MZ> for AnnotatedPeak<A> {
    #[inline]
    fn coordinate_mut(&mut self) -> &mut f64 {
        &mut self.mz.value
    }
}

impl<A> IntensityMeasurement for AnnotatedPeak<A> {
    #[inline]
    fn intensity(&self) -> f32 {
        self.intensity
    }
}

impl<A> IntensityMeasurementMut for AnnotatedPeak<A> {
    #[inline]
    fn intensity_mut(&mut self) -> &mut f32 {
        &mut self.intensity
    }
}

impl<A> IndexedCoordinate<mzdata::mzpeaks::MZ> for AnnotatedPeak<A> {
    #[inline]
    fn get_index(&self) -> mzdata::mzpeaks::IndexType {
        self.index
    }

    #[inline]
    fn set_index(&mut self, index: mzdata::mzpeaks::IndexType) {
        self.index = index;
    }
}

impl<A> From<AnnotatedPeak<A>> for mzdata::mzpeaks::CentroidPeak {
    fn from(peak: AnnotatedPeak<A>) -> Self {
        peak.as_centroid()
    }
}

impl<A> From<AnnotatedPeak<A>> for mzdata::mzpeaks::peak::MZPoint {
    fn from(peak: AnnotatedPeak<A>) -> Self {
        Self {
            mz: peak.coordinate(),
            intensity: peak.intensity(),
        }
    }
}

impl<A> From<mzdata::mzpeaks::CentroidPeak> for AnnotatedPeak<A> {
    fn from(peak: mzdata::mzpeaks::CentroidPeak) -> Self {
        let mut inst = Self {
            mz: MassOverCharge::new::<mzcore::system::thomson>(peak.coordinate()),
            intensity: peak.intensity(),
            ..Self::default()
        };
        inst.set_index(peak.index);
        inst
    }
}

impl<A: Display> Display for AnnotatedPeak<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}", self.mz.value, self.intensity)?;

        if !self.annotations.is_empty() {
            f.write_str("\t")?;
            for (i, a) in self.annotations.iter().enumerate() {
                if i == 0 {
                    write!(f, "{a}")?;
                } else {
                    write!(f, ",{a}")?;
                }
            }
        }
        if !self.aggregations.is_empty() {
            f.write_str("\t")?;
            write!(f, "{}", self.aggregations.join(","))?;
        }
        Ok(())
    }
}
