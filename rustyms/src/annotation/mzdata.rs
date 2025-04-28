use mzdata::{prelude::*, spectrum::RefPeakDataLevel};

use crate::{
    annotation::{AnnotatableSpectrum, AnnotatedPeak, AnnotatedSpectrum},
    sequence::CompoundPeptidoformIon,
    spectrum::RawPeak,
    system::MassOverCharge,
};

impl<S: SpectrumLike> AnnotatableSpectrum for S {
    type Tolerance = Tolerance;

    fn empty_annotated(&self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum {
        AnnotatedSpectrum {
            title: self.description().id.clone(),
            num_scans: self.description().acquisition.scans.len() as u64,
            rt: None,
            charge: None,
            mass: None,
            peptide,
            spectrum: match self.peaks() {
                RefPeakDataLevel::Missing | RefPeakDataLevel::RawData(_) => Vec::new(),
                RefPeakDataLevel::Centroid(data) => data
                    .iter()
                    .map(|p| {
                        AnnotatedPeak::background(&RawPeak {
                            mz: MassOverCharge::new::<crate::system::mz>(p.mz),
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
                RefPeakDataLevel::Deconvoluted(data) => data
                    .iter()
                    .map(|p| {
                        AnnotatedPeak::background(&RawPeak {
                            mz: MassOverCharge::new::<crate::system::mz>(p.neutral_mass), // TODO: This is M (not MH+) which is not very well supported in the current matching
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
            },
        }
    }

    fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize> {
        self.peaks().search(query.value, tolerance)
    }
}

impl From<crate::quantities::Tolerance<MassOverCharge>> for Tolerance {
    fn from(value: crate::quantities::Tolerance<MassOverCharge>) -> Self {
        match value {
            crate::quantities::Tolerance::Absolute(value) => {
                Self::Da(value.get::<crate::system::mz>()) // This is in Thompson (verified with crate author)
            }
            crate::quantities::Tolerance::Relative(value) => {
                Self::PPM(value.get::<crate::system::ratio::ppm>())
            }
        }
    }
}
