use mzdata::{prelude::*, spectrum::RefPeakDataLevel};

use crate::{
    annotation::{AnnotatableSpectrum, AnnotatedPeak, AnnotatedSpectrum},
    spectrum::RawPeak,
};
use mzcore::{
    sequence::CompoundPeptidoformIon,
    system::{MassOverCharge, ratio::ppm},
};

impl<S: SpectrumLike> AnnotatableSpectrum for S {
    type Tolerance = mzcore::quantities::Tolerance<MassOverCharge>;

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
                            mz: MassOverCharge::new::<mzcore::system::thomson>(p.mz),
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
                RefPeakDataLevel::Deconvoluted(data) => data
                    .iter()
                    .map(|p| {
                        AnnotatedPeak::background(&RawPeak {
                            mz: MassOverCharge::new::<mzcore::system::thomson>(p.neutral_mass), // TODO: This is M (not MH+) which is not very well supported in the current matching
                            intensity: ordered_float::OrderedFloat(f64::from(p.intensity)),
                        })
                    })
                    .collect(),
            },
        }
    }

    fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize> {
        self.peaks().search(
            query.value,
            match tolerance {
                mzcore::quantities::Tolerance::Absolute(mz) => Tolerance::Da(mz.value),
                mzcore::quantities::Tolerance::Relative(ratio) => {
                    Tolerance::PPM(ratio.get::<ppm>())
                }
            },
        )
    }
}

// TODO: fix somehow
// impl From<mzcore::quantities::Tolerance<MassOverCharge>> for Tolerance {
//     fn from(value: mzcore::quantities::Tolerance<MassOverCharge>) -> Self {
//         match value {
//             mzcore::quantities::Tolerance::Absolute(value) => {
//                 Self::Da(value.get::<mzcore::system::thomson>()) // This is in thomson (verified with crate author)
//             }
//             mzcore::quantities::Tolerance::Relative(value) => {
//                 Self::PPM(value.get::<mzcore::system::ratio::ppm>())
//             }
//         }
//     }
// }
