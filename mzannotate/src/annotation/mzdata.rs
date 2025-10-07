use crate::{annotation::AnnotatableSpectrum, mzspeclib::Analyte, spectrum::AnnotatedSpectrum};
use mzcore::{
    prelude::*,
    system::{MassOverCharge, ratio::ppm},
};
use mzdata::prelude::*;

// impl<S: SpectrumLike> AnnotatableSpectrum for S {
//     type Tolerance = mzcore::quantities::Tolerance<MassOverCharge>;

//     fn empty_annotated(&self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum {
//         AnnotatedSpectrum {
//             title: self.description().id.clone(),
//             num_scans: self.description().acquisition.scans.len() as u64,
//             rt: None,
//             charge: None,
//             mass: None,
//             peptide,
//             spectrum: match self.peaks() {
//                 RefPeakDataLevel::Missing
//                 | RefPeakDataLevel::RawData(_)
//                 | RefPeakDataLevel::Deconvoluted(_) => Vec::new(), // TODO: handle deconvoluted data better
//                 RefPeakDataLevel::Centroid(data) => data.iter().map(|p| p.clone().into()).collect(),
//             },
//         }
//     }

//     fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize> {
//         self.peaks().search(
//             query.value,
//             match tolerance {
//                 mzcore::quantities::Tolerance::Absolute(mz) => Tolerance::Da(mz.value),
//                 mzcore::quantities::Tolerance::Relative(ratio) => {
//                     Tolerance::PPM(ratio.get::<ppm>())
//                 }
//             },
//         )
//     }
// }

impl AnnotatableSpectrum for mzdata::Spectrum {
    type Tolerance = mzcore::quantities::Tolerance<MassOverCharge>;

    fn empty_annotated(&self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum {
        // TODO: Check if this can be rewritten to use `SpectrumLike`
        let mut spectrum = AnnotatedSpectrum::from(self.clone());
        for (index, peptidoform_ion) in peptide.into_peptidoform_ions().into_iter().enumerate() {
            let charge = peptidoform_ion
                .get_charge_carriers()
                .map(MolecularCharge::charge);
            spectrum.add_analyte(Analyte::new(
                index as u32,
                Some(peptidoform_ion),
                charge,
                Vec::new(),
            ));
        }

        spectrum
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
