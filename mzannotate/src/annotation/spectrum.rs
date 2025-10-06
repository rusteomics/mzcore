use mzcore::{
    prelude::CompoundPeptidoformIon,
    quantities::{Tolerance, WithinTolerance},
    system::MassOverCharge,
};

use crate::{
    annotation::{AnnotatableSpectrum, AnnotatedPeak, AnnotatedSpectrum},
    spectrum::RawSpectrum,
};

impl AnnotatableSpectrum for RawSpectrum {
    type Tolerance = Tolerance<MassOverCharge>;

    fn empty_annotated(&self, peptide: CompoundPeptidoformIon) -> AnnotatedSpectrum {
        AnnotatedSpectrum {
            title: self.title.clone(),
            num_scans: self.num_scans,
            rt: self.rt,
            charge: self.charge,
            mass: self.mass,
            peptide,
            spectrum: self
                .spectrum
                .iter()
                .map(AnnotatedPeak::background)
                .collect(),
        }
    }

    fn search(&self, query: MassOverCharge, tolerance: Self::Tolerance) -> Option<usize> {
        let index = self
            .spectrum
            .binary_search_by(|p| p.mz.value.total_cmp(&query.value))
            .unwrap_or_else(|i| i);

        // Check index-1, index and index+1 (if existing) to find the one with the lowest ppm
        let mut closest = (0, f64::INFINITY);
        for i in if index == 0 { 0 } else { index - 1 }..=(index + 1).min(self.spectrum.len() - 1) {
            let ppm = self.spectrum[i].ppm(query).value;
            if ppm < closest.1 {
                closest = (i, ppm);
            }
        }

        tolerance
            .within(&self.spectrum[closest.0].mz, &query)
            .then_some(closest.0)
    }
}
