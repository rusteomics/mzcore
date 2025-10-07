use mzcore::{
    chemistry::AmbiguousLabel,
    prelude::{AminoAcid, MassMode, SequencePosition},
};

use crate::{
    annotation::{AnnotatedSpectrum, Recovered, model::MatchingParameters},
    fragment::Fragment,
};

impl AnnotatedSpectrum {
    /// Get the spectrum scores for this annotated spectrum.
    /// The returned tuple has the scores for all peptides combined as first item
    /// and as second item a vector with for each peptide its individual scores.
    pub fn ambigous_statistics(
        &self,
        fragments: &[Fragment],
        parameters: &MatchingParameters,
        mass_mode: MassMode,
    ) -> AmbiguousStatistics {
        let fragments = fragments
            .iter()
            .filter(|f| {
                f.mz(mass_mode)
                    .is_some_and(|mz| parameters.mz_range.contains(&mz))
            })
            .collect::<Vec<_>>();
        let total_intensity: f32 = self.spectrum.iter().map(|p| p.intensity).sum();

        let theoretical = |label: &AmbiguousLabel| {
            fragments
                .iter()
                .filter(|f| {
                    f.formula
                        .as_ref()
                        .is_some_and(|f| f.labels().contains(label))
                })
                .count() as u32
        };
        let annotated = |label: &AmbiguousLabel| {
            self.spectrum.iter().fold((0_u32, 0.0), |acc, p| {
                let count = p
                    .annotations
                    .iter()
                    .filter(|f| {
                        f.formula
                            .as_ref()
                            .is_some_and(|f| f.labels().contains(label))
                    })
                    .count();
                if count == 0 {
                    acc
                } else {
                    (acc.0 + count as u32, acc.1 + p.intensity)
                }
            })
        };
        let statistics = |label: &AmbiguousLabel| {
            let th_count = theoretical(label);
            let (an_count, an_intensity) = annotated(label);
            (
                Recovered {
                    found: an_count,
                    total: th_count,
                },
                Recovered {
                    found: an_intensity,
                    total: total_intensity,
                },
            )
        };

        let mut result = AmbiguousStatistics {
            aminoacids: Vec::new(),
            modifications: Vec::new(),
        };
        for (peptidoform_ion_index, peptidoform_ion) in
            self.peptide.peptidoform_ions().iter().enumerate()
        {
            for (peptidoform_index, peptidoform) in
                peptidoform_ion.peptidoforms().iter().enumerate()
            {
                for (sequence_index, aa) in peptidoform.sequence().iter().enumerate() {
                    match aa.aminoacid.aminoacid() {
                        AminoAcid::AmbiguousLeucine => {
                            let isoleucine = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::Isoleucine,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            let leucine = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::Leucine,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            result.aminoacids.push(AmbiguousAminoAcid {
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                                optiona_a: (AminoAcid::Isoleucine, isoleucine.0, isoleucine.1),
                                optiona_b: (AminoAcid::Leucine, leucine.0, leucine.1),
                            });
                        }
                        AminoAcid::AmbiguousAsparagine => {
                            let aspartic_acid = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::AsparticAcid,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            let asparagine = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::Asparagine,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            result.aminoacids.push(AmbiguousAminoAcid {
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                                optiona_a: (
                                    AminoAcid::AsparticAcid,
                                    aspartic_acid.0,
                                    aspartic_acid.1,
                                ),
                                optiona_b: (AminoAcid::Asparagine, asparagine.0, asparagine.1),
                            });
                        }
                        AminoAcid::AmbiguousGlutamine => {
                            let glutamic_acid = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::GlutamicAcid,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            let glutamine = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::Glutamine,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            result.aminoacids.push(AmbiguousAminoAcid {
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                                optiona_a: (
                                    AminoAcid::GlutamicAcid,
                                    glutamic_acid.0,
                                    glutamic_acid.1,
                                ),
                                optiona_b: (AminoAcid::Glutamine, glutamine.0, glutamine.1),
                            });
                        }
                        _ => (),
                    }
                }

                for (id, locations) in peptidoform.get_ambiguous_modifications().iter().enumerate()
                {
                    result.modifications.push(AmbiguousModification {
                        id,
                        peptidoform_index,
                        peptidoform_ion_index,
                        options: locations
                            .iter()
                            .copied()
                            .map(|sequence_index| {
                                let stats = statistics(&AmbiguousLabel::Modification {
                                    id,
                                    sequence_index,
                                    peptidoform_index,
                                    peptidoform_ion_index,
                                });
                                (sequence_index, stats.0, stats.1)
                            })
                            .collect(),
                    });
                }
            }
        }
        result
    }
}

/// All statistics for ambiguous parts of the peptidoform definition
#[derive(Clone, Debug)]
pub struct AmbiguousStatistics {
    /// All ambiguous amino acids
    pub aminoacids: Vec<AmbiguousAminoAcid>,
    /// All ambiguous modifications (modifications of unknown position)
    pub modifications: Vec<AmbiguousModification>,
}

/// The statistics on an ambiguous amino acid and the support for both options
#[derive(Clone, Copy, Debug)]
pub struct AmbiguousAminoAcid {
    /// What location in the sequence are we talking about
    pub sequence_index: usize,
    /// Peptidoform index
    pub peptidoform_index: usize,
    /// Peptidoform ion index
    pub peptidoform_ion_index: usize,
    /// First option, the amino acid, the fraction of theoretical fragments found, and the fraction of TIC that is annotated
    pub optiona_a: (AminoAcid, Recovered<u32>, Recovered<f32>),
    /// Second option, the amino acid, the fraction of theoretical fragments found, and the fraction of TIC that is annotated
    pub optiona_b: (AminoAcid, Recovered<u32>, Recovered<f32>),
}

/// The statistics on an ambiguous modification and the support for each of the possible locations
#[derive(Clone, Debug)]
pub struct AmbiguousModification {
    /// Which ambiguous modification
    pub id: usize,
    /// Peptidoform index
    pub peptidoform_index: usize,
    /// Peptidoform ion index
    pub peptidoform_ion_index: usize,
    /// All options, with the location, the fraction of theoretical fragments found, and the fraction of TIC that is annotated
    pub options: Vec<(SequencePosition, Recovered<u32>, Recovered<f32>)>,
}
