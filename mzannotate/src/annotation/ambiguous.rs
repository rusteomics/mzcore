use mzcore::{
    chemistry::{AmbiguousLabel, MassOutputMode, MassOutputType},
    prelude::{AminoAcid, MassMode, SequencePosition},
    sequence::{AtLeast, CheckedAminoAcid, SemiAmbiguous, SequenceElement, UnAmbiguous},
    system::f32::Ratio,
};

use crate::{
    annotation::{Recovered, model::MatchingParameters},
    fragment::Fragment,
    mzspeclib::AnalyteTarget,
    spectrum::AnnotatedSpectrum,
};

impl<Mode: MassOutputMode> AnnotatedSpectrum<Mode> {
    /// Get the spectrum scores for this annotated spectrum.
    /// The returned tuple has the scores for all peptides combined as first item
    /// and as second item a vector with for each peptide its individual scores.
    pub fn ambigous_statistics(
        &self,
        fragments: &[Fragment<Mode>],
        parameters: &MatchingParameters,
        mass_mode: MassMode,
    ) -> AmbiguousStatistics {
        let fragments = fragments
            .iter()
            .filter(|f| f.mz(mass_mode).is_some_and(|mz| parameters.mz_range.contains(&mz)))
            .collect::<Vec<_>>();
        let total_intensity: f32 = self.peaks.iter().map(|p| p.intensity).sum();

        let theoretical = |label: &AmbiguousLabel| {
            fragments
                .iter()
                .filter(|f| f.formula.as_ref().is_some_and(|f| f.labels().contains(label)))
                .count() as u32
        };
        let annotated = |label: &AmbiguousLabel| {
            self.peaks.iter().fold((0_u32, 0.0), |acc, p| {
                let count = p
                    .annotations
                    .iter()
                    .filter(|f| f.formula.as_ref().is_some_and(|f| f.labels().contains(label)))
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
            deamidation: Vec::new(),
        };
        for (peptidoform_ion_index, analyte) in self.analytes.iter().enumerate() {
            let AnalyteTarget::PeptidoformIon(peptidoform_ion) = &analyte.target else {
                continue;
            };
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
                        AminoAcid::AsparticAcid
                            if aa.modifications.iter().any(|m| {
                                m.description().is_some_and(|d| {
                                    d.id()
                                        == mzcv::curie!(MOD:02117|L-aspartic acid isomer).accession
                                })
                            }) =>
                        {
                            let asp = statistics(&AmbiguousLabel::AminoAcid {
                                option: AminoAcid::AsparticAcid,
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });
                            let iso_asp = statistics(&AmbiguousLabel::IsoAsparticAcid {
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            });

                            result.deamidation.push(AmbiguousAsparticAcid {
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                                theoretical_fragments_asp: asp.0,
                                theoretical_fragments_iso_asp: iso_asp.0,
                                tic_asp: asp.1,
                                tic_iso_asp: iso_asp.1,
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
    /// All ambiguous aspartic acids
    pub deamidation: Vec<AmbiguousAsparticAcid>,
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
    /// First option, the amino acid, the fraction of theoretical fragments found, and the fraction
    /// of TIC that is annotated
    pub optiona_a: (AminoAcid, Recovered<u32>, Recovered<f32>),
    /// Second option, the amino acid, the fraction of theoretical fragments found, and the
    /// fraction of TIC that is annotated
    pub optiona_b: (AminoAcid, Recovered<u32>, Recovered<f32>),
}

/// The statistics on an ambiguous aspartic acid and the support for both options
#[derive(Clone, Copy, Debug)]
pub struct AmbiguousAsparticAcid {
    /// What location in the sequence are we talking about
    pub sequence_index: usize,
    /// Peptidoform index
    pub peptidoform_index: usize,
    /// Peptidoform ion index
    pub peptidoform_ion_index: usize,
    /// The fraction of theoretical fragments found for Asp
    pub theoretical_fragments_asp: Recovered<u32>,
    /// The fraction of theoretical fragments found for IsoAsp
    pub theoretical_fragments_iso_asp: Recovered<u32>,
    /// The fraction of TIC that could be explained by these fragments for Asp
    pub tic_asp: Recovered<f32>,
    /// The fraction of TIC that could be explained by these fragments for IsoAsp
    pub tic_iso_asp: Recovered<f32>,
}

impl AmbiguousAsparticAcid {
    /// The ratio between the fragments (TIC IsoAsp / Asp), this is None if Asp is not found
    pub fn ratio(&self) -> Option<Ratio> {
        (self.tic_asp.found > 0.0).then(|| {
            Ratio::new::<mzcore::system::ratio::fraction>(
                self.tic_iso_asp.found / self.tic_asp.found,
            )
        })
    }

    /// Resolve this location, based on the ration threshold below which this is `D`, above this it
    /// is `D[L-aspartic acid isomer]` (it returns the base as that is assumed to already be in this
    /// state) and the upper threshold above which it is `N[deamidated L-asparagine]`. It retrieves
    /// the deamidation modification from the ontologies, if the modification does not exist
    /// regardless of if it is used in this instance or not this function returns None.
    ///
    /// # Panics
    /// It has debug assertions for the lower threshold being smaller or equal to the upper
    /// threshold, the lower threshold being 0.0 or more, and the upper threshold being 1.0 or
    /// below.
    #[allow(clippy::zero_prefixed_literal)] // Looks more recognisable that way
    pub fn resolve<Complexity>(
        &self,
        lower_threshold: f32,
        upper_threshold: f32,
        ontologies: &mzcore::ontology::Ontologies,
        base: &SequenceElement<Complexity>,
    ) -> Option<SequenceElement<Complexity>>
    where
        CheckedAminoAcid<Complexity>: From<CheckedAminoAcid<UnAmbiguous>>,
    {
        debug_assert!(lower_threshold <= upper_threshold);
        debug_assert!(upper_threshold <= 1.0);
        debug_assert!(lower_threshold >= 0.0);

        ontologies
            .get_by_index(
                mzcore::ontology::Ontology::Psimod,
                &mzcv::curie!(MOD:00684|deamidated L-asparagine).accession,
            )
            .map(|deamidation| {
                self.ratio()
                    .map(|r| r.get::<mzcore::system::ratio::fraction>())
                    .map_or_else(
                        || base.clone(),
                        |r| {
                            if r <= lower_threshold {
                                SequenceElement::new_with_modifications(
                                    CheckedAminoAcid::AsparticAcid.into(),
                                    base.modifications
                                        .iter()
                                        .filter(|m| {
                                            m.description().is_none_or(|d| {
                                                d.id()
                                                != mzcv::curie!(MOD:02117|L-aspartic acid isomer)
                                                    .accession
                                            })
                                        })
                                        .cloned()
                                        .collect(),
                                    base.ambiguous,
                                )
                            } else if r >= upper_threshold {
                                SequenceElement::new_with_modifications(
                                    CheckedAminoAcid::Asparagine.into(),
                                    base.modifications
                                        .iter()
                                        .filter(|m| {
                                            m.description().is_none_or(|d| {
                                                d.id()
                                                != mzcv::curie!(MOD:02117|L-aspartic acid isomer)
                                                    .accession
                                            })
                                        })
                                        .cloned()
                                        .chain(std::iter::once(deamidation.into()))
                                        .collect(),
                                    base.ambiguous,
                                )
                            } else {
                                base.clone()
                            }
                        },
                    )
            })
    }
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
    /// All options, with the location, the fraction of theoretical fragments found, and the
    /// fraction of TIC that is annotated
    pub options: Vec<(SequencePosition, Recovered<u32>, Recovered<f32>)>,
}
