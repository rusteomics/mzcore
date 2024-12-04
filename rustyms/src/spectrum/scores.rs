//! Scoring of annotated spectra

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    fragment::{Fragment, FragmentKind},
    peptide::UnAmbiguous,
    AnnotatedSpectrum, LinearPeptide, MassMode, Model,
};

impl AnnotatedSpectrum {
    /// Get the spectrum scores for this annotated spectrum.
    /// The returned tuple has the scores for all peptides combined as first item
    /// and as second item a vector with for each peptide its individual scores.
    pub fn scores(
        &self,
        fragments: &[Fragment],
        model: &Model,
        mass_mode: MassMode,
    ) -> (Scores, Vec<Vec<Scores>>) {
        let fragments = fragments
            .iter()
            .filter(|f| {
                f.mz(mass_mode)
                    .is_some_and(|mz| model.mz_range.contains(&mz))
            })
            .collect_vec();
        let total_intensity: f64 = self.spectrum.iter().map(|p| *p.intensity).sum();
        let individual_peptides = self
            .peptide
            .peptidoforms()
            .iter()
            .enumerate()
            .map(|(peptidoform_index, peptidoform)| {
                peptidoform
                    .peptides()
                    .iter()
                    .enumerate()
                    .map(|(peptide_index, peptide)| {
                        let (recovered_fragments, peaks, intensity_annotated) = self
                            .filtered_base_score(
                                &fragments,
                                Some(peptidoform_index),
                                Some(peptide_index),
                                None,
                            );
                        let (positions, expected_positions) = self.score_positions(
                            &fragments,
                            peptidoform_index,
                            peptide_index,
                            None,
                        );
                        Scores {
                            score: Score::Position {
                                fragments: recovered_fragments,
                                peaks,
                                intensity: Recovered::new(intensity_annotated, total_intensity),
                                theoretical_positions: Recovered::new(
                                    positions,
                                    peptide.len() as u32,
                                ),
                                expected_positions: Recovered::new(positions, expected_positions),
                            },
                            ions: self.score_individual_ions(
                                &fragments,
                                Some((peptidoform_index, peptide_index, peptide)),
                                total_intensity,
                            ),
                        }
                    })
                    .collect()
            })
            .collect();
        // Get the statistics for the combined peptides
        let (recovered_fragments, peaks, intensity_annotated) =
            self.filtered_base_score(&fragments, None, None, None);
        let unique_formulas = self.score_unique_formulas(&fragments, None, None);
        (
            Scores {
                score: Score::UniqueFormulas {
                    fragments: recovered_fragments,
                    peaks,
                    intensity: Recovered::new(intensity_annotated, total_intensity),
                    unique_formulas,
                },
                ions: self.score_individual_ions::<UnAmbiguous>(&fragments, None, total_intensity),
            },
            individual_peptides,
        )
    }

    /// Get the base score of this spectrum
    /// (Fragments, peaks, intensity)
    fn filtered_base_score(
        &self,
        fragments: &[&Fragment],
        peptidoform_index: Option<usize>,
        peptide_index: Option<usize>,
        ion: Option<FragmentKind>,
    ) -> (Recovered<u32>, Recovered<u32>, f64) {
        let (peaks_annotated, fragments_found, intensity_annotated) = self
            .spectrum
            .iter()
            .filter_map(|p| {
                let number = p
                    .annotation
                    .iter()
                    .filter(|a| {
                        peptidoform_index.map_or(true, |i| a.peptidoform_index == Some(i))
                            && peptide_index.map_or(true, |i| a.peptide_index == Some(i))
                            && ion.map_or(true, |kind| a.ion.kind() == kind)
                    })
                    .count() as u32;
                if number == 0 {
                    None
                } else {
                    Some((number, *p.intensity))
                }
            })
            .fold((0u32, 0u32, 0.0), |(n, f, intensity), p| {
                (n + 1, f + p.0, intensity + p.1)
            });
        let total_fragments = fragments
            .iter()
            .filter(|f| {
                peptidoform_index.map_or(true, |i| f.peptidoform_index == Some(i))
                    && peptide_index.map_or(true, |i| f.peptide_index == Some(i))
                    && ion.map_or(true, |kind| f.ion.kind() == kind)
            })
            .count() as u32;
        (
            Recovered::new(fragments_found, total_fragments),
            Recovered::new(peaks_annotated, self.spectrum.len() as u32),
            intensity_annotated,
        )
    }

    /// Get the total number of positions covered
    fn score_positions(
        &self,
        fragments: &[&Fragment],
        peptidoform_index: usize,
        peptide_index: usize,
        ion: Option<FragmentKind>,
    ) -> (u32, u32) {
        (
            self.spectrum
                .iter()
                .flat_map(|p| {
                    p.annotation
                        .iter()
                        .filter(|a| {
                            a.peptidoform_index == Some(peptidoform_index)
                                && a.peptide_index == Some(peptide_index)
                                && ion.map_or(true, |kind| a.ion.kind() == kind)
                        })
                        .filter_map(|a| a.ion.position())
                })
                .map(|pos| pos.sequence_index)
                .unique()
                .count() as u32,
            fragments
                .iter()
                .filter(|f| {
                    f.peptidoform_index == Some(peptidoform_index)
                        && f.peptide_index == Some(peptide_index)
                        && ion.map_or(true, |i| i == f.ion.kind())
                })
                .filter_map(|f| f.ion.position().map(|p| p.sequence_index))
                .unique()
                .count() as u32,
        )
    }

    /// Get the amount of unique formulas recovered
    fn score_unique_formulas(
        &self,
        fragments: &[&Fragment],
        peptide_index: Option<usize>,
        ion: Option<FragmentKind>,
    ) -> Recovered<u32> {
        let num_annotated = self
            .spectrum
            .iter()
            .flat_map(|p| {
                p.annotation.iter().filter(|a| {
                    peptide_index.map_or(true, |i| a.peptide_index == Some(i))
                        && ion.map_or(true, |kind| a.ion.kind() == kind)
                })
            })
            .map(|f| f.formula.clone())
            .unique()
            .count() as u32;
        let total_fragments = fragments
            .iter()
            .filter(|f| {
                peptide_index.map_or(true, |i| f.peptide_index == Some(i))
                    && ion.map_or(true, |kind| f.ion.kind() == kind)
            })
            .map(|f| f.formula.clone())
            .unique()
            .count() as u32;
        Recovered::new(num_annotated, total_fragments)
    }

    /// Get the scores for the individual ion series
    fn score_individual_ions<T>(
        &self,
        fragments: &[&Fragment],
        peptide: Option<(usize, usize, &LinearPeptide<T>)>,
        total_intensity: f64,
    ) -> Vec<(FragmentKind, Score)> {
        [
            FragmentKind::a,
            FragmentKind::b,
            FragmentKind::c,
            FragmentKind::d,
            FragmentKind::v,
            FragmentKind::w,
            FragmentKind::x,
            FragmentKind::y,
            FragmentKind::z,
        ]
        .iter()
        .copied()
        .filter_map(|ion| {
            let (recovered_fragments, peaks, intensity_annotated) = self.filtered_base_score(
                fragments,
                peptide.as_ref().map(|p| p.0),
                peptide.as_ref().map(|p| p.1),
                Some(ion),
            );
            if let Some((peptidoform_index, peptide_index, peptide)) = peptide {
                if recovered_fragments.total > 0 {
                    let (positions, expected_positions) = self.score_positions(
                        fragments,
                        peptidoform_index,
                        peptide_index,
                        Some(ion),
                    );
                    Some((
                        ion,
                        Score::Position {
                            fragments: recovered_fragments,
                            peaks,
                            intensity: Recovered::new(intensity_annotated, total_intensity),
                            theoretical_positions: Recovered::new(positions, peptide.len() as u32),
                            expected_positions: Recovered::new(positions, expected_positions),
                        },
                    ))
                } else {
                    None
                }
            } else if recovered_fragments.total > 0 {
                let unique_formulas = self.score_unique_formulas(fragments, None, Some(ion));
                Some((
                    ion,
                    Score::UniqueFormulas {
                        fragments: recovered_fragments,
                        peaks,
                        intensity: Recovered::new(intensity_annotated, total_intensity),
                        unique_formulas,
                    },
                ))
            } else {
                None
            }
        })
        .chain(
            [
                FragmentKind::Y,
                FragmentKind::Oxonium,
                FragmentKind::immonium,
                FragmentKind::precursor_side_chain_loss,
                FragmentKind::diagnostic,
                FragmentKind::precursor,
            ]
            .iter()
            .copied()
            .filter_map(|ion| {
                let (recovered_fragments, peaks, intensity_annotated) = self.filtered_base_score(
                    fragments,
                    peptide.as_ref().map(|p| p.0),
                    peptide.as_ref().map(|p| p.1),
                    Some(ion),
                );
                if recovered_fragments.total > 0 {
                    let unique_formulas = self.score_unique_formulas(
                        fragments,
                        peptide.as_ref().map(|p| p.0),
                        Some(ion),
                    );
                    Some((
                        ion,
                        Score::UniqueFormulas {
                            fragments: recovered_fragments,
                            peaks,
                            intensity: Recovered::new(intensity_annotated, total_intensity),
                            unique_formulas,
                        },
                    ))
                } else {
                    None
                }
            }),
        )
        .collect()
    }
}

/// The scores for an annotated spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Scores {
    /// The scores, based on unique formulas for all peptides combined or based on positions for single peptides.
    pub score: Score,
    /// The scores per [`FragmentKind`], based on unique formulas for all peptides combined or any fragment kind that is not an ion series, or based on positions in the other case.
    pub ions: Vec<(FragmentKind, Score)>,
}

/// The scores for a single fragment series for a single peptide in an annotated spectrum
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub enum Score {
    /// A score for a something that has peptide position coverage
    Position {
        /// The fraction of the total fragments that could be annotated
        fragments: Recovered<u32>,
        /// The fraction of the total peaks that could be annotated
        peaks: Recovered<u32>,
        /// The fraction of the total intensity that could be annotated
        intensity: Recovered<f64>,
        /// The fraction of the total positions (all positions on the peptide) that has at least one fragment found
        theoretical_positions: Recovered<u32>,
        /// The fraction of the total positions (all positions with fragments) that has at least one fragment found
        expected_positions: Recovered<u32>,
    },
    /// A score for something that does not have position coverage, but instead is scored on the number of unique formulas
    UniqueFormulas {
        /// The fraction of the total fragments that could be annotated
        fragments: Recovered<u32>,
        /// The fraction of the total peaks that could be annotated
        peaks: Recovered<u32>,
        /// The fraction of the total intensity that could be annotated
        intensity: Recovered<f64>,
        /// The fraction of with unique formulas that has been found
        unique_formulas: Recovered<u32>,
    },
}
/// A single statistic that has a total number and a subset of that found
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
#[non_exhaustive]
pub struct Recovered<T> {
    /// The number actually found
    pub found: T,
    /// The total number
    pub total: T,
}

impl<T> Recovered<T> {
    /// Create a new recovered statistic
    fn new(found: impl Into<T>, total: impl Into<T>) -> Self {
        Self {
            found: found.into(),
            total: total.into(),
        }
    }
}

impl<T> Recovered<T>
where
    f64: From<T>,
    T: Copy,
{
    /// Get the recovered amount as fraction
    pub fn fraction(&self) -> f64 {
        f64::from(self.found) / f64::from(self.total)
    }
}
