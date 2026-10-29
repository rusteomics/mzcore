use mzcore::{
    chemistry::MassMode,
    quantities::Multi,
    sequence::{HasPeptidoform, Linear},
    system::Mass,
};
#[cfg(feature = "rayon")]
use rayon::prelude::*;

use crate::{
    AlignScoring, AlignType, Alignment,
    diagonal_array::DiagonalArray,
    mass_alignment::{align_cached, calculate_masses},
};

/// An alignment index to store a collection of sequences to do alignment against with one or more
/// other sequences. This improves performance because it can cache some intermediate results.
#[derive(Clone, Debug)]
pub struct AlignIndex<const STEPS: u16, Source> {
    pub(super) sequences: Vec<(Source, DiagonalArray<Multi<Mass>, STEPS>)>,
    distance_matrix: Option<(
        DiagonalArray<f64, { u16::MAX }>,
        AlignType,
        AlignScoring<'static>,
    )>,
    pub(super) mode: MassMode,
}

impl<const S: u16, P> Default for AlignIndex<S, P> {
    fn default() -> Self {
        Self {
            sequences: Vec::new(),
            distance_matrix: None,
            mode: MassMode::Monoisotopic,
        }
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<Linear> + Clone> AlignIndex<STEPS, Source> {
    /// Check if this index is empty (contains no sequences)
    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    /// Get the number of sequences in the index
    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    /// Create a new index for alignments
    pub fn new(sequences: impl IntoIterator<Item = Source>, mode: MassMode) -> Self {
        let sequences = sequences
            .into_iter()
            .map(|p| {
                let masses = calculate_masses::<STEPS>(p.cast_peptidoform(), mode);
                (p, masses)
            })
            .collect::<Vec<_>>();

        Self {
            sequences,
            distance_matrix: None,
            mode,
        }
    }

    /// Add a new sequence to the index
    pub fn add_sequence(&mut self, sequence: Source) {
        let masses = calculate_masses(sequence.cast_peptidoform(), self.mode);
        if let Some((matrix, align_type, scoring)) = self.distance_matrix.as_mut() {
            matrix.push_row(
                self.sequences
                    .iter()
                    .map(|s| {
                        align_cached(
                            s.0.cast_peptidoform(),
                            &s.1,
                            sequence.cast_peptidoform(),
                            &masses,
                            *scoring,
                            *align_type,
                        )
                        .distance()
                    })
                    .chain([0.0])
                    .collect(),
            );
        }
        self.sequences.push((sequence, masses));
    }

    /// Retains only the elements specified by the predicate. Works using the same logic as
    /// [`Vec::retain`]
    pub fn retain(&mut self, f: impl Fn(&Source) -> bool) {
        let mut delete = Vec::new();
        let mut index = 0;
        self.sequences.retain(|(s, _)| {
            let retain = f(s);
            if !retain {
                delete.push(index);
            }
            index += 1;
            retain
        });
        if let Some((matrix, ..)) = self.distance_matrix.as_mut() {
            matrix.delete_rows(&delete);
        }
    }

    /// Align a single peptidoform to the index
    pub fn align_one<'a, Other: HasPeptidoform<Linear> + Clone + 'a>(
        &'a self,
        other: Other,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ExactSizeIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences.iter().map(move |(template, template_masses)| {
            align_cached::<STEPS, Source, Other>(
                template.clone(),
                template_masses,
                other.clone(),
                &other_masses,
                scoring,
                align_type,
            )
        })
    }

    /// Align a single peptidoform to the index
    pub fn align_one_specific<'a, Other: HasPeptidoform<Linear> + Clone + 'a>(
        &'a self,
        other: Other,
        scoring: AlignScoring<'a>,
        align_type: impl Fn(&Source) -> AlignType,
    ) -> impl ExactSizeIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences.iter().map(move |(template, template_masses)| {
            align_cached::<STEPS, Source, Other>(
                template.clone(),
                template_masses,
                other.clone(),
                &other_masses,
                scoring,
                align_type(template),
            )
        })
    }

    /// Align a single peptidoform to the index but only for the source sequences that match the
    /// given filter.
    pub fn align_one_filtered<'a, Other: HasPeptidoform<Linear> + Clone + 'a>(
        &'a self,
        other: Other,
        filter: impl Fn(&Source) -> bool + 'a,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl Iterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences.iter().filter(move |s| filter(&s.0)).map(
            move |(template, template_masses)| {
                align_cached::<STEPS, Source, Other>(
                    template.clone(),
                    template_masses,
                    other.clone(),
                    &other_masses,
                    scoring,
                    align_type,
                )
            },
        )
    }

    /// Align multiple peptides to this index
    pub fn align<'a, Other: HasPeptidoform<Linear> + Clone + 'a>(
        &'a self,
        others: impl IntoIterator<Item = Other>,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl Iterator<Item = impl ExactSizeIterator<Item = Alignment<Source, Other>>> {
        others.into_iter().map(move |o| self.align_one(o, scoring, align_type))
    }

    /// Get matrix with the distances between all sequences in this index. It uses
    /// [`crate::Alignment::distance`] as metric.
    /// # Panics
    /// If more than [`u16::MAX`] sequences are given.
    pub fn distance_matrix(
        &mut self,
        scoring: AlignScoring<'static>,
        align_type: AlignType,
    ) -> &DiagonalArray<f64, { u16::MAX }> {
        debug_assert!(self.sequences.len() <= u16::MAX.into());
        if self.distance_matrix.is_none()
            || self
                .distance_matrix
                .as_ref()
                .is_some_and(|(_, t, s)| *t != align_type || *s != scoring)
        {
            let mut matrix = DiagonalArray::<f64, { u16::MAX }>::new(self.sequences.len());
            for (i, s) in self.sequences.iter().enumerate() {
                matrix[[i, i]] = 0.0;
                for (o, so) in self.sequences.iter().enumerate().take(i) {
                    matrix[[i, o]] = align_cached(
                        s.0.cast_peptidoform(),
                        &s.1,
                        so.0.cast_peptidoform(),
                        &so.1,
                        scoring,
                        align_type,
                    )
                    .distance();
                }
            }
            self.distance_matrix = Some((matrix, align_type, scoring));
        }
        &self.distance_matrix.as_ref().unwrap().0
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<Linear> + Clone + Sync + Send> Extend<Source>
    for AlignIndex<STEPS, Source>
{
    fn extend<T: IntoIterator<Item = Source>>(&mut self, iter: T) {
        for v in iter {
            self.add_sequence(v);
        }
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<Linear> + Clone + Sync + Send>
    AlignIndex<STEPS, Source>
{
    #[cfg(feature = "rayon")]
    /// Align a single peptidoform to the index in parallel
    pub fn par_align_one<'a, Other: HasPeptidoform<Linear> + Clone + 'a + Send + Sync>(
        &'a self,
        other: Other,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences.par_iter().map(move |(template, template_masses)| {
            align_cached::<STEPS, Source, Other>(
                template.clone(),
                template_masses,
                other.clone(),
                &other_masses,
                scoring,
                align_type,
            )
        })
    }

    #[cfg(feature = "rayon")]
    /// Align a single peptidoform to the index in parallel but only for the source sequences that
    /// match the given filter.
    pub fn par_align_one_filtered<'a, Other: HasPeptidoform<Linear> + Clone + 'a + Send + Sync>(
        &'a self,
        other: Other,
        filter: impl Fn(&Source) -> bool + Sync + Send + 'a,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences.par_iter().filter(move |s| filter(&s.0)).map(
            move |(template, template_masses)| {
                align_cached::<STEPS, Source, Other>(
                    template.clone(),
                    template_masses,
                    other.clone(),
                    &other_masses,
                    scoring,
                    align_type,
                )
            },
        )
    }

    /// Create a new index for alignments in parallel
    #[cfg(feature = "rayon")]
    pub fn par_new(sequences: impl IntoParallelIterator<Item = Source>, mode: MassMode) -> Self {
        let sequences = sequences
            .into_par_iter()
            .map(|p| {
                let masses = calculate_masses::<STEPS>(p.cast_peptidoform(), mode);
                (p, masses)
            })
            .collect::<Vec<_>>();

        Self {
            sequences,
            distance_matrix: None,
            mode,
        }
    }

    /// Get matrix with the distances between all sequences in this index. It uses
    /// [`crate::Alignment::distance`] as metric.
    /// # Panics
    /// If more than [`u16::MAX`] sequences are given.
    #[cfg(feature = "rayon")]
    pub fn par_distance_matrix(
        &mut self,
        scoring: AlignScoring<'static>,
        align_type: AlignType,
    ) -> &DiagonalArray<f64, { u16::MAX }> {
        debug_assert!(self.sequences.len() <= u16::MAX.into());
        if self.distance_matrix.is_none()
            || self
                .distance_matrix
                .as_ref()
                .is_some_and(|(_, t, s)| *t != align_type || *s != scoring)
        {
            let mut matrix = DiagonalArray::<f64, { u16::MAX }>::new(self.sequences.len());
            for results in (&self)
                .sequences
                .par_iter()
                .enumerate()
                .flat_map(|(i, s)| {
                    (&self)
                        .sequences
                        .par_iter()
                        .enumerate()
                        .take(i)
                        .map(move |(o, so)| (i, s, o, so))
                })
                .map(|(i, s, o, so)| {
                    (
                        i,
                        o,
                        align_cached(
                            s.0.cast_peptidoform(),
                            &s.1,
                            so.0.cast_peptidoform(),
                            &so.1,
                            scoring,
                            align_type,
                        )
                        .distance(),
                    )
                })
                .collect_vec_list()
            {
                for (i, o, d) in results {
                    matrix[[i, o]] = d;
                }
            }
            for i in 0..self.sequences.len() {
                matrix[[i, i]] = 0.0;
            }
            self.distance_matrix = Some((matrix, align_type, scoring));
        }
        &self.distance_matrix.as_ref().unwrap().0
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<Linear> + Clone + Sync> AlignIndex<STEPS, Source> {
    #[cfg(feature = "rayon")]
    /// Align multiple peptides to this index in parallel
    pub fn par_align<'a, Other: HasPeptidoform<Linear> + Clone + 'a + Sync + Send>(
        &'a self,
        others: impl IntoParallelIterator<Item = Other>,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = impl ExactSizeIterator<Item = Alignment<Source, Other>>> {
        others.into_par_iter().map(move |o| self.align_one(o, scoring, align_type))
    }
}
