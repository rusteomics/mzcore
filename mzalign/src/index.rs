#[cfg(feature = "rayon")]
use rayon::prelude::*;

use crate::{
    AlignScoring, AlignType, Alignment,
    diagonal_array::DiagonalArray,
    mass_alignment::{align_cached, calculate_masses},
};
use mzcore::{
    chemistry::MassMode,
    quantities::Multi,
    sequence::{HasPeptidoform, SimpleLinear},
    system::Mass,
};

/// An alignment index to store a collection of sequences to do alignment against with one or more
/// other sequences. This improves performance because it can cache some intermediate results.
#[derive(Debug)]
pub struct AlignIndex<const STEPS: u16, Source> {
    pub(super) sequences: Vec<(Source, DiagonalArray<Multi<Mass>, STEPS>)>,
    pub(super) mode: MassMode,
}

impl<const STEPS: u16, Source: HasPeptidoform<SimpleLinear> + Clone> AlignIndex<STEPS, Source> {
    /// Create a new index for alignments
    pub fn new(sequences: impl IntoIterator<Item = Source>, mode: MassMode) -> Self {
        let sequences = sequences
            .into_iter()
            .map(|p| {
                let masses = calculate_masses::<STEPS>(p.cast_peptidoform(), mode);
                (p, masses)
            })
            .collect::<Vec<_>>();

        Self { sequences, mode }
    }

    /// Align a single peptidoform to the index
    pub fn align_one<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a>(
        &'a self,
        other: Other,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ExactSizeIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences
            .iter()
            .map(move |(template, template_masses)| {
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

    /// Align a single peptidoform to the index but only for the source sequences that match the given filter.
    pub fn align_one_filtered<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a>(
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
    pub fn align<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a>(
        &'a self,
        others: impl IntoIterator<Item = Other>,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl Iterator<Item = impl ExactSizeIterator<Item = Alignment<Source, Other>>> {
        others
            .into_iter()
            .map(move |o| self.align_one(o, scoring, align_type))
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<SimpleLinear> + Clone + Sync + Send>
    AlignIndex<STEPS, Source>
{
    #[cfg(feature = "rayon")]
    /// Align a single peptidoform to the index in parallel
    pub fn par_align_one<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a + Send + Sync>(
        &'a self,
        other: Other,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences
            .par_iter()
            .map(move |(template, template_masses)| {
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
    /// Align a single peptidoform to the index in parallel but only for the source sequences that match the given filter.
    pub fn par_align_one_filtered<
        'a,
        Other: HasPeptidoform<SimpleLinear> + Clone + 'a + Send + Sync,
    >(
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
}

impl<const STEPS: u16, Source: HasPeptidoform<SimpleLinear> + Clone + Sync>
    AlignIndex<STEPS, Source>
{
    #[cfg(feature = "rayon")]
    /// Align multiple peptides to this index in parallel
    pub fn par_align<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a + Sync + Send>(
        &'a self,
        others: impl IntoParallelIterator<Item = Other>,
        scoring: AlignScoring<'a>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = impl ExactSizeIterator<Item = Alignment<Source, Other>>> {
        others
            .into_par_iter()
            .map(move |o| self.align_one(o, scoring, align_type))
    }
}
