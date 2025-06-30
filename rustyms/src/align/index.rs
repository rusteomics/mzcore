#[cfg(feature = "rayon")]
use rayon::prelude::*;

use crate::{
    align::{
        AlignScoring, AlignType, Alignment,
        diagonal_array::DiagonalArray,
        mass_alignment::{align_cached, calculate_masses},
    },
    prelude::MassMode,
    quantities::Multi,
    sequence::{HasPeptidoform, SimpleLinear},
    system::Mass,
};

/// An alignment index to store a collection of sequences to do alignment against with one or more
/// other sequences. This improves performance because it can cache some intermediate results.
#[derive(Debug)]
pub struct AlignIndex<'lifetime, const STEPS: u16, Source> {
    sequences: Vec<(Source, DiagonalArray<Multi<Mass>>)>,
    mode: MassMode,
    scoring: AlignScoring<'lifetime>,
    align_type: AlignType,
}

impl<'lifetime, const STEPS: u16, Source: HasPeptidoform<SimpleLinear> + Clone>
    AlignIndex<'lifetime, STEPS, Source>
{
    /// Create a new index for alignments
    pub fn new(
        sequences: impl IntoIterator<Item = Source>,
        mode: MassMode,
        scoring: AlignScoring<'lifetime>,
        align_type: AlignType,
    ) -> Self {
        let sequences = sequences
            .into_iter()
            .map(|p| {
                let masses = calculate_masses::<STEPS>(p.cast_peptidoform(), mode);
                (p, masses)
            })
            .collect::<Vec<_>>();

        Self {
            sequences,
            mode,
            scoring,
            align_type,
        }
    }

    /// Align a single peptidoform to the index
    pub fn align_one<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a>(
        &'a self,
        other: Other,
    ) -> impl Iterator<Item = Alignment<Source, Other>> {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.mode);

        self.sequences
            .iter()
            .map(move |(template, template_masses)| {
                align_cached::<STEPS, Source, Other>(
                    template.clone(),
                    template_masses,
                    other.clone(),
                    &other_masses,
                    self.scoring,
                    self.align_type,
                )
            })
    }

    /// Align multiple peptides to this index
    pub fn align<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a>(
        &'a self,
        others: impl IntoIterator<Item = Other>,
    ) -> impl Iterator<Item = impl Iterator<Item = Alignment<Source, Other>>> {
        others.into_iter().map(|o| self.align_one(o))
    }
}

impl<const STEPS: u16, Source: HasPeptidoform<SimpleLinear> + Clone + Sync>
    AlignIndex<'_, STEPS, Source>
{
    #[cfg(feature = "rayon")]
    /// Align multiple peptides to this index in parallel
    pub fn par_align<'a, Other: HasPeptidoform<SimpleLinear> + Clone + 'a + Sync + Send>(
        &'a self,
        others: impl IntoParallelIterator<Item = Other>,
    ) -> impl ParallelIterator<Item = impl Iterator<Item = Alignment<Source, Other>>> {
        others.into_par_iter().map(|o| self.align_one(o))
    }
}
