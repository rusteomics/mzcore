use crate::*;
use imgt::*;
use mzcore::sequence::{
    AnnotatedPeptide, AtMax, HasPeptidoform, Linear, Peptidoform, Region, SimpleLinear, UnAmbiguous,
};
use mzcv::CVIndex;
use std::collections::HashSet;

use itertools::Itertools;

/// A consecutive alignment, which align one sequence to multiple sequences.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConsecutiveAlignment<'lifetime, A> {
    /// All underlying alignments, per gene there is a vector containing all options for that gene.
    pub alignments: Vec<
        Vec<(
            Allele<'lifetime>,
            Alignment<&'lifetime Peptidoform<UnAmbiguous>, Peptidoform<A>>,
        )>,
    >,
}

impl<'lifetime, A> ConsecutiveAlignment<'lifetime, A> {
    /// Get the main alignment, the alignment taking the best alignment for each gene.
    pub fn main_alignment(
        &self,
    ) -> Vec<&(
        Allele<'lifetime>,
        Alignment<&'lifetime Peptidoform<UnAmbiguous>, Peptidoform<A>>,
    )> {
        self.alignments.iter().filter_map(|a| a.first()).collect()
    }
}

impl<A: AtMax<Linear>> ConsecutiveAlignment<'_, A> {
    /// Break up in the main alignment into the regions as annotated in the alleles.
    #[expect(clippy::missing_panics_doc)]
    pub fn regions(&self) -> Vec<(Peptidoform<A>, Region)> {
        let mut b_offset = 0;
        self.alignments
            .iter()
            .filter_map(|a| {
                a.first().map(|(allele, alignment)| {
                    let mut index_a = alignment.start_a;
                    let mut start_region_b = alignment.start_b;
                    let mut index_b = alignment.start_b;
                    let mut region = allele
                        .get_region(index_a)
                        .map_or(&Region::Framework(1), |(r, _)| r);
                    let mut ranges = Vec::new();

                    for step in &alignment.path {
                        let new_region = allele
                            .get_region(index_a + step.step_a as usize)
                            .map_or(&Region::Framework(1), |(r, _)| r);
                        if region != new_region {
                            ranges.push((b_offset + start_region_b, b_offset + index_b, region));
                            start_region_b = index_b + 1;
                            region = new_region;
                        }
                        index_a += step.step_a as usize;
                        index_b += step.step_b as usize;
                    }
                    b_offset += alignment.len_b() + alignment.start_b;

                    ranges
                })
            })
            .flatten()
            .chunk_by(|(_, _, r)| *r)
            .into_iter()
            .map(|(region, mut chunk)| {
                let first = chunk.next().unwrap();
                let last = chunk.last().unwrap_or(first);
                (
                    self.alignments[0][0]
                        .1
                        .seq_b()
                        .cast_peptidoform()
                        .sub_peptide(first.0..=last.1),
                    region.clone(),
                )
            })
            .collect()
    }
}

/// Only available if features `align` and `imgt` are turned on.
/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
/// # Panics
/// If there are not two or more genes listed. If the return number is 0.
#[expect(clippy::needless_pass_by_value)]
pub fn consecutive_align<'imgt, const STEPS: u16, A: HasPeptidoform<SimpleLinear> + Eq + Clone>(
    sequence: A,
    genes: &[(GeneType, AlignType)],
    species: Option<
        HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default + 'imgt>,
    >,
    chains: Option<
        HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default + 'imgt>,
    >,
    allele: AlleleSelection,
    scoring: AlignScoring<'_>,
    return_number: usize,
    imgt: &'imgt CVIndex<IMGT>,
) -> ConsecutiveAlignment<'imgt, SimpleLinear> {
    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<
        Vec<(
            Allele<'imgt>,
            Alignment<&'imgt Peptidoform<UnAmbiguous>, Peptidoform<SimpleLinear>>,
        )>,
    > = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || {
                    (
                        sequence.cast_peptidoform().clone(),
                        species.clone(),
                        chains.clone(),
                    )
                },
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.cast_peptidoform().sub_peptide(prev..),
                        Some(std::iter::once(last.0.species).collect()),
                        Some(std::iter::once(last.0.gene.chain).collect()),
                    )
                },
            );

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: use_chains,
                allele,
                genes: Some([gene.0].into()),
            }
            .germlines(imgt)
            .map(|seq| {
                let alignment = align::<
                    STEPS,
                    &'imgt Peptidoform<UnAmbiguous>,
                    Peptidoform<SimpleLinear>,
                >(
                    seq.sequence, left_sequence.clone(), scoring, gene.1
                );
                (seq, alignment)
            })
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    ConsecutiveAlignment { alignments: output }
}

/// Only available with if features `align`, `rayon`, and `imgt` are turned on.
/// Align one sequence to multiple consecutive genes. Each gene can be controlled to be global to the left or free to allow unmatched residues between it and the previous gene.
/// If the sequence is too short to cover all genes only the genes that could be matched are returned.
/// # Panics
/// If there are not two or more genes listed. If the return number is 0.
#[cfg(feature = "rayon")]
#[expect(clippy::needless_pass_by_value)]
pub fn par_consecutive_align<
    'imgt,
    const STEPS: u16,
    A: HasPeptidoform<SimpleLinear> + Send + Sync + Eq + Clone,
>(
    sequence: A,
    genes: &[(GeneType, AlignType)],
    species: Option<
        HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default + 'imgt>,
    >,
    chains: Option<
        HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default + 'imgt>,
    >,
    allele: AlleleSelection,
    scoring: AlignScoring<'_>,
    return_number: usize,
    imgt: &'imgt CVIndex<IMGT>,
) -> ConsecutiveAlignment<'imgt, SimpleLinear> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<
        Vec<(
            Allele<'imgt>,
            Alignment<&'imgt Peptidoform<UnAmbiguous>, Peptidoform<SimpleLinear>>,
        )>,
    > = Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || {
                    (
                        sequence.cast_peptidoform().clone(),
                        species.clone(),
                        chains.clone(),
                    )
                },
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.cast_peptidoform().sub_peptide(prev..),
                        Some(std::iter::once(last.0.species).collect()),
                        Some(std::iter::once(last.0.gene.chain).collect()),
                    )
                },
            );

        if left_sequence.is_empty() {
            break;
        }

        output.push(
            Selection {
                species: use_species,
                chains: use_chains,
                allele,
                genes: Some([gene.0].into()),
            }
            .par_germlines(imgt)
            .map(|seq| {
                let alignment = align::<
                    STEPS,
                    &'imgt Peptidoform<UnAmbiguous>,
                    Peptidoform<SimpleLinear>,
                >(
                    seq.sequence, left_sequence.clone(), scoring, gene.1
                );
                (seq, alignment)
            })
            .collect::<Vec<_>>()
            .into_iter()
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    ConsecutiveAlignment { alignments: output }
}
