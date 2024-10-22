use crate::{
    align::AlignScoring,
    align::*,
    imgt::*,
    peptide::{AtMax, SimpleLinear, UnAmbiguous},
    *,
};
use std::collections::HashSet;

use itertools::Itertools;

/// A consecutive alignment, which align one sequence to multiple sequences.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConsecutiveAlignment<'lifetime, A> {
    /// All underlying alignments, per gene there is a vector containing all options for that gene.
    pub alignments: Vec<Vec<(Allele<'lifetime>, Alignment<'lifetime, UnAmbiguous, A>)>>,
}

impl<'lifetime, A> ConsecutiveAlignment<'lifetime, A> {
    /// Get the main alignment, the alignment taking the best alignment for each gene.
    pub fn main_alignment(
        &self,
    ) -> Vec<&(Allele<'lifetime>, Alignment<'lifetime, UnAmbiguous, A>)> {
        self.alignments.iter().filter_map(|a| a.first()).collect()
    }
}

impl<'lifetime, A: AtMax<Linear>> ConsecutiveAlignment<'lifetime, A> {
    /// Break up in the main alignment into the regions as annotated in the alleles.
    #[allow(clippy::missing_panics_doc)]
    pub fn regions(&self) -> Vec<(LinearPeptide<A>, Region)> {
        let mut b_offset = 0;
        self.alignments
            .iter()
            .filter_map(|a| {
                a.first().map(|(allele, alignment)| {
                    let mut index_a = alignment.start_a;
                    let mut start_region_b = alignment.start_b;
                    let mut index_b = alignment.start_b;
                    let mut region = allele.region(index_a).map_or(Region::FR1, |(r, _)| r);
                    let mut ranges = Vec::new();

                    for step in &alignment.path {
                        let new_region = allele
                            .region(index_a + step.step_a as usize)
                            .map_or(Region::FR1, |(r, _)| r);
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
                        .sub_peptide(first.0..=last.1),
                    region,
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
#[allow(clippy::needless_pass_by_value)]
pub fn consecutive_align<const STEPS: u16, A: AtMax<SimpleLinear> + AtMax<Linear>>(
    sequence: &LinearPeptide<A>,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    chains: Option<HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    allele: AlleleSelection,
    scoring: AlignScoring<'_>,
    return_number: usize,
) -> ConsecutiveAlignment<'static, A> {
    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> =
        Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || (sequence.clone(), species.clone(), chains.clone()),
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.sub_peptide(prev..),
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
            .germlines()
            .map(|seq| {
                let alignment =
                    align::<STEPS, UnAmbiguous, A>(seq.sequence, &left_sequence, scoring, gene.1)
                        .to_owned();
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
#[cfg(all(feature = "rayon"))]
#[allow(clippy::needless_pass_by_value)]
pub fn par_consecutive_align<
    const STEPS: u16,
    A: AtMax<SimpleLinear> + AtMax<Linear> + Send + Sync,
>(
    sequence: &LinearPeptide<A>,
    genes: &[(GeneType, AlignType)],
    species: Option<HashSet<Species, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    chains: Option<HashSet<ChainType, impl std::hash::BuildHasher + Clone + Send + Sync + Default>>,
    allele: AlleleSelection,
    scoring: AlignScoring<'_>,
    return_number: usize,
) -> ConsecutiveAlignment<'static, A> {
    use rayon::iter::ParallelIterator;

    assert!(genes.len() >= 2);
    assert!(return_number != 0);

    let mut output: Vec<Vec<(Allele<'static>, Alignment<'static, UnAmbiguous, A>)>> =
        Vec::with_capacity(genes.len());

    let mut prev = 0;
    for gene in genes {
        let (left_sequence, use_species, use_chains) =
            output.last().and_then(|v| v.first()).map_or_else(
                || (sequence.clone(), species.clone(), chains.clone()),
                |last| {
                    prev += last.1.start_b() + last.1.len_b();
                    (
                        sequence.sub_peptide(prev..),
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
            .par_germlines()
            .map(|seq| {
                let alignment =
                    align::<STEPS, UnAmbiguous, A>(seq.sequence, &left_sequence, scoring, gene.1);
                (seq, alignment.to_owned())
            })
            .collect::<Vec<_>>()
            .into_iter()
            .k_largest_by(return_number, |a, b| a.1.cmp(&b.1))
            .collect_vec(),
        );
    }
    ConsecutiveAlignment { alignments: output }
}
