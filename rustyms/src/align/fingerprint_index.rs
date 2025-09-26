#[cfg(feature = "rayon")]
use rayon::prelude::*;

use std::cmp;
use std::hash::{DefaultHasher, Hasher};
use std::{collections::HashMap, hash::Hash};

use crate::align::scoring::AlignScoring;
use crate::align::AlignType;
use crate::align::Alignment;
use crate::quantities::Tolerance;
use crate::system::OrderedMass;
use crate::{align::{diagonal_array::DiagonalArray, mass_alignment::calculate_masses}, prelude::Peptidoform, quantities::Multi, sequence::{HasPeptidoform, SimpleLinear}, system::Mass};

/// This trait generates a number of fingerprints for a given peptidoform.
/// Fingerprints have a property that if two peptidoforms have good alignment, they share at least one fingerprint.
/// This is used as a pre-filter for batch alignment: only peptidoforms sharing a fingerprint with target are selected.
pub trait FingerprintGenerator<S: HasPeptidoform<SimpleLinear>> {
    /// Generate a fingerprint library for the given peptidoform.
    /// If None is returned, the peptidoform will always be aligned againts all queries.
    fn gen_fingerprint_library(&self, peptidoform: &S, scoring: &AlignScoring) -> Option<Vec<u64>>;
    /// Generate a fingerprint query for the given peptidoform.
    /// If None is returned, the peptidoform will always be aligned againts all templates in library.
    fn gen_fingerprint_query(&self, peptidoform: &S, scoring: &AlignScoring) -> Option<Vec<u64>>;
}

/// An implementation of `FingerprintGenerator` that generates fingerprints based on mass
/// of segments of the peptidoform.
/// For each segment of length between `min_length` and `max_length`, all possible splits
/// into at most `max_segment_count` subsegments with at most `max_gaps` single amino acid
/// gaps are considered. For every subsegment, the mass is calculated in quantized into
/// bin index based on the tolerance. The fingerprint is the hash of the vector.
/// Example: consier peptidoform "ABCD" with `min_length=2`, max_length=3, max_segment_count=2, max_gaps=1.
/// The following fingerprints are possible (not an exhaustive list):
/// [A][B] -- just two consecutive segments of length 1.
/// [AB][C] -- two segments, one of length 2.
/// [A]x[CD] -- two segments, one of length 1, one of length 2, with a gap in between.
/// Note, that [A]xx[D], [A][B][C], and [ABCD] are not allowed due to the constraints.
#[derive(Debug)]
pub struct MassFingerprintGenerator<const STEPS: u16> {
    options: MassFingerprintGeneratorOptions,
}

impl<const STEPS: u16> MassFingerprintGenerator<STEPS> {
    /// Create a new `MassFingerprintGenerator` with the given options.
    pub fn new(options: MassFingerprintGeneratorOptions) -> Self {
        Self { options }
    }

    fn generate(
        &self,
        peptidoform: &Peptidoform<SimpleLinear>,
        index: usize,
        masses: &DiagonalArray<Multi<Mass>, STEPS>,
        current: &mut Vec<i64>,
        current_length: usize,
        current_gaps: usize,
        current_segment_count: usize,
        generate_neighbor_bins: bool,
        result: &mut Vec<u64>,
        scoring: &AlignScoring)
    {
        if current_length >= self.options.min_length {
            result.push(Self::hash_vec(current));
        }

        // Try to add another segment.
        if current_segment_count < self.options.max_segment_count && current_length < self.options.max_length && index > 0 {
            let mut max_length = cmp::min(index, self.options.max_length - current_length);
            max_length = cmp::min(max_length, STEPS as usize);
            for length in 1..=max_length {
                for mass in masses[[index - 1, length - 1]].iter() {
                    let bin_index = get_bin_index(*mass, scoring.tolerance);
                    let min_bin_index = if generate_neighbor_bins {
                        bin_index.saturating_sub(1)
                    } else {
                        bin_index
                    };
                    let max_bin_index = if generate_neighbor_bins {
                        bin_index.saturating_add(1)
                    } else {
                        bin_index
                    };
                    for bin in min_bin_index..=max_bin_index {
                        current.push(bin);
                        self.generate(
                            peptidoform,
                            index - length,
                            masses,
                            current,
                            current_length + length,
                            current_gaps,
                            current_segment_count + 1,
                            generate_neighbor_bins,
                            result,
                            scoring,
                        );
                        current.pop();
                    }
                }
            }
        }

        // Try to add a gap.
        if current_gaps < self.options.max_gaps && current_length < self.options.max_length && index > 0 {
            self.generate(
                peptidoform,
                index - 1,
                masses,
                current,
                current_length + 1,
                current_gaps + 1,
                current_segment_count,
                generate_neighbor_bins,
                result,
                scoring,
            );
        }
    }

    fn hash_vec(vec: &Vec<i64>) -> u64 {
        let mut hasher = DefaultHasher::new();
        vec.hash(&mut hasher);
        hasher.finish()
    }
}

impl<const STEPS: u16> FingerprintGenerator<Peptidoform<SimpleLinear>> for MassFingerprintGenerator<STEPS> {
    fn gen_fingerprint_library(&self, peptidoform: &Peptidoform<SimpleLinear>, scoring: &AlignScoring) -> Option<Vec<u64>> {
        let masses = calculate_masses::<STEPS>(peptidoform, scoring.mass_mode);

        let mut result = Vec::new();
        let mut current = Vec::new();

        for index in 1..=peptidoform.len() {
            self.generate(
                peptidoform,
                index,
                &masses,
                &mut current,
                0,
                0,
                0,
                false,
                &mut result,
                scoring,
            );
        }

        Some(result)
    }

    fn gen_fingerprint_query(&self, peptidoform: &Peptidoform<SimpleLinear>, scoring: &AlignScoring) -> Option<Vec<u64>> {
        let masses = calculate_masses::<STEPS>(peptidoform, scoring.mass_mode);

        let mut result = Vec::new();
        let mut current = Vec::new();

        for index in 1..=peptidoform.len() {
            self.generate(
                peptidoform,
                index,
                &masses,
                &mut current,
                0,
                0,
                0,
                true,
                &mut result,
                scoring,
            );
        }

        Some(result)
    }
}

#[derive(Clone, Copy, Debug)]
/// Options for the mass fingerprint generator.
pub struct MassFingerprintGeneratorOptions {
    /// Minimum length of the peptidoform segment to consider.
    pub min_length: usize,
    /// Maximum length of the peptidoform segment to consider.
    pub max_length: usize,
    /// Maximum number of segments to split into.
    pub max_segment_count: usize,
    /// Maximum number of allowed single amino acid gaps.
    pub max_gaps: usize,
}

impl Default for MassFingerprintGeneratorOptions {
    fn default() -> Self {
        Self {
            min_length: 6,
            max_length: 7,
            max_segment_count: 3,
            max_gaps: 1,
        }
    }
}

fn get_bin_index(mass: Mass, tolerance: Tolerance<OrderedMass>) -> i64 {
    match tolerance {
        Tolerance::Absolute(eps) => (mass.value / eps.value).round() as i64,
        Tolerance::Relative(ppm) => {
            // To get to a bin index for relative tolerance, we use a logarithmic scale.
            // Original range is [value * (1 - ppm / 10^6), value * (1 + ppm / 10^6)]
            // Log range is [ln(value) + ln(1 - ppm / 10^6), ln(value) + ln(1 + ppm / 10^6)]
            // We may assume conservatively that the range is
            // [ln(value) + ln(1 - ppm / 10^6), ln(value) - ln(1 - ppm / 10^6)]
            let log_mass = mass.value.ln();
            let log_tolerance = (1.0 - ppm.value / 1_000_000.0).ln().abs();
            (log_mass / log_tolerance).round() as i64
        }
    }
}

/// The index that uses fingerprint generation for fast pre-filtering of candidate peptidoforms for alignment.
#[derive(Debug)]
pub struct FingerprintIndex<'a, const STEPS: u16, S: HasPeptidoform<SimpleLinear>, G: FingerprintGenerator<S>> {
    sequences: Vec<(S, DiagonalArray<Multi<Mass>, STEPS>)>,
    generator: G,
    scoring: AlignScoring<'a>,

    index: HashMap<u64, Vec<usize>>,

    naughty_sequences: Vec<usize>,
}

impl<'a, const STEPS: u16, S: HasPeptidoform<SimpleLinear> + Clone + 'a, G: FingerprintGenerator<S>> FingerprintIndex<'a, STEPS, S, G> {
    /// Constructs a new index based on the given sequences, fingerprint generator, and scoring.
    pub fn new(sequences: Vec<S>, generator: G, scoring: AlignScoring<'a>) -> Self {
        let sequences = sequences
            .into_iter()
            .map(|p| {
                let masses = calculate_masses::<STEPS>(p.cast_peptidoform(), scoring.mass_mode);
                (p, masses)
            })
            .collect::<Vec<_>>();

        let mut index: HashMap<u64, Vec<usize>> = HashMap::new();
        let mut naughty_sequences = Vec::new();

        for (i, seq) in sequences.iter().enumerate() {
            if let Some(fingerprints) = generator.gen_fingerprint_library(&seq.0, &scoring) {
                for fingerprint in fingerprints {
                    index.entry(fingerprint).or_default().push(i);
                }
            } else {
                naughty_sequences.push(i);
            }
        }

        Self {
            sequences,
            generator,
            scoring,
            index,
            naughty_sequences,
        }
    }

    /// Align one peptidoform to the index, returning an iterator over alignments.
    pub fn align_one(
        &'a self,
        other: S,
        align_type: AlignType,
    ) -> impl ExactSizeIterator<Item = Alignment<S, S>> + 'a {
        let other_masses = calculate_masses::<STEPS>(other.cast_peptidoform(), self.scoring.mass_mode);

        let query_indices = self.get_query_indices(&other);

        query_indices.into_iter().map(move |i| {
            let (template, template_masses) = &self.sequences[i];
            crate::align::mass_alignment::align_cached::<STEPS, S, S>(
                template.clone(),
                template_masses,
                other.clone(),
                &other_masses,
                self.scoring,
                align_type,
            )
        })
    }

    /// Align multiple peptidoforms to the index.
    pub fn align(
        &'a self,
        others: impl IntoIterator<Item = S>,
        align_type: AlignType,
    ) -> impl Iterator<Item = impl ExactSizeIterator<Item = Alignment<S, S>>> {
        others
            .into_iter()
            .map(move |o| self.align_one(o, align_type))
    }

    fn get_query_indices(&self, peptidoform: &S) -> Vec<usize> {
        if let Some(fingerprints) = self.generator.gen_fingerprint_query(peptidoform, &self.scoring) {
            let mut indices = Vec::new();
            for fingerprint in fingerprints {
                if let Some(idxs) = self.index.get(&fingerprint) {
                    indices.extend_from_slice(idxs);
                }
            }
            indices.sort_unstable();
            indices.dedup();
            indices.extend_from_slice(&self.naughty_sequences);
            indices
        } else {
            (0..self.sequences.len()).collect()
        }
    }
}

impl<'a, const STEPS: u16, S: HasPeptidoform<SimpleLinear> + Clone + Send + Sync + 'a, G: FingerprintGenerator<S> + Send + Sync> FingerprintIndex<'a, STEPS, S, G> {
    #[cfg(feature = "rayon")]
    /// Align multiple peptides to this index in parallel
    pub fn par_align(
        &'a self,
        others: impl IntoParallelIterator<Item = S>,
        align_type: AlignType,
    ) -> impl ParallelIterator<Item = impl ExactSizeIterator<Item = Alignment<S, S>>> {
        others
            .into_par_iter()
            .map(move |o| self.align_one(o, align_type))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn create_peptide(sequence: &str) -> Peptidoform<SimpleLinear> {
        Peptidoform::pro_forma(sequence, None)
            .unwrap()
            .into_simple_linear()
            .unwrap()
    }

    fn default_scoring() -> AlignScoring<'static> {
        AlignScoring::default()
    }

    fn default_options() -> MassFingerprintGeneratorOptions {
        MassFingerprintGeneratorOptions {
            min_length: 4,
            max_length: 6,
            max_segment_count: 4,
            max_gaps: 1,
        }
    }

    #[test]
    fn fingerprint_intersection_simple() {
        let generator = MassFingerprintGenerator::<4>::new(default_options());
        let scoring = default_scoring();

        let intersects = |a: &str, b: &str| -> bool {
            let peptide1 = create_peptide(a);
            let peptide2 = create_peptide(b);

            let fp1 = generator.gen_fingerprint_library(&peptide1, &scoring).unwrap();
            let fp2 = generator.gen_fingerprint_query(&peptide2, &scoring).unwrap();

            let set1: HashSet<_> = fp1.into_iter().collect();
            let set2: HashSet<_> = fp2.into_iter().collect();

            let intersection: Vec<_> = set1.intersection(&set2).collect();

            !intersection.is_empty()
        };

        assert!(intersects("ANNA", "ANNA"));
        assert!(intersects("ATNA", "ANNA"));
        assert!(intersects("ANNA", "AAAAANNAAAAA"));
        assert!(intersects("IIII", "LLLL"));
        assert!(intersects("ANGGNA", "ANNNA"));
        assert!(intersects("ANSNA", "ANFNA"));
        assert!(intersects("ANTNA", "ANNA"));
        assert!(intersects("STNQ", "TNQS"));

        assert!(!intersects("STNQ", "AAAA"));
        assert!(!intersects("WWWW", "AAAA"));
    }

    #[test]
    fn test_fingerprint_index_query() {
        let sequences = vec![
            create_peptide("ANNA"),
            create_peptide("AGLA"),
            create_peptide("LLLLL"),
        ];

        let generator = MassFingerprintGenerator::<4>::new(default_options());
        let scoring = default_scoring();

        let index: FingerprintIndex<4, _, _> = FingerprintIndex::new(sequences, generator, scoring);

        let query = create_peptide("ATNA");
        let results: Vec<_> = index.align_one(query, AlignType::GLOBAL).collect();

        assert!(results.len() == 1); // ATNA should match only ANNA
        assert_eq!(results[0].seq_a.to_string(), "ANNA");
        assert_eq!(results[0].seq_b.to_string(), "ATNA");
    }
}
