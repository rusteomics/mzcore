use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{
    AlignIndex, AlignScoring, AlignType, MatchType, Piece, Score,
    align_matrix::Matrix,
    diagonal_array::DiagonalArray,
    mass_alignment::{
        align_cached, mass_range, mass_range_expanded, score, score_pair, score_pair_mass_mismatch,
    },
};
use mzcore::{
    prelude::*,
    quantities::Multi,
    sequence::{HasPeptidoform, Linear},
    system::Mass,
};

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignment<Sequence> {
    lines: Vec<MultiAlignmentLine<Sequence>>,
    score: Score,
    maximal_step: u16,
    align_type: MultiAlignType,
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignment<Sequence> {
    fn debug_display(&self, mut w: impl std::fmt::Write) {
        writeln!(
            w,
            "Multi score: {} ({}/{}) max_step: {}",
            self.score.normalised, self.score.absolute, self.score.max, self.maximal_step
        )
        .unwrap();
        for line in &self.lines {
            line.debug_display(&mut w);
        }
    }

    fn new<const STEPS: u16>(
        sequences: Vec<MultiAlignmentLineTemp<'_, Sequence, STEPS>>,
        align_type: MultiAlignType,
    ) -> Self {
        Self {
            lines: sequences
                .into_iter()
                .map(|temp| MultiAlignmentLine {
                    original_index: temp.original_index,
                    sequence: temp.sequence,
                    path: temp.path,
                    start: temp.start,
                })
                .collect(),
            score: Score {
                // TODO: build a scoring function that gets a score from a MultiAlignment
                normalised: 0.0.into(),
                absolute: 0,
                max: 0,
            },
            align_type,
            maximal_step: STEPS,
        }
    }

    /// Get the sequence variance for this MMSA
    pub fn variance(&self) -> SequenceVariance {
        let mut variance = Vec::new();

        for aligned_index in 0..self
            .lines
            .iter()
            .map(|l| l.aligned_length())
            .max()
            .unwrap_or_default()
        {
            // Break when the last item is reached
            let mut element = BTreeMap::new();
            for line in &self.lines {
                if let Some(item) = line.get_item(aligned_index) {
                    let values: &mut (usize, f64) = element.entry(item).or_default();
                    values.0 += 1;
                }
            }
            variance.push(element);
        }

        variance
    }
}

/// For each location the possible sequences and their length, their depth of coverage (how often seen) and their average local confidence.
pub type SequenceVariance = Vec<BTreeMap<(SequenceElement<Linear>, u16), (usize, f64)>>;

impl<Sequence: HasPeptidoform<Linear>> std::fmt::Display for MultiAlignment<Sequence> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.debug_display(f);
        Ok(())
    }
}

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiAlignmentLine<Sequence> {
    original_index: usize,
    sequence: Sequence,
    path: Vec<MultiPiece>,
    start: usize,
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignmentLine<Sequence> {
    // Get the sequence element at the aligned index. Only returns something if the start of a step is selected.
    fn get_item(&self, index: usize) -> Option<(SequenceElement<Linear>, u16)> {
        let mut path_index = 0;
        let mut aligned_index = 0;
        let mut sequence_index = 0;
        loop {
            sequence_index += self.path[path_index].sequence_length as usize;
            aligned_index += self.path[path_index].aligned_length as usize;
            if aligned_index >= index + self.path[path_index].aligned_length as usize {
                break;
            }
            path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
        }
        let seq_index = sequence_index
            .min(self.sequence.cast_peptidoform().len())
            .saturating_sub(1);
        (aligned_index == index + self.path[path_index].aligned_length as usize
            && self.path[path_index].match_type != MatchType::Gap)
            .then(|| {
                (
                    self.sequence.cast_peptidoform().sequence()[seq_index].clone(),
                    self.path[path_index].aligned_length,
                )
            })
    }

    fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let mut seq_index = self.start;
        for piece in &self.path {
            if piece.sequence_length == 0 {
                write!(w, "{}", "-".repeat(piece.aligned_length as usize)).unwrap();
            } else {
                let subseq = &sequence[seq_index..seq_index + piece.sequence_length as usize];
                // Obviously misses mods now
                let display = subseq
                    .iter()
                    .map(|s| s.aminoacid.one_letter_code().unwrap_or('X'))
                    .collect::<String>();
                write!(
                    w,
                    "{display}{}",
                    "·".repeat((piece.aligned_length - piece.sequence_length) as usize)
                )
                .unwrap();
                seq_index += piece.sequence_length as usize;
            }
        }
        writeln!(w).unwrap();
    }

    fn aligned_length(&self) -> usize {
        let mut aligned_index = 0;
        for piece in &self.path {
            aligned_index += piece.aligned_length as usize;
        }
        aligned_index
    }
}

#[derive(Clone, Debug, PartialEq)]
struct MultiAlignmentLineTemp<'a, Sequence, const STEPS: u16> {
    original_index: usize,
    sequence: Sequence,
    masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    path: Vec<MultiPiece>,
    start: usize,
}

impl<'a, Sequence: HasPeptidoform<Linear>, const STEPS: u16>
    MultiAlignmentLineTemp<'a, Sequence, STEPS>
{
    fn single(
        original_index: usize,
        sequence: Sequence,
        masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    ) -> Self {
        let len = sequence.cast_peptidoform().len();
        Self {
            original_index,
            sequence,
            masses,
            path: vec![
                MultiPiece {
                    match_type: MatchType::FullIdentity,
                    aligned_length: 1,
                    sequence_length: 1,
                };
                len
            ],
            start: 0,
        }
    }

    /// Update this line with the given merged alignment, this does not add end gaps
    fn update(mut self, path: &[Piece], is_a: bool, offset: u16) -> Self {
        // println!(
        //     "S: {}, P: {}, a: {is_a} offset: {offset}",
        //     self.sequence.cast_peptidoform(),
        //     Piece::cigar(path, 0, 0),
        // );
        let mut path_index = 0;
        let mut aligned_index = 0;
        let mut goal = 0;
        let mut internal_index = 0;
        if offset > 0 {
            if self.path[0].match_type == MatchType::Gap {
                self.path[0].aligned_length += offset;
                internal_index = offset;
            } else {
                self.path.insert(
                    path_index,
                    MultiPiece {
                        match_type: MatchType::Gap,
                        aligned_length: offset,
                        sequence_length: 0,
                    },
                );
                path_index += 1;
            }
        }
        for piece in path {
            let aligned_step = if is_a { piece.step_b } else { piece.step_a };
            let sequence_step = if is_a { piece.step_a } else { piece.step_b };
            if sequence_step == 0 {
                self.path.insert(
                    path_index,
                    MultiPiece {
                        match_type: MatchType::Gap,
                        aligned_length: aligned_step,
                        sequence_length: sequence_step,
                    },
                );
                aligned_index += aligned_step;
                path_index += 1;
            } else {
                self.path[path_index].aligned_length += aligned_step.saturating_sub(sequence_step);
                goal += aligned_step.max(sequence_step);
                while aligned_index + internal_index < goal {
                    aligned_index += self.path[path_index].aligned_length - internal_index;
                    internal_index = 0;
                    path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
                    if path_index == self.path.len() - 1 {
                        break;
                    }
                }
            }
        }
        self
    }

    fn aligned_length(&self) -> usize {
        let mut aligned_index = 0;
        for piece in &self.path {
            aligned_index += piece.aligned_length as usize;
        }
        aligned_index
    }

    /// Pad to length
    fn pad(&mut self, goal_length: usize) {
        let mut aligned_index = 0;
        for piece in &self.path {
            aligned_index += piece.aligned_length as usize;
        }
        if aligned_index < goal_length {
            let last = self.path.len() - 1;
            let dif = goal_length - aligned_index;
            if self.path[last].match_type == MatchType::Gap {
                self.path[last].aligned_length += dif as u16;
            } else {
                self.path.push(MultiPiece {
                    match_type: MatchType::Gap,
                    aligned_length: dif as u16,
                    sequence_length: 0,
                });
            }
        }
    }

    // Get the sequence index for an aligned index
    fn get_sequence_index(&self, index: usize) -> usize {
        let mut path_index = 0;
        let mut aligned_index = 0;
        let mut sequence_index = 0;
        while aligned_index < index {
            sequence_index += self.path[path_index].sequence_length as usize;
            aligned_index += self.path[path_index].aligned_length as usize;
            path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
        }
        sequence_index.min(self.sequence.cast_peptidoform().len())
    }

    // Sequence index
    fn get_aa_at(&self, sequence_index: usize) -> (&SequenceElement<Linear>, &Multi<Mass>) {
        let c = self.sequence.cast_peptidoform();
        let i = (sequence_index).min(c.len()).saturating_sub(1);
        unsafe {
            (
                &c.sequence().get_unchecked(i),
                &self.masses.get_unchecked([i, 0]),
            )
        }
    }

    // Sequence index, sequence len
    fn get_block_at(
        &self,
        sequence_index: usize,
        len: usize,
    ) -> Option<(&[SequenceElement<Linear>], &Multi<Mass>)> {
        if sequence_index >= len {
            Some(unsafe {
                (
                    &self
                        .sequence
                        .cast_peptidoform()
                        .sequence()
                        .get_unchecked(sequence_index - len..sequence_index),
                    &self.masses.get_unchecked([sequence_index - 1, len - 1]),
                )
            })
        } else {
            None
        }
    }

    fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let mut seq_index = self.start;
        for piece in &self.path {
            if piece.sequence_length == 0 {
                write!(w, "{}", "-".repeat(piece.aligned_length as usize)).unwrap();
            } else {
                let subseq = &sequence[seq_index..seq_index + piece.sequence_length as usize];
                // Obviously misses mods now
                let display = subseq
                    .iter()
                    .map(|s| s.aminoacid.one_letter_code().unwrap_or('X'))
                    .collect::<String>();
                write!(
                    w,
                    "{display}{}",
                    "·".repeat((piece.aligned_length - piece.sequence_length) as usize)
                )
                .unwrap();
                seq_index += piece.sequence_length as usize;
            }
        }
        writeln!(w).unwrap();
    }
}

#[derive(
    Clone, Copy, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize,
)]
struct MultiPiece {
    match_type: MatchType,
    /// aligned_length is required to always be at least sequence_length
    aligned_length: u16,
    sequence_length: u16,
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignType {
    pub left: MultiAlignSide,
    pub right: MultiAlignSide,
}

impl MultiAlignType {
    pub const EITHER_GLOBAL: Self = Self {
        left: MultiAlignSide::EitherGlobal,
        right: MultiAlignSide::EitherGlobal,
    };
    pub const GLOBAL: Self = Self {
        left: MultiAlignSide::Global,
        right: MultiAlignSide::Global,
    };
}

#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum MultiAlignSide {
    Global,
    EitherGlobal,
}

impl From<MultiAlignType> for AlignType {
    fn from(value: MultiAlignType) -> Self {
        AlignType {
            left: match value.left {
                MultiAlignSide::Global => crate::Side::Specified { a: true, b: true },
                MultiAlignSide::EitherGlobal => crate::Side::EitherGlobal,
            },
            right: match value.right {
                MultiAlignSide::Global => crate::Side::Specified { a: true, b: true },
                MultiAlignSide::EitherGlobal => crate::Side::EitherGlobal,
            },
        }
    }
}

impl<const STEPS: u16, Sequence: HasPeptidoform<Linear> + Clone> AlignIndex<STEPS, Sequence> {
    /// Get matrix with the distances between all sequences in this index. It uses [`Alignment::distance`] as metric.
    /// # Panics
    /// If more than [`u16::MAX`] sequences are given.
    pub fn distance_matrix(
        &self,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> DiagonalArray<f64, { u16::MAX }> {
        debug_assert!(self.sequences.len() <= u16::MAX.into());
        let align_type = align_type.into();
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
        matrix
    }

    pub fn multi_align(
        &self,
        maximal_distance: Option<f64>, // Maximal distance to join clusters, or join all if set to None
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Vec<MultiAlignment<Sequence>> {
        let distance_matrix = self.distance_matrix(scoring, align_type);
        self.multi_align_inner(maximal_distance, distance_matrix, scoring, align_type)
    }

    fn multi_align_inner(
        &self,
        maximal_distance: Option<f64>, // Maximal distance to join clusters, or join all if set to None
        mut distance_matrix: DiagonalArray<f64, { u16::MAX }>,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Vec<MultiAlignment<Sequence>> {
        assert!(self.sequences.len() > 1);

        let mut nodes = Vec::with_capacity(self.sequences.len());
        for (i, s) in self.sequences.iter().enumerate() {
            nodes.push(vec![MultiAlignmentLineTemp::single(i, s.0.clone(), &s.1)]);
        }

        // Progressively merge (this creates a phylogenetic tree that is not stored in any way)
        while nodes.len() > 1 {
            let (distance, [far_index, close_index]) = distance_matrix.min::<true>();
            if maximal_distance.is_some_and(|m| distance > m) {
                // Stop if this distance is more than the threshold, or if the outgroup would be merged
                break;
            }
            // Merge nodes (this is where the actual alignment happens)
            let close_node = std::mem::take(&mut nodes[close_index]);
            let far_node = nodes.remove(far_index);
            nodes[close_index] =
                multi_align_cached::<STEPS, Sequence>(close_node, far_node, scoring, align_type);

            // Update matrix to merge these columns
            distance_matrix = distance_matrix.merge_columns(close_index, far_index, distance);
        }

        // Create final data structures from the temporary ones
        nodes
            .into_iter()
            .map(|seqs| MultiAlignment::new(seqs, align_type))
            .collect()
    }
}

impl<const STEPS: u16, Sequence: HasPeptidoform<Linear> + Clone + Send + Sync>
    AlignIndex<STEPS, Sequence>
{
    /// Get matrix with the distances between all sequences in this index. It uses [`Alignment::distance`] as metric.
    /// # Panics
    /// If more than [`u16::MAX`] sequences are given.
    #[cfg(feature = "rayon")]
    pub fn par_distance_matrix(
        &self,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> DiagonalArray<f64, { u16::MAX }> {
        use rayon::prelude::*;
        debug_assert!(self.sequences.len() <= u16::MAX.into());
        let align_type = align_type.into();
        let mut matrix = DiagonalArray::<f64, { u16::MAX }>::new(self.sequences.len());
        for results in self
            .sequences
            .par_iter()
            .enumerate()
            .flat_map(move |(i, s)| {
                self.sequences
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

        matrix
    }

    /// This parallelizes the distance matrix calculation, the multiple sequence merging is still done sequentially
    #[cfg(feature = "rayon")]
    pub fn par_multi_align(
        &self,
        maximal_distance: Option<f64>, // Maximal distance to join clusters, or join all if set to None
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Vec<MultiAlignment<Sequence>> {
        let distance_matrix = self.par_distance_matrix(scoring, align_type);
        self.multi_align_inner(maximal_distance, distance_matrix, scoring, align_type)
    }
}

/// Do mass based alignment of two MMSA clusters, but with precomputed masses
#[expect(clippy::too_many_lines)]
#[allow(clippy::similar_names)]
pub(super) fn multi_align_cached<'a, const STEPS: u16, Sequence: HasPeptidoform<Linear>>(
    a: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    b: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    scoring: AlignScoring<'_>,
    align_type: MultiAlignType,
) -> Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>> {
    let len_a = a
        .iter()
        .fold(0, |acc, i| acc.max(i.sequence.cast_peptidoform().len()));
    let len_b = b
        .iter()
        .fold(0, |acc, i| acc.max(i.sequence.cast_peptidoform().len()));
    let mut matrix = Matrix::new(len_a, len_b);
    let mut global_highest = (0, 0, 0);

    if align_type.left == MultiAlignSide::Global {
        matrix.global_start(true, scoring);
        matrix.global_start(false, scoring);
    }

    // These two arrays serve for optimizational purposes.
    // In `score` routine equality of masses is checked
    // between masses_a[index_a - len_a .. index_a] and masses_b[index_b - len_b .. index_b]
    // ranges which is expensive.
    // However, in most of the cases masses are very different there and the check fails.
    // We precompute two following arrays:
    // - ranges_a[index_a - len_a .. index_a] = (min_mass, max_mass)
    // min_mass is the minimal mass that matches at least one mass in masses_a[index_a - len_a .. index_a].
    // max_mass is the same for maximal mass.
    // - ranges_b[index_b - len_b .. index_b] = (min_mass, max_mass)
    // min_mass in the minimal mass across masses_b[index_b - len_b .. index_b].
    // max_mass is the same for maximal mass.
    // Then, before calling `score` we can perform a quick check whether
    // the ranges overlap for given slices. If they do not overlap, we skip calling `score`
    // as it is guaranteed that there are no two masses that are within tolerance.
    let mass_ranges_a = a
        .iter()
        .map(|a| {
            let mut ranges: DiagonalArray<(Mass, Mass), STEPS> =
                DiagonalArray::new(a.sequence.cast_peptidoform().len());
            for i in 0..a.sequence.cast_peptidoform().len() {
                for j in 0..=i.min(STEPS as usize) {
                    let (min, max) = mass_range_expanded(
                        unsafe { a.masses.get_unchecked([i, j]) },
                        scoring.tolerance,
                    );
                    ranges[[i, j]] = (min, max);
                }
            }
            ranges
        })
        .collect::<Vec<_>>();

    let mass_ranges_b = b
        .iter()
        .map(|b| {
            let mut ranges: DiagonalArray<(Mass, Mass), STEPS> =
                DiagonalArray::new(b.sequence.cast_peptidoform().len());
            for i in 0..b.sequence.cast_peptidoform().len() {
                for j in 0..=i.min(STEPS as usize) {
                    let (min, max) = mass_range(unsafe { b.masses.get_unchecked([i, j]) });
                    ranges[[i, j]] = (min, max);
                }
            }
            ranges
        })
        .collect::<Vec<_>>();

    // First index_b then the sequence at that index
    let sequence_indices_b: Vec<Vec<usize>> = (1..=len_b)
        .map(|index_b| {
            b.iter()
                .map(|s| s.get_sequence_index(index_b))
                .collect::<Vec<_>>()
        })
        .collect();

    let sorted_pair_masses_b: Vec<Vec<((Mass, Mass), usize)>> = sequence_indices_b
        .iter()
        .map(|seq_b_indices| {
            let mut sorted_masses_b: Vec<((Mass, Mass), usize)> = seq_b_indices
                .iter()
                .enumerate()
                .map(|(ib, is)| {
                    (
                        mass_range(unsafe {
                            b[ib].masses.get_unchecked([is.saturating_sub(1), 0])
                        }),
                        ib,
                    )
                })
                .collect();
            sorted_masses_b.sort_by(|a, b| a.0.0.value.total_cmp(&b.0.0.value));
            sorted_masses_b
        })
        .collect();

    for index_a in 1..=len_a {
        // Precalculate sequence indices and sorted masses
        let sequence_indices_a = a
            .iter()
            .map(|s| s.get_sequence_index(index_a))
            .collect::<Vec<_>>();
        let mut sorted_pair_masses_a: Vec<((Mass, Mass), usize)> = sequence_indices_a
            .iter()
            .enumerate()
            .map(|(ia, is)| {
                (
                    mass_range_expanded(
                        unsafe { a[ia].masses.get_unchecked([is.saturating_sub(1), 0]) },
                        scoring.tolerance,
                    ),
                    ia,
                )
            })
            .collect();
        sorted_pair_masses_a.sort_by(|a, b| a.0.0.value.total_cmp(&b.0.0.value));

        for index_b in 1..=len_b {
            // Returns the score for a gap transition.
            // gap_a controls whether the gap is in a or b.
            let score_gap = |gap_a: bool| {
                let prev = if gap_a {
                    unsafe { matrix.get_unchecked([index_a - 1, index_b]) }
                } else {
                    unsafe { matrix.get_unchecked([index_a, index_b - 1]) }
                };

                let is_first_step = prev.step_a == 0 && prev.step_b == 0;
                let is_previous_gap = prev.step_a == 0 && !gap_a || prev.step_b == 0 && gap_a;
                let is_gap_start = is_first_step || !is_previous_gap;
                // First check the score to be used for affine gaps
                let score = scoring.gap_extend as isize
                    + scoring.gap_start as isize * isize::from(is_gap_start);

                let len_a = u16::from(gap_a);
                let len_b = u16::from(!gap_a);
                Piece {
                    score: prev.score + score,
                    local_score: score,
                    match_type: MatchType::Gap,
                    step_a: len_a,
                    step_b: len_b,
                }
            };

            // First try all gap possibilities.
            let gap_score_a = score_gap(true);
            let gap_score_b = score_gap(false);

            let mut highest = if gap_score_a.score >= gap_score_b.score {
                gap_score_a
            } else {
                gap_score_b
            };

            // Now try matching single aminoacids.
            let prev = unsafe { matrix.get_unchecked([index_a - 1, index_b - 1]) };

            let mut sorted_masses_b_index = 0;
            let sorted_masses_b_local = &sorted_pair_masses_b[index_b - 1];
            let total_masses_b = sorted_masses_b_local.len();

            // Go over all masses at this position in A, this allows finding perfect matches quickly
            'perfect_loop: for ((low_mass_a, high_mass_a), ia) in &sorted_pair_masses_a {
                // Skip all masses in B that are too low
                while sorted_masses_b_local[sorted_masses_b_index].0.1 < *low_mass_a {
                    sorted_masses_b_index += 1;
                    if sorted_masses_b_index == total_masses_b {
                        break 'perfect_loop; // No more masses in B so all masses left in A can be disregarded
                    }
                }
                // Take all masses in B that are within the range
                let mut offset = 0;
                while sorted_masses_b_index + offset < total_masses_b
                    && sorted_masses_b_local[sorted_masses_b_index + offset].0.0 <= *high_mass_a
                {
                    let ib = sorted_masses_b_local[sorted_masses_b_index + offset].1;
                    let pair_score = score_pair(
                        a[*ia].get_aa_at(sequence_indices_a[*ia]),
                        b[ib].get_aa_at(sequence_indices_b[index_b - 1][ib]),
                        scoring,
                        prev.score,
                    );
                    if pair_score.score > highest.score {
                        highest = pair_score;
                    }

                    offset += 1;
                }
            }

            if highest.match_type != MatchType::FullIdentity {
                // Check if an identity but mass mismatch is a possible alignment
                for line_a in 0..a.len() {
                    let aa_a = unsafe {
                        &a[line_a]
                            .sequence
                            .cast_peptidoform()
                            .sequence()
                            .get_unchecked(sequence_indices_a[line_a].saturating_sub(1))
                    };
                    for line_b in 0..b.len() {
                        let pair_score = score_pair_mass_mismatch(
                            aa_a,
                            unsafe {
                                &b[line_b]
                                    .sequence
                                    .cast_peptidoform()
                                    .sequence()
                                    .get_unchecked(
                                        sequence_indices_b[index_b - 1][line_b].saturating_sub(1),
                                    )
                            },
                            scoring,
                            prev.score,
                        );
                        if pair_score.score > highest.score {
                            highest = pair_score;
                        }
                    }
                }

                // Check all mass based steps
                for line_a in 0..a.len() {
                    let sequence_index_a = sequence_indices_a[line_a];
                    for line_b in 0..b.len() {
                        let sequence_index_b = sequence_indices_b[index_b - 1][line_b];
                        for len_a in 1..=sequence_index_a.min(STEPS as usize) {
                            let range_a = unsafe {
                                mass_ranges_a[line_a]
                                    .get_unchecked([sequence_index_a.saturating_sub(1), len_a - 1])
                            };

                            let min_len_b = if len_a == 1 { 2 } else { 1 };

                            for len_b in min_len_b..=sequence_index_b.min(STEPS as usize) {
                                let range_b = unsafe {
                                    mass_ranges_b[line_b].get_unchecked([
                                        sequence_index_b.saturating_sub(1),
                                        len_b - 1,
                                    ])
                                };
                                // Note that ranges are already expanded by tolerance, so
                                // exact comparison is fine here.
                                if range_a.0 > range_b.1 || range_b.0 > range_a.1 {
                                    continue;
                                }
                                // len_a and b are always <= STEPS
                                let match_score = {
                                    let prev = unsafe {
                                        matrix.get_unchecked([index_a - len_a, index_b - len_b])
                                    };
                                    let base_score = prev.score;
                                    a[line_a]
                                        .get_block_at(sequence_index_a, len_a) // Think about if this len should be sequence or aligned index
                                        .and_then(|block_a| {
                                            b[line_b]
                                                .get_block_at(sequence_index_b, len_b)
                                                .map(|block_b| (block_a, block_b))
                                        })
                                        .and_then(|(block_a, block_b)| {
                                            score(block_a, block_b, scoring, base_score)
                                        })
                                };
                                if let Some(p) = match_score
                                    && p.score > highest.score
                                {
                                    highest = p;
                                }
                            }
                        }
                    }
                }
            }
            if highest.score >= global_highest.0 {
                global_highest = (highest.score, index_a, index_b);
            }
            unsafe {
                *matrix.get_unchecked_mut([index_a, index_b]) = highest;
            }
        }
    }

    // Finish up by tracing the path and updating all enclosed sequences to this path
    let (start_a, start_b, path) = matrix.trace_path(align_type.into(), global_highest);

    let mut sequences = Vec::with_capacity(a.len() + b.len());

    for line in a {
        sequences.push(line.update(&path, true, start_b as u16));
    }

    for line in b {
        sequences.push(line.update(&path, false, start_a as u16));
    }

    // Fix end gaps
    let goal_length = sequences
        .iter()
        .fold(0, |acc, s| acc.max(s.aligned_length()));
    for s in &mut sequences {
        s.pad(goal_length);
    }

    sequences
}

#[allow(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::mass_alignment::calculate_masses;
    use itertools::Itertools;
    use mzcore::ontology::STATIC_ONTOLOGIES;

    fn seq(def: &str) -> Peptidoform<Linear> {
        Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_linear()
            .unwrap()
    }

    #[test]
    fn simple() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let sequences = vec![seq("AGGWHD"), seq("ANWHN[Deamidated]"), seq("AHYDH")];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "AGGWHD\nAN·WHN\nAH·YDH\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn many() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let sequences = vec![
            seq("WRGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("WRNDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWNDGFYW[U:Oxidation]DYWGQG"),
            seq("RWGGN[U:Deamidated]GFYW[U:Oxidation]DYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "WRGGDGFYAMDYWGQG\nWRN·DGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWN·DGFY·WDYWGQG\nRWGGNGFY·WDYWGQG\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn gaps() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let sequences = vec![
            seq("WRGGDGFYAMDYWGQG"),
            seq("RWGDGFYAMYDWGQG"),
            seq("RWNDGYAMEDYWGQG"),
            seq("RWGGDGFMDYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected =
            "WRGGDGFYAM-DYWGQG\nRWGGDGF--M-DYWGQG\nRWG-DGFYAM-YDWGQG\nRWN·DG-YAMEDYWGQG\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn either_global() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let sequences = vec![
            seq("WRGGDGFYAMDYWGQG"),
            seq("RWGDGFYAMYD"),
            seq("RWNDGYAMDYWGQG"),
            seq("GDGFMDYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment =
            index.multi_align(None, AlignScoring::default(), MultiAlignType::EITHER_GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "WRGGDGFYAMDYWGQG\n---GDGF--MDYWGQG\nRWN·DG-YAMDYWGQG\nRWG-DGFYAMYD----\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn either_global_2() {
        let sequences = vec![
            seq("CSRWRGGDGF"),
            seq("WRNDGFYAM"),
            seq("RWGGDGFYAMDYWG"),
            seq("YW[U:Oxidation]DYWGQG"),
            seq("RWGGNGFYW[U:Oxidation]DYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let scoring = AlignScoring::<'_> {
            tolerance: mzcore::quantities::Tolerance::Relative(
                mzcore::system::Ratio::new::<mzcore::system::ratio::ppm>(20.0).into(), // AM W[Oxidation] are 16ppm apart
            ),
            ..Default::default()
        };
        let alignment = index.multi_align(None, scoring, MultiAlignType::EITHER_GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "CSRWRGGDGF---------\n---WRN·DGFYAM------\n---RWGGDGFYAMDYWG--\n----------YW·DYWGQG\n---RWGGNGFYW·DYWGQG\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn crash() {
        let sequences = vec![
            seq("HYTTPPTFGQGT"),
            seq("WGG"),
            seq("VTC[U:Carboxymethyl]QGLSSPKSL"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let scoring = AlignScoring::<'_> {
            tolerance: mzcore::quantities::Tolerance::Relative(
                mzcore::system::Ratio::new::<mzcore::system::ratio::ppm>(20.0).into(),
            ),
            ..Default::default()
        };
        let alignment = index.multi_align(None, scoring, MultiAlignType::EITHER_GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "----HYTTPPTFGQGT\n-----------WG-G-\nVTCQGLSSPKSL----\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn properties() {
        let sequence = seq("AGGWHD");
        let masses = calculate_masses::<4>(&sequence, MassMode::Monoisotopic);
        let mut line = MultiAlignmentLineTemp::single(0, &sequence, &masses);
        assert_eq!(line.get_sequence_index(4), 4);
        assert_eq!(line.get_sequence_index(5), 5);
        assert_eq!(
            line.get_aa_at(5),
            (
                &SequenceElement::new(AminoAcid::Histidine.into(), None).into(),
                &Multi::from(CheckedAminoAcid::Histidine.formula().monoisotopic_mass())
            )
        );
        assert_eq!(
            line.get_block_at(5, 2),
            Some((
                [
                    SequenceElement::new(AminoAcid::Tryptophan.into(), None).into(),
                    SequenceElement::new(AminoAcid::Histidine.into(), None).into(),
                ]
                .as_slice(),
                &Multi::from(
                    (CheckedAminoAcid::Tryptophan.formula()
                        + CheckedAminoAcid::Histidine.formula())
                    .monoisotopic_mass()
                )
            ))
        );
        let updated = line.clone().update(
            &[
                Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // A
                Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
                Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
                Piece::new(0, 0, MatchType::Isobaric, 1, 2),     // W
                Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // H
                Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
            ],
            true,
            0,
        );
        line.path[3].aligned_length = 2; // The W is bigger
        assert_eq!(updated, line);
        assert_eq!(line.get_sequence_index(4), 4);
        assert_eq!(line.get_sequence_index(5), 4);
        assert_eq!(line.get_sequence_index(6), 5);
        assert_eq!(
            line.get_aa_at(5),
            (
                &SequenceElement::new(AminoAcid::Histidine.into(), None).into(),
                &Multi::from(CheckedAminoAcid::Histidine.formula().monoisotopic_mass())
            )
        );
        assert_eq!(
            line.get_block_at(5, 2),
            Some((
                [
                    SequenceElement::new(AminoAcid::Tryptophan.into(), None).into(),
                    SequenceElement::new(AminoAcid::Histidine.into(), None).into(),
                ]
                .as_slice(),
                &Multi::from(
                    (CheckedAminoAcid::Tryptophan.formula()
                        + CheckedAminoAcid::Histidine.formula())
                    .monoisotopic_mass()
                )
            ))
        );
    }

    #[test]
    fn properties_gaps() {
        let sequence = seq("WRGGDGFYAMDYWGQG");
        let sequence2 = seq("RWNDGYAMEDYWGQG");
        let path = &[
            Piece::new(0, 0, MatchType::Rotation, 2, 2),     // WR
            Piece::new(0, 0, MatchType::Isobaric, 2, 1),     // GG N
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
            Piece::new(0, 0, MatchType::Gap, 1, 0),          // F
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Y
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // A
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // M
            Piece::new(0, 0, MatchType::Gap, 0, 1),          // E
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Y
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // W
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Q
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
        ];
        let masses = calculate_masses::<4>(&sequence, MassMode::Monoisotopic);
        let masses2 = calculate_masses::<4>(&sequence2, MassMode::Monoisotopic);
        let line = MultiAlignmentLineTemp::single(0, &sequence, &masses);
        let line2 = MultiAlignmentLineTemp::single(0, &sequence2, &masses2);
        let updated = line.update(path, true, 0);
        let updated2 = line2.update(path, false, 0);
        let mut buf1 = String::new();
        updated.debug_display(&mut buf1);
        let mut buf2 = String::new();
        updated2.debug_display(&mut buf2);
        assert_eq!(buf1, "WRGGDGFYAM-DYWGQG\n");
        assert_eq!(buf2, "RWN·DG-YAMEDYWGQG\n");
    }
}
