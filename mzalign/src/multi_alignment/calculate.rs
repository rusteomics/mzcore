use std::borrow::Cow;

use mzcore::{
    prelude::*,
    quantities::Multi,
    sequence::{HasPeptidoform, Linear},
    system::Mass,
};

use crate::{
    AlignIndex, AlignScoring, MatchType, Piece,
    align_matrix::Matrix,
    diagonal_array::DiagonalArray,
    mass_alignment::{
        align_cached, mass_range, mass_range_expanded, score, score_pair, score_pair_mass_mismatch,
    },
    multi_alignment::structs::{MultiAlignSide, MultiAlignType, MultiAlignment, MultiPiece},
};

#[derive(Clone, Debug, PartialEq)]
pub(super) struct MultiAlignmentLineTemp<'a, Sequence, const STEPS: u16> {
    pub original_index: usize,
    pub sequence: Sequence,
    pub masses: Cow<'a, DiagonalArray<Multi<Mass>, STEPS>>,
    pub path: Vec<MultiPiece>,
}

impl<'a, Sequence: HasPeptidoform<Linear>, const STEPS: u16>
    MultiAlignmentLineTemp<'a, Sequence, STEPS>
{
    pub(super) fn single(
        original_index: usize,
        sequence: Sequence,
        masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    ) -> Self {
        let len = sequence.cast_peptidoform().len();
        Self {
            original_index,
            sequence,
            masses: Cow::Borrowed(masses),
            path: vec![
                MultiPiece {
                    match_type: MatchType::FullIdentity,
                    aligned_length: 1,
                    sequence_length: 1,
                };
                len
            ],
        }
    }

    /// Update this line with the given merged alignment, this does not add end gaps
    pub(super) fn update(mut self, path: &[Piece], is_a: bool, offset: u16) -> Self {
        let mut path_index = 0;
        let mut aligned_index = 0;
        let mut goal = 0;
        let mut internal_index = 0;
        if offset > 0 {
            if self.path[0].match_type == MatchType::Gap {
                self.path[0].aligned_length += offset;
                internal_index = offset;
            } else {
                self.path.insert(path_index, MultiPiece {
                    match_type: MatchType::Gap,
                    aligned_length: offset,
                    sequence_length: 0,
                });
                path_index += 1;
            }
        }
        for piece in path {
            let aligned_step = if is_a { piece.step_b } else { piece.step_a };
            let sequence_step = if is_a { piece.step_a } else { piece.step_b };
            if sequence_step == 0 {
                self.path.insert(path_index, MultiPiece {
                    match_type: MatchType::Gap,
                    aligned_length: aligned_step,
                    sequence_length: sequence_step,
                });
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
    pub(super) fn get_sequence_index(&self, index: usize) -> (usize, u16) {
        let mut path_index = 0;
        let mut aligned_index = 0;
        let mut sequence_index = 0;
        while aligned_index < index {
            sequence_index += self.path[path_index].sequence_length as usize;
            aligned_index += self.path[path_index].aligned_length as usize;
            path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
        }
        (
            sequence_index.min(self.sequence.cast_peptidoform().len()),
            (aligned_index - index) as u16,
        )
    }

    // Sequence index
    pub(super) fn get_aa_at(
        &self,
        sequence_index: usize,
    ) -> (&SequenceElement<Linear>, &Multi<Mass>) {
        let c = self.sequence.cast_peptidoform();
        let i = (sequence_index).min(c.len()).saturating_sub(1);
        unsafe {
            (
                c.sequence().get_unchecked(i),
                self.masses.get_unchecked([i, 0]),
            )
        }
    }

    // Sequence index, sequence len
    pub(super) fn get_block_at(
        &self,
        sequence_index: usize,
        len: usize,
    ) -> Option<(&[SequenceElement<Linear>], &Multi<Mass>)> {
        if sequence_index >= len {
            Some(unsafe {
                (
                    self.sequence
                        .cast_peptidoform()
                        .sequence()
                        .get_unchecked(sequence_index - len..sequence_index),
                    self.masses.get_unchecked([sequence_index - 1, len - 1]),
                )
            })
        } else {
            None
        }
    }

    #[allow(dead_code)] // Used for debugging purposes
    pub(super) fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let mut seq_index = 0;
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

    #[allow(dead_code)]
    fn short(&self) -> String {
        MultiPiece::short(&self.path)
    }
}

impl<const STEPS: u16, Sequence: HasPeptidoform<Linear> + Clone> AlignIndex<STEPS, Sequence> {
    /// Get matrix with the distances between all sequences in this index. It uses
    /// [`mzalign::Alignment::distance`] as metric.
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

    /// Create a MMSA out of the sequences in this index. The maximal distance indicates the
    /// maximal distance between two clusters that can still be joined. If set to None this will
    /// merge all clusters. This can be used to return a set of motifs in otherwise unrelated
    /// sequences.
    pub fn multi_align(
        &self,
        maximal_distance: Option<f64>,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Vec<MultiAlignment<Sequence>> {
        let distance_matrix = self.distance_matrix(scoring, align_type);
        self.multi_align_inner(maximal_distance, distance_matrix, scoring, align_type)
    }

    fn multi_align_inner(
        &self,
        maximal_distance: Option<f64>, /* Maximal distance to join clusters, or join all if set
                                        * to None */
        mut distance_matrix: DiagonalArray<f64, { u16::MAX }>,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Vec<MultiAlignment<Sequence>> {
        if self.sequences.is_empty() {
            return Vec::new();
        } else if self.sequences.len() == 1 {
            return vec![MultiAlignment::new(
                vec![MultiAlignmentLineTemp::single(
                    0,
                    self.sequences[0].0.clone(),
                    &self.sequences[0].1,
                )],
                align_type,
            )];
        }

        let mut nodes = Vec::with_capacity(self.sequences.len());
        for (i, s) in self.sequences.iter().enumerate() {
            nodes.push(vec![MultiAlignmentLineTemp::single(i, s.0.clone(), &s.1)]);
        }

        // Progressively merge (this creates a phylogenetic tree that is not stored in any way)
        while nodes.len() > 1 {
            let (distance, [far_index, close_index]) = distance_matrix.min::<true>();
            if maximal_distance.is_some_and(|m| distance > m) {
                // Stop if this distance is more than the threshold
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
    /// Get matrix with the distances between all sequences in this index. It uses
    /// [`Alignment::distance`] as metric.
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

    /// This parallelizes the distance matrix calculation, the multiple sequence merging is still
    /// done sequentially
    #[cfg(feature = "rayon")]
    pub fn par_multi_align(
        &self,
        maximal_distance: Option<f64>, /* Maximal distance to join clusters, or join all if set
                                        * to None */
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
    let len_a = a.first().map_or(0, MultiAlignmentLineTemp::aligned_length);
    let len_b = b.first().map_or(0, MultiAlignmentLineTemp::aligned_length);
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
    // min_mass is the minimal mass that matches at least one mass in masses_a[index_a - len_a ..
    // index_a]. max_mass is the same for maximal mass.
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
    let sequence_indices_b: Vec<Vec<(usize, u16)>> = (1..=len_b)
        .map(|index_b| b.iter().map(|s| s.get_sequence_index(index_b)).collect::<Vec<_>>())
        .collect();

    let sorted_pair_masses_b: Vec<Vec<((Mass, Mass), usize)>> = sequence_indices_b
        .iter()
        .map(|seq_b_indices| {
            let mut sorted_masses_b: Vec<((Mass, Mass), usize)> = seq_b_indices
                .iter()
                .enumerate()
                .map(|(ib, (is, _))| {
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
        let sequence_indices_a =
            a.iter().map(|s| s.get_sequence_index(index_a)).collect::<Vec<_>>();
        let mut sorted_pair_masses_a: Vec<((Mass, Mass), usize)> = sequence_indices_a
            .iter()
            .enumerate()
            .map(|(ia, (is, _))| {
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
                    // TODO: Remember how much further A or B goes and update the piece
                    let pair_score = score_pair(
                        a[*ia].get_aa_at(sequence_indices_a[*ia].0),
                        b[ib].get_aa_at(sequence_indices_b[index_b - 1][ib].0),
                        scoring,
                        prev.score,
                    );
                    if pair_score.score > highest.score {
                        // if sequence_indices_b[index_b - 1][ib].1 == 0 {
                        //     pair_score.step_a += sequence_indices_a[*ia].1;
                        // }
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
                            .get_unchecked(sequence_indices_a[line_a].0.saturating_sub(1))
                    };
                    for line_b in 0..b.len() {
                        let pair_score = score_pair_mass_mismatch(
                            aa_a,
                            unsafe {
                                b[line_b].sequence.cast_peptidoform().sequence().get_unchecked(
                                    sequence_indices_b[index_b - 1][line_b].0.saturating_sub(1),
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
                    let sequence_index_a = sequence_indices_a[line_a].0;
                    for line_b in 0..b.len() {
                        let sequence_index_b = sequence_indices_b[index_b - 1][line_b].0;
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
                                        .get_block_at(sequence_index_a, len_a)
                                        .zip(b[line_b].get_block_at(sequence_index_b, len_b))
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
    let goal_length = sequences.iter().fold(0, |acc, s| acc.max(s.aligned_length()));
    for s in &mut sequences {
        s.pad(goal_length);
    }

    sequences
}

#[allow(dead_code)]
fn debug_mmsa_path_short(path: &[Piece]) -> String {
    use std::fmt::Write;
    let mut output = String::new();
    let mut print = |piece: Piece, count: usize| {
        if count > 1 {
            let _ = write!(&mut output, "{count}*");
        }
        let _ = write!(&mut output, "{}:{}", piece.step_a, piece.step_b);
        let _ = output.write_char(match piece.match_type {
            MatchType::FullIdentity => '=',
            MatchType::IdentityMassMismatch => 'm',
            MatchType::Mismatch => 'X',
            MatchType::Isobaric => 'i',
            MatchType::Rotation => 'r',
            MatchType::Gap => 'D',
        });
    };

    let mut last: Option<(Piece, usize)> = None;

    for piece in path {
        if let Some((last_piece, count)) = &mut last {
            if last_piece.match_type == piece.match_type
                && last_piece.step_a == piece.step_a
                && last_piece.step_b == piece.step_b
            {
                *count += 1;
            } else {
                print(*last_piece, *count);
                last = Some((*piece, 1));
            }
        } else {
            last = Some((*piece, 1));
        }
    }
    if let Some((piece, count)) = last {
        print(piece, count);
    }

    output
}
