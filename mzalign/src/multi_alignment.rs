use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::{
    AlignIndex, AlignScoring, AlignType, MatchType, Piece, Score,
    align_matrix::Matrix,
    diagonal_array::DiagonalArray,
    mass_alignment::{align_cached, calculate_masses, score, score_pair},
};
use mzcore::{
    prelude::*,
    quantities::{Multi, Tolerance},
    sequence::{HasPeptidoform, Linear},
    system::{Mass, OrderedMass, dalton},
};

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignment<Sequence> {
    lines: Vec<MultiAlignmentLine<Sequence>>,
    score: Score,
    maximal_step: u16,
    align_type: AlignType,
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
        align_type: AlignType,
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

        for ref_index in 0..self
            .lines
            .iter()
            .map(|l| l.ref_length())
            .max()
            .unwrap_or_default()
        {
            // Break when the last item is reached
            let mut element = BTreeMap::new();
            for line in &self.lines {
                if let Some(item) = line.get_item(ref_index) {
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
    // Outer 'ref' index
    fn get_item(&self, index: usize) -> Option<(SequenceElement<Linear>, u16)> {
        let mut path_index = 0;
        let mut ref_index = 0;
        let mut seq_index = 0;
        loop {
            seq_index += self.path[path_index].seq_step as usize;
            ref_index += self.path[path_index].ref_step as usize;
            if ref_index >= index + self.path[path_index].ref_step as usize {
                break;
            }
            path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
        }
        let seq_index = seq_index
            .min(self.sequence.cast_peptidoform().len())
            .saturating_sub(1);
        (ref_index == index + self.path[path_index].ref_step as usize
            && self.path[path_index].match_type != MatchType::Gap)
            .then(|| {
                (
                    self.sequence.cast_peptidoform().sequence()[seq_index].clone(),
                    self.path[path_index].ref_step,
                )
            })
    }

    fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let mut seq_index = self.start;
        for piece in &self.path {
            if piece.seq_step == 0 {
                write!(w, "{}", "-".repeat(piece.ref_step as usize)).unwrap();
            } else {
                let subseq = &sequence[seq_index..seq_index + piece.seq_step as usize];
                // Obviously misses mods now
                let display = subseq
                    .iter()
                    .map(|s| s.aminoacid.one_letter_code().unwrap_or('X'))
                    .collect::<String>();
                write!(
                    w,
                    "{display}{}",
                    "·".repeat((piece.ref_step - piece.seq_step) as usize)
                )
                .unwrap();
                seq_index += piece.seq_step as usize;
            }
        }
        writeln!(w).unwrap();
    }

    fn ref_length(&self) -> usize {
        let mut ref_index = 0;
        for piece in &self.path {
            ref_index += piece.ref_step as usize;
        }
        ref_index
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
                    ref_step: 1,
                    seq_step: 1,
                };
                len
            ],
            start: 0,
        }
    }

    /// Update this line with the given merged alignment, this does not add end gaps
    fn update(mut self, path: &[Piece], is_a: bool, start_self: u16, start_other: u16) -> Self {
        // println!(
        //     "S: {}, P: {}, a: {is_a} ss: {start_self}, so: {start_other}",
        //     self.sequence.cast_peptidoform(),
        //     Piece::cigar(path, 0, 0),
        // );
        let mut path_index = 0;
        let mut ref_index = 0;
        let mut goal = 0;
        let mut internal_index = 0;
        if start_other > 0 {
            if self.path[0].match_type == MatchType::Gap {
                self.path[0].ref_step += start_other;
                internal_index = start_other;
            } else {
                self.path.insert(
                    path_index,
                    MultiPiece {
                        match_type: MatchType::Gap,
                        ref_step: start_other,
                        seq_step: 0,
                    },
                );
                path_index += 1;
            }
        }
        for piece in path {
            let ref_step = if is_a { piece.step_b } else { piece.step_a };
            let seq_step = if is_a { piece.step_a } else { piece.step_b };
            // println!(
            //     "r{ref_step} s{seq_step} ri{ref_index} pi{path_index} ii{internal_index} {:?}",
            //     piece.match_type
            // );
            if seq_step == 0 {
                self.path.insert(
                    path_index,
                    MultiPiece {
                        match_type: MatchType::Gap,
                        ref_step,
                        seq_step,
                    },
                );
                ref_index += ref_step;
                path_index += 1;
            } else {
                self.path[path_index].ref_step += ref_step.saturating_sub(seq_step);
                goal += ref_step.max(seq_step);
                while ref_index + internal_index < goal {
                    ref_index += self.path[path_index].ref_step - internal_index;
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

    fn ref_length(&self) -> usize {
        let mut ref_index = 0;
        for piece in &self.path {
            ref_index += piece.ref_step as usize;
        }
        ref_index
    }

    /// Pad to length
    fn pad(&mut self, goal_length: usize) {
        let mut ref_index = 0;
        for piece in &self.path {
            ref_index += piece.ref_step as usize;
        }
        if ref_index < goal_length {
            let last = self.path.len() - 1;
            let dif = goal_length - ref_index;
            if self.path[last].match_type == MatchType::Gap {
                self.path[last].ref_step += dif as u16;
            } else {
                self.path.push(MultiPiece {
                    match_type: MatchType::Gap,
                    ref_step: dif as u16,
                    seq_step: 0,
                });
            }
        }
    }

    // Outer 'ref' index
    fn get_seq_index(&self, index: usize) -> usize {
        let mut path_index = 0;
        let mut ref_index = 0;
        let mut seq_index = 0;
        while ref_index < index {
            seq_index += self.path[path_index].seq_step as usize;
            ref_index += self.path[path_index].ref_step as usize;
            path_index += usize::from(path_index != self.path.len() - 1); // just saturate for now
        }
        seq_index.min(self.sequence.cast_peptidoform().len())
    }

    // Seq index
    fn get_aa_at(&self, seq_index: usize) -> (&SequenceElement<Linear>, &Multi<Mass>) {
        let c = self.sequence.cast_peptidoform();
        let i = (seq_index).min(c.len()).saturating_sub(1);
        unsafe {
            (
                &c.sequence().get_unchecked(i),
                &self.masses.get_unchecked([i, 0]),
            )
        }
    }

    // Seq index, seq len
    fn get_block_at(
        &self,
        seq_index: usize,
        len: usize,
    ) -> Option<(&[SequenceElement<Linear>], &Multi<Mass>)> {
        if seq_index >= len {
            Some(unsafe {
                (
                    &self
                        .sequence
                        .cast_peptidoform()
                        .sequence()
                        .get_unchecked(seq_index - len..seq_index),
                    &self.masses.get_unchecked([seq_index - 1, len - 1]),
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
            if piece.seq_step == 0 {
                write!(w, "{}", "-".repeat(piece.ref_step as usize)).unwrap();
            } else {
                let subseq = &sequence[seq_index..seq_index + piece.seq_step as usize];
                // Obviously misses mods now
                let display = subseq
                    .iter()
                    .map(|s| s.aminoacid.one_letter_code().unwrap_or('X'))
                    .collect::<String>();
                write!(
                    w,
                    "{display}{}",
                    "·".repeat((piece.ref_step - piece.seq_step) as usize)
                )
                .unwrap();
                seq_index += piece.seq_step as usize;
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
    /// Ref seq is required to always be at least `seq_step`
    ref_step: u16,
    seq_step: u16,
}

impl<const STEPS: u16, Sequence: HasPeptidoform<Linear> + Clone> AlignIndex<STEPS, Sequence> {
    pub fn multi_align(
        &self,
        maximal_distance: Option<f64>, // Maximal distance to join clusters, or join all if set to None
        scoring: AlignScoring<'_>,
        align_type: AlignType, // TODO: figure out what to do for align_type
    ) -> Vec<MultiAlignment<Sequence>> {
        assert!(self.sequences.len() > 1);
        // Create an outgroup of 'random' sequence with the length of the average of all used sequences

        // Make distance matrix
        let mut matrix = DiagonalArray::<f64, { u16::MAX }>::new(self.sequences.len());
        let mut nodes = Vec::with_capacity(self.sequences.len());
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
            nodes.push(vec![MultiAlignmentLineTemp::single(i, s.0.clone(), &s.1)]);
        }
        // Remember that the last node is not made (the outgroup node) so do not try to merge that one

        // Create tree
        // It might work out to now rely on the fact that the outgroup is always the furthest distance from all other nodes that we can just immediately start merging the nodes to produce the multi alignment

        // Progressively merge
        while nodes.len() > 1 {
            // println!("{}", matrix.to_csv());

            let (distance, [far_index, close_index]) = matrix.min::<true>();
            if maximal_distance.is_some_and(|m| distance > m) {
                // Stop if this distance is more than the threshold, or if the outgroup would be merged
                break;
            }
            // Merge nodes (this is where the actual alignment happens)
            let close_node = std::mem::take(&mut nodes[close_index]);
            let far_node = nodes.remove(far_index);
            nodes[close_index] =
                multi_align_cached::<STEPS, Sequence>(close_node, far_node, scoring, align_type);

            // for l in &nodes[close_index] {
            //     let mut buf = String::new();
            //     l.debug_display(&mut buf);
            //     print!("{buf}");
            // }

            for line in &nodes[close_index] {
                let mut buf = String::new();
                line.debug_display(&mut buf);
                // print!("{buf}");
            }

            // Update matrix to merge these columns
            matrix = matrix.merge_columns(close_index, far_index, distance);
            // println!(
            //     "Merged {close_index} with {far_index} (left: nodes {}, matrix {})",
            //     nodes.len(),
            //     matrix.len(),
            // );
        }

        nodes
            .into_iter()
            .map(|seqs| MultiAlignment::new(seqs, align_type))
            .collect()
    }
}

fn mass_range_expanded(masses: &Multi<Mass>, tolerance: Tolerance<OrderedMass>) -> (Mass, Mass) {
    masses.iter().fold(
        (
            Mass::new::<dalton>(f64::INFINITY),
            Mass::new::<dalton>(f64::NEG_INFINITY),
        ),
        |(min, max), &m| {
            let range = tolerance.bounds(m);
            (min.min(range.0), max.max(range.1))
        },
    )
}

fn mass_range(masses: &Multi<Mass>) -> (Mass, Mass) {
    masses.iter().fold(
        (
            Mass::new::<dalton>(f64::INFINITY),
            Mass::new::<dalton>(f64::NEG_INFINITY),
        ),
        |(min, max), &m| (min.min(m), max.max(m)),
    )
}

/// Do mass based alignment, but with precomputed masses, see [align] for the normal entry point variant.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger than [`isize::MAX`].
#[expect(clippy::too_many_lines)]
#[allow(clippy::similar_names)]
pub(super) fn multi_align_cached<'a, const STEPS: u16, Sequence: HasPeptidoform<Linear>>(
    mut a: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    mut b: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
) -> Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>> {
    // println!("Max steps: {STEPS}");
    let len_a = a
        .iter()
        .fold(0, |acc, i| acc.max(i.sequence.cast_peptidoform().len()));
    let len_b = b
        .iter()
        .fold(0, |acc, i| acc.max(i.sequence.cast_peptidoform().len()));
    let mut matrix = Matrix::new(len_a, len_b);
    let mut global_highest = (0, 0, 0);

    if align_type.left.global_a() {
        matrix.global_start(true, scoring);
    }
    if align_type.left.global_b() {
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
    let ranges_a = a
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

    let ranges_b = b
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

    for index_a in 1..=len_a {
        let seq_a_indices = a
            .iter()
            .map(|s| s.get_seq_index(index_a))
            .collect::<Vec<_>>();
        for index_b in 1..=len_b {
            let seq_b_indices = b
                .iter()
                .map(|s| s.get_seq_index(index_b))
                .collect::<Vec<_>>();
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

            for ia in 0..a.len() {
                for ib in 0..b.len() {
                    let pair_score = score_pair(
                        a[ia].get_aa_at(seq_a_indices[ia]),
                        b[ib].get_aa_at(seq_b_indices[ib]),
                        scoring,
                        prev.score,
                    );
                    if pair_score.score > highest.score {
                        highest = pair_score;
                    }
                }
            }

            if highest.match_type != MatchType::FullIdentity {
                for ia in 0..a.len() {
                    let seq_index_a = seq_a_indices[ia];
                    for ib in 0..b.len() {
                        let seq_index_b = seq_b_indices[ib];
                        // Now try matching longer sequences.
                        for len_a in 1..=index_a.min(STEPS as usize) {
                            let range_a = unsafe {
                                ranges_a[ia]
                                    .get_unchecked([seq_index_a.saturating_sub(1), len_a - 1])
                            };

                            let min_len_b = if len_a == 1 { 2 } else { 1 };

                            for len_b in min_len_b..=index_b.min(STEPS as usize) {
                                let range_b = unsafe {
                                    ranges_b[ib]
                                        .get_unchecked([seq_index_b.saturating_sub(1), len_b - 1])
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
                                    a[ia]
                                        .get_block_at(seq_index_a, len_a)
                                        .and_then(|a| {
                                            b[ib].get_block_at(seq_index_b, len_b).map(|b| (a, b))
                                        })
                                        .and_then(|(a, b)| score(a, b, scoring, base_score))
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
            if align_type.left.global() || highest.score > 0 {
                unsafe {
                    *matrix.get_unchecked_mut([index_a, index_b]) = highest;
                }
            }
        }
    }
    let (start_a, start_b, path) = matrix.trace_path(align_type, global_highest);

    let mut sequences = Vec::with_capacity(a.len() + b.len());

    for seq in a {
        sequences.push(seq.update(&path, true, start_a as u16, start_b as u16));
    }

    for seq in b {
        sequences.push(seq.update(&path, false, start_b as u16, start_a as u16));
    }

    // Fix end gaps
    let goal_length = sequences.iter().fold(0, |acc, s| acc.max(s.ref_length()));
    for s in &mut sequences {
        s.pad(goal_length);
    }

    sequences
}

#[allow(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use mzcore::ontology::STATIC_ONTOLOGIES;

    #[test]
    fn simple() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap()
        };
        let sequences = vec![seq("AGGWHD"), seq("ANWHN[Deamidated]"), seq("AHYDH")];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), AlignType::GLOBAL);
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
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap()
        };
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
        let alignment = index.multi_align(None, AlignScoring::default(), AlignType::GLOBAL);
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
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap()
        };
        let sequences = vec![
            seq("WRGGDGFYAMDYWGQG"),
            seq("RWGDGFYAMYDWGQG"),
            seq("RWNDGYAMEDYWGQG"),
            seq("RWGGDGFMDYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), AlignType::GLOBAL);
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
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap()
        };
        let sequences = vec![
            seq("WRGGDGFYAMDYWGQG"),
            seq("RWGDGFYAMYD"),
            seq("RWNDGYAMDYWGQG"),
            seq("GDGFMDYWGQG"),
        ];
        let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), AlignType::EITHER_GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "WRGGDGFYAMDYWGQG\n---GDGF--MDYWGQG\nRWN·DG-YAMDYWGQG\nRWG-DGFYAMYD----\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn either_global_2() {
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap()
        };
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
        let alignment = index.multi_align(None, scoring, AlignType::EITHER_GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        println!("{buf}");
        buf = buf.split('\n').skip(1).join("\n");
        let expected = "CSRWRGGDGF---------\n---WRN·DGFYAM------\n---RWGGDGFYAMDYWG--\n----------YW·DYWGQG\n---RWGGNGFYW·DYWGQG\n";
        assert_eq!(buf, expected, "Expected:\n{expected}");
    }

    #[test]
    fn properties() {
        let sequence = Peptidoform::pro_forma("AGGWHD", &STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_linear()
            .unwrap();
        let masses = calculate_masses::<4>(&sequence, MassMode::Monoisotopic);
        let mut line = MultiAlignmentLineTemp::single(0, &sequence, &masses);
        assert_eq!(line.get_seq_index(4), 4);
        assert_eq!(line.get_seq_index(5), 5);
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
            0,
        );
        line.path[3].ref_step = 2; // The W is bigger
        assert_eq!(updated, line);
        assert_eq!(line.get_seq_index(4), 4);
        assert_eq!(line.get_seq_index(5), 4);
        assert_eq!(line.get_seq_index(6), 5);
        assert_eq!(
            line.get_aa_at(6),
            (
                &SequenceElement::new(AminoAcid::Histidine.into(), None).into(),
                &Multi::from(CheckedAminoAcid::Histidine.formula().monoisotopic_mass())
            )
        );
        assert_eq!(
            line.get_block_at(6, 2),
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
        let sequence = Peptidoform::pro_forma("WRGGDGFYAMDYWGQG", &STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_linear()
            .unwrap();
        let sequence2 = Peptidoform::pro_forma("RWNDGYAMEDYWGQG", &STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_linear()
            .unwrap();
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
        let updated = line.update(path, true, 0, 0);
        let updated2 = line2.update(path, false, 0, 0);
        let mut buf1 = String::new();
        updated.debug_display(&mut buf1);
        let mut buf2 = String::new();
        updated2.debug_display(&mut buf2);
        assert_eq!(buf1, "WRGGDGFYAM-DYWGQG\n");
        assert_eq!(buf2, "RWN·DG-YAMEDYWGQG\n");
    }
}
