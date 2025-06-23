//! Functions to generate alignments of peptides based on homology, while taking mass spectrometry errors into account.

use std::borrow::Cow;

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::Deserialize;
use serde::Serialize;

use super::align_type::*;
use super::piece::*;
use super::scoring::*;

use crate::sequence::HasPeptidoform;
use crate::sequence::SemiAmbiguous;
use crate::{
    align::mass_alignment::{determine_final_score, score_pair},
    annotation::model::GlycanModel,
    chemistry::MolecularFormula,
    helper_functions::next_num,
    quantities::Multi,
    sequence::{AtMax, Linear, Peptidoform, SequencePosition, SimpleLinear},
    system::{Mass, Ratio},
};

/// An alignment of two reads. It has either a reference to the two sequences to prevent overzealous use of memory, or if needed use [`Self::to_owned`] to get a variant that clones the sequences and so can be used in more places.
#[derive(Debug, Deserialize, Serialize)]
pub struct Alignment<A, B> {
    /// The first sequence
    pub(super) seq_a: A,
    /// The second sequence
    pub(super) seq_b: B,
    /// The scores of this alignment
    pub(super) score: Score,
    /// The path or steps taken for the alignment
    pub(super) path: Vec<Piece>,
    /// The position in the first sequence where the alignment starts
    pub(super) start_a: usize,
    /// The position in the second sequence where the alignment starts
    pub(super) start_b: usize,
    /// The alignment type
    pub(super) align_type: AlignType,
    /// The maximal step size (the const generic STEPS)
    pub(super) maximal_step: u16,
}

impl<A: Clone, B: Clone> Clone for Alignment<A, B> {
    fn clone(&self) -> Self {
        Self {
            seq_a: self.seq_a.clone(),
            seq_b: self.seq_b.clone(),
            score: self.score,
            path: self.path.clone(),
            start_a: self.start_a,
            start_b: self.start_b,
            align_type: self.align_type,
            maximal_step: self.maximal_step,
        }
    }
}

impl<A: PartialEq, B: PartialEq> PartialEq for Alignment<A, B> {
    fn eq(&self, other: &Self) -> bool {
        self.seq_a == other.seq_a
            && self.seq_b == other.seq_b
            && self.score == other.score
            && self.path == other.path
            && self.start_a == other.start_a
            && self.start_b == other.start_b
            && self.align_type == other.align_type
            && self.maximal_step == other.maximal_step
    }
}

impl<A: std::hash::Hash, B: std::hash::Hash> std::hash::Hash for Alignment<A, B> {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.seq_a.hash(state);
        self.seq_b.hash(state);
        self.score.hash(state);
        self.path.hash(state);
        self.start_a.hash(state);
        self.start_b.hash(state);
        self.align_type.hash(state);
        self.maximal_step.hash(state);
    }
}

impl<A: PartialEq + Eq, B: PartialEq + Eq> Eq for Alignment<A, B> {}

impl<A: Clone, B: Clone> Alignment<Cow<'_, A>, Cow<'_, B>> {
    /// Clone the referenced sequences to make an alignment that owns the sequences.
    /// This can be necessary in some context where the references cannot be guaranteed to stay as long as you need the alignment.
    #[must_use]
    pub fn to_owned(&self) -> Alignment<Cow<'static, A>, Cow<'static, B>> {
        Alignment {
            seq_a: Cow::Owned(self.seq_a.clone().into_owned()),
            seq_b: Cow::Owned(self.seq_b.clone().into_owned()),
            ..self.clone()
        }
    }
}

impl<A: Eq, B: Eq> PartialOrd for Alignment<A, B> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<A: Eq, B: Eq> Ord for Alignment<A, B> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.score.normalised.cmp(&other.score.normalised)
    }
}

impl<'lifetime, A: AtMax<SimpleLinear>, B: AtMax<SimpleLinear>>
    Alignment<&'lifetime Peptidoform<A>, &'lifetime Peptidoform<B>>
{
    /// Recreate an alignment from a path, the path is [`Self::short`].
    #[expect(clippy::missing_panics_doc, clippy::too_many_arguments)]
    pub fn create_from_path(
        seq_a: &'lifetime Peptidoform<A>,
        seq_b: &'lifetime Peptidoform<B>,
        start_a: usize,
        start_b: usize,
        path: &str,
        scoring: AlignScoring,
        align_type: AlignType,
        maximal_step: u16,
    ) -> Option<Self> {
        let mut index = 0;
        let mut steps = Vec::new();
        while index < path.len() {
            if let Some((offset, num)) = next_num(path.as_bytes(), index, false) {
                let num = u16::try_from(num).unwrap();
                index += offset + 1;
                match path.as_bytes()[index - 1] {
                    b'I' => steps.push((MatchType::Gap, 0, num)),
                    b'D' => steps.push((MatchType::Gap, num, 0)),
                    b'X' => steps.push((MatchType::Mismatch, num, num)),
                    b'm' => steps.push((MatchType::IdentityMassMismatch, num, num)),
                    b'=' => steps.push((MatchType::FullIdentity, num, num)),
                    b'r' => steps.push((MatchType::Rotation, num, num)),
                    b'i' => steps.push((MatchType::Isobaric, num, num)),
                    b':' => {
                        if let Some((offset, num2)) = next_num(path.as_bytes(), index, false) {
                            let num2 = u16::try_from(num2).unwrap();
                            index += offset + 1;
                            match path.as_bytes()[index - 1] {
                                b'i' => steps.push((MatchType::Isobaric, num, num2)),
                                _ => return None,
                            }
                        } else {
                            return None;
                        }
                    }
                    _ => return None,
                }
            } else {
                return None;
            }
        }

        let mut index_a = start_a;
        let mut index_b = start_b;
        let mut score = 0;
        let path = steps
            .into_iter()
            .flat_map(|(step, a, b)| match step {
                MatchType::Gap => {
                    let step_a = u16::from(a != 0);
                    let step_b = u16::from(b != 0);
                    index_a += a as usize;
                    index_b += b as usize;

                    (0..a.max(b))
                        .map(|i| {
                            // The first gap will have the gap_start score
                            if i == 0 {
                                let local_score =
                                    scoring.gap_start as isize + scoring.gap_extend as isize;
                                score += local_score;
                                Piece {
                                    score,
                                    local_score,
                                    match_type: MatchType::Gap,
                                    step_a,
                                    step_b,
                                }
                            } else {
                                let local_score = scoring.gap_extend as isize;
                                score += scoring.gap_extend as isize;
                                Piece {
                                    score,
                                    local_score,
                                    match_type: MatchType::Gap,
                                    step_a,
                                    step_b,
                                }
                            }
                        })
                        .collect_vec()
                }
                MatchType::FullIdentity | MatchType::IdentityMassMismatch | MatchType::Mismatch => {
                    (0..a)
                        .map(|_| {
                            let mass_a = seq_a[index_a]
                                .formulas_all(
                                    &[],
                                    &[],
                                    &mut Vec::new(),
                                    false,
                                    SequencePosition::Index(index_a),
                                    0,
                                    0,
                                    &GlycanModel::DISALLOW,
                                )
                                .0
                                .iter()
                                .map(MolecularFormula::monoisotopic_mass)
                                .collect();
                            let mass_b = seq_b[index_b]
                                .formulas_all(
                                    &[],
                                    &[],
                                    &mut Vec::new(),
                                    false,
                                    SequencePosition::Index(index_b),
                                    0,
                                    0,
                                    &GlycanModel::DISALLOW,
                                )
                                .0
                                .iter()
                                .map(MolecularFormula::monoisotopic_mass)
                                .collect();
                            let piece = score_pair(
                                (&seq_a[index_a], &mass_a),
                                (&seq_b[index_b], &mass_b),
                                scoring,
                                score,
                            );
                            index_a += 1;
                            index_b += 1;
                            score += piece.local_score;
                            piece
                        })
                        .collect_vec()
                }
                MatchType::Rotation => {
                    let local_score = scoring.mass_base as isize
                        + scoring.rotated as isize
                            * isize::try_from(a).expect("Could not fit number in isize");
                    score += local_score;
                    index_a += a as usize;
                    index_b += b as usize;
                    vec![Piece {
                        score,
                        local_score,
                        match_type: MatchType::Rotation,
                        step_a: a,
                        step_b: b,
                    }]
                }
                MatchType::Isobaric => {
                    let local_score = scoring.mass_base as isize
                        + scoring.isobaric as isize
                            * isize::try_from(a + b).expect("Could not fit number in isize")
                            / 2;
                    score += local_score;
                    index_a += a as usize;
                    index_b += b as usize;
                    vec![Piece {
                        score,
                        local_score,
                        match_type: MatchType::Isobaric,
                        step_a: a,
                        step_b: b,
                    }]
                }
            })
            .collect_vec();

        Some(Self {
            seq_a,
            seq_b,
            score: determine_final_score(seq_a, seq_b, start_a, start_b, &path, scoring),
            path,
            start_a,
            start_b,
            align_type,
            maximal_step,
        })
    }
}

impl<A, B> Alignment<A, B> {
    /// The first sequence
    pub const fn seq_a(&self) -> &A {
        &self.seq_a
    }
    /// The second sequence
    pub const fn seq_b(&self) -> &B {
        &self.seq_b
    }

    /// The normalised score, normalised for the alignment length and for the used alphabet.
    /// The normalisation is calculated as follows `absolute_score / max_score`.
    pub const fn normalised_score(&self) -> f64 {
        self.score.normalised.0
    }

    /// All three scores for this alignment.
    pub const fn score(&self) -> Score {
        self.score
    }

    /// The path or steps taken for the alignment
    pub fn path(&self) -> &[Piece] {
        &self.path
    }

    /// The position in the sequences where the alignment starts (a, b)
    pub const fn start(&self) -> (usize, usize) {
        (self.start_a, self.start_b)
    }

    /// The position in the first sequence where the alignment starts
    pub const fn start_a(&self) -> usize {
        self.start_a
    }

    /// The position in the second sequence where the alignment starts
    pub const fn start_b(&self) -> usize {
        self.start_b
    }

    /// The alignment type
    pub const fn align_type(&self) -> AlignType {
        self.align_type
    }

    /// The maximal step size (the const generic STEPS)
    pub const fn max_step(&self) -> u16 {
        self.maximal_step
    }

    /// The total number of residues matched on the first sequence
    pub fn len_a(&self) -> usize {
        self.path().iter().map(|p| p.step_a as usize).sum()
    }

    /// The total number of residues matched on the second sequence
    pub fn len_b(&self) -> usize {
        self.path().iter().map(|p| p.step_b as usize).sum()
    }

    /// Returns statistics for this match.
    pub fn stats(&self) -> Stats {
        let (identical, mass_similar, similar, gaps, length) =
            self.path().iter().fold((0, 0, 0, 0, 0), |acc, p| {
                let m = p.match_type;
                (
                    acc.0
                        + usize::from(
                            m == MatchType::IdentityMassMismatch || m == MatchType::FullIdentity,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.1
                        + usize::from(
                            m == MatchType::FullIdentity
                                || m == MatchType::Isobaric
                                || m == MatchType::Rotation,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.2
                        + usize::from(
                            (m == MatchType::IdentityMassMismatch
                                || m == MatchType::FullIdentity
                                || m == MatchType::Mismatch)
                                && p.local_score >= 0,
                        ) * p.step_a.max(p.step_b) as usize,
                    acc.3 + usize::from(m == MatchType::Gap),
                    acc.4 + p.step_a.max(p.step_b) as usize,
                )
            });
        Stats {
            identical,
            mass_similar,
            similar,
            gaps,
            length,
        }
    }
}

impl<A: HasPeptidoform<Linear>, B: HasPeptidoform<Linear>> Alignment<A, B> {
    /// The mass(es) for the matched portion of the first sequence TODO: this assumes no terminal mods
    pub fn mass_a(&self) -> Multi<MolecularFormula> {
        let seq_a = self.seq_a().peptidoform();
        if self.align_type().left.global_a() && self.align_type().right.global_a() {
            seq_a.bare_formulas()
        } else {
            let mut placed_a = vec![false; seq_a.number_of_ambiguous_modifications()];
            seq_a[self.start_a()..self.start_a() + self.len_a()]
                .iter()
                .enumerate()
                .fold(Multi::default(), |acc, (index, s)| {
                    acc * s
                        .formulas_greedy(
                            &mut placed_a,
                            &[],
                            &[],
                            &mut Vec::new(),
                            false,
                            SequencePosition::Index(index),
                            0,
                            0,
                            &GlycanModel::DISALLOW,
                        )
                        .0
                })
        }
    }

    /// The mass(es) for the matched portion of the second sequence
    pub fn mass_b(&self) -> Multi<MolecularFormula> {
        let seq_b = self.seq_b().peptidoform();
        if self.align_type().left.global_b() && self.align_type().right.global_b() {
            seq_b.bare_formulas()
        } else {
            let mut placed_b = vec![false; seq_b.number_of_ambiguous_modifications()];
            seq_b[self.start_b()..self.start_b() + self.len_b()]
                .iter()
                .enumerate()
                .fold(Multi::default(), |acc, (index, s)| {
                    acc * s
                        .formulas_greedy(
                            &mut placed_b,
                            &[],
                            &[],
                            &mut Vec::new(),
                            false,
                            SequencePosition::Index(index),
                            0,
                            0,
                            &GlycanModel::DISALLOW,
                        )
                        .0
                })
        }
    }

    /// Get the mass delta for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[expect(clippy::missing_panics_doc)]
    pub fn mass_difference(&self) -> Mass {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass() - b.monoisotopic_mass())
            .min_by(|a, b| a.abs().value.total_cmp(&b.abs().value))
            .expect("An empty Multi<MolecularFormula>  was detected")
    }

    /// Get the error in ppm for this match, if it is a (partial) local match it will only take the matched amino acids into account.
    /// If there are multiple possible masses for any of the stretches it returns the smallest difference.
    #[expect(clippy::missing_panics_doc)]
    pub fn ppm(&self) -> Ratio {
        self.mass_a()
            .iter()
            .cartesian_product(self.mass_b().iter())
            .map(|(a, b)| a.monoisotopic_mass().ppm(b.monoisotopic_mass()))
            .min_by(|a, b| a.value.total_cmp(&b.value))
            .expect("An empty Multi<MolecularFormula> was detected")
    }

    /// Get a short representation of the alignment in CIGAR like format.
    /// It has three additional classes `{a}(:{b})?(r|i)` and `{a}m` denoting any special step with the given a and b step size, if b is not given it is the same as a.
    /// `r` is rotation, `i` is isobaric, and `m` is identity but mass mismatch (modification).
    pub fn short(&self) -> String {
        #[derive(Eq, PartialEq)]
        enum StepType {
            Insertion,
            Deletion,
            Match,
            Mismatch,
            Massmismatch,
            Special(MatchType, u16, u16),
        }
        impl std::fmt::Display for StepType {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(
                    f,
                    "{}",
                    match self {
                        Self::Insertion => String::from("I"),
                        Self::Deletion => String::from("D"),
                        Self::Match => String::from("="),
                        Self::Mismatch => String::from("X"),
                        Self::Massmismatch => String::from("m"),
                        Self::Special(MatchType::Rotation, a, _) => format!("{a}r"),
                        Self::Special(MatchType::Isobaric, a, b) if a == b => format!("{a}i"),
                        Self::Special(MatchType::Isobaric, a, b) => format!("{a}:{b}i"),
                        Self::Special(..) => panic!("A special match cannot be of this match type"),
                    }
                )
            }
        }
        let (_, _, output, last) = self.path().iter().fold(
            (self.start_a(), self.start_b(), String::new(), None),
            |(a, b, output, last), step| {
                let current_type = match (step.match_type, step.step_a, step.step_b) {
                    (MatchType::Isobaric, a, b) => StepType::Special(MatchType::Isobaric, a, b), // Catch any 1/1 isobaric sets before they are counted as Match/Mismatch
                    (_, 0, 1) => StepType::Insertion,
                    (_, 1, 0) => StepType::Deletion,
                    (MatchType::IdentityMassMismatch, 1, 1) => StepType::Massmismatch,
                    (MatchType::FullIdentity, 1, 1) => StepType::Match,
                    (MatchType::Mismatch, 1, 1) => StepType::Mismatch,
                    (m, a, b) => StepType::Special(m, a, b),
                };
                let (str, last) = match last {
                    Some((t @ StepType::Special(..), _)) => {
                        (format!("{output}{t}"), Some((current_type, 1)))
                    }
                    Some((t, n)) if t == current_type => (output, Some((t, n + 1))),
                    Some((t, n)) => (format!("{output}{n}{t}"), Some((current_type, 1))),
                    None => (output, Some((current_type, 1))),
                };
                (
                    a + step.step_a as usize,
                    b + step.step_b as usize,
                    str,
                    last,
                )
            },
        );
        match last {
            Some((t @ StepType::Special(..), _)) => format!("{output}{t}"),
            Some((t, n)) => format!("{output}{n}{t}"),
            _ => output,
        }
    }
}

/// Statistics for an alignment with some helper functions to easily retrieve the number of interest.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Stats {
    /// The total number of identical positions
    pub identical: usize,
    /// The total number of mass similar positions, so including isobaric and rotations
    pub mass_similar: usize,
    /// The total number of similar positions, where the scoring matrix scores above 0, and not isobaric
    pub similar: usize,
    /// The total number of gap positions
    pub gaps: usize,
    /// The length of the alignment, the sum of the max of the step for A and B for each position.
    pub length: usize,
}

impl Stats {
    /// Get the identity as fraction.
    pub fn identity(self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.identical as f64 / self.length as f64
        }
    }

    /// Get the mass similarity as fraction.
    pub fn mass_similarity(self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.mass_similar as f64 / self.length as f64
        }
    }

    /// Get the similarity as fraction.
    pub fn similarity(self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.similar as f64 / self.length as f64
        }
    }

    /// Get the gaps as fraction.
    pub fn gaps_fraction(self) -> f64 {
        if self.length == 0 {
            0.0
        } else {
            self.gaps as f64 / self.length as f64
        }
    }
}

/// The score of an alignment
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Score {
    /// The normalised score (absolute / max)
    pub normalised: OrderedFloat<f64>,
    /// The absolute score
    pub absolute: isize,
    /// The maximal possible score, the average score of the sequence slices on sequence a and b if they were aligned to themself, rounded down.
    ///    Think of it like this: `align(sequence_a.sequence[start_a..len_a], sequence_a.sequence[start_a..len_a])`.
    pub max: isize,
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        align::{AlignScoring, AlignType, align},
        chemistry::MultiChemical,
        sequence::{AminoAcid, Peptidoform, SimpleLinear},
    };

    #[test]
    fn mass_difference() {
        // Test if the mass difference calculation is correct for some harder alignments.
        // A has an ambiguous AA, B, and C have the two options, while D has a sub peptide of A.
        let a = Peptidoform::pro_forma("AABAA", None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let b = Peptidoform::pro_forma("AANAA", None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let c = Peptidoform::pro_forma("AADAA", None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let d = Peptidoform::pro_forma("ADA", None)
            .unwrap()
            .into_simple_linear()
            .unwrap();

        assert!(
            align::<1, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
                &a,
                &b,
                AlignScoring::default(),
                AlignType::GLOBAL
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        assert!(
            align::<1, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
                &a,
                &c,
                AlignScoring::default(),
                AlignType::GLOBAL
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        assert!(
            align::<1, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
                &a,
                &d,
                AlignScoring::default(),
                AlignType::GLOBAL_B
            )
            .mass_difference()
            .value
            .abs()
                < f64::EPSILON
        );
        let mass_diff_nd = (AminoAcid::Asparagine.formulas()[0].monoisotopic_mass()
            - AminoAcid::AsparticAcid.formulas()[0].monoisotopic_mass())
        .value
        .abs();
        let mass_diff_bc = align::<1, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
            &b,
            &c,
            AlignScoring::default(),
            AlignType::GLOBAL_B,
        )
        .mass_difference()
        .value
        .abs();
        assert!(
            (mass_diff_bc - mass_diff_nd).abs() < 1E-10,
            "{mass_diff_bc} (peptides) should be equal to {mass_diff_nd} (ND)"
        );
    }
}
