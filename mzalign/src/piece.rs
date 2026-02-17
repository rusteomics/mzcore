use super::scoring::MatchType;
use serde::{Deserialize, Serialize};

/// A piece in an alignment, determining what step was taken in the alignment and how this impacted the score
#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub struct Piece {
    /// The total score of the path up till now
    pub score: isize,
    /// The local contribution to the score of this piece
    pub local_score: isize,
    /// The type of the match
    pub match_type: MatchType,
    /// The number of steps on the first sequence
    pub step_a: u16,
    /// The number of steps on the second sequence
    pub step_b: u16,
}

impl Piece {
    /// Create a new alignment piece
    pub const fn new(
        score: isize,
        local_score: isize,
        match_type: MatchType,
        step_a: u16,
        step_b: u16,
    ) -> Self {
        Self {
            score,
            local_score,
            match_type,
            step_a,
            step_b,
        }
    }

    /// Create a CIGAR like string from a path
    pub(crate) fn cigar(path: &[Self], start_a: usize, start_b: usize) -> String {
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

        let (_, _, output, last) = path.iter().fold(
            (start_a, start_b, String::new(), None),
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
