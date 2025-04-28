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
}
