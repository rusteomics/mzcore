#![allow(dead_code)]
use std::borrow::Cow;

use crate::Peptidoform;

use super::{AlignType, MatchType, Score};

use serde::{Deserialize, Serialize};

type MultiAlignment<'lifetime, Complexity> = Vec<MultiAlignmentLine<'lifetime, Complexity>>;

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
struct MultiAlignmentLine<'lifetime, Complexity> {
    sequence: Cow<'lifetime, Peptidoform<Complexity>>,
    path: Vec<MultiPiece>,
    score: Score,
    start: usize,
    align_type: AlignType,
    maximal_step: u16,
}

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
struct MultiPiece {
    score: isize,
    local_score: isize,
    match_type: MatchType,
    step: u16,
}

impl<Complexity> MultiAlignmentLine<'_, Complexity> {
    fn debug_display(&self) {
        for piece in self
            .path
            .iter()
            .zip(self.sequence.sequence().iter().skip(self.start))
        {
            print!(
                "{}{}",
                piece.1.aminoacid,
                "·".repeat(piece.0.step as usize - 1)
            );
        }
        println!();
    }
}
