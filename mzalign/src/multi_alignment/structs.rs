use std::{borrow::Cow, collections::BTreeMap};

use mzcore::{
    prelude::*,
    sequence::{HasPeptidoform, Linear},
};
use serde::{Deserialize, Serialize};

use crate::{
    AlignScoring, AlignType, MatchType, Score,
    helper_functions::next_num,
    multi_alignment::calculate::{MultiAlignmentLineTemp, multi_align_cached},
};

/// A mass-based multiple sequence alignment (MMSA).
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignment<Sequence> {
    lines: Vec<MultiAlignmentLine<Sequence>>,
    score: Score,
    maximal_step: u16,
    align_type: MultiAlignType,
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignment<Sequence> {
    pub(super) fn debug_display(&self, mut w: impl std::fmt::Write) {
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

    pub(super) fn new<const STEPS: u16>(
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
            .map(MultiAlignmentLine::aligned_length)
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

    /// Get the lines
    pub fn iter(&self) -> std::slice::Iter<'_, MultiAlignmentLine<Sequence>> {
        self.lines.iter()
    }

    /// Get the score (TODO: calculate)
    pub const fn score(&self) -> Score {
        self.score
    }

    /// Get the maximal step
    pub const fn maximal_step(&self) -> u16 {
        self.maximal_step
    }

    /// Get the align type
    pub const fn align_type(&self) -> MultiAlignType {
        self.align_type
    }

    /// Combine two `MultiAlignments` into one
    #[must_use]
    pub fn combine<const STEPS: u16>(
        self,
        other: Self,
        mass_mode: MassMode,
        scoring: AlignScoring<'_>,
        align_type: MultiAlignType,
    ) -> Self {
        let temp_self = self.lines.into_iter().map(|l| l.into_temp::<STEPS>(mass_mode)).collect();
        let temp_other = other.lines.into_iter().map(|l| l.into_temp::<STEPS>(mass_mode)).collect();

        Self::new::<STEPS>(
            multi_align_cached(temp_self, temp_other, scoring, align_type),
            align_type,
        )
    }
}

impl<'a, Sequence: HasPeptidoform<Linear>> IntoIterator for &'a MultiAlignment<Sequence> {
    type IntoIter = std::slice::Iter<'a, MultiAlignmentLine<Sequence>>;
    type Item = &'a MultiAlignmentLine<Sequence>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// For each location the possible sequences and their length, their depth of coverage (how often
/// seen) and their average local confidence.
pub type SequenceVariance = Vec<BTreeMap<(SequenceElement<Linear>, u16), (usize, f64)>>;

impl<Sequence: HasPeptidoform<Linear>> std::fmt::Display for MultiAlignment<Sequence> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.debug_display(f);
        Ok(())
    }
}

/// The alignment of a single peptidoform in an MMSA.
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignmentLine<Sequence> {
    original_index: usize,
    sequence: Sequence,
    path: Vec<MultiPiece>,
}

impl<Sequence> MultiAlignmentLine<Sequence> {
    /// Get the index into the original list of the MMSA
    pub const fn original_index(&self) -> usize {
        self.original_index
    }

    /// Get the sequence
    pub const fn sequence(&self) -> &Sequence {
        &self.sequence
    }

    /// Get the path
    pub fn path(&self) -> &[MultiPiece] {
        &self.path
    }

    /// Get a short representation of the alignment in CIGAR like format.
    /// It contains the aligned length, and a character for the match
    /// type (same as for the normal alignments). If it is an Isobaric or Rotation step the sequence
    /// length is written as `:n` after the aligned length, for all other steps the sequence length
    /// is constant and can be seen in the table below. It also does a simple run length encoding by
    /// grouping any identical steps as `n*(step)`, if this is not written this step happened once.
    ///
    /// | Match type | Symbol | Sequence length |
    /// | --- | --- | --- |
    /// | [`MatchType::FullIdentity`] | `=` | 1 |
    /// | [`MatchType::IdentityMassMismatch`] | `m` | 1 |
    /// | [`MatchType::Mismatch`] | `X` | 1|
    /// | [`MatchType::Isobaric`] | `i` | `:n` |
    /// | [`MatchType::Rotation`] | `r` | `:n` |
    /// | [`MatchType::Gap`] | `D` | 0 |
    ///
    /// Example line: `2*1=1:2i1=6D`
    pub fn short(&self) -> String {
        MultiPiece::short(&self.path)
    }

    /// Get the aligned length (pad start, length, pad end)
    pub fn placement(&self) -> (u16, usize, u16) {
        let mut total_length = 0;
        for piece in &self.path {
            total_length += piece.aligned_length as usize;
        }

        let start_pad = if self.path.len() > 1 && self.path[0].match_type == MatchType::Gap {
            self.path[0].aligned_length
        } else {
            0
        };
        let end_pad =
            if self.path.len() > 1 && self.path[self.path.len() - 1].match_type == MatchType::Gap {
                self.path[self.path.len() - 1].aligned_length
            } else {
                0
            };
        (
            start_pad,
            total_length - start_pad as usize - end_pad as usize,
            end_pad,
        )
    }
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignmentLine<Sequence> {
    // Get the sequence element at the aligned index. Only returns something if the start of a step
    // is selected.
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

    fn aligned_length(&self) -> usize {
        let mut aligned_index = 0;
        for piece in &self.path {
            aligned_index += piece.aligned_length as usize;
        }
        aligned_index
    }

    /// Create a temporary multi alignment line again (by calculating the masses)
    fn into_temp<const STEPS: u16>(
        self,
        mass_mode: MassMode,
    ) -> MultiAlignmentLineTemp<'static, Sequence, STEPS> {
        let masses =
            crate::mass_alignment::calculate_masses(self.sequence.cast_peptidoform(), mass_mode);
        MultiAlignmentLineTemp {
            original_index: self.original_index,
            sequence: self.sequence,
            path: self.path,
            masses: Cow::Owned(masses),
        }
    }

    /// Recreate an alignment from a path, the path is [`Self::short`]. It returns `None` when the
    /// path is invalid, or the length of the path does not match the given sequence.
    pub fn create_from_path(sequence: Sequence, path: &str, original_index: usize) -> Option<Self> {
        let mut index = 0;
        let mut parsed_path = Vec::new();
        while index < path.len() {
            let mut count = None;
            let aligned_length;
            let mut sequence_length = None;
            let (offset, num) = next_num(path.as_bytes(), index)?;
            index += offset;
            if path[index..].starts_with('*') {
                count = Some(num);
                let (offset, num) = next_num(path.as_bytes(), index)?;
                aligned_length = num;
                index += offset;
            } else {
                aligned_length = num;
            }
            if path[index..].starts_with(':') {
                let (offset, num) = next_num(path.as_bytes(), index)?;
                sequence_length = Some(num);
                index += offset;
            }
            let match_type = match path.as_bytes()[index] {
                b'=' => MatchType::FullIdentity,
                b'm' => MatchType::IdentityMassMismatch,
                b'X' => MatchType::Mismatch,
                b'i' => MatchType::Isobaric,
                b'r' => MatchType::Rotation,
                b'D' => MatchType::Gap,
                _ => return None,
            };
            index += 1;
            parsed_path.extend(std::iter::repeat_n(
                MultiPiece {
                    match_type,
                    aligned_length,
                    sequence_length: sequence_length.or_else(|| {
                        if match_type == MatchType::Isobaric || match_type == MatchType::Rotation {
                            None
                        } else {
                            Some(u16::from(match_type != MatchType::Gap))
                        }
                    })?,
                },
                count.map_or(1, |c| c as usize),
            ));
        }

        if sequence.cast_peptidoform().len()
            == parsed_path.iter().map(|m| m.sequence_length as usize).sum::<usize>()
        {
            Some(Self {
                sequence,
                original_index,
                path: parsed_path,
            })
        } else {
            None
        }
    }
}

/// How a single piece of sequence if aligned in an MMSA, analogous to [Piece] from a pairwise
/// alignment.
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize,
)]
pub struct MultiPiece {
    /// How this piece was matched, for now can only be [`MatchType::FullIdentity`] or
    /// [`MatchType::Gap`]
    pub match_type: MatchType,
    /// How long this piece of sequence is stretched to, this is required to always be at least
    /// `sequence_length`
    pub aligned_length: u16,
    /// The number of sequence elements in this step
    pub sequence_length: u16,
}

impl MultiPiece {
    pub(super) fn short(path: &[Self]) -> String {
        use std::fmt::Write;
        let mut output = String::new();
        let mut print = |piece: Self, count: usize| {
            if count > 1 {
                let _ = write!(&mut output, "{count}*");
            }
            let _ = write!(&mut output, "{}", piece.aligned_length);
            if piece.match_type == MatchType::Isobaric || piece.match_type == MatchType::Rotation {
                let _ = write!(&mut output, ":{}", piece.sequence_length);
            }
            let _ = output.write_char(match piece.match_type {
                MatchType::FullIdentity => '=',
                MatchType::IdentityMassMismatch => 'm',
                MatchType::Mismatch => 'X',
                MatchType::Isobaric => 'i',
                MatchType::Rotation => 'r',
                MatchType::Gap => 'D',
            });
        };

        let mut last: Option<(Self, usize)> = None;

        for piece in path {
            if let Some((last_piece, count)) = &mut last {
                if *last_piece == *piece {
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
}

/// How to align sequences in an MMSA
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MultiAlignType {
    /// What are the requirements for the left side
    pub left: MultiAlignSide,
    /// What are the requirements for the right side
    pub right: MultiAlignSide,
}

impl MultiAlignType {
    /// Both side either global
    pub const EITHER_GLOBAL: Self = Self {
        left: MultiAlignSide::EitherGlobal,
        right: MultiAlignSide::EitherGlobal,
    };
    /// Both side global
    pub const GLOBAL: Self = Self {
        left: MultiAlignSide::Global,
        right: MultiAlignSide::Global,
    };
}

/// In an MMSA the sides can either be global, meaning that all peptidoforms have to end at the
/// same time or either global meaning that ragged edges are allowed, but diverging ends are not
/// allowed.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum MultiAlignSide {
    /// A global alignment, 'straight' edge
    Global,
    /// An either global alignment, 'ragged' edge
    EitherGlobal,
}

impl From<MultiAlignType> for AlignType {
    fn from(value: MultiAlignType) -> Self {
        Self {
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
