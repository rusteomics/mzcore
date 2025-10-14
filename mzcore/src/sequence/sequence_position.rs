use serde::{Deserialize, Serialize};

use crate::parse_json::{ParseJson, use_serde};

/// A position on a sequence
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, serde::Serialize, serde::Deserialize,
)]
pub enum SequencePosition {
    /// N-terminal
    NTerm,
    /// An amino acid at the given index
    Index(usize),
    /// C-terminal
    CTerm,
}

/// Add to the index, the onus of making sure the index is still valid for the peptide is on the caller.
impl std::ops::Add<u8> for SequencePosition {
    type Output = Self;
    fn add(self, rhs: u8) -> Self::Output {
        match self {
            Self::Index(i) => Self::Index(i.saturating_add(rhs as usize)),
            n => n,
        }
    }
}

/// Subtract from the index, the onus of making sure the index is still valid for the peptide is on the caller.
impl std::ops::Sub<u8> for SequencePosition {
    type Output = Self;
    fn sub(self, rhs: u8) -> Self::Output {
        match self {
            Self::Index(i) => Self::Index(i.saturating_sub(rhs as usize)),
            n => n,
        }
    }
}

impl Default for SequencePosition {
    fn default() -> Self {
        Self::Index(0)
    }
}

impl std::fmt::Display for SequencePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NTerm => write!(f, "N-terminal"),
            Self::Index(index) => write!(f, "{index}"),
            Self::CTerm => write!(f, "C-terminal"),
        }
    }
}

impl SequencePosition {
    /// Reverse this position, if the peptide would be reversed what would this location be in that reversed peptide.
    #[must_use]
    pub const fn reverse(self, peptide_length: usize) -> Self {
        match self {
            Self::NTerm => Self::CTerm,
            Self::Index(i) => Self::Index(peptide_length - i),
            Self::CTerm => Self::NTerm,
        }
    }

    /// Convert an index to a sequence position.
    /// The index is defined as follows:
    /// * `0` is N term
    /// * `1..=peptide_length` is in the sequence
    /// * `peptide_length + 1` is C term
    ///
    /// # Panics
    /// Anything outside of this range will panic.
    pub fn from_index(index: usize, peptide_length: usize) -> Self {
        match index {
            0 => Self::NTerm,
            c if c == peptide_length + 1 => Self::CTerm,
            i => {
                if i <= peptide_length {
                    Self::Index(i - 1)
                } else {
                    panic!(
                        "Index {index} it outside of range for a peptide of length {peptide_length}"
                    )
                }
            }
        }
    }
}

impl ParseJson for SequencePosition {
    fn from_json_value(
        value: serde_json::Value,
    ) -> Result<Self, context_error::BoxedError<'static, context_error::BasicKind>> {
        use_serde(value)
    }
}

/// The definition of the position of an ion
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
pub struct PeptidePosition {
    /// The sequence index (0 based into the peptide sequence)
    pub sequence_index: SequencePosition,
    /// The series number (1 based from the ion series terminal)
    pub series_number: usize,
    /// The length of the whole sequence
    pub sequence_length: usize,
}

impl PeptidePosition {
    /// Generate a position for N terminal ion series
    pub const fn n(sequence_index: SequencePosition, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: match sequence_index {
                SequencePosition::NTerm => 0,
                SequencePosition::Index(i) => i + 1,
                SequencePosition::CTerm => length,
            },
            sequence_length: length,
        }
    }
    /// Generate a position for C terminal ion series
    pub const fn c(sequence_index: SequencePosition, length: usize) -> Self {
        Self {
            sequence_index,
            series_number: match sequence_index {
                SequencePosition::NTerm => length,
                SequencePosition::Index(i) => length.saturating_sub(i),
                SequencePosition::CTerm => 0,
            },
            sequence_length: length,
        }
    }
    /// Check if this position is on the N terminus
    pub fn is_n_terminal(&self) -> bool {
        self.sequence_index == SequencePosition::NTerm
    }
    /// Check if this position is on the C terminus
    pub fn is_c_terminal(&self) -> bool {
        self.sequence_index == SequencePosition::CTerm
    }
    /// Flip to the other series (N->C and C->N)
    #[must_use]
    pub const fn flip_terminal(self) -> Self {
        Self {
            sequence_index: self.sequence_index,
            series_number: self.sequence_length + 1 - self.series_number,
            sequence_length: self.sequence_length,
        }
    }
}

/// A flanking sequence
// Impossible to get the Sequence option smaller (size of a pointer plus alignment of a pointer so the discriminator is 8 bytes as well)
#[allow(variant_size_differences)]
#[derive(Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub enum FlankingSequence {
    /// If the flanking sequence is unknown (in _de novo_ for example)
    #[default]
    Unknown,
    /// If this is the terminus
    Terminal,
    /// If only a single amino acid is known (added to prevent overhead of needing to create a sequence)
    AminoAcid(crate::sequence::AminoAcid),
    /// If a (small part of the) sequence is known, always written in N to C direction
    Sequence(Box<crate::sequence::Peptidoform<crate::sequence::SemiAmbiguous>>),
}

impl FlankingSequence {
    /// A static reference to an unknown flanking sequence
    pub const UNKNOWN: &Self = &Self::Unknown;
}
