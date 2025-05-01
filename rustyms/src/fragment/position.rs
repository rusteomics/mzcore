use std::fmt::{Debug, Display};

use serde::{Deserialize, Serialize};

use crate::{
    glycan::{GlycanPosition, MonoSaccharide},
    sequence::{AminoAcid, Modification, SequencePosition},
};

// /// An isotope annotation.
// #[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
// pub struct MatchedIsotopeDistribution {
//     /// The index of the matched peak in the spectrum, if found
//     pub peak_index: Option<usize>,
//     /// The isotope offset in whole daltons from the monoisotopic peak
//     pub isotope_offset: usize,
//     /// The theoretical abundance of this isotope (normalised to 1 for the whole distribution)
//     pub theoretical_isotope_abundance: OrderedFloat<f64>,
// }

/// The definition of the position of an ion
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
#[non_exhaustive]
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
                SequencePosition::Index(i) => length - i,
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

/// Any position on a glycan or a peptide
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum DiagnosticPosition {
    /// A position on a glycan
    Glycan(GlycanPosition, MonoSaccharide),
    /// A position on a compositional glycan (attachment AA + sequence index + the sugar)
    GlycanCompositional(MonoSaccharide, Option<(AminoAcid, SequencePosition)>),
    /// A position on a peptide
    Peptide(PeptidePosition, AminoAcid),
    /// Labile modification
    Labile(Modification),
    /// Reporter ion
    Reporter,
}

/// A label for a satellite ion, none for most amino acids but a or b for Thr and Ile
#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize, Default,
)]
pub enum SatelliteLabel {
    /// No label needed
    #[default]
    None,
    /// Heaviest of the two options
    A,
    /// Lightest of the two options
    B,
}

impl Display for SatelliteLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::None => "",
                Self::A => "a",
                Self::B => "b",
            }
        )
    }
}

/// All positions where a glycan can break
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum GlycanBreakPos {
    /// No breaks just until the end of a chain
    End(GlycanPosition),
    /// Break at a Y position
    Y(GlycanPosition),
    /// Break at a B position
    B(GlycanPosition),
}

impl GlycanBreakPos {
    /// Get the position of this breaking position
    pub const fn position(&self) -> &GlycanPosition {
        match self {
            Self::B(p) | Self::End(p) | Self::Y(p) => p,
        }
    }

    /// Get the label for this breaking position
    pub const fn label(&self) -> &str {
        match self {
            Self::End(_) => "End",
            Self::Y(_) => "Y",
            Self::B(_) => "B",
        }
    }
}

impl Display for GlycanBreakPos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.label(), self.position().label())
    }
}
