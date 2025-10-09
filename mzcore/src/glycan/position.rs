use serde::{Deserialize, Serialize};

use crate::sequence::{AminoAcid, SequencePosition};

use super::{GlycanBranchIndex, GlycanBranchMassIndex};

/// The definition of the position of an ion inside a glycan
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct GlycanPosition {
    /// The depth starting at the amino acid
    pub inner_depth: usize,
    /// The series number (from the ion series terminal)
    pub series_number: usize,
    /// The branch naming
    pub branch: Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>,
    /// The aminoacid index where this glycan is attached
    pub attachment: Option<(AminoAcid, SequencePosition)>,
}

impl GlycanPosition {
    /// Get the branch names
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn branch_names(&self) -> String {
        self.branch
            .iter()
            .enumerate()
            .map(|(i, (_, b))| {
                if i == 0 {
                    char::from_u32(
                        (0x03B1..=0x03C9)
                            .chain(0x0391..=0x03A9)
                            .nth(*b)
                            .expect("Too many branches in glycan, out of greek letters"),
                    )
                    .unwrap()
                    .to_string()
                } else if i == 1 {
                    "\'".repeat(*b)
                } else {
                    format!(",{b}")
                }
            })
            .collect::<String>()
    }
    /// Generate the label for this glycan position, example: `1α'`
    /// # Panics
    /// Panics if the first branch number is outside the range of the greek alphabet (small and caps together).
    pub fn label(&self) -> String {
        format!("{}{}", self.series_number, self.branch_names())
    }
    /// Generate the label for this glycan attachment eg N1 (1 based numbering) or an empty string if the attachment is unknown
    pub fn attachment(&self) -> String {
        self.attachment
            .map(|(aa, pos)| format!("{aa}{pos}"))
            .unwrap_or_default()
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

impl std::fmt::Display for GlycanBreakPos {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}{}", self.label(), self.position().label())
    }
}

/// The selected (part) of a glycan to render, using [`Self::FULL`] is a shortcut to get the full glycan.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GlycanSelection<'a> {
    /// A subtree of the glycan, with potentially a break of the root of the subtree and breaks in the branches.
    /// If no breaks are specified the full glycan is shown. The root is the first monosaccharide to be included
    /// in the rendering. The fragment will not include the indicated glycan positions for the branch breaks.
    Subtree(Option<&'a GlycanPosition>, &'a [GlycanPosition]),
    /// A single sugar, all it branches will be shown as broken.
    SingleSugar(&'a GlycanPosition),
}

impl GlycanSelection<'static> {
    /// A shorthand for a full glycan.
    pub const FULL: Self = Self::Subtree(None, &[]);
}
