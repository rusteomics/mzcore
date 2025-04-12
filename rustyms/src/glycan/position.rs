use serde::{Deserialize, Serialize};

use crate::{AminoAcid, SequencePosition};

use super::{GlycanBranchIndex, GlycanBranchMassIndex};

/// The definition of the position of an ion inside a glycan
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
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
