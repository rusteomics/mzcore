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
}
