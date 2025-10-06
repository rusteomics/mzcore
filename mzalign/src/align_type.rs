use std::str::FromStr;

use serde::{Deserialize, Serialize};

/// The type of alignment to perform
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct AlignType {
    /// The settings for the left side
    pub left: Side,
    /// The settings for the right side
    pub right: Side,
}

impl AlignType {
    /// Global alignment, which tries to find the best alignment to link both sequences fully to each other, like the Needleman Wunsch algorithm
    pub const GLOBAL: Self = Self::new(Some((true, true)), Some((true, true)));
    /// Local alignment, which tries to find the best patch of both sequences to align to each other, this could lead to trailing ends on both sides of both sequences, like the Smith Waterman
    pub const LOCAL: Self = Self::new(Some((false, false)), Some((false, false)));
    /// Hybrid alignment, the first sequence will be fully aligned to the second sequence, this could lead to trailing ends on the second sequence but not on the first.
    pub const GLOBAL_A: Self = Self::new(Some((true, false)), Some((true, false)));
    /// Hybrid alignment, the second sequence will be fully aligned to the first sequence, this could lead to trailing ends on the first sequence but not on the second.
    pub const GLOBAL_B: Self = Self::new(Some((false, true)), Some((false, true)));
    /// Hybrid alignment, globally align the left (start) side of the alignment
    pub const GLOBAL_LEFT: Self = Self::new(Some((true, true)), Some((false, false)));
    /// Hybrid alignment, globally align the right (end) side of the alignment
    pub const GLOBAL_RIGHT: Self = Self::new(Some((false, false)), Some((true, true)));
    /// Hybrid alignment, extend sequence a with sequence b (▀▀██▄▄)
    pub const EXTEND_A: Self = Self::new(Some((false, true)), Some((true, false)));
    /// Hybrid alignment, extend sequence b with sequence a (▄▄██▀▀)
    pub const EXTEND_B: Self = Self::new(Some((true, false)), Some((false, true)));
    /// Alignment where on both side either of the sequences has to be aligned globally, see [`Side::EitherGlobal`].
    pub const EITHER_GLOBAL: Self = Self::new(None, None);

    /// Create a new alignment type, specifying None on any side indicates [`Side::EitherGlobal`].
    pub(crate) const fn new(left: Option<(bool, bool)>, right: Option<(bool, bool)>) -> Self {
        Self {
            left: if let Some((a, b)) = left {
                Side::Specified { a, b }
            } else {
                Side::EitherGlobal
            },
            right: if let Some((a, b)) = right {
                Side::Specified { a, b }
            } else {
                Side::EitherGlobal
            },
        }
    }

    /// Get a descriptive name for this alignment type
    pub fn description(self) -> String {
        let left = self.left.description();
        let right = self.right.description();

        if left == right {
            left.to_string()
        } else {
            format!("{left} left, {right} right")
        }
    }

    /// Get a concise symbolic representation for this alignment type
    pub fn symbol(self) -> String {
        format!("{}{}", self.left.symbol_left(), self.right.symbol_right())
    }
}

impl std::fmt::Display for AlignType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{} {}", self.symbol(), self.description())
    }
}

impl FromStr for AlignType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(stripped) = s.strip_prefix('-') {
            Ok(Self {
                left: Side::EitherGlobal,
                right: Side::from_str(stripped)?,
            })
        } else if s.len() > 2 {
            Ok(Self {
                left: Side::from_str(&s[0..=1])?,
                right: Side::from_str(&s[2..])?,
            })
        } else {
            Err(())
        }
    }
}

impl Default for AlignType {
    /// Defaults to global alignment
    fn default() -> Self {
        Self::GLOBAL
    }
}

/// The alignment specification for a single side
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum Side {
    /// Align with the specified rules for the two peptides
    Specified {
        /// Rules for peptide A
        a: bool,
        /// Rules for peptide B
        b: bool,
    },
    /// Align either A or B globally, meaning that either one of the peptides has to end before the other,
    /// the other peptide could then have unmatched sequence at the end.
    EitherGlobal,
}

impl Side {
    /// Check if A is global, if this is [`Self::EitherGlobal`] it returns false.
    pub const fn global_a(self) -> bool {
        match self {
            Self::EitherGlobal => false,
            Self::Specified { a, .. } => a,
        }
    }
    /// Check if B is global, if this is [`Self::EitherGlobal`] it returns false.
    pub const fn global_b(self) -> bool {
        match self {
            Self::EitherGlobal => false,
            Self::Specified { b, .. } => b,
        }
    }
    /// Check if any of the peptides is global, or if this is [`Self::EitherGlobal`].
    pub const fn global(self) -> bool {
        match self {
            Self::EitherGlobal => true,
            Self::Specified { a, b } => a || b,
        }
    }
    /// Get a text description of the alignment type for this side.
    pub const fn description(self) -> &'static str {
        match self {
            Self::EitherGlobal => "either global",
            Self::Specified { a, b } => match (a, b) {
                (true, true) => "global",
                (false, false) => "local",
                (true, false) => "global A",
                (false, true) => "global B",
            },
        }
    }
    /// Get a symbolic representation of the alignment type for this side.
    pub const fn symbol_left(self) -> &'static str {
        match self {
            Self::EitherGlobal => "-=",
            Self::Specified { a, b } => match (a, b) {
                (true, true) => "=",
                (false, false) => "⁐=",
                (true, false) => "‿=",
                (false, true) => "⁀=",
            },
        }
    }
    /// Get a symbolic representation of the alignment type for this side.
    pub const fn symbol_right(self) -> &'static str {
        match self {
            Self::EitherGlobal => "=-",
            Self::Specified { a, b } => match (a, b) {
                (true, true) => "=",
                (false, false) => "=⁐",
                (true, false) => "=‿",
                (false, true) => "=⁀",
            },
        }
    }
}

impl FromStr for Side {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "-" => Ok(Self::EitherGlobal),
            "00" => Ok(Self::Specified { a: false, b: false }),
            "01" => Ok(Self::Specified { a: false, b: true }),
            "10" => Ok(Self::Specified { a: true, b: false }),
            "11" => Ok(Self::Specified { a: true, b: true }),
            _ => Err(()),
        }
    }
}
