use std::{ops::RangeInclusive, path::PathBuf};

use serde::{Deserialize, Serialize};

use crate::system::OrderedTime;

/// Multiple spectrum identifiers
#[derive(Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub enum SpectrumIds {
    /// When no spectra references are known at all
    #[default]
    None,
    /// When the source file is now known
    FileNotKnown(Vec<SpectrumId>),
    /// When the source file is known, grouped per file
    FileKnown(Vec<(PathBuf, Vec<SpectrumId>)>),
}

/// A spectrum identifier
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum SpectrumId {
    /// A native id, the format differs between vendors
    Native(String),
    /// A spectrum index
    Index(usize),
    /// A scan number, unless there is a better alternative should be interpreted as index+1
    Number(usize),
    /// Time range, assumes all MS2 spectra within this range are selected
    RetentionTime(RangeInclusive<OrderedTime>),
}

impl Default for SpectrumId {
    fn default() -> Self {
        Self::Index(0)
    }
}

impl std::fmt::Display for SpectrumId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Index(i) => write!(f, "{i}"),
            Self::Number(i) => write!(f, "\x23{i}"),
            Self::Native(n) => write!(f, "{n}"),
            Self::RetentionTime(n) => {
                write!(f, "{:.3} — {:.3} min", n.start().value, n.end().value)
            }
        }
    }
}

impl SpectrumId {
    /// Get the index if this is an index or scan number
    pub fn index(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::Number(i) => Some(i - 1),
            Self::Native(_) | Self::RetentionTime(_) => None,
        }
    }

    /// Get the native ID if this is a native ID
    pub fn native(&self) -> Option<&str> {
        match self {
            Self::Native(n) => Some(n),
            Self::Index(_) | Self::RetentionTime(_) | Self::Number(_) => None,
        }
    }
}
