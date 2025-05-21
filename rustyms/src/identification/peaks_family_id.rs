use std::fmt::Display;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

/// The scans identifier for a peaks identification
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct PeaksFamilyId {
    /// The file, if defined
    pub file: Option<usize>,
    /// The scan numbers
    pub scans: Vec<usize>,
}

impl Display for PeaksFamilyId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self.file.map_or(String::new(), |f| format!("F{f}:")),
            self.scans.iter().join(",")
        )
    }
}

impl std::str::FromStr for PeaksFamilyId {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((start, end)) = s.split_once(':') {
            if start.is_empty() || end.is_empty() {
                Err(())
            } else {
                Ok(Self {
                    file: Some(start[1..].parse().map_err(|_| ())?),
                    scans: end
                        .split(' ')
                        .map(str::parse)
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|_| ())?,
                })
            }
        } else {
            Ok(Self {
                file: None,
                scans: s
                    .split(' ')
                    .map(str::parse)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| ())?,
            })
        }
    }
}
