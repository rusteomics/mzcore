use itertools::Itertools;
use serde::{Deserialize, Serialize};

use std::{cmp::Ordering, collections::BTreeSet};

use crate::{chemistry::MolecularFormula, fragment::*, sequence::PlacementRule};

/// Indicate the cross-link side, it contains a set of all placement rules that apply for the placed
/// location to find all possible ways of breaking and/or neutral losses. These numbers are the
/// index into the [`LinkerSpecificity`] rules.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum CrossLinkSide {
    /// The cross-link is symmetric, or if asymmetric it can be placed in both orientations
    Symmetric(BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'left' side
    Left(BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'right' side
    Right(BTreeSet<usize>),
}

impl PartialOrd for CrossLinkSide {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CrossLinkSide {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Symmetric(_), Self::Symmetric(_)) | (Self::Left(_), Self::Left(_)) => {
                Ordering::Equal
            }
            (Self::Symmetric(_), _) => Ordering::Greater,
            (_, Self::Symmetric(_)) => Ordering::Less,
            (Self::Left(_), _) => Ordering::Greater,
            (_, Self::Left(_)) => Ordering::Less,
            (Self::Right(_), Self::Right(_)) => Ordering::Equal,
        }
    }
}

impl std::hash::Hash for CrossLinkSide {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let (i, r) = match self {
            Self::Symmetric(r) => (0, r),
            Self::Left(r) => (1, r),
            Self::Right(r) => (2, r),
        };
        state.write_u8(i);
        state.write(
            &r.iter()
                .sorted()
                .flat_map(|r| r.to_ne_bytes())
                .collect_vec(),
        );
    }
}

/// The name of a cross-link
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum CrossLinkName {
    /// A branch
    Branch,
    /// A cross-link
    Name(String),
}

/// The linker position specificities for a linker
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum LinkerSpecificity {
    /// A symmetric specificity where both ends have the same specificity.
    Symmetric {
        /// The placement rules for both ends.
        rules: Vec<PlacementRule>,
        /// All stubs that can be left after cleaving or breaking of the cross-link.
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        /// All possible neutral losses from the intact cross-linker.
        neutral_losses: Vec<NeutralLoss>,
        /// All diagnostic ions from the cross-linker
        diagnostic: Vec<DiagnosticIon>,
    },
    /// An asymmetric specificity where both ends have a different specificity..
    Asymmetric {
        /// The placement rules for both ends, these can be asymmetric thus are provided for 'right' and 'left' separately.
        rules: (Vec<PlacementRule>, Vec<PlacementRule>),
        /// All stubs that can be left after cleaving or breaking of the cross-link. The stubs are specific for right and left in the same orientation as the rules.
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        /// All possible neutral losses from the intact cross-linker.
        neutral_losses: Vec<NeutralLoss>,
        /// All diagnostic ions from the cross-linker
        diagnostic: Vec<DiagnosticIon>,
    },
}
