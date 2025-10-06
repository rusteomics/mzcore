use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Bound, RangeBounds},
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::ChargeRange,
    fragment::{FragmentKind, NeutralLoss},
    glycan::MonoSaccharide,
    sequence::AminoAcid,
};

use super::built_in::GLYCAN_LOSSES;

/// The settings for glycan fragmentation
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct GlycanModel {
    /// Allows fragments from glycans with defined structures (i.e. GNO modifications)
    pub allow_structural: bool,
    /// Allows fragments from glycans where only the composition is known (i.e. `Glycan:Hex1`).
    /// This allows any fragment containing any number of monosaccharides within this range.
    pub compositional_range: (Option<usize>, Option<usize>),
    /// The allowed neutral losses
    pub neutral_losses: Vec<NeutralLoss>,
    /// Allowed neutral losses on diagnostic ions based on monosaccharides with a flag to indicate precise isomeric state matching
    pub specific_neutral_losses: Vec<(MonoSaccharide, bool, Vec<NeutralLoss>)>,
    /// Glycan fragmentation on peptide fragments when no rules apply
    pub default_peptide_fragment: GlycanPeptideFragment,
    /// Peptide fragment rules
    pub peptide_fragment_rules: Vec<(Vec<AminoAcid>, Vec<FragmentKind>, GlycanPeptideFragment)>,
    /// The allowed charges for oxonium ions (B, internal fragments etc)
    pub oxonium_charge_range: ChargeRange,
    /// The allowed charges for other glycan fragments (Y)
    pub other_charge_range: ChargeRange,
}

impl GlycanModel {
    /// Sets the status of glycan fragments from structural modifications
    #[must_use]
    pub fn allow_structural(self, allow_structural: bool) -> Self {
        Self {
            allow_structural,
            ..self
        }
    }

    /// Set the range of monosaccharides that can result in composition fragments, see [`Self::compositional_range`].
    #[must_use]
    pub fn compositional_range(self, compositional_range: (Option<usize>, Option<usize>)) -> Self {
        Self {
            compositional_range,
            ..self
        }
    }

    /// Set the range of monosaccharides that can result in composition fragments, see [`Self::compositional_range`].
    #[must_use]
    pub fn compositional_range_generic(self, compositional_range: impl RangeBounds<usize>) -> Self {
        Self {
            compositional_range: (
                match compositional_range.start_bound() {
                    Bound::Unbounded => None,
                    Bound::Included(n) => Some(*n),
                    Bound::Excluded(n) => Some(n + 1),
                },
                match compositional_range.end_bound() {
                    Bound::Unbounded => None,
                    Bound::Included(n) => Some(*n),
                    Bound::Excluded(n) => Some(n - 1),
                },
            ),
            ..self
        }
    }

    /// Replace the neutral losses
    #[must_use]
    pub fn neutral_losses(self, neutral_losses: Vec<NeutralLoss>) -> Self {
        Self {
            neutral_losses,
            ..self
        }
    }

    /// Replace the charge range for oxonium ions (B, internal fragments etc)
    #[must_use]
    pub fn oxonium_charge_range(self, oxonium_charge_range: ChargeRange) -> Self {
        Self {
            oxonium_charge_range,
            ..self
        }
    }

    /// Replace the charge range for other glycan ions (Y etc)
    #[must_use]
    pub fn other_charge_range(self, other_charge_range: ChargeRange) -> Self {
        Self {
            other_charge_range,
            ..self
        }
    }

    /// Replace the default rules for glycans on peptide fragments
    #[must_use]
    pub fn default_peptide_fragment(self, default_peptide_fragment: GlycanPeptideFragment) -> Self {
        Self {
            default_peptide_fragment,
            ..self
        }
    }

    /// Replace the specific rules for glycans on peptide fragments
    #[must_use]
    pub fn peptide_fragment_rules(
        self,
        peptide_fragment_rules: Vec<(Vec<AminoAcid>, Vec<FragmentKind>, GlycanPeptideFragment)>,
    ) -> Self {
        Self {
            peptide_fragment_rules,
            ..self
        }
    }

    /// Default set for models that allow glycan fragmentation
    pub fn default_allow() -> Self {
        Self {
            allow_structural: true,
            compositional_range: (None, None),
            neutral_losses: Vec::new(),
            specific_neutral_losses: GLYCAN_LOSSES.clone(),
            default_peptide_fragment: GlycanPeptideFragment::CORE_AND_FREE,
            peptide_fragment_rules: Vec::new(),
            oxonium_charge_range: ChargeRange::ONE,
            other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
        }
    }

    /// Default set for models that disallow glycan fragmentation
    pub const DISALLOW: Self = Self {
        allow_structural: false,
        compositional_range: (None, Some(0)),
        neutral_losses: Vec::new(),
        specific_neutral_losses: Vec::new(),
        default_peptide_fragment: GlycanPeptideFragment::FULL,
        peptide_fragment_rules: Vec::new(),
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };

    /// Get the possible glycan peptide fragments based on this attachment location.
    /// This simplifies the rules somewhat to mostly contain unique rules in the fragment specific
    /// rules. But it is not guaranteed to be fully unique.
    pub fn get_peptide_fragments(
        &self,
        attachment: Option<AminoAcid>,
    ) -> (
        GlycanPeptideFragment,
        HashMap<FragmentKind, GlycanPeptideFragment>,
    ) {
        let base = self
            .peptide_fragment_rules
            .iter()
            .find(|rule| attachment.is_some_and(|a| rule.0.contains(&a)) && rule.1.is_empty())
            .map_or(self.default_peptide_fragment, |rule| rule.2);
        (
            base,
            self.peptide_fragment_rules
                .iter()
                .filter(|rule| attachment.is_some_and(|a| rule.0.contains(&a)))
                .flat_map(|rule| rule.1.iter().map(|f| (f, rule.2)))
                .into_group_map()
                .into_iter()
                .map(|(f, settings)| (*f, settings.iter().fold(settings[0], |acc, v| acc + v)))
                .collect(),
        )
    }
}

/// Rules to determine the glycan fragmentation for glycans on other fragments.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct GlycanPeptideFragment {
    /// The full glycan stays attached
    pub(super) full: bool,
    /// The glycan fragments and at any number of monosaccharides within the range (min, max, inclusive) stay attached (any fucoses on these fragments are always included and do not count towards the limit)
    pub(super) core: Option<(Option<usize>, Option<usize>)>,
}

impl std::ops::Add<&Self> for GlycanPeptideFragment {
    type Output = Self;
    fn add(self, rhs: &Self) -> Self::Output {
        Self {
            full: self.full || rhs.full,
            core: self
                .core
                .and_then(|c| rhs.core.map(|r| (c, r)))
                .map(|(c, r)| (c.0.min(r.0), c.1.max(r.1)))
                .or(self.core)
                .or(rhs.core),
        }
    }
}

impl GlycanPeptideFragment {
    /// A default model that only allows full fragments
    pub const FULL: Self = Self {
        full: true,
        core: None,
    };
    /// A default model that only allows core and free (0..=1)
    pub const CORE_AND_FREE: Self = Self {
        full: false,
        core: Some((None, Some(1))),
    };
    /// A default model that only allows core (1..=1)
    pub const CORE: Self = Self {
        full: false,
        core: Some((Some(1), Some(1))),
    };
    /// A default model that only allows free (0..=0)
    pub const FREE: Self = Self {
        full: false,
        core: Some((None, Some(0))),
    };

    /// Get if this models allows full glycans on peptide fragments
    pub const fn full(self) -> bool {
        self.full
    }

    /// Get if this models allows core fragments on peptide fragments and the depth range of those fragments
    pub const fn core(self) -> Option<(Option<usize>, Option<usize>)> {
        self.core
    }
}

// impl GlycanPeptideFragment {
//     /// Check if this fragment contains anything not contained in the other fragment. None if this
//     /// fragment is a subset of other, Some with the exclusive options if this fragment is disjoint
//     /// or not a true subset of other. If the range of other fully fits inside the self range (so
//     /// subtraction would result in two disjoint ranges for self) the full range of self is returned.
//     fn simplify(self, other: Self) -> Option<Self> {
//         match (self, other) {
//             (Self::Full, Self::Full) => None,
//             (Self::Full, Self::Core(_, _)) => Some(self),
//             (Self::Core(_, _), Self::Full) => Some(self),
//             (Self::Core(mins, maxs), Self::Core(mino, maxo)) => {
//                 if mins >= mino || maxs <= maxo {
//                     None
//                 } else if mins < mino && maxs <= maxo {
//                     Some(Self::Core(mins, mino))
//                 } else if mins >= mino && maxs > maxo {
//                     Some(Self::Core(maxo, maxs))
//                 } else {
//                     Some(self)
//                 }
//             }
//         }
//     }
// }
