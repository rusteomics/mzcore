use std::ops::RangeInclusive;

use serde::{Deserialize, Serialize};

use crate::{
    fragment::FragmentKind, glycan::MonoSaccharide, model::ChargeRange, AminoAcid, NeutralLoss,
};

use super::built_in::glycan_losses;

/// The settings for glycan fragmentation
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct GlycanModel {
    /// Allows fragments from glycans with defined structures (i.e. GNO modifications)
    pub allow_structural: bool,
    /// Allows fragments from glycans where only the composition is known (i.e. `Glycan:Hex1`).
    /// This allows any fragment containing any number of monosaccharides within this range.
    pub compositional_range: RangeInclusive<usize>,
    /// The allowed neutral losses
    pub neutral_losses: Vec<NeutralLoss>,
    /// Allowed neutral losses based on monosaccharides with a flag to indicate precise isomeric state matching
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
    pub fn compositional_range(self, compositional_range: RangeInclusive<usize>) -> Self {
        Self {
            compositional_range,
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
    /// Default set for models that allow glycan fragmentation
    pub fn default_allow() -> Self {
        Self {
            allow_structural: true,
            compositional_range: 1..=10,
            neutral_losses: Vec::new(),
            specific_neutral_losses: glycan_losses().clone(),
            default_peptide_fragment: GlycanPeptideFragment::Core(0, 1),
            peptide_fragment_rules: Vec::new(),
            oxonium_charge_range: ChargeRange::ONE,
            other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
        }
    }
    /// Default set for models that disallow glycan fragmentation
    pub const DISALLOW: Self = Self {
        allow_structural: false,
        compositional_range: 0..=0,
        neutral_losses: Vec::new(),
        specific_neutral_losses: Vec::new(),
        default_peptide_fragment: GlycanPeptideFragment::Full,
        peptide_fragment_rules: Vec::new(),
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };

    /// Get the possible glycan peptide fragments based on this attachment location.
    /// This simplifies the rules somewhat to mostly contain unique rules in the fragment specific
    /// rules. But it is not guarenteed to be fully unique.
    pub fn get_peptide_fragments(
        &self,
        attachment: Option<AminoAcid>,
    ) -> (
        GlycanPeptideFragment,
        Vec<(&[FragmentKind], GlycanPeptideFragment)>,
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
                .filter_map(|rule| rule.2.simplify(base).map(|f| (rule.1.as_slice(), f)))
                .collect(),
        )
    }
}

/// Rules to determine the glycan fragmentation for glycans on other fragments.
#[derive(Clone, Copy, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub enum GlycanPeptideFragment {
    /// The full glycan stays attached
    Full,
    /// The glycan fragments and at any number of monosaccharides within the range (min, max) stay attached (any fucoses on these fragments are always included and do not count towards the limit)
    Core(u8, u8),
}

impl GlycanPeptideFragment {
    /// Check if this fragment contains anything not contained in the other fragment. None if this
    /// fragment is a subset of other, Some with the exclusive options if this fragment is disjoint
    /// or not a true subset of other. If the range of other fully fits inside the self range (so
    /// subtraction would result in two disjoint ranges for self) the full range of self is returned.
    fn simplify(self, other: Self) -> Option<Self> {
        match (self, other) {
            (Self::Full, Self::Full) => None,
            (Self::Full, Self::Core(_, _)) => Some(self),
            (Self::Core(_, _), Self::Full) => Some(self),
            (Self::Core(mins, maxs), Self::Core(mino, maxo)) => {
                if mins >= mino || maxs <= maxo {
                    None
                } else if mins < mino && maxs <= maxo {
                    Some(Self::Core(mins, mino))
                } else if mins >= mino && maxs > maxo {
                    Some(Self::Core(maxo, maxs))
                } else {
                    Some(self)
                }
            }
        }
    }
}
