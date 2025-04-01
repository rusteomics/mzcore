use std::ops::RangeInclusive;

use serde::{Deserialize, Serialize};

use crate::{
    fragment::FragmentKind, glycan::MonoSaccharide, model::ChargeRange, AminoAcid, NeutralLoss,
};

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
    /// Allowed neutral losses based on monosaccharides
    pub specific_neutral_losses: Vec<(MonoSaccharide, Vec<NeutralLoss>)>,
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
    pub const ALLOW: Self = Self {
        allow_structural: true,
        compositional_range: 1..=10,
        neutral_losses: Vec::new(),
        specific_neutral_losses: Vec::new(),
        default_peptide_fragment: GlycanPeptideFragment::CORE_AND_FREE,
        peptide_fragment_rules: Vec::new(),
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };
    /// Default set for models that disallow glycan fragmentation
    pub const DISALLOW: Self = Self {
        allow_structural: false,
        compositional_range: 0..=0,
        neutral_losses: Vec::new(),
        specific_neutral_losses: Vec::new(),
        default_peptide_fragment: GlycanPeptideFragment::FREE,
        peptide_fragment_rules: Vec::new(),
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };
}

/// Rules to determine the glycan fragmentation for glycans on other fragments.
///
/// * free indicates that the entire glycan is removed
/// * core indicates that the core hexnac with potentially any fucoses stays attached
/// * full indicates that the entire glycan stays put
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct GlycanPeptideFragment {
    free: bool,
    core: bool,
    full: bool,
}

impl GlycanPeptideFragment {
    /// Create a new set of rules, at least one rule needs to be true, if this is not the case None is returned.
    pub fn new(free: bool, core: bool, full: bool) -> Option<Self> {
        (free || core || full).then_some(Self { free, core, full })
    }

    /// A constant with only free
    pub const FREE: Self = Self {
        free: true,
        core: false,
        full: false,
    };
    /// A constant with free and core
    pub const CORE_AND_FREE: Self = Self {
        free: true,
        core: true,
        full: false,
    };
}
