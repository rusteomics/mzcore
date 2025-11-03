use std::{
    collections::HashMap,
    hash::Hash,
    ops::{Bound, RangeBounds},
};

use itertools::Itertools;
use mzcore::{
    chemistry::{ChargeRange, NeutralLoss},
    glycan::{BackboneFragmentKind, GlycanAttachement, GlycanPeptideFragment, MonoSaccharide},
    prelude::AminoAcid,
};
use serde::{Deserialize, Serialize};

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
    pub peptide_fragment_rules: Vec<(
        Vec<AminoAcid>,
        Vec<BackboneFragmentKind>,
        GlycanPeptideFragment,
    )>,
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
        peptide_fragment_rules: Vec<(
            Vec<AminoAcid>,
            Vec<BackboneFragmentKind>,
            GlycanPeptideFragment,
        )>,
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

    /// Default set for models that use electron based dissociation with some additional collision induced dissociation.
    /// This sets the peptide fragment rules correctly.
    pub fn default_exd_allow() -> Self {
        Self::default_allow()
            .default_peptide_fragment(GlycanPeptideFragment::FULL)
            .peptide_fragment_rules(vec![
                (
                    vec![
                        AminoAcid::Asparagine,
                        AminoAcid::Tryptophan,
                        AminoAcid::Serine,
                        AminoAcid::Threonine,
                    ],
                    vec![BackboneFragmentKind::c, BackboneFragmentKind::z],
                    GlycanPeptideFragment::FULL,
                ),
                (
                    vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                    vec![BackboneFragmentKind::b, BackboneFragmentKind::y],
                    GlycanPeptideFragment::CORE,
                ),
                (
                    vec![AminoAcid::Serine, AminoAcid::Threonine],
                    vec![BackboneFragmentKind::b, BackboneFragmentKind::y],
                    GlycanPeptideFragment::FREE,
                ),
            ])
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
}

impl GlycanAttachement for GlycanModel {
    fn get_default_fragments(&self, attachment: Option<AminoAcid>) -> GlycanPeptideFragment {
        self.peptide_fragment_rules
            .iter()
            .find(|rule| attachment.is_some_and(|a| rule.0.contains(&a)) && rule.1.is_empty())
            .map_or(self.default_peptide_fragment, |rule| rule.2)
    }

    fn get_specific_fragments(
        &self,
        attachment: Option<AminoAcid>,
    ) -> HashMap<BackboneFragmentKind, GlycanPeptideFragment> {
        self.peptide_fragment_rules
            .iter()
            .filter(|rule| attachment.is_some_and(|a| rule.0.contains(&a)))
            .flat_map(|rule| rule.1.iter().map(|f| (f, rule.2)))
            .into_group_map()
            .into_iter()
            .map(|(f, settings)| (*f, settings.iter().fold(settings[0], |acc, v| acc + v)))
            .collect()
    }
}
