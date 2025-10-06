use std::collections::HashMap;

use context_error::{BasicKind, BoxedError};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    parse_json::{ParseJson, use_serde},
    prelude::AminoAcid,
};

/// Just keep the full glycan on the structure. This is the behaviour needed for most use cases.
#[derive(Debug, Copy, Clone)]
pub struct FullGlycan {}

impl GlycanAttachement for FullGlycan {
    fn get_default_fragments(&self, _attachment: Option<AminoAcid>) -> GlycanPeptideFragment {
        GlycanPeptideFragment::FULL
    }

    fn get_specific_fragments(
        &self,
        _attachment: Option<AminoAcid>,
    ) -> HashMap<BackboneFragmentKind, GlycanPeptideFragment> {
        HashMap::new()
    }
}

/// Logic to get the rules on how glycans should be calculated. This is built-in to easily support
/// fragmentation of glycans as is seen in CID fragmentation without code duplication.
pub trait GlycanAttachement {
    /// Get the default rules
    fn get_default_fragments(&self, attachment: Option<AminoAcid>) -> GlycanPeptideFragment;

    /// Get the possible glycan peptide fragments based on this attachment location.
    /// This simplifies the rules somewhat to mostly contain unique rules in the fragment specific
    /// rules. But it is not guaranteed to be fully unique.
    fn get_specific_fragments(
        &self,
        attachment: Option<AminoAcid>,
    ) -> HashMap<BackboneFragmentKind, GlycanPeptideFragment>;
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

impl ParseJson for GlycanPeptideFragment {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

/// The possible kinds of fragments, same options as [`FragmentType`] but without any additional data
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[expect(non_camel_case_types)]
pub enum BackboneFragmentKind {
    /// a
    a,
    /// b
    b,
    /// c
    c,
    /// d
    d,
    /// v
    v,
    /// w
    w,
    /// x
    x,
    /// y
    y,
    /// z
    z,
}

impl std::fmt::Display for BackboneFragmentKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a => "a",
                Self::b => "b",
                Self::c => "c",
                Self::d => "d",
                Self::x => "x",
                Self::y => "y",
                Self::v => "v",
                Self::w => "w",
                Self::z => "z",
            }
        )
    }
}

impl ParseJson for BackboneFragmentKind {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}
