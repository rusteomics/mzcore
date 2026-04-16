use std::ops::RangeInclusive;

use mzcore::{
    quantities::Tolerance,
    system::{MassOverCharge, thomson},
};
use serde::{Deserialize, Serialize};

/// Parameters for the matching, allowing control over when a match is allowed.
#[non_exhaustive]
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct MatchingParameters {
    /// The matching tolerance
    pub tolerance: Tolerance<MassOverCharge>,
    /// The range in which fragments fall, can be used to limit the theoretical fragments to a known window
    pub mz_range: RangeInclusive<MassOverCharge>,
    /// Indicate of isotopes should be matched
    pub match_isotopes: bool,
    /// The matching tolerance for isotope peaks
    pub isotope_tolerance: Tolerance<MassOverCharge>,
    /// The minimum cosine similarity for an isotope envelop to be considered a true match
    pub isotope_filter: f64,
}

impl MatchingParameters {
    /// Set the tolerance
    #[must_use]
    pub fn tolerance(self, tolerance: impl Into<Tolerance<MassOverCharge>>) -> Self {
        Self {
            tolerance: tolerance.into(),
            ..self
        }
    }
    /// Set the mz range
    #[must_use]
    pub fn mz_range(self, mz_range: RangeInclusive<MassOverCharge>) -> Self {
        Self { mz_range, ..self }
    }
}

impl Default for MatchingParameters {
    fn default() -> Self {
        Self {
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<thomson>(0.0)
                ..=MassOverCharge::new::<thomson>(f64::MAX),
            match_isotopes: false,
            isotope_tolerance: Tolerance::new_ppm(1.0),
            isotope_filter: 0.5,
        }
    }
}
