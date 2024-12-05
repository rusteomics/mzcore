//! Handle model instantiation.

use std::ops::RangeInclusive;

use serde::{Deserialize, Serialize};

use crate::{
    fragment::PeptidePosition,
    system::{e, f64::MassOverCharge, isize::Charge, mz},
    NeutralLoss, Tolerance,
};

/// Control what charges are allowed for an ion series. Defined as an inclusive range.
/// Any charge above the precursor charge will result in the quotient time the precursor
/// charge carriers + all options for the remainder within the limits of the precursor
/// charge carriers.
#[non_exhaustive]
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord, Serialize, Deserialize)]
pub struct ChargeRange {
    /// Start point
    start: ChargePoint,
    /// End point (inclusive)
    end: ChargePoint,
}

impl ChargeRange {
    /// Get all possible charges for the given precursor charge.
    pub fn charges(&self, precursor: Charge) -> RangeInclusive<Charge> {
        Charge::new::<e>(self.start.to_absolute(precursor).value.max(1))
            ..=self.end.to_absolute(precursor)
    }

    /// Get all possible charges for the given precursor charge.
    pub fn charges_iter(
        &self,
        precursor: Charge,
    ) -> impl DoubleEndedIterator<Item = Charge> + Clone {
        (self.start.to_absolute(precursor).value.max(1)..=self.end.to_absolute(precursor).value)
            .map(Charge::new::<e>)
    }

    /// Solely single charged
    pub const ONE: Self = Self {
        start: ChargePoint::Absolute(1),
        end: ChargePoint::Absolute(1),
    };
    /// Only the exact precursor charge
    pub const PRECURSOR: Self = Self {
        start: ChargePoint::Relative(0),
        end: ChargePoint::Relative(0),
    };
    /// Range from 1 to the precursor
    pub const ONE_TO_PRECURSOR: Self = Self {
        start: ChargePoint::Absolute(1),
        end: ChargePoint::Relative(0),
    };
}

/// A reference point for charge range definition.
#[derive(Copy, Clone, PartialEq, Eq, Hash, Debug, PartialOrd, Ord, Serialize, Deserialize)]
pub enum ChargePoint {
    /// Relative to the precursor, with the given offset.
    Relative(isize),
    /// Absolute charge.
    Absolute(isize),
}

impl ChargePoint {
    /// Get the absolute charge of this charge point given a precursor charge
    fn to_absolute(self, precursor: Charge) -> Charge {
        match self {
            Self::Absolute(a) => Charge::new::<e>(a),
            Self::Relative(r) => Charge::new::<e>(precursor.value + r),
        }
    }
}
/// A model for the fragmentation, allowing control over what theoretical fragments to generate.
#[non_exhaustive]
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[allow(clippy::struct_excessive_bools)]
pub struct Model {
    /// a series ions
    pub a: PrimaryIonSeries,
    /// b series ions
    pub b: PrimaryIonSeries,
    /// c series ions
    pub c: PrimaryIonSeries,
    /// d series ions (side chain fragmentation from a)
    pub d: PrimaryIonSeries,
    /// v series ions (full side chain broken off)
    pub v: PrimaryIonSeries,
    /// w series ions (side chain fragmentation from z)
    pub w: PrimaryIonSeries,
    /// x series ions
    pub x: PrimaryIonSeries,
    /// y series ions
    pub y: PrimaryIonSeries,
    /// z series ions
    pub z: PrimaryIonSeries,
    /// precursor ions
    pub precursor: (Vec<NeutralLoss>, ChargeRange),
    /// immonium ions
    pub immonium: (bool, ChargeRange),
    /// m ions, loss of the amino acid side chain from the precursor (follows precursor charge)
    pub m: bool,
    /// If the neutral losses specific for modifications should be generated
    pub modification_specific_neutral_losses: bool,
    /// If the diagnostic ions specific for modifications should be generated with the allowed charge range
    pub modification_specific_diagnostic_ions: (bool, ChargeRange),
    /// Glycan fragmentation
    pub glycan: GlycanModel,
    /// Allow any MS cleavable cross-link to be cleaved
    pub allow_cross_link_cleavage: bool,
    /// The matching tolerance
    pub tolerance: Tolerance<MassOverCharge>,
    /// The range in which fragments fall, can be used to limit the theoretical fragments to a known window
    pub mz_range: RangeInclusive<MassOverCharge>,
}

/// The settings for any primary ion series
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct PrimaryIonSeries {
    /// Which locations are assumed to lead to fragmentation
    pub location: Location,
    /// The allowed neutral losses
    pub neutral_losses: Vec<NeutralLoss>,
    /// The allowed charges
    pub charge_range: ChargeRange,
}

impl PrimaryIonSeries {
    /// Replace the location
    #[must_use]
    pub fn location(self, location: Location) -> Self {
        Self { location, ..self }
    }
    /// Replace the neutral losses
    #[must_use]
    pub fn neutral_losses(self, neutral_losses: Vec<NeutralLoss>) -> Self {
        Self {
            neutral_losses,
            ..self
        }
    }
    /// Replace the charge range
    #[must_use]
    pub fn charge_range(self, charge_range: ChargeRange) -> Self {
        Self {
            charge_range,
            ..self
        }
    }
}

impl std::default::Default for PrimaryIonSeries {
    fn default() -> Self {
        Self {
            location: Location::All,
            neutral_losses: Vec::new(),
            charge_range: ChargeRange::ONE_TO_PRECURSOR,
        }
    }
}

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
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };
    /// Default set for models that disallow glycan fragmentation
    pub const DISALLOW: Self = Self {
        allow_structural: false,
        compositional_range: 0..=0,
        neutral_losses: Vec::new(),
        oxonium_charge_range: ChargeRange::ONE,
        other_charge_range: ChargeRange::ONE_TO_PRECURSOR,
    };
}

/// A struct to handle all possible fragments that could be generated on a single location
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[non_exhaustive]
pub struct PossibleIons<'a> {
    /// a series ions
    pub a: (bool, &'a [NeutralLoss], ChargeRange),
    /// b series ions
    pub b: (bool, &'a [NeutralLoss], ChargeRange),
    /// c series ions
    pub c: (bool, &'a [NeutralLoss], ChargeRange),
    /// d series ions (side chain fragmentation from a)
    pub d: (bool, &'a [NeutralLoss], ChargeRange),
    /// v series ions (full side chain broken off)
    pub v: (bool, &'a [NeutralLoss], ChargeRange),
    /// w series ions (side chain fragmentation from z)
    pub w: (bool, &'a [NeutralLoss], ChargeRange),
    /// x series ions
    pub x: (bool, &'a [NeutralLoss], ChargeRange),
    /// y series ions
    pub y: (bool, &'a [NeutralLoss], ChargeRange),
    /// z series ions
    pub z: (bool, &'a [NeutralLoss], ChargeRange),
    /// precursor ions
    pub precursor: (&'a [NeutralLoss], ChargeRange),
    /// immonium
    pub immonium: (bool, ChargeRange),
}

impl<'a> PossibleIons<'a> {
    /// Give an upper bound for the number of theoretical fragment for these possible ions
    pub fn size_upper_bound(&self) -> usize {
        usize::from(self.a.0) * (self.a.1.len() + 1)
            + usize::from(self.b.0) * (self.b.1.len() + 1)
            + usize::from(self.c.0) * (self.c.1.len() + 1)
            + usize::from(self.d.0) * 2 * (self.d.1.len() + 1)
            + usize::from(self.v.0) * (self.v.1.len() + 1)
            + usize::from(self.w.0) * 2 * (self.w.1.len() + 1)
            + usize::from(self.x.0) * (self.x.1.len() + 1)
            + usize::from(self.y.0) * (self.y.1.len() + 1)
            + usize::from(self.z.0) * 2 * (self.z.1.len() + 1)
            + self.precursor.0.len()
            + 1
    }
}

/// Builder style methods
impl Model {
    /// Set a
    #[must_use]
    pub fn a(self, a: PrimaryIonSeries) -> Self {
        Self { a, ..self }
    }
    /// Set b
    #[must_use]
    pub fn b(self, b: PrimaryIonSeries) -> Self {
        Self { b, ..self }
    }
    /// Set c
    #[must_use]
    pub fn c(self, c: PrimaryIonSeries) -> Self {
        Self { c, ..self }
    }
    /// Set d
    #[must_use]
    pub fn d(self, d: PrimaryIonSeries) -> Self {
        Self { d, ..self }
    }
    /// Set v
    #[must_use]
    pub fn v(self, v: PrimaryIonSeries) -> Self {
        Self { v, ..self }
    }
    /// Set w
    #[must_use]
    pub fn w(self, w: PrimaryIonSeries) -> Self {
        Self { w, ..self }
    }
    /// Set x
    #[must_use]
    pub fn x(self, x: PrimaryIonSeries) -> Self {
        Self { x, ..self }
    }
    /// Set y
    #[must_use]
    pub fn y(self, y: PrimaryIonSeries) -> Self {
        Self { y, ..self }
    }
    /// Set z
    #[must_use]
    pub fn z(self, z: PrimaryIonSeries) -> Self {
        Self { z, ..self }
    }
    /// Set glycan
    #[must_use]
    pub fn glycan(self, glycan: GlycanModel) -> Self {
        Self { glycan, ..self }
    }
    /// Overwrite the precursor neutral losses
    #[must_use]
    pub fn precursor(self, neutral_loss: Vec<NeutralLoss>, charges: ChargeRange) -> Self {
        Self {
            precursor: (neutral_loss, charges),
            ..self
        }
    }
    /// Set immonium
    #[must_use]
    pub fn immonium(self, state: (bool, ChargeRange)) -> Self {
        Self {
            immonium: state,
            ..self
        }
    }
    /// Set m
    #[must_use]
    pub fn m(self, state: bool) -> Self {
        Self { m: state, ..self }
    }
    /// Set modification specific neutral losses
    #[must_use]
    pub fn modification_specific_neutral_losses(self, state: bool) -> Self {
        Self {
            modification_specific_neutral_losses: state,
            ..self
        }
    }
    /// Set modification specific diagnostic ions
    #[must_use]
    pub fn modification_specific_diagnostic_ions(self, state: (bool, ChargeRange)) -> Self {
        Self {
            modification_specific_diagnostic_ions: state,
            ..self
        }
    }
    /// Set the tolerance
    #[must_use]
    pub fn allow_cross_link_cleavage(self, state: bool) -> Self {
        Self {
            allow_cross_link_cleavage: state,
            ..self
        }
    }
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

impl Model {
    /// Give all possible ions for the given N position
    pub fn ions(&self, position: PeptidePosition) -> PossibleIons {
        let c_position = position.flip_terminal();
        PossibleIons {
            a: (
                self.a.location.possible(position),
                self.a.neutral_losses.as_slice(),
                self.a.charge_range,
            ),
            b: (
                self.b.location.possible(position),
                self.b.neutral_losses.as_slice(),
                self.b.charge_range,
            ),
            c: (
                self.c.location.possible(position),
                self.c.neutral_losses.as_slice(),
                self.c.charge_range,
            ),
            d: (
                self.d.location.possible(position),
                self.d.neutral_losses.as_slice(),
                self.d.charge_range,
            ),
            v: (
                self.v.location.possible(c_position),
                self.v.neutral_losses.as_slice(),
                self.v.charge_range,
            ),
            w: (
                self.w.location.possible(c_position),
                self.w.neutral_losses.as_slice(),
                self.w.charge_range,
            ),
            x: (
                self.x.location.possible(c_position),
                self.x.neutral_losses.as_slice(),
                self.x.charge_range,
            ),
            y: (
                self.y.location.possible(c_position),
                self.y.neutral_losses.as_slice(),
                self.y.charge_range,
            ),
            z: (
                self.z.location.possible(c_position),
                self.z.neutral_losses.as_slice(),
                self.z.charge_range,
            ),
            precursor: (self.precursor.0.as_slice(), self.precursor.1),
            immonium: self.immonium,
        }
    }

    /// Generate all possible fragments
    pub fn all() -> Self {
        Self {
            a: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            b: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            d: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            w: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                ChargeRange::PRECURSOR,
            ),
            immonium: (true, ChargeRange::ONE),
            m: true,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: (true, ChargeRange::ONE),
            glycan: GlycanModel::ALLOW
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// Generate no fragments (except for precursor)
    pub fn none() -> Self {
        Self {
            a: PrimaryIonSeries::default().location(Location::None),
            b: PrimaryIonSeries::default().location(Location::None),
            c: PrimaryIonSeries::default().location(Location::None),
            d: PrimaryIonSeries::default().location(Location::None),
            v: PrimaryIonSeries::default().location(Location::None),
            w: PrimaryIonSeries::default().location(Location::None),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default().location(Location::None),
            z: PrimaryIonSeries::default().location(Location::None),
            precursor: (vec![], ChargeRange::PRECURSOR),
            immonium: (false, ChargeRange::ONE),
            m: false,
            modification_specific_neutral_losses: false,
            modification_specific_diagnostic_ions: (false, ChargeRange::ONE),
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: false,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// electron-transfer/higher-energy collisional dissociation
    pub fn ethcd() -> Self {
        Self {
            a: PrimaryIonSeries::default().location(Location::TakeN { skip: 0, take: 1 }),
            b: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            d: PrimaryIonSeries::default().location(Location::TakeN { skip: 0, take: 1 }),
            v: PrimaryIonSeries::default().location(Location::None),
            w: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                ChargeRange::ONE_TO_PRECURSOR,
            ),
            immonium: (true, ChargeRange::ONE),
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: (true, ChargeRange::ONE),
            glycan: GlycanModel::ALLOW
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// hot EACID
    pub fn hot_eacid() -> Self {
        Self {
            a: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            b: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            d: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: PrimaryIonSeries::default().location(Location::None),
            w: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                ChargeRange::ONE_TO_PRECURSOR,
            ),
            immonium: (true, ChargeRange::ONE),
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: (true, ChargeRange::ONE),
            glycan: GlycanModel::ALLOW
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// CID Hcd
    pub fn cid_hcd() -> Self {
        Self {
            a: PrimaryIonSeries::default()
                .location(Location::TakeN { skip: 0, take: 1 })
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            b: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            c: PrimaryIonSeries::default().location(Location::None),
            d: PrimaryIonSeries::default()
                .location(Location::TakeN { skip: 0, take: 1 })
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: PrimaryIonSeries::default().location(Location::None),
            w: PrimaryIonSeries::default().location(Location::None),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default().location(Location::None),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                ChargeRange::PRECURSOR,
            ),
            immonium: (false, ChargeRange::ONE),
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: (true, ChargeRange::ONE),
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// ETD
    pub fn etd() -> Self {
        Self {
            a: PrimaryIonSeries::default().location(Location::None),
            b: PrimaryIonSeries::default().location(Location::None),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            d: PrimaryIonSeries::default().location(Location::None),
            v: PrimaryIonSeries::default().location(Location::None),
            w: PrimaryIonSeries::default().location(Location::None), // TODO: Are w ions also formed here?
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                    NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2)),
                    NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2)),
                ],
                ChargeRange {
                    start: ChargePoint::Relative(-2),
                    end: ChargePoint::Relative(0),
                },
            ),
            immonium: (false, ChargeRange::ONE),
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: (true, ChargeRange::ONE),
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }
}

/// A location, or range of locations where an ion can be generated
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize)]
pub enum Location {
    /// Skip the given number from the N terminal side
    SkipN(usize),
    /// Skip the given number of aminoacids from the N terminal and C terminal side respectively, only using the positions between these two
    SkipNC(usize, usize),
    /// Skip a certain number and then take a certain number of aminoacids
    TakeN {
        /// Skip this number of aminoacids
        skip: usize,
        /// Take this number of aminoacids
        take: usize,
    },
    /// Skip a given number from the C terminal side
    SkipC(usize),
    /// Take a given number of aminoacids from the C terminal side
    TakeC(usize),
    /// All positions (including 0 and len-1)
    All,
    /// Do not allow it anywhere
    #[default]
    None,
}

impl Location {
    /// Determine if an ion is possible on this location
    /// # Panics
    /// If the peptide position is a terminal position
    pub const fn possible(&self, position: PeptidePosition) -> bool {
        let crate::SequencePosition::Index(sequence_index) = position.sequence_index else {
            panic!("Not allowed to call possible with a terminal PeptidePosition")
        };
        match self {
            Self::SkipN(n) => sequence_index >= *n,
            Self::SkipNC(n, c) => {
                sequence_index >= *n && position.sequence_length - sequence_index > *c
            }
            Self::TakeN { skip, take } => sequence_index >= *skip && sequence_index < *skip + *take,
            Self::SkipC(n) => position.sequence_length - sequence_index > *n,
            Self::TakeC(n) => position.sequence_length - sequence_index <= *n,
            Self::All => position.series_number != position.sequence_length,
            Self::None => false,
        }
    }
}

#[test]
#[allow(clippy::missing_panics_doc, clippy::similar_names)]
fn location_all() {
    let all = Model::all();
    let ions_n0 = all.ions(PeptidePosition::n(crate::SequencePosition::default(), 2));
    let ions_c0 = all.ions(PeptidePosition::c(crate::SequencePosition::default(), 2));
    assert!(ions_n0.a.0);
    assert!(!ions_n0.x.0);
    assert!(!ions_c0.a.0);
    assert!(ions_c0.x.0);
}
