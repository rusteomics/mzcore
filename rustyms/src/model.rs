//! Handle model instantiation.

use std::ops::RangeInclusive;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    fragment::PeptidePosition,
    system::{e, f64::MassOverCharge, isize::Charge, mz},
    AminoAcid, MultiChemical, NeutralLoss, Peptidoform, SequenceElement, Tolerance,
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
    /// Get the number of possible charges for the given precursor charge.
    pub fn len(&self, precursor: Charge) -> usize {
        (self.end.to_absolute(precursor).value - self.start.to_absolute(precursor).value.max(1))
            .unsigned_abs()
    }

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
pub struct Model {
    /// a series ions
    pub a: PrimaryIonSeries,
    /// b series ions
    pub b: PrimaryIonSeries,
    /// c series ions
    pub c: PrimaryIonSeries,
    /// d series ions (side chain fragmentation from a)
    pub d: SatelliteIonSeries,
    /// v series ions (full side chain broken off from y)
    pub v: SatelliteIonSeries,
    /// w series ions (side chain fragmentation from z)
    pub w: SatelliteIonSeries,
    /// x series ions
    pub x: PrimaryIonSeries,
    /// y series ions
    pub y: PrimaryIonSeries,
    /// z series ions
    pub z: PrimaryIonSeries,
    /// precursor ions, standard neutral losses, amino acid specific neutral losses, and charge range
    pub precursor: (
        Vec<NeutralLoss>,
        Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
        (u8, Option<Vec<AminoAcid>>),
        ChargeRange,
    ),
    /// immonium ions
    pub immonium: Option<ChargeRange>,
    /// If the neutral losses specific for modifications should be generated
    pub modification_specific_neutral_losses: bool,
    /// If the diagnostic ions specific for modifications should be generated with the allowed charge range
    pub modification_specific_diagnostic_ions: Option<ChargeRange>,
    /// Glycan fragmentation
    pub glycan: GlycanModel,
    /// Allow any MS cleavable cross-link to be cleaved
    pub allow_cross_link_cleavage: bool,
    /// The matching tolerance
    pub tolerance: Tolerance<MassOverCharge>,
    /// The range in which fragments fall, can be used to limit the theoretical fragments to a known window
    pub mz_range: RangeInclusive<MassOverCharge>,
}

/// The settings for any satellite ion series
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct SatelliteIonSeries {
    /// Which locations are assumed to lead to fragmentation
    pub location: SatelliteLocation,
    /// The allowed neutral losses
    pub neutral_losses: Vec<NeutralLoss>,
    /// The allowed amino acid specific neutral losses
    pub amino_acid_neutral_losses: Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
    /// The maximal number of side chain lost and the amino acids from which the side chains can be lost (or all if no selection is given)
    pub amino_acid_side_chain_losses: (u8, Option<Vec<AminoAcid>>),
    /// The allowed charges
    pub charge_range: ChargeRange,
    /// The allowed ion variants (e.g. w vs w+1 vs w-1)
    pub allowed_variants: Vec<i8>,
}

impl SatelliteIonSeries {
    /// A default value where the base is `Some(0)`
    #[must_use]
    pub fn base() -> Self {
        Self {
            location: SatelliteLocation {
                rules: Vec::new(),
                base: Some(0),
            },
            ..Self::default()
        }
    }

    /// Replace the location
    #[must_use]
    pub fn location(self, location: SatelliteLocation) -> Self {
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
    /// Replace the amino acid specific neutral losses
    #[must_use]
    pub fn amino_acid_neutral_losses(
        self,
        amino_acid_neutral_losses: Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
    ) -> Self {
        Self {
            amino_acid_neutral_losses,
            ..self
        }
    }
    /// Replace the amino acid side chain losses
    #[must_use]
    pub fn amino_acid_side_chain_losses(
        self,
        amino_acid_side_chain_losses: (u8, Option<Vec<AminoAcid>>),
    ) -> Self {
        Self {
            amino_acid_side_chain_losses,
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
    /// Replace the allowed variants, e.g. a vs a+1 vs a+2
    #[must_use]
    pub fn variants(self, allowed_variants: Vec<i8>) -> Self {
        Self {
            allowed_variants,
            ..self
        }
    }
}

impl std::default::Default for SatelliteIonSeries {
    fn default() -> Self {
        Self {
            location: SatelliteLocation::default(),
            neutral_losses: Vec::new(),
            amino_acid_neutral_losses: Vec::new(),
            amino_acid_side_chain_losses: (0, None),
            charge_range: ChargeRange::ONE_TO_PRECURSOR,
            allowed_variants: vec![0],
        }
    }
}

/// The settings for any primary ion series
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
pub struct PrimaryIonSeries {
    /// Which locations are assumed to lead to fragmentation
    pub location: Location,
    /// The allowed neutral losses
    pub neutral_losses: Vec<NeutralLoss>,
    /// The allowed amino acid specific neutral losses
    pub amino_acid_neutral_losses: Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
    /// The maximal number of side chain lost and the amino acids from which the side chains can be lost (or all if no selection is given)
    pub amino_acid_side_chain_losses: (u8, Option<Vec<AminoAcid>>),
    /// The allowed charges
    pub charge_range: ChargeRange,
    /// The allowed ion variants (e.g. a vs a+1 vs a+2)
    pub allowed_variants: Vec<i8>,
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
    /// Replace the amino acid specific neutral losses
    #[must_use]
    pub fn amino_acid_neutral_losses(
        self,
        amino_acid_neutral_losses: Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
    ) -> Self {
        Self {
            amino_acid_neutral_losses,
            ..self
        }
    }
    /// Replace the amino acid side chain losses
    #[must_use]
    pub fn amino_acid_side_chain_losses(
        self,
        amino_acid_side_chain_losses: (u8, Option<Vec<AminoAcid>>),
    ) -> Self {
        Self {
            amino_acid_side_chain_losses,
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
    /// Replace the allowed variants, e.g. a vs a+1 vs a+2
    #[must_use]
    pub fn variants(self, allowed_variants: Vec<i8>) -> Self {
        Self {
            allowed_variants,
            ..self
        }
    }
}

impl std::default::Default for PrimaryIonSeries {
    fn default() -> Self {
        Self {
            location: Location::All,
            neutral_losses: Vec::new(),
            amino_acid_neutral_losses: Vec::new(),
            amino_acid_side_chain_losses: (0, None),
            charge_range: ChargeRange::ONE_TO_PRECURSOR,
            allowed_variants: vec![0],
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

/// The possibilities for primary ions, a list of all allowed neutral losses, all charge options, and all allowed variant ions
pub type PossiblePrimaryIons<'a> = (Vec<NeutralLoss>, ChargeRange, &'a [i8]);
/// The possibilities for satellite ions, a list of all satellite ions with their amino acid and
/// distance from the parent backbone cleavage, as well as all ion settings as for primary series.
pub type PossibleSatelliteIons<'a> = (Vec<(AminoAcid, u8)>, PossiblePrimaryIons<'a>);

/// A struct to handle all possible fragments that could be generated on a single location
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
#[non_exhaustive]
pub struct PossibleIons<'a> {
    /// a series ions
    pub a: Option<PossiblePrimaryIons<'a>>,
    /// b series ions
    pub b: Option<PossiblePrimaryIons<'a>>,
    /// c series ions
    pub c: Option<PossiblePrimaryIons<'a>>,
    /// d series ions (side chain fragmentation from a)
    pub d: PossibleSatelliteIons<'a>,
    /// v series ions (full side chain broken off from y)
    pub v: PossibleSatelliteIons<'a>,
    /// w series ions (side chain fragmentation from z)
    pub w: PossibleSatelliteIons<'a>,
    /// x series ions
    pub x: Option<PossiblePrimaryIons<'a>>,
    /// y series ions
    pub y: Option<PossiblePrimaryIons<'a>>,
    /// z series ions
    pub z: Option<PossiblePrimaryIons<'a>>,
    /// immonium
    pub immonium: Option<ChargeRange>,
}

impl PossibleIons<'_> {
    /// Give an upper bound for the number of theoretical fragment for these possible ions
    pub fn size_upper_bound(&self) -> usize {
        self.a
            .as_ref()
            .map(|o| (o.0.len() + 1) * o.2.len())
            .unwrap_or_default()
            + self
                .b
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .c
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self.d.0.len() * 2 * (self.d.1 .0.len() + 1) * self.d.1 .2.len()
            + self.v.0.len() * (self.v.1 .0.len() + 1) * self.v.1 .2.len()
            + self.w.0.len() * 2 * (self.w.1 .0.len() + 1) * self.w.1 .2.len()
            + self
                .x
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .y
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .z
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + usize::from(self.immonium.is_some())
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
    pub fn d(self, d: SatelliteIonSeries) -> Self {
        Self { d, ..self }
    }
    /// Set v
    #[must_use]
    pub fn v(self, v: SatelliteIonSeries) -> Self {
        Self { v, ..self }
    }
    /// Set w
    #[must_use]
    pub fn w(self, w: SatelliteIonSeries) -> Self {
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
    pub fn precursor(
        self,
        neutral_loss: Vec<NeutralLoss>,
        amino_acid_specific_neutral_losses: Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
        amino_acid_side_chain_losses: (u8, Option<Vec<AminoAcid>>),
        charges: ChargeRange,
    ) -> Self {
        Self {
            precursor: (
                neutral_loss,
                amino_acid_specific_neutral_losses,
                amino_acid_side_chain_losses,
                charges,
            ),
            ..self
        }
    }
    /// Set immonium
    #[must_use]
    pub fn immonium(self, state: Option<ChargeRange>) -> Self {
        Self {
            immonium: state,
            ..self
        }
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
    pub fn modification_specific_diagnostic_ions(self, state: Option<ChargeRange>) -> Self {
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
    /// Give all possible ions for the given position.
    /// # Panics
    /// If the position is a terminal position.
    pub fn ions<Complexity>(
        &self,
        position: PeptidePosition,
        peptidoform: &Peptidoform<Complexity>,
    ) -> PossibleIons {
        let crate::SequencePosition::Index(sequence_index) = position.sequence_index else {
            panic!("Not allowed to call possible with a terminal PeptidePosition")
        };

        let get_neutral_losses = |neutral_losses: &Vec<NeutralLoss>,
                                  amino_acid_specific: &Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
                                  amino_acid_side_chains: &(u8, Option<Vec<AminoAcid>>),
                                  c_terminal| {
            let peptide_slice = &peptidoform[if c_terminal {
                sequence_index..peptidoform.len()
            } else {
                0..sequence_index
            }];
            neutral_losses
                .iter()
                .chain(
                    peptide_slice
                        .iter()
                        .flat_map(|seq| {
                            amino_acid_specific.iter().filter_map(|(rule, loss)| {
                                rule.contains(&seq.aminoacid.aminoacid()).then_some(loss)
                            })
                        })
                        .flatten(),
                )
                .cloned()
                .chain(get_all_sidechain_losses(
                    peptide_slice,
                    &amino_acid_side_chains,
                ))
                .collect()
        };

        let c_position = position.flip_terminal();

        PossibleIons {
            a: self.a.location.possible(position).then_some((
                get_neutral_losses(
                    &self.a.neutral_losses,
                    &self.a.amino_acid_neutral_losses,
                    &self.a.amino_acid_side_chain_losses,
                    false,
                ),
                self.a.charge_range,
                self.a.allowed_variants.as_slice(),
            )),
            b: self.b.location.possible(position).then_some((
                get_neutral_losses(
                    &self.b.neutral_losses,
                    &self.b.amino_acid_neutral_losses,
                    &self.b.amino_acid_side_chain_losses,
                    false,
                ),
                self.b.charge_range,
                self.b.allowed_variants.as_slice(),
            )),
            c: self.c.location.possible(position).then_some((
                get_neutral_losses(
                    &self.c.neutral_losses,
                    &self.c.amino_acid_neutral_losses,
                    &self.c.amino_acid_side_chain_losses,
                    false,
                ),
                self.c.charge_range,
                self.c.allowed_variants.as_slice(),
            )),
            d: (
                self.d.location.possible(position, peptidoform, false),
                (
                    get_neutral_losses(
                        &self.d.neutral_losses,
                        &self.d.amino_acid_neutral_losses,
                        &self.d.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.d.charge_range,
                    self.d.allowed_variants.as_slice(),
                ),
            ),
            v: (
                self.v.location.possible(c_position, peptidoform, true),
                (
                    get_neutral_losses(
                        &self.v.neutral_losses,
                        &self.v.amino_acid_neutral_losses,
                        &self.v.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.v.charge_range,
                    self.v.allowed_variants.as_slice(),
                ),
            ),
            w: (
                self.w.location.possible(c_position, peptidoform, true),
                (
                    get_neutral_losses(
                        &self.w.neutral_losses,
                        &self.w.amino_acid_neutral_losses,
                        &self.w.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.w.charge_range,
                    self.w.allowed_variants.as_slice(),
                ),
            ),
            x: self.x.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.x.neutral_losses,
                    &self.x.amino_acid_neutral_losses,
                    &self.x.amino_acid_side_chain_losses,
                    false,
                ),
                self.x.charge_range,
                self.x.allowed_variants.as_slice(),
            )),
            y: self.y.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.y.neutral_losses,
                    &self.y.amino_acid_neutral_losses,
                    &self.y.amino_acid_side_chain_losses,
                    false,
                ),
                self.y.charge_range,
                self.y.allowed_variants.as_slice(),
            )),
            z: self.z.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.z.neutral_losses,
                    &self.z.amino_acid_neutral_losses,
                    &self.z.amino_acid_side_chain_losses,
                    false,
                ),
                self.z.charge_range,
                self.z.allowed_variants.as_slice(),
            )),
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
            d: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            w: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (1, None),
                ChargeRange::PRECURSOR,
            ),
            immonium: Some(ChargeRange::ONE),
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
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
            d: SatelliteIonSeries::default(),
            v: SatelliteIonSeries::default(),
            w: SatelliteIonSeries::default(),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default().location(Location::None),
            z: PrimaryIonSeries::default().location(Location::None),
            precursor: (Vec::new(), Vec::new(), (0, None), ChargeRange::PRECURSOR),
            immonium: None,
            modification_specific_neutral_losses: false,
            modification_specific_diagnostic_ions: None,
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: false,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// UVPD
    /// 10.1021/acs.chemrev.9b00440 and 10.1021/jacs.6b05147
    pub fn uvpd() -> Self {
        Self {
            a: PrimaryIonSeries::default().variants(vec![0, 1, 2]),
            b: PrimaryIonSeries::default().variants(vec![0, 2]),
            c: PrimaryIonSeries::default(),
            d: SatelliteIonSeries::base(),
            v: SatelliteIonSeries::base(),
            w: SatelliteIonSeries::base(),
            x: PrimaryIonSeries::default().variants(vec![-1, 0, 1, 2]),
            y: PrimaryIonSeries::default().variants(vec![-2, -1, 0]),
            z: PrimaryIonSeries::default(),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (0, None),
                ChargeRange::PRECURSOR,
            ),
            immonium: Some(ChargeRange::ONE),
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
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
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
                .variants(vec![0, 1]),
            d: SatelliteIonSeries::base(),
            v: SatelliteIonSeries::default(),
            w: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (0, None),
                ChargeRange::ONE_TO_PRECURSOR,
            ),
            immonium: None,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
            glycan: GlycanModel::ALLOW
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// EAD
    pub fn ead() -> Self {
        Self {
            a: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            b: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            d: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            w: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (0, None),
                ChargeRange::ONE_TO_PRECURSOR,
            ),
            immonium: Some(ChargeRange::ONE),
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
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
            d: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            w: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            x: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (0, None),
                ChargeRange::ONE_TO_PRECURSOR,
            ),
            immonium: None,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
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
            d: SatelliteIonSeries::base()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            v: SatelliteIonSeries::default(),
            w: SatelliteIonSeries::default(),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default().location(Location::None),
            precursor: (
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
                Vec::new(),
                (0, None),
                ChargeRange::PRECURSOR,
            ),
            immonium: None,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// ETD 10.1002/jms.3919
    pub fn etd() -> Self {
        Self {
            a: PrimaryIonSeries::default().location(Location::None),
            b: PrimaryIonSeries::default().location(Location::None),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
                .variants(vec![-1, 0, 2]),
            d: SatelliteIonSeries::default(),
            v: SatelliteIonSeries::base(),
            w: SatelliteIonSeries::default()
                .location(SatelliteLocation {
                    rules: vec![
                        (vec![AminoAcid::Methionine], 5),
                        (vec![AminoAcid::Leucine, AminoAcid::GlutamicAcid], 2),
                        (vec![AminoAcid::Isoleucine, AminoAcid::AsparticAcid], 1),
                        (
                            vec![
                                AminoAcid::Valine,
                                AminoAcid::Asparagine,
                                AminoAcid::Threonine,
                                AminoAcid::Serine,
                                AminoAcid::Tryptophan,
                                AminoAcid::Histidine,
                                AminoAcid::Phenylalanine,
                                AminoAcid::Tyrosine,
                            ],
                            0,
                        ),
                    ],
                    base: None,
                })
                .variants(vec![-1, 0, 1]),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
                .variants(vec![-1, 0, 1, 2]),
            precursor: (
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                ],
                vec![
                    (
                        vec![AminoAcid::AsparticAcid],
                        vec![NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2))],
                    ),
                    (
                        vec![AminoAcid::GlutamicAcid],
                        vec![NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2))],
                    ),
                ],
                (2, None),
                ChargeRange {
                    start: ChargePoint::Relative(-2),
                    end: ChargePoint::Relative(0),
                },
            ),
            immonium: None,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
            glycan: GlycanModel::DISALLOW,
            allow_cross_link_cleavage: true,
            tolerance: Tolerance::new_ppm(20.0),
            mz_range: MassOverCharge::new::<mz>(0.0)..=MassOverCharge::new::<mz>(f64::MAX),
        }
    }

    /// Top Down ETD
    pub fn td_etd() -> Self {
        Self {
            a: PrimaryIonSeries::default().location(Location::None),
            b: PrimaryIonSeries::default().location(Location::None),
            c: PrimaryIonSeries::default()
                .neutral_losses(vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                ])
                .variants(vec![0, 1, 2]),
            d: SatelliteIonSeries::default(),
            v: SatelliteIonSeries::default(),
            w: SatelliteIonSeries::default(),
            x: PrimaryIonSeries::default().location(Location::None),
            y: PrimaryIonSeries::default().location(Location::None),
            z: PrimaryIonSeries::default()
                .neutral_losses(vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                ])
                .variants(vec![-1, 0, 1, 2]),
            precursor: (
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                    NeutralLoss::Gain(molecular_formula!(H 1)),
                    NeutralLoss::Gain(molecular_formula!(H 2)),
                    NeutralLoss::Gain(molecular_formula!(H 3)),
                ],
                vec![
                    (
                        vec![AminoAcid::AsparticAcid],
                        vec![NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2))],
                    ),
                    (
                        vec![AminoAcid::GlutamicAcid],
                        vec![NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2))],
                    ),
                ],
                (0, None),
                ChargeRange::PRECURSOR,
            ),
            immonium: None,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
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

/// The locations for a satellite ion. These are defined to form for any location where the parent
/// ion forms. And the maximal distance from the original backbone cleavage can be defined.
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize)]
pub struct SatelliteLocation {
    /// A set of rules each indicate a set of aminoacids and the (0 based) distance they can occur
    /// away from the original cleavage.
    pub rules: Vec<(Vec<AminoAcid>, u8)>,
    /// The base distance that this satellite ion can occur for any amino acid undefined in the
    /// specific rules.
    pub base: Option<u8>,
}

impl SatelliteLocation {
    /// Determine if a satellite ion is possible on this location
    /// # Panics
    /// If the peptide position is a terminal position
    pub fn possible<Complexity>(
        &self,
        position: PeptidePosition,
        peptidoform: &Peptidoform<Complexity>,
        c_terminal: bool,
    ) -> Vec<(AminoAcid, u8)> {
        let crate::SequencePosition::Index(sequence_index) = position.sequence_index else {
            panic!("Not allowed to call possible with a terminal PeptidePosition")
        };
        let mut output = Vec::new();
        let max_distance = match (self.base, self.rules.iter().map(|r| r.1).max()) {
            (Some(b), Some(r)) => b.max(r) + 1,
            (Some(b), None) => b + 1,
            (None, Some(r)) => r + 1,
            (None, None) => return Vec::new(),
        };
        let range = if c_terminal {
            sequence_index
                ..sequence_index
                    .saturating_add(max_distance as usize)
                    .min(peptidoform.len() - 1)
        } else {
            sequence_index.saturating_sub(max_distance as usize)..sequence_index
        };
        for (index, seq) in peptidoform[range].iter().enumerate() {
            if let Ok(distance) = u8::try_from(if c_terminal {
                index
            } else {
                max_distance as usize - index
            }) {
                let mut allowed = None;
                for rule in &self.rules {
                    if rule.0.contains(&seq.aminoacid.aminoacid()) {
                        allowed = Some(distance <= rule.1);
                        break;
                    }
                }
                if allowed
                    .or_else(|| self.base.map(|b| distance <= b))
                    .is_some_and(|a| a)
                {
                    output.push((seq.aminoacid.aminoacid(), distance));
                }
            }
        }
        output
    }
}

/// Get all possible side chain losses for a given stretch of amino acids. The number indicates
/// the maximum total side chains lost, the selection restricts the set of amino acids that can
/// lose their side chain, if this is None all amino acids can loose their side chain.
///
/// This might generate non sensible options for satellite ions (e.g. double side chain loss for v)
pub(crate) fn get_all_sidechain_losses<Complexity>(
    slice: &[SequenceElement<Complexity>],
    settings: &(u8, Option<Vec<AminoAcid>>),
) -> Vec<NeutralLoss> {
    if settings.0 == 0 {
        Vec::new()
    } else {
        let options: Vec<NeutralLoss> = slice
            .iter()
            .filter_map(|seq| {
                settings
                    .1
                    .as_ref()
                    .is_none_or(|aa| aa.contains(&seq.aminoacid.aminoacid()))
                    .then_some(seq.aminoacid.aminoacid())
            })
            .flat_map(|aa| {
                aa.formulas()
                    .iter()
                    .map(|f| {
                        NeutralLoss::SideChainLoss(f - molecular_formula!(H 3 C 2 N 1 O 1), aa)
                    })
                    .collect::<Vec<_>>()
            })
            .collect();
        (1..=settings.0)
            .flat_map(|k| options.iter().combinations(k as usize))
            .flatten()
            .cloned()
            .collect()
    }
}
