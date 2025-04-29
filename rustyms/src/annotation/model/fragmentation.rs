//! Handle model instantiation.

use std::sync::LazyLock;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::{ChargeRange, GlycanModel},
    chemistry::MultiChemical,
    fragment::{NeutralLoss, PeptidePosition},
    sequence::{AminoAcid, BACKBONE, Peptidoform, SequenceElement, SequencePosition},
};

/// A model for the fragmentation, allowing control over what theoretical fragments to generate.
#[non_exhaustive]
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct FragmentationModel {
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
    pub immonium: Option<(ChargeRange, Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>)>,
    /// If the neutral losses specific for modifications should be generated
    pub modification_specific_neutral_losses: bool,
    /// If the diagnostic ions specific for modifications should be generated with the allowed charge range
    pub modification_specific_diagnostic_ions: Option<ChargeRange>,
    /// Glycan fragmentation
    pub glycan: GlycanModel,
    /// Allow any MS cleavable cross-link to be cleaved
    pub allow_cross_link_cleavage: bool,
}

/// The settings for any satellite ion series
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
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

impl Default for SatelliteIonSeries {
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
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
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
    /// Generate a new empty series
    pub fn none() -> Self {
        Self {
            location: Location::None,
            ..Default::default()
        }
    }

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

impl Default for PrimaryIonSeries {
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

/// Builder style methods
impl FragmentationModel {
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
    pub fn immonium(
        self,
        state: Option<(ChargeRange, Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>)>,
    ) -> Self {
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
}

/// A location, or range of locations where an ion can be generated
#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
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
    pub const fn possible(self, position: PeptidePosition) -> bool {
        let SequencePosition::Index(sequence_index) = position.sequence_index else {
            panic!("Not allowed to call possible with a terminal PeptidePosition")
        };
        match self {
            Self::SkipN(n) => sequence_index >= n,
            Self::SkipNC(n, c) => {
                sequence_index >= n && position.sequence_length - sequence_index > c
            }
            Self::TakeN { skip, take } => sequence_index >= skip && sequence_index < skip + take,
            Self::SkipC(n) => position.sequence_length - sequence_index > n,
            Self::TakeC(n) => position.sequence_length - sequence_index <= n,
            Self::All => position.series_number != position.sequence_length,
            Self::None => false,
        }
    }
}

/// The locations for a satellite ion. These are defined to form for any location where the parent
/// ion forms. And the maximal distance from the original backbone cleavage can be defined.
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
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
        let SequencePosition::Index(sequence_index) = position.sequence_index else {
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
) -> Vec<Vec<NeutralLoss>> {
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
            .unique()
            .flat_map(|aa| {
                aa.formulas()
                    .iter()
                    .map(|f| NeutralLoss::SideChainLoss(f - LazyLock::force(&BACKBONE), aa))
                    .filter(|l| !l.is_empty())
                    .collect::<Vec<_>>()
            })
            .collect();
        (1..=settings.0)
            .flat_map(|k| options.iter().combinations(k as usize))
            .map(|o| o.into_iter().cloned().collect_vec())
            .collect()
    }
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn side_chain_losses() {
    let peptide = Peptidoform::pro_forma("FGGGTKLELKR", None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    assert_eq!(
        0,
        get_all_sidechain_losses(peptide.sequence(), &(0, None)).len()
    );
    assert_eq!(
        1,
        get_all_sidechain_losses(
            peptide.sequence(),
            &(1, Some(vec![AminoAcid::Phenylalanine]))
        )
        .len()
    );
    assert_eq!(
        0,
        get_all_sidechain_losses(peptide.sequence(), &(1, Some(vec![AminoAcid::Glycine]))).len()
    );
    assert_eq!(
        1,
        get_all_sidechain_losses(peptide.sequence(), &(1, Some(vec![AminoAcid::Leucine]))).len()
    );
    assert_eq!(
        6,
        get_all_sidechain_losses(peptide.sequence(), &(1, None)).len()
    );
    assert_eq!(
        3,
        dbg!(get_all_sidechain_losses(
            peptide.sequence(),
            &(2, Some(vec![AminoAcid::Phenylalanine, AminoAcid::Leucine]))
        ))
        .len()
    );
}
