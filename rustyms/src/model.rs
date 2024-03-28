//! Handle model instantiation.

use serde::{Deserialize, Serialize};

use crate::{
    system::{f64::MassOverCharge, mass_over_charge::mz},
    Element, MolecularFormula, NeutralLoss,
};

/// A model for the fragmentation, allowing control over what theoretical fragments to generate.
#[non_exhaustive]
#[derive(Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct Model {
    /// a series ions
    pub a: (Location, Vec<NeutralLoss>),
    /// b series ions
    pub b: (Location, Vec<NeutralLoss>),
    /// c series ions
    pub c: (Location, Vec<NeutralLoss>),
    /// d series ions (side chain fragmentation from a)
    pub d: (Location, Vec<NeutralLoss>),
    /// v series ions (full side chain broken off)
    pub v: (Location, Vec<NeutralLoss>),
    /// w series ions (side chain fragmentation from z)
    pub w: (Location, Vec<NeutralLoss>),
    /// x series ions
    pub x: (Location, Vec<NeutralLoss>),
    /// y series ions
    pub y: (Location, Vec<NeutralLoss>),
    /// z series ions
    pub z: (Location, Vec<NeutralLoss>),
    /// precursor ions
    pub precursor: Vec<NeutralLoss>,
    /// immonium ions
    pub immonium: bool,
    /// m ions, loss of the amino acid side chain from the precursor
    pub m: bool,
    /// If the neutral losses specific for modifications should be generated
    pub modification_specific_neutral_losses: bool,
    /// If the diagnostic ions specific for modifications should be generated
    pub modification_specific_diagnostic_ions: bool,
    /// If some search for glycan (B/Y/Internal) with the given neutral losses
    pub glycan: Option<Vec<NeutralLoss>>,
    /// The matching tolerance
    pub ppm: MassOverCharge,
}

/// A struct to handle all possible fragments that could be generated on a single location
#[allow(clippy::struct_excessive_bools)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Hash)]
pub struct PossibleIons<'a> {
    /// a series ions
    pub a: (bool, &'a [NeutralLoss]),
    /// b series ions
    pub b: (bool, &'a [NeutralLoss]),
    /// c series ions
    pub c: (bool, &'a [NeutralLoss]),
    /// d series ions (side chain fragmentation from a)
    pub d: (bool, &'a [NeutralLoss]),
    /// v series ions (full side chain broken off)
    pub v: (bool, &'a [NeutralLoss]),
    /// w series ions (side chain fragmentation from z)
    pub w: (bool, &'a [NeutralLoss]),
    /// x series ions
    pub x: (bool, &'a [NeutralLoss]),
    /// y series ions
    pub y: (bool, &'a [NeutralLoss]),
    /// z series ions
    pub z: (bool, &'a [NeutralLoss]),
    /// precursor ions
    pub precursor: &'a [NeutralLoss],
    /// immonium
    pub immonium: bool,
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
            + self.precursor.len()
            + 1
    }
}

impl Model {
    /// Give all possible ions for the given position
    pub fn ions(&self, index: usize, length: usize) -> PossibleIons {
        PossibleIons {
            a: (self.a.0.possible(index, length), self.a.1.as_slice()),
            b: (self.b.0.possible(index, length), self.b.1.as_slice()),
            c: (self.c.0.possible(index, length), self.c.1.as_slice()),
            d: (self.d.0.possible(index, length), self.d.1.as_slice()),
            v: (self.v.0.possible(index, length), self.v.1.as_slice()),
            w: (self.w.0.possible(index, length), self.w.1.as_slice()),
            x: (self.x.0.possible(index, length), self.x.1.as_slice()),
            y: (self.y.0.possible(index, length), self.y.1.as_slice()),
            z: (self.z.0.possible(index, length), self.z.1.as_slice()),
            precursor: self.precursor.as_slice(),
            immonium: self.immonium,
        }
    }

    /// Build a new model
    #[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
    pub fn new(
        a: (Location, Vec<NeutralLoss>),
        b: (Location, Vec<NeutralLoss>),
        c: (Location, Vec<NeutralLoss>),
        d: (Location, Vec<NeutralLoss>),
        v: (Location, Vec<NeutralLoss>),
        w: (Location, Vec<NeutralLoss>),
        x: (Location, Vec<NeutralLoss>),
        y: (Location, Vec<NeutralLoss>),
        z: (Location, Vec<NeutralLoss>),
        precursor: Vec<NeutralLoss>,
        immonium: bool,
        m: bool,
        modification_specific_diagnostic_ions: bool,
        modification_specific_neutral_losses: bool,
        glycan: Option<Vec<NeutralLoss>>,
        ppm: MassOverCharge,
    ) -> Self {
        Self {
            a,
            b,
            c,
            d,
            v,
            w,
            x,
            y,
            z,
            precursor,
            immonium,
            m,
            modification_specific_neutral_losses,
            modification_specific_diagnostic_ions,
            glycan,
            ppm,
        }
    }

    /// Generate all possible fragments
    pub fn all() -> Self {
        Self {
            a: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            b: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            c: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            d: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            v: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            w: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            x: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            z: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            immonium: true,
            m: true,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: Some(vec![
                NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap()),
                NeutralLoss::Loss(molecular_formula!(H 4 O 2).unwrap()),
            ]),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    /// Generate no fragments (except for precursor)
    pub fn none() -> Self {
        Self {
            a: (Location::None, vec![]),
            b: (Location::None, vec![]),
            c: (Location::None, vec![]),
            d: (Location::None, vec![]),
            v: (Location::None, vec![]),
            w: (Location::None, vec![]),
            x: (Location::None, vec![]),
            y: (Location::None, vec![]),
            z: (Location::None, vec![]),
            precursor: vec![],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: false,
            modification_specific_diagnostic_ions: false,
            glycan: None,
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    /// electron-transfer/higher-energy collisional dissociation
    pub fn ethcd() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            c: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            z: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: Some(vec![
                NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap()),
                NeutralLoss::Loss(molecular_formula!(H 4 O 2).unwrap()),
            ]),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    /// CID Hcd
    pub fn cid_hcd() -> Self {
        Self {
            a: (
                Location::TakeN { skip: 1, take: 1 },
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            b: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            c: (Location::None, Vec::new()),
            d: (
                Location::TakeN { skip: 1, take: 1 },
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()),
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            z: (Location::None, Vec::new()),
            precursor: vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: None,
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    /// ETD
    pub fn etd() -> Self {
        Self {
            a: (Location::None, Vec::new()),
            b: (Location::None, Vec::new()),
            c: (
                Location::SkipC(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            d: (Location::None, Vec::new()),
            v: (Location::None, Vec::new()),
            w: (Location::None, Vec::new()),
            x: (Location::None, Vec::new()),
            y: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            z: (
                Location::SkipN(1),
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap())],
            ),
            precursor: vec![
                NeutralLoss::Loss(molecular_formula!(H 2 O 1).unwrap()),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1).unwrap()),
            ],
            immonium: false,
            m: false,
            modification_specific_neutral_losses: true,
            modification_specific_diagnostic_ions: true,
            glycan: None,
            ppm: MassOverCharge::new::<mz>(20.0),
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
    pub const fn possible(&self, index: usize, length: usize) -> bool {
        match self {
            Self::SkipN(n) => index >= *n,
            Self::SkipNC(n, c) => index >= *n && length - index > *c,
            Self::TakeN { skip, take } => index >= *skip && index < *skip + *take,
            Self::SkipC(n) => length - index > *n,
            Self::TakeC(n) => length - index <= *n,
            Self::All => true,
            Self::None => false,
        }
    }
}
