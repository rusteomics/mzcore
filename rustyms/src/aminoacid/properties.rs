//! All amino acid property classes according to IMGT.
//! > IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties
//! >
//! > Christelle Pommié, Séverine Levadoux, Robert Sabatier, Gérard Lefranc, Marie-Paule Lefranc
//! >
//! > <https://doi.org/10.1002/jmr.647>
//!
//! Adapted to include J,B,Z,U,O, and X.
#![allow(missing_docs)]

use serde::{Deserialize, Serialize};

/// All amino acid property classes according to IMGT.
/// > IMGT standardized criteria for statistical analysis of immunoglobulin V-REGION amino acid properties
/// >
/// > Christelle Pommié, Séverine Levadoux, Robert Sabatier, Gérard Lefranc, Marie-Paule Lefranc
/// >
/// > <https://doi.org/10.1002/jmr.647>
///
/// Adapted to include J,B,Z,U,O, and X.
impl crate::AminoAcid {
    pub const fn physiochemical_class(self) -> PhysiochemicalClass {
        match self {
            Self::Alanine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Valine => PhysiochemicalClass::Aliphatic,
            Self::Proline => PhysiochemicalClass::Proline,
            Self::Tyrosine => PhysiochemicalClass::Tyrosine,
            Self::Tryptophan => PhysiochemicalClass::Tryptophan,
            Self::Glycine => PhysiochemicalClass::Glycine,
            Self::Phenylalanine => PhysiochemicalClass::Phenylalanine,
            Self::Cysteine | Self::Methionine => PhysiochemicalClass::Sulphur,
            Self::Serine | Self::Threonine => PhysiochemicalClass::Hydroxyl,
            Self::AsparticAcid | Self::GlutamicAcid => PhysiochemicalClass::Acidic,
            Self::Asparagine | Self::Glutamine => PhysiochemicalClass::Amide,
            Self::Arginine | Self::Histidine | Self::Lysine => PhysiochemicalClass::Basic,
            Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Selenocysteine
            | Self::Unknown => PhysiochemicalClass::Unknown,
        }
    }
    pub const fn chemical_class(self) -> ChemicalClass {
        match self {
            Self::Alanine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Proline
            | Self::Glycine
            | Self::Valine => ChemicalClass::Aliphatic,
            Self::Tyrosine | Self::Tryptophan | Self::Phenylalanine => ChemicalClass::Aromatic,
            Self::Cysteine | Self::Methionine => ChemicalClass::Sulphur,
            Self::Serine | Self::Threonine => ChemicalClass::Hydroxyl,
            Self::AsparticAcid | Self::GlutamicAcid => ChemicalClass::Acidic,
            Self::Asparagine | Self::Glutamine => ChemicalClass::Amide,
            Self::Arginine | Self::Histidine | Self::Lysine => ChemicalClass::Basic,
            Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Selenocysteine
            | Self::Unknown => ChemicalClass::Unknown,
        }
    }
    pub const fn volume_class(self) -> VolumeClass {
        match self {
            Self::Alanine | Self::Glycine | Self::Serine => VolumeClass::VerySmall,
            Self::Proline
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::AmbiguousAsparagine
            | Self::Threonine
            | Self::Cysteine
            | Self::Selenocysteine => VolumeClass::Small,
            Self::GlutamicAcid
            | Self::Glutamine
            | Self::AmbiguousGlutamine
            | Self::Histidine
            | Self::Valine => VolumeClass::Medium,
            Self::Arginine
            | Self::Lysine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine => VolumeClass::Large,
            Self::Tyrosine | Self::Tryptophan | Self::Phenylalanine | Self::Pyrrolysine => {
                VolumeClass::VeryLarge
            }
            Self::Unknown => VolumeClass::Unknown,
        }
    }
    pub const fn hydropathy_class(self) -> HydropathyClass {
        match self {
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine => HydropathyClass::Hydrophobic,
            Self::Glycine
            | Self::Serine
            | Self::Histidine
            | Self::Proline
            | Self::Threonine
            | Self::Tyrosine => HydropathyClass::Neutral,
            Self::Lysine
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::AmbiguousAsparagine
            | Self::GlutamicAcid
            | Self::Glutamine
            | Self::AmbiguousGlutamine
            | Self::Pyrrolysine
            | Self::Arginine => HydropathyClass::Hydrophilic,
            Self::Unknown => HydropathyClass::Unknown,
        }
    }
    pub const fn charge_class(self) -> ChargeClass {
        match self {
            Self::Arginine | Self::Histidine | Self::Lysine => ChargeClass::Positive,
            Self::AsparticAcid | Self::GlutamicAcid => ChargeClass::Negative,
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine
            | Self::Glycine
            | Self::Serine
            | Self::Proline
            | Self::Threonine
            | Self::Tyrosine
            | Self::Asparagine
            | Self::Glutamine
            | Self::Pyrrolysine => ChargeClass::Uncharged,
            Self::Unknown | Self::AmbiguousAsparagine | Self::AmbiguousGlutamine => {
                ChargeClass::Unknown
            }
        }
    }
    pub const fn polarity_class(self) -> PolarityClass {
        match self {
            Self::Arginine
            | Self::Histidine
            | Self::Lysine
            | Self::Pyrrolysine
            | Self::AmbiguousAsparagine
            | Self::AmbiguousGlutamine
            | Self::AsparticAcid
            | Self::Asparagine
            | Self::Glutamine
            | Self::Threonine
            | Self::Tyrosine
            | Self::Serine
            | Self::GlutamicAcid => PolarityClass::Polar,
            Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Tryptophan
            | Self::Valine
            | Self::Glycine
            | Self::Proline => PolarityClass::Nonpolar,
            Self::Unknown => PolarityClass::Unknown,
        }
    }
    pub const fn hydrogen_bond_class(self) -> HydrogenBondClass {
        match self {
            Self::Arginine | Self::Lysine | Self::Tryptophan => HydrogenBondClass::Donor,
            Self::AsparticAcid | Self::GlutamicAcid => HydrogenBondClass::Acceptor,
            Self::Histidine
            | Self::Pyrrolysine
            | Self::Asparagine
            | Self::Threonine
            | Self::Tyrosine
            | Self::Serine => HydrogenBondClass::Both,
            Self::Glutamine
            | Self::Alanine
            | Self::Cysteine
            | Self::Selenocysteine
            | Self::Isoleucine
            | Self::AmbiguousLeucine
            | Self::Leucine
            | Self::Methionine
            | Self::Phenylalanine
            | Self::Valine
            | Self::Glycine
            | Self::Proline => HydrogenBondClass::None,
            Self::AmbiguousAsparagine | Self::AmbiguousGlutamine | Self::Unknown => {
                HydrogenBondClass::Unknown
            }
        }
    }
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum PhysiochemicalClass {
    Aliphatic,
    Proline,
    Tyrosine,
    Tryptophan,
    Glycine,
    Phenylalanine,
    Sulphur,
    Hydroxyl,
    Acidic,
    Amide,
    Basic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum ChemicalClass {
    Aliphatic,
    Sulphur,
    Hydroxyl,
    Acidic,
    Amide,
    Basic,
    Aromatic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum VolumeClass {
    VerySmall,
    Small,
    Medium,
    Large,
    VeryLarge,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum HydropathyClass {
    Hydrophobic,
    Neutral,
    Hydrophilic,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum ChargeClass {
    Uncharged,
    Positive,
    Negative,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum PolarityClass {
    Polar,
    Nonpolar,
    Unknown,
}

#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum HydrogenBondClass {
    Donor,
    Acceptor,
    Both,
    None,
    Unknown,
}
