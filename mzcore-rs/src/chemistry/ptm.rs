use std::fmt::Formatter;
use std::str::FromStr;
use anyhow::*;
use serde::{Deserialize, Serialize};

// TODO: implement Eq&Hash
#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
struct AminoAcidPtm {
    pub id: i64,
    pub name: String,
    pub formula: String,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub position_constraint: PtmLocation, // any N-term, any C-term, protein N-term, protein C-term
    pub residue_constraint: u8
}

impl AminoAcidPtm {
    pub fn new(id: i64, name: &str, formula: &str, mono_mass: f64, average_mass: f64, position_constraint: PtmLocation, residue_constraint: u8) -> anyhow::Result<AminoAcidPtm> {
        if name.is_empty() { bail!("name is empty") }
        if formula.is_empty() { bail!("formula is empty") }

        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(AminoAcidPtm {
            id: id,
            name: name.to_string(),
            formula: formula.to_string(),
            mono_mass: mono_mass,
            average_mass: average_mass,
            position_constraint: position_constraint,
            residue_constraint: residue_constraint,
        })
    }
}

// --- PtmLocation definition --- //

/// Enumeration that represents different PTM locations.
/// The values of this enum are serialized to strings matching the Unimod specification.
/// See the `position_t` type in the
/// [Unimod 2.0 Schema](https://www.unimod.org/xmlns/schema/unimod_2/unimod_2.xsd).
///
/// ```text
/// <xs:simpleType name="position_t">
/// <xs:restriction base="xs:string">
///     <xs:enumeration value="Anywhere"/>
///     <xs:enumeration value="Any N-term"/>
///     <xs:enumeration value="Any C-term"/>
///     <xs:enumeration value="Protein N-term"/>
///     <xs:enumeration value="Protein C-term"/>
/// </xs:restriction>
/// </xs:simpleType>
/// ```
///
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum PtmLocation {
    ProteinNTerm, // = Value("Protein N-term")
    ProteinCTerm, // = Value("Protein C-term")
    AnyNTerm, // = Value("Any N-term")
    AnyCTerm, // = Value("Any C-term")
    Anywhere, // = Value("Anywhere")
}

const PROTEIN_NTERM: &'static str = "Protein N-term";
const PROTEIN_CTERM: &'static str = "Protein C-term";
const ANY_NTERM: &'static str = "Any N-term";
const ANY_CTERM: &'static str = "Any C-term";
const ANYWHERE: &'static str = "Anywhere";

impl PtmLocation {
    pub fn to_str(&self) -> &'static str {
        (*self).into()
    }
}

impl std::fmt::Display for PtmLocation {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_str())
    }
}

impl From<PtmLocation> for &'static str {
    fn from(e: PtmLocation) -> Self {
        match e {
            PtmLocation::ProteinNTerm => PROTEIN_NTERM,
            PtmLocation::ProteinCTerm => PROTEIN_CTERM,
            PtmLocation::AnyNTerm => ANY_NTERM,
            PtmLocation::AnyCTerm => ANY_CTERM,
            PtmLocation::Anywhere => ANYWHERE,
        }
    }
}

impl FromStr for PtmLocation {
    type Err = Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            PROTEIN_NTERM => Ok(Self::ProteinNTerm),
            PROTEIN_CTERM => Ok(Self::ProteinCTerm),
            ANY_NTERM => Ok(Self::AnyNTerm),
            ANY_CTERM => Ok(Self::AnyCTerm),
            ANYWHERE => Ok(Self::Anywhere),
            _ => Err(anyhow!("Unknown element {}", s)),
        }
    }
}