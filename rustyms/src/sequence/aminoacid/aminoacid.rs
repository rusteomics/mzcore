//! Module used define the implementations for the [`IsAminoAcid`] trait

use std::{borrow::Cow, collections::HashMap, sync::LazyLock};

use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::*,
    chemistry::{AmbiguousLabel, CachedCharge, MolecularFormula, MultiChemical},
    fragment::{Fragment, FragmentKind, FragmentType, PeptidePosition, SatelliteLabel},
    molecular_formula,
    quantities::Multi,
    sequence::SequencePosition,
};

use super::is_amino_acid::IsAminoAcid;

impl std::fmt::Display for dyn IsAminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pro_forma_definition())
    }
}

/// An amino acid, alongside the standard ones some [ambiguous (B/J/Z/X) and non-standard (U/O)](https://www.insdc.org/submitting-standards/feature-table/#7.4.3) are included.
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum AminoAcid {
    /// Ala, A
    #[default]
    Alanine = 0,
    /// Arg, R
    Arginine,
    /// Asn, N
    Asparagine,
    /// Asp, D
    AsparticAcid,
    /// Cys, C
    Cysteine,
    /// Gln, Q
    Glutamine,
    /// Glu, E
    GlutamicAcid,
    /// Gly, G
    Glycine,
    /// His, H
    Histidine,
    /// Ile, I
    Isoleucine,
    /// Leu, L
    Leucine,
    /// Lys, K
    Lysine,
    /// Met, M
    Methionine,
    /// Phe, F
    Phenylalanine,
    /// Pro, P
    Proline,
    /// Ser, S
    Serine,
    /// Thr, T
    Threonine,
    /// Trp, W
    Tryptophan,
    /// Tyr, Y
    Tyrosine,
    /// Val, V
    Valine,
    /// Asx, B
    AmbiguousAsparagine,
    /// Xle, J
    AmbiguousLeucine,
    /// Glx, Z
    AmbiguousGlutamine,
    /// Sec, U
    Selenocysteine,
    /// Pyl, O
    Pyrrolysine,
    /// Xxx, X
    Unknown,
}

/// The error that a given sequence is not a valid codon
#[derive(Clone, Copy, Debug)]
pub struct NotACodon;

impl std::fmt::Display for NotACodon {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Not a valid codon")
    }
}

impl std::error::Error for NotACodon {}

#[allow(dead_code)]
impl AminoAcid {
    /// The total number of amino acids
    pub const TOTAL_NUMBER: usize = Self::Unknown as usize + 1;
    /// Translate the dna codon into the corresponding amino acid according to the standard DNA codon table.
    /// It returns None for a stop codon.
    /// <https://en.wikipedia.org/wiki/DNA_and_RNA_codon_tables>
    /// # Errors
    /// It returns `Err(NotACodon)` when the given codon is not a valid dna codon.
    #[allow(dead_code)]
    pub fn from_dna(dna: &str) -> Result<Option<Self>, NotACodon> {
        match dna.to_lowercase().as_str() {
            "ttt" | "ttc" => Ok(Some(Self::Phenylalanine)),
            "tta" | "ttg" | "ctt" | "ctc" | "cta" | "ctg" => Ok(Some(Self::Leucine)),
            "att" | "atc" | "ata" => Ok(Some(Self::Isoleucine)),
            "atg" => Ok(Some(Self::Methionine)),
            "gtt" | "gtc" | "gta" | "gtg" => Ok(Some(Self::Valine)),
            "tct" | "tcc" | "tca" | "tcg" | "agt" | "agc" => Ok(Some(Self::Serine)),
            "cct" | "ccc" | "cca" | "ccg" => Ok(Some(Self::Proline)),
            "act" | "acc" | "aca" | "acg" => Ok(Some(Self::Threonine)),
            "gct" | "gcc" | "gca" | "gcg" => Ok(Some(Self::Alanine)),
            "tat" | "tac" => Ok(Some(Self::Tyrosine)),
            "taa" | "tag" | "tga" => Ok(None),
            "cat" | "cac" => Ok(Some(Self::Histidine)),
            "caa" | "cag" => Ok(Some(Self::Glutamine)),
            "aat" | "aac" => Ok(Some(Self::Asparagine)),
            "cgt" | "cgc" | "cga" | "cgg" | "aga" | "agg" => Ok(Some(Self::Arginine)),
            "aaa" | "aag" => Ok(Some(Self::Lysine)),
            "gat" | "gac" => Ok(Some(Self::AsparticAcid)),
            "gaa" | "gag" => Ok(Some(Self::GlutamicAcid)),
            "tgt" | "tgc" => Ok(Some(Self::Cysteine)),
            "tgg" => Ok(Some(Self::Tryptophan)),
            "ggt" | "ggc" | "gga" | "ggg" => Ok(Some(Self::Glycine)),
            _ => Err(NotACodon),
        }
    }
}

impl std::str::FromStr for AminoAcid {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::try_from(s)
    }
}

impl TryFrom<&str> for AminoAcid {
    type Error = ();
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        if value.is_ascii() && value.len() == 1 {
            let ch = value.chars().next().unwrap();
            ch.try_into()
        } else {
            Err(())
        }
    }
}

impl TryFrom<String> for AminoAcid {
    type Error = ();
    fn try_from(value: String) -> Result<Self, Self::Error> {
        value.as_str().try_into()
    }
}

impl TryFrom<&String> for AminoAcid {
    type Error = ();
    fn try_from(value: &String) -> Result<Self, Self::Error> {
        value.as_str().try_into()
    }
}

impl TryFrom<char> for AminoAcid {
    type Error = ();
    fn try_from(value: char) -> Result<Self, Self::Error> {
        if value.is_ascii() {
            let num = value as u8;
            num.try_into()
        } else {
            Err(())
        }
    }
}

impl TryFrom<&u8> for AminoAcid {
    type Error = ();
    fn try_from(value: &u8) -> Result<Self, Self::Error> {
        match value {
            b'A' | b'a' => Ok(Self::Alanine),
            b'B' | b'b' => Ok(Self::AmbiguousAsparagine),
            b'C' | b'c' => Ok(Self::Cysteine),
            b'D' | b'd' => Ok(Self::AsparticAcid),
            b'E' | b'e' => Ok(Self::GlutamicAcid),
            b'F' | b'f' => Ok(Self::Phenylalanine),
            b'G' | b'g' => Ok(Self::Glycine),
            b'H' | b'h' => Ok(Self::Histidine),
            b'I' | b'i' => Ok(Self::Isoleucine),
            b'J' | b'j' => Ok(Self::AmbiguousLeucine),
            b'K' | b'k' => Ok(Self::Lysine),
            b'L' | b'l' => Ok(Self::Leucine),
            b'M' | b'm' => Ok(Self::Methionine),
            b'N' | b'n' => Ok(Self::Asparagine),
            b'O' | b'o' => Ok(Self::Pyrrolysine),
            b'P' | b'p' => Ok(Self::Proline),
            b'Q' | b'q' => Ok(Self::Glutamine),
            b'R' | b'r' => Ok(Self::Arginine),
            b'S' | b's' => Ok(Self::Serine),
            b'T' | b't' => Ok(Self::Threonine),
            b'U' | b'u' => Ok(Self::Selenocysteine),
            b'V' | b'v' => Ok(Self::Valine),
            b'W' | b'w' => Ok(Self::Tryptophan),
            b'X' | b'x' => Ok(Self::Unknown),
            b'Y' | b'y' => Ok(Self::Tyrosine),
            b'Z' | b'z' => Ok(Self::AmbiguousGlutamine),
            _ => Err(()),
        }
    }
}

impl TryFrom<u8> for AminoAcid {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        Self::try_from(&value)
    }
}

impl MultiChemical for AminoAcid {
    /// Get all possible formulas for an amino acid (has one for all except B/Z has two for these)
    /// # Panics
    /// Is the sequence index is a terminal index
    fn formulas_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Multi<MolecularFormula> {
        let SequencePosition::Index(sequence_index) = sequence_index else {
            panic!("Not allowed to call amino acid formulas with a terminal sequence index")
        };
        match self {
            Self::Alanine => molecular_formula!(H 5 C 3 O 1 N 1).into(),
            Self::Arginine => molecular_formula!(H 12 C 6 O 1 N 4).into(), // One of the H's counts as the charge carrier and is added later
            Self::Asparagine => molecular_formula!(H 6 C 4 O 2 N 2).into(),
            Self::AsparticAcid => molecular_formula!(H 5 C 4 O 3 N 1).into(),
            Self::AmbiguousAsparagine => vec![
                molecular_formula!(H 6 C 4 O 2 N 2 (AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                molecular_formula!(H 5 C 4 O 3 N 1 (AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
            ]
            .into(),
            Self::Cysteine => molecular_formula!(H 5 C 3 O 1 N 1 S 1).into(),
            Self::Glutamine => molecular_formula!(H 8 C 5 O 2 N 2).into(),
            Self::GlutamicAcid => molecular_formula!(H 7 C 5 O 3 N 1).into(),
            Self::AmbiguousGlutamine => vec![
                molecular_formula!(H 8 C 5 O 2 N 2 (AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                molecular_formula!(H 7 C 5 O 3 N 1 (AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
            ]
            .into(),
            Self::Glycine => molecular_formula!(H 3 C 2 O 1 N 1).into(),
            Self::Histidine => molecular_formula!(H 7 C 6 O 1 N 3).into(),
            Self::AmbiguousLeucine | Self::Isoleucine | Self::Leucine => {
                molecular_formula!(H 11 C 6 O 1 N 1).into()
            }
            Self::Lysine => molecular_formula!(H 12 C 6 O 1 N 2).into(),
            Self::Methionine => molecular_formula!(H 9 C 5 O 1 N 1 S 1).into(),
            Self::Phenylalanine => molecular_formula!(H 9 C 9 O 1 N 1).into(),
            Self::Proline => molecular_formula!(H 7 C 5 O 1 N 1).into(),
            Self::Pyrrolysine => molecular_formula!(H 19 C 11 O 2 N 3).into(),
            Self::Selenocysteine => molecular_formula!(H 5 C 3 O 1 N 1 Se 1).into(),
            Self::Serine => molecular_formula!(H 5 C 3 O 2 N 1).into(),
            Self::Threonine => molecular_formula!(H 7 C 4 O 2 N 1).into(),
            Self::Tryptophan => molecular_formula!(H 10 C 11 O 1 N 2).into(),
            Self::Tyrosine => molecular_formula!(H 9 C 9 O 2 N 1).into(),
            Self::Valine => molecular_formula!(H 9 C 5 O 1 N 1).into(),
            Self::Unknown => molecular_formula!().into(),
        }
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pro_forma_definition())
    }
}

impl IsAminoAcid for AminoAcid {
    /// Get the single letter representation of the amino acid
    fn one_letter_code(&self) -> Option<char> {
        Some(match self {
            Self::Alanine => 'A',
            Self::AmbiguousAsparagine => 'B',
            Self::Cysteine => 'C',
            Self::AsparticAcid => 'D',
            Self::GlutamicAcid => 'E',
            Self::Phenylalanine => 'F',
            Self::Glycine => 'G',
            Self::Histidine => 'H',
            Self::Isoleucine => 'I',
            Self::AmbiguousLeucine => 'J',
            Self::Lysine => 'K',
            Self::Leucine => 'L',
            Self::Methionine => 'M',
            Self::Asparagine => 'N',
            Self::Pyrrolysine => 'O',
            Self::Proline => 'P',
            Self::Glutamine => 'Q',
            Self::Arginine => 'R',
            Self::Serine => 'S',
            Self::Threonine => 'T',
            Self::Selenocysteine => 'U',
            Self::Valine => 'V',
            Self::Tryptophan => 'W',
            Self::Unknown => 'X',
            Self::Tyrosine => 'Y',
            Self::AmbiguousGlutamine => 'Z',
        })
    }

    fn pro_forma_definition(&self) -> Cow<'_, str> {
        Cow::Borrowed(match self {
            Self::Alanine => "A",
            Self::AmbiguousAsparagine => "B",
            Self::Cysteine => "C",
            Self::AsparticAcid => "D",
            Self::GlutamicAcid => "E",
            Self::Phenylalanine => "F",
            Self::Glycine => "G",
            Self::Histidine => "H",
            Self::Isoleucine => "I",
            Self::AmbiguousLeucine => "J",
            Self::Lysine => "K",
            Self::Leucine => "L",
            Self::Methionine => "M",
            Self::Asparagine => "N",
            Self::Pyrrolysine => "O",
            Self::Proline => "P",
            Self::Glutamine => "Q",
            Self::Arginine => "R",
            Self::Serine => "S",
            Self::Threonine => "T",
            Self::Selenocysteine => "U",
            Self::Valine => "V",
            Self::Tryptophan => "W",
            Self::Unknown => "X",
            Self::Tyrosine => "Y",
            Self::AmbiguousGlutamine => "Z",
        })
    }

    /// Get the 3 letter code for the amino acid
    fn three_letter_code(&self) -> Option<Cow<'_, str>> {
        Some(Cow::Borrowed(match self {
            Self::Alanine => "Ala",
            Self::AmbiguousAsparagine => "Asx",
            Self::Cysteine => "Cys",
            Self::AsparticAcid => "Asp",
            Self::GlutamicAcid => "Glu",
            Self::Phenylalanine => "Phe",
            Self::Glycine => "Gly",
            Self::Histidine => "His",
            Self::Isoleucine => "Ile",
            Self::AmbiguousLeucine => "Xle",
            Self::Lysine => "Lys",
            Self::Leucine => "Leu",
            Self::Methionine => "Met",
            Self::Asparagine => "Asn",
            Self::Pyrrolysine => "Pyl",
            Self::Proline => "Pro",
            Self::Glutamine => "Gln",
            Self::Arginine => "Arg",
            Self::Serine => "Ser",
            Self::Threonine => "Thr",
            Self::Selenocysteine => "Sec",
            Self::Valine => "Val",
            Self::Tryptophan => "Trp",
            Self::Unknown => "Xaa",
            Self::Tyrosine => "Tyr",
            Self::AmbiguousGlutamine => "Glx",
        }))
    }

    /// Get the full name for the amino acid
    fn name(&self) -> Cow<'_, str> {
        Cow::Borrowed(match self {
            Self::Alanine => "Alanine",
            Self::AmbiguousAsparagine => "AmbiguousAsparagine",
            Self::Cysteine => "Cysteine",
            Self::AsparticAcid => "AsparticAcid",
            Self::GlutamicAcid => "GlutamicAcid",
            Self::Phenylalanine => "Phenylalanine",
            Self::Glycine => "Glycine",
            Self::Histidine => "Histidine",
            Self::Isoleucine => "Isoleucine",
            Self::AmbiguousLeucine => "AmbiguousLeucine",
            Self::Lysine => "Lysine",
            Self::Leucine => "Leucine",
            Self::Methionine => "Methionine",
            Self::Asparagine => "Asparagine",
            Self::Pyrrolysine => "Pyrrolysine",
            Self::Proline => "Proline",
            Self::Glutamine => "Glutamine",
            Self::Arginine => "Arginine",
            Self::Serine => "Serine",
            Self::Threonine => "Threonine",
            Self::Selenocysteine => "Selenocysteine",
            Self::Valine => "Valine",
            Self::Tryptophan => "Tryptophan",
            Self::Unknown => "Unknown",
            Self::Tyrosine => "Tyrosine",
            Self::AmbiguousGlutamine => "AmbiguousGlutamine",
        })
    }

    fn side_chain(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Cow<'_, Multi<MolecularFormula>> {
        let SequencePosition::Index(sequence_index) = sequence_index else {
            return Cow::Owned(Multi::default());
        };
        Cow::Owned(match self {
            Self::Alanine => molecular_formula!(H 3 C 1).into(),
            Self::Arginine => molecular_formula!(H 10 C 4 N 3).into(), // One of the H's counts as the charge carrier and is added later
            Self::Asparagine => molecular_formula!(H 4 C 2 O 1 N 1).into(),
            Self::AsparticAcid => molecular_formula!(H 3 C 2 O 2).into(),
            Self::AmbiguousAsparagine => vec![
                molecular_formula!(H 4 C 2 O 1 N 1 (AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                molecular_formula!(H 3 C 2 O 2 (AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
            ]
            .into(),
            Self::Cysteine => molecular_formula!(H 3 C 1 S 1).into(),
            Self::Glutamine => molecular_formula!(H 6 C 3 O 1 N 1).into(),
            Self::GlutamicAcid => molecular_formula!(H 5 C 3 O 2).into(),
            Self::AmbiguousGlutamine => vec![
                molecular_formula!(H 6 C 3 O 1 N 1 (AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                molecular_formula!(H 5 C 3 O 2 (AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
            ]
            .into(),
            Self::Glycine => molecular_formula!(H 1).into(),
            Self::Histidine => molecular_formula!(H 5 C 4 N 2).into(),
            Self::AmbiguousLeucine | Self::Isoleucine | Self::Leucine => {
                molecular_formula!(H 9 C 4).into()
            }
            Self::Lysine => molecular_formula!(H 10 C 4 N 1).into(),
            Self::Methionine => molecular_formula!(H 7 C 3 S 1).into(),
            Self::Phenylalanine => molecular_formula!(H 7 C 7).into(),
            Self::Proline => molecular_formula!(H 5 C 3).into(),
            Self::Pyrrolysine => molecular_formula!(H 17 C 9 O 1 N 2).into(),
            Self::Selenocysteine => molecular_formula!(H 3 C 1 Se 1).into(),
            Self::Serine => molecular_formula!(H 3 C 1 O 1).into(),
            Self::Threonine => molecular_formula!(H 5 C 2 O 1).into(),
            Self::Tryptophan => molecular_formula!(H 8 C 9 N 1).into(),
            Self::Tyrosine => molecular_formula!(H 7 C 7 O 1).into(),
            Self::Valine => molecular_formula!(H 7 C 3).into(),
            Self::Unknown => molecular_formula!().into(),
        })
    }

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    fn satellite_ion_fragments(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Option<Cow<'_, Vec<(SatelliteLabel, MolecularFormula)>>> {
        let SequencePosition::Index(sequence_index) = sequence_index else {
            return None;
        };

        match self {
            Self::Alanine
            | Self::Glycine
            | Self::Histidine
            | Self::Phenylalanine
            | Self::Proline
            | Self::Tryptophan
            | Self::Tyrosine
            | Self::Unknown => None,
            Self::Arginine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 9 C 2 N 2),
            )])),
            Self::Asparagine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 2 C 1 N 1 O 1),
            )])),
            Self::AsparticAcid => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 1 C 1 O 2),
            )])),
            Self::AmbiguousAsparagine => Some(Cow::Owned(vec![
                (
                    SatelliteLabel::None,
                    molecular_formula!(H 2 C 1 N 1 O 1 (AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
                (
                    SatelliteLabel::None,
                    molecular_formula!(H 1 C 1 O 2 (AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
            ])),
            Self::Cysteine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 1 S 1),
            )])),
            Self::Glutamine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 4 C 2 N 1 O 1),
            )])),
            Self::GlutamicAcid => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 3 C 2 O 2),
            )])),
            Self::AmbiguousGlutamine => Some(Cow::Owned(vec![
                (
                    SatelliteLabel::None,
                    molecular_formula!(H 4 C 2 N 1 O 1 (AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
                (
                    SatelliteLabel::None,
                    molecular_formula!(H 3 C 2 O 2 (AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
            ])),
            Self::Isoleucine => Some(Cow::Owned(vec![
                (SatelliteLabel::A, molecular_formula!(H 3 C 1)),
                (SatelliteLabel::B, molecular_formula!(H 5 C 2)),
            ])),
            Self::Leucine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 7 C 3),
            )])),
            Self::AmbiguousLeucine => Some(Cow::Owned(vec![
                (
                    SatelliteLabel::A,
                    molecular_formula!(H 3 C 1 (AmbiguousLabel::AminoAcid{option: Self::Isoleucine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
                (
                    SatelliteLabel::B,
                    molecular_formula!(H 5 C 2 (AmbiguousLabel::AminoAcid{option: Self::Isoleucine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
                (
                    SatelliteLabel::None,
                    molecular_formula!(H 7 C 3 (AmbiguousLabel::AminoAcid{option: Self::Leucine, sequence_index, peptidoform_index, peptidoform_ion_index})),
                ),
            ])),
            Self::Lysine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 8 C 3 N 1),
            )])),
            Self::Methionine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 5 C 2 S 1),
            )])),
            Self::Pyrrolysine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 15 C 9 N 2 O 1),
            )])),
            Self::Selenocysteine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(Se 1),
            )])),
            Self::Serine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 1 O 1),
            )])),
            Self::Threonine => Some(Cow::Owned(vec![
                (SatelliteLabel::None, molecular_formula!(H 1 O 1)),
                (SatelliteLabel::None, molecular_formula!(H 3 C 1)),
            ])),
            Self::Valine => Some(Cow::Owned(vec![(
                SatelliteLabel::None,
                molecular_formula!(H 3 C 1),
            )])), // Technically two options, but both have the same mass
        }
    }
}

/// The mass of the backbone of an amino acid
pub static BACKBONE: LazyLock<MolecularFormula> =
    LazyLock::new(|| molecular_formula!(H 3 C 2 N 1 O 1));

impl AminoAcid {
    /// All amino acids with a unique mass (no I/L in favour of J, no B, no Z, and no X)
    pub const UNIQUE_MASS_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::AmbiguousLeucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
        Self::Selenocysteine,
        Self::Pyrrolysine,
    ];

    /// All 20 canonical amino acids
    pub const CANONICAL_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::Leucine,
        Self::Isoleucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
    ];

    /// All amino acids (including I/L/J/B/Z but excluding X)
    pub const ALL_AMINO_ACIDS: &'static [Self] = &[
        Self::Alanine,
        Self::AmbiguousAsparagine,
        Self::AmbiguousGlutamine,
        Self::AmbiguousLeucine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::GlutamicAcid,
        Self::Glutamine,
        Self::Glycine,
        Self::Histidine,
        Self::Isoleucine,
        Self::Leucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Pyrrolysine,
        Self::Selenocysteine,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
    ];

    // TODO: generalise over used storage type, so using molecularformula, monoisotopic mass, or average mass, also make sure that AAs can return these numbers in a const fashion
    #[expect(clippy::too_many_lines, clippy::too_many_arguments)]
    pub(crate) fn fragments(
        self,
        n_term: &(
            Multi<MolecularFormula>,
            HashMap<FragmentKind, Multi<MolecularFormula>>,
        ),
        c_term: &(
            Multi<MolecularFormula>,
            HashMap<FragmentKind, Multi<MolecularFormula>>,
        ),
        modifications: &(
            Multi<MolecularFormula>,
            HashMap<FragmentKind, Multi<MolecularFormula>>,
        ),
        charge_carriers: &mut CachedCharge,
        sequence_index: SequencePosition,
        sequence_length: usize,
        ions: &PossibleIons,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        allow_terminal: (bool, bool),
    ) -> Vec<Fragment> {
        let mut base_fragments = Vec::with_capacity(ions.size_upper_bound());
        let n_pos = PeptidePosition::n(sequence_index, sequence_length);
        let c_pos = PeptidePosition::c(sequence_index, sequence_length);

        if allow_terminal.0 {
            if let Some(settings) = &ions.a {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::a)
                        .unwrap_or(&modifications.0)
                        - molecular_formula!(H 1 C 1 O 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::a(n_pos, 0),
                    n_term.1.get(&FragmentKind::a).unwrap_or(&n_term.0),
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.b {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::b)
                        .unwrap_or(&modifications.0)
                        - molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::b(n_pos, 0),
                    n_term.1.get(&FragmentKind::b).unwrap_or(&n_term.0),
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.c {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::c)
                        .unwrap_or(&modifications.0)
                        + molecular_formula!(H 2 N 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::c(n_pos, 0),
                    n_term.1.get(&FragmentKind::c).unwrap_or(&n_term.0),
                    charge_carriers,
                    settings,
                ));
            }
            for (aa, distance) in &ions.d.0 {
                if let Some(satellite_fragments) = aa.satellite_ion_fragments(
                    sequence_index - *distance,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) {
                    for (label, formula) in satellite_fragments.iter() {
                        base_fragments.extend(Fragment::generate_series(
                            &(modifications
                                .1
                                .get(&FragmentKind::d)
                                .unwrap_or(&modifications.0)
                                * self.formulas_inner(
                                    sequence_index,
                                    peptidoform_index,
                                    peptidoform_ion_index,
                                )
                                + molecular_formula!(H 1 C 1 O 1)
                                - formula),
                            peptidoform_ion_index,
                            peptidoform_index,
                            &FragmentType::d(n_pos, *aa, *distance, 0, *label),
                            n_term.1.get(&FragmentKind::d).unwrap_or(&n_term.0),
                            charge_carriers,
                            &ions.d.1,
                        ));
                    }
                }
            }
        }
        if allow_terminal.1 {
            for (aa, distance) in &ions.v.0 {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * -aa.formulas_inner(
                        sequence_index + *distance,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) + LazyLock::force(&BACKBONE)),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::v(c_pos, *aa, *distance, 0),
                    c_term.1.get(&FragmentKind::v).unwrap_or(&c_term.0),
                    charge_carriers,
                    &ions.v.1,
                ));
            }
            for (aa, distance) in &ions.w.0 {
                if let Some(satellite_fragments) = aa.satellite_ion_fragments(
                    sequence_index - *distance,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) {
                    for (label, formula) in satellite_fragments.iter() {
                        base_fragments.extend(Fragment::generate_series(
                            &(modifications
                                .1
                                .get(&FragmentKind::w)
                                .unwrap_or(&modifications.0)
                                * self.formulas_inner(
                                    sequence_index,
                                    peptidoform_index,
                                    peptidoform_ion_index,
                                )
                                + molecular_formula!(H 2 N 1)
                                - formula),
                            peptidoform_ion_index,
                            peptidoform_index,
                            &FragmentType::w(c_pos, *aa, *distance, 0, *label),
                            c_term.1.get(&FragmentKind::w).unwrap_or(&c_term.0),
                            charge_carriers,
                            &ions.w.1,
                        ));
                    }
                }
            }
            if let Some(settings) = &ions.x {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::x)
                        .unwrap_or(&modifications.0)
                        + molecular_formula!(C 1 O 1)
                        - molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::x(c_pos, 0),
                    c_term.1.get(&FragmentKind::x).unwrap_or(&c_term.0),
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.y {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::y)
                        .unwrap_or(&modifications.0)
                        + molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::y(c_pos, 0),
                    c_term.1.get(&FragmentKind::y).unwrap_or(&c_term.0),
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.z {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::z)
                        .unwrap_or(&modifications.0)
                        - molecular_formula!(H 2 N 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::z(c_pos, 0),
                    c_term.1.get(&FragmentKind::z).unwrap_or(&c_term.0),
                    charge_carriers,
                    settings,
                ));
            }
        }

        if allow_terminal.0 && allow_terminal.1 {
            if let Some((charge, losses)) = &ions.immonium {
                base_fragments.extend(Fragment::generate_all(
                    &(self.formulas_inner(
                        sequence_index,
                        peptidoform_index,
                        peptidoform_ion_index,
                    ) * (modifications
                        .1
                        .get(&FragmentKind::immonium)
                        .unwrap_or(&modifications.0)
                        - molecular_formula!(C 1 O 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::Immonium(n_pos, self.into()), // TODO: get the actual sequence element here
                    &Multi::default(),
                    &losses
                        .iter()
                        .filter(|(aa, _)| aa.contains(&self))
                        .flat_map(|(_, l)| l.iter())
                        .map(|l| vec![l.clone()])
                        .collect::<Vec<_>>(),
                    charge_carriers,
                    *charge,
                ));
            }
        }
        base_fragments
    }

    /// Check if two amino acids are considered identical. X is identical to anything, J to IL, B to ND, Z to EQ.
    pub(crate) fn canonical_identical(self, rhs: Self) -> bool {
        match (self, rhs) {
            (a, b) if a == b => true,
            (Self::Unknown, _)
            | (_, Self::Unknown)
            | (Self::AmbiguousLeucine, Self::Leucine | Self::Isoleucine)
            | (Self::Leucine | Self::Isoleucine, Self::AmbiguousLeucine)
            | (Self::AmbiguousAsparagine, Self::Asparagine | Self::AsparticAcid)
            | (Self::Asparagine | Self::AsparticAcid, Self::AmbiguousAsparagine)
            | (Self::AmbiguousGlutamine, Self::Glutamine | Self::GlutamicAcid)
            | (Self::Glutamine | Self::GlutamicAcid, Self::AmbiguousGlutamine) => true,
            _ => false,
        }
    }
}

#[cfg(test)]
#[expect(clippy::unreadable_literal, clippy::missing_panics_doc)]
mod tests {
    use super::*;

    #[test]
    fn mass() {
        let weight_ala = AminoAcid::Alanine.formulas()[0].average_weight();
        let mass_ala = AminoAcid::Alanine.formulas()[0].monoisotopic_mass();
        assert_ne!(weight_ala, mass_ala);
        assert!((weight_ala.value - 71.07793).abs() < 1e-5);
        assert!((mass_ala.value - 71.037113783).abs() < 1e-5);
    }

    #[test]
    fn mass_lysine() {
        let weight_lys = AminoAcid::Lysine.formulas()[0].average_weight();
        let mass_lys = AminoAcid::Lysine.formulas()[0].monoisotopic_mass();
        assert_ne!(weight_lys, mass_lys);
        assert!((weight_lys.value - 128.17240999999999).abs() < 1e-5);
        assert!((mass_lys.value - 128.094963010536).abs() < 1e-5);
    }

    #[test]
    fn masses() {
        let known = &[
            ('A', 71.03711, 71.08),
            ('R', 156.10111, 156.2),
            ('N', 114.04293, 114.1),
            ('D', 115.02694, 115.1),
            ('C', 103.00919, 103.1),
            ('E', 129.04259, 129.1),
            ('Q', 128.05858, 128.1),
            ('G', 57.02146, 57.05),
            ('H', 137.05891, 137.1),
            ('I', 113.08406, 113.2),
            ('L', 113.08406, 113.2),
            ('K', 128.09496, 128.2),
            ('M', 131.04049, 131.2),
            ('F', 147.06841, 147.2),
            ('P', 97.05276, 97.12),
            ('S', 87.03203, 87.08),
            ('T', 101.04768, 101.1),
            ('W', 186.07931, 186.2),
            ('Y', 163.06333, 163.2),
            ('V', 99.06841, 99.13),
        ];

        for (aa, mono_mass, average_weight) in known {
            let aa = AminoAcid::try_from(*aa).unwrap();
            let (mono, weight) = (
                aa.formulas()[0].monoisotopic_mass().value,
                aa.formulas()[0].average_weight().value,
            );
            println!(
                "{}: {} {} {} {}",
                aa.pro_forma_definition(),
                mono,
                mono_mass,
                weight,
                average_weight
            );
            assert!((mono - *mono_mass).abs() < 1e-5);
            assert!((weight - *average_weight).abs() < 1e-1);
        }
    }

    #[test]
    fn read_aa() {
        assert_eq!(
            AminoAcid::try_from('B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(
            AminoAcid::try_from(b'B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(AminoAcid::try_from('c'), Ok(AminoAcid::Cysteine));
        assert_eq!(AminoAcid::try_from('🦀'), Err(()));
    }
}
