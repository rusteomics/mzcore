use std::borrow::Cow;
use std::fmt::Formatter;
use std::str::FromStr;

use anyhow::*;
use serde::{Deserialize, Serialize};

use crate::chemistry::api::*;
use crate::chemistry::table::proteinogenic_amino_acid_table;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[repr(u8)]
pub enum AA {
    // The 20 standard AAs
    A = b'A',
    C = b'C',
    D = b'D',
    E = b'E',
    F = b'F',
    G = b'G',
    H = b'H',
    I = b'I',
    K = b'K',
    L = b'L',
    M = b'M',
    N = b'N',
    P = b'P',
    Q = b'Q',
    R = b'R',
    S = b'S',
    T = b'T',
    V = b'V',
    W = b'W',
    Y = b'Y',

    // Proteinogenic AAs
    U = b'U',
    O = b'O',

    // Ambiguous AAs
    B = b'B',
    J = b'J',
    X = b'X',
    Z = b'Z',
}

impl AA {
    pub fn to_str(&self) -> &'static str {
        (*self).into()
    }

    pub fn definition(&self) -> &AminoAcidDefinition {
        let aa_as_byte = *self as u8;

        proteinogenic_amino_acid_table().aa_from_byte(&aa_as_byte).unwrap()
    }
}

impl std::fmt::Display for AA {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_str())
    }
}

impl From<AA> for &'static str {
    fn from(aa: AA) -> Self {
        match aa {
            AA::A => "A",
            AA::C => "C",
            AA::D => "D",
            AA::E => "E",
            AA::F => "F",
            AA::G => "G",
            AA::H => "H",
            AA::I => "I",
            AA::K => "K",
            AA::L => "L",
            AA::M => "M",
            AA::N => "N",
            AA::P => "P",
            AA::Q => "Q",
            AA::R => "R",
            AA::S => "S",
            AA::T => "T",
            AA::V => "V",
            AA::W => "W",
            AA::Y => "Y",
            AA::U => "U",
            AA::O => "O",
            AA::B => "B",
            AA::J => "J",
            AA::X => "X",
            AA::Z => "Z",
        }
    }
}

impl FromStr for AA {
    type Err = Error;

    #[allow(clippy::too_many_lines)]
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value {
            "A" => Ok(AA::A),
            "C" => Ok(AA::C),
            "D" => Ok(AA::D),
            "E" => Ok(AA::E),
            "F" => Ok(AA::F),
            "G" => Ok(AA::G),
            "H" => Ok(AA::H),
            "I" => Ok(AA::I),
            "K" => Ok(AA::K),
            "L" => Ok(AA::L),
            "M" => Ok(AA::M),
            "N" => Ok(AA::N),
            "P" => Ok(AA::P),
            "Q" => Ok(AA::Q),
            "R" => Ok(AA::R),
            "S" => Ok(AA::S),
            "T" => Ok(AA::T),
            "V" => Ok(AA::V),
            "W" => Ok(AA::W),
            "Y" => Ok(AA::Y),
            "U" => Ok(AA::U),
            "O" => Ok(AA::O),
            "B" => Ok(AA::B),
            "J" => Ok(AA::J),
            "X" => Ok(AA::X),
            "Z" => Ok(AA::Z),
            _ => Err(anyhow::anyhow!("Unknown AA {}", value)),
        }
    }
}

impl TryFrom<u8> for AA {
    type Error = Error;

    #[allow(clippy::too_many_lines)]
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'A' => Ok(AA::A),
            b'C' => Ok(AA::C),
            b'D' => Ok(AA::D),
            b'E' => Ok(AA::E),
            b'F' => Ok(AA::F),
            b'G' => Ok(AA::G),
            b'H' => Ok(AA::H),
            b'I' => Ok(AA::I),
            b'K' => Ok(AA::K),
            b'L' => Ok(AA::L),
            b'M' => Ok(AA::M),
            b'N' => Ok(AA::N),
            b'P' => Ok(AA::P),
            b'Q' => Ok(AA::Q),
            b'R' => Ok(AA::R),
            b'S' => Ok(AA::S),
            b'T' => Ok(AA::T),
            b'V' => Ok(AA::V),
            b'W' => Ok(AA::W),
            b'Y' => Ok(AA::Y),
            b'U' => Ok(AA::U),
            b'O' => Ok(AA::O),
            b'B' => Ok(AA::B),
            b'J' => Ok(AA::J),
            b'X' => Ok(AA::X),
            b'Z' => Ok(AA::Z),
            _ => Err(anyhow::anyhow!("Unknown AA {}", value)),
        }
    }
}

// Kind of dangerous because it implicitly means that u8 <-> amino acid (it might not be desired)
impl HasMass for u8 {
    fn mono_mass(&self) -> f64 {
        let aa: AA = (*self).try_into().unwrap();
        aa.definition().mono_mass
    }
    fn average_mass(&self) -> Option<f64> {
        let aa: AA = (*self).try_into().unwrap();
        Some(aa.definition().average_mass)
    }
}

impl HasNameAndSymbol for u8 {
    fn name(&self) -> Cow<str> {
        let aa: AA = (*self).try_into().unwrap();
        Cow::from(aa.definition().name().to_string())
    }
    fn symbol(&self) -> Cow<str> {
        let aa: AA = (*self).try_into().unwrap();
        Cow::from(aa.definition().symbol().to_string())
    }
}

impl IsAminoAcid for u8 {
    fn is_valid(&self) -> bool {
        let aa_res: Result<AA> = (*self).try_into();
        aa_res.is_ok()
    }
    fn single_letter_code(&self) -> u8 {
        let aa: AA = (*self).try_into().unwrap();
        aa.definition().single_letter_code()
    }
    fn three_letter_code(&self) -> Cow<str> {
        let aa: AA = (*self).try_into().unwrap();
        Cow::from(aa.definition().three_letter_code().to_string())
    }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct AminoAcidDefinition {
    pub code1: u8, // FIXME: shall we use AA enum type instead?
    pub code3: String,
    pub name: String,
    pub formula: Option<String>,
    pub mono_mass: f64,
    pub average_mass: f64,
    pub occurrence: f32, // occurrence in human proteins
    pub pka1: f32, // C-term pKa
    pub pka2: f32, // N-term pKa
    pub pka3: f32, // side chain pKa
    pub pi: f32, // isoelectric point
    pub codons: Vec<String>
}

const ACUG_STR: &'static str = "ACUG";

impl AminoAcidDefinition {
    pub fn new(
        code1: u8,
        code3: &str,
        name: &str,
        formula: &str,
        mono_mass: f64,
        average_mass: f64,
        codons: Vec<String>
    ) -> Result<AminoAcidDefinition> {

        if code3.len() < 3 { bail!("code3 must contain three characters") }
        if name.is_empty() { bail!("name is empty") }
        if formula.is_empty() { bail!("formula is empty") }

        for codon in &*codons {
            if codon.len() != 3 { bail!("a codon must contain three characters") }
            let only_acug_chars = codon.chars().all(|c| ACUG_STR.contains(c));
            if !only_acug_chars { bail!("a codon must only contain ACUG letters") }
        }

        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(AminoAcidDefinition {
            code1: code1,
            code3: code3.to_string(),
            name: name.to_string(),
            formula: Some(formula.to_string()),
            mono_mass: mono_mass,
            average_mass: average_mass,
            occurrence: 0.0,
            pka1: 0.0,
            pka2: 0.0,
            pka3: 0.0,
            pi: 0.0,
            codons: codons //.iter().map(|s| *s.to_string()).collect()
        })
    }
}

impl HasMass for AminoAcidDefinition {
    fn mono_mass(&self) -> f64 {
        self.mono_mass
    }
    fn average_mass(&self) -> Option<f64> {
        Some(self.average_mass)
    }
}

impl HasNameAndSymbol for AminoAcidDefinition {
    fn name(&self) -> Cow<str> {
        Cow::Borrowed(&self.name)
    }
    fn symbol(&self) -> Cow<str> {
        Cow::Borrowed(&self.code3)
    }
}

impl IsAminoAcid for AminoAcidDefinition {
    fn is_valid(&self) -> bool {
        true
    }
    fn single_letter_code(&self) -> u8 {
        self.code1
    }
    fn three_letter_code(&self) -> Cow<str> {
        Cow::Borrowed(&self.code3)
    }
}


