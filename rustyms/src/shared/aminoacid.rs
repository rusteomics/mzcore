/// An amino acid, alongside the standard ones some [ambiguous (B/J/Z/X) and non-standard (U/O)](https://www.insdc.org/submitting-standards/feature-table/#7.4.3) are included.
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[expect(missing_docs)]
pub enum AminoAcid {
    #[default]
    Alanine = 0,
    Arginine,
    Asparagine,
    AsparticAcid,
    Cysteine,
    Glutamine,
    GlutamicAcid,
    Glycine,
    Histidine,
    Isoleucine,
    Leucine,
    Lysine,
    Methionine,
    Phenylalanine,
    Proline,
    Serine,
    Threonine,
    Tryptophan,
    Tyrosine,
    Valine,
    AmbiguousAsparagine,
    AmbiguousLeucine,
    AmbiguousGlutamine,
    Selenocysteine,
    Pyrrolysine,
    Unknown,
}
//ARNDCQEGHILKMFPSTWYVBJZUOX

/// The error that a given sequence is not a valid codon
#[derive(Debug)]
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
        sequence_index: crate::SequencePosition,
        peptidoform_index: usize,
    ) -> Multi<MolecularFormula> {
        let crate::SequencePosition::Index(sequence_index) = sequence_index else {
            panic!("Not allowed to call amino acid formulas with a terminal sequence index")
        };
        match self {
            Self::Alanine => molecular_formula!(H 5 C 3 O 1 N 1).into(),
            Self::Arginine => molecular_formula!(H 12 C 6 O 1 N 4).into(), // One of the H's counts as the charge carrier and is added later
            Self::Asparagine => molecular_formula!(H 6 C 4 O 2 N 2).into(),
            Self::AsparticAcid => molecular_formula!(H 5 C 4 O 3 N 1).into(),
            Self::AmbiguousAsparagine => vec![
                molecular_formula!(H 6 C 4 O 2 N 2 (crate::AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index})),
                molecular_formula!(H 5 C 4 O 3 N 1 (crate::AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Cysteine => molecular_formula!(H 5 C 3 O 1 N 1 S 1).into(),
            Self::Glutamine => molecular_formula!(H 8 C 5 O 2 N 2).into(),
            Self::GlutamicAcid => molecular_formula!(H 7 C 5 O 3 N 1).into(),
            Self::AmbiguousGlutamine => vec![
                molecular_formula!(H 8 C 5 O 2 N 2 (crate::AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index})),
                molecular_formula!(H 7 C 5 O 3 N 1 (crate::AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index})),
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
