use std::num::NonZeroU8;

/// Controlled vocabularies used in mass spectrometry data files
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[repr(u8)]
#[allow(clippy::upper_case_acronyms)]
pub enum ControlledVocabulary {
    /// The PSI-MS Controlled Vocabulary [https://www.ebi.ac.uk/ols4/ontologies/ms](https://www.ebi.ac.uk/ols4/ontologies/ms)
    MS = 1,
    /// The Unit Ontology [https://www.ebi.ac.uk/ols4/ontologies/uo](https://www.ebi.ac.uk/ols4/ontologies/uo)
    UO,
    /// The Experimental Factor Ontology <https://www.ebi.ac.uk/ols4/ontologies/efo>
    EFO,
    /// The Ontology for Biomedical Investigations <https://www.ebi.ac.uk/ols4/ontologies/obi>
    OBI,
    /// The Human Ancestry Ontology <https://www.ebi.ac.uk/ols4/ontologies/hancestro>
    HANCESTRO,
    /// The Basic Formal Ontology <https://www.ebi.ac.uk/ols4/ontologies/bfo>
    BFO,
    /// The NCI Thesaurus OBO Edition <https://www.ebi.ac.uk/ols4/ontologies/ncit>
    NCIT,
    /// The BRENDA Tissue Ontology <https://www.ebi.ac.uk/ols4/ontologies/bto>
    BTO,
    /// The PRIDE Controlled Vocabulary <https://www.ebi.ac.uk/ols4/ontologies/pride>
    PRIDE,
    GNO,
    XLMOD,
    RESID,
    /// PSI-MOD
    MOD,
    Unknown,
}

impl std::str::FromStr for ControlledVocabulary {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "MS" | "PSI-MS" => Ok(Self::MS),
            "UO" => Ok(Self::UO),
            "EFO" => Ok(Self::EFO),
            "OBI" => Ok(Self::OBI),
            "HANCESTRO" => Ok(Self::HANCESTRO),
            "BFO" => Ok(Self::BFO),
            "NCIT" => Ok(Self::NCIT),
            "BTO" => Ok(Self::BTO),
            "PRIDE" => Ok(Self::PRIDE),
            "GNO" | "GNOME" | "G" => Ok(Self::GNO),
            "XLMOD" => Ok(Self::XLMOD),
            "RESID" => Ok(Self::RESID),
            "MOD" | "PSI-MOD" => Ok(Self::MOD),
            _ => Ok(Self::Unknown),
        }
    }
}

impl std::fmt::Display for ControlledVocabulary {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::MS => "MS",
                Self::UO => "UO",
                Self::EFO => "EFO",
                Self::OBI => "OBI",
                Self::HANCESTRO => "HANCESTRO",
                Self::BFO => "BFO",
                Self::NCIT => "NCIT",
                Self::BTO => "BTO",
                Self::PRIDE => "PRIDE",
                Self::GNO => "GNO",
                Self::XLMOD => "XLMOD",
                Self::RESID => "RESID",
                Self::MOD => "MOD",
                Self::Unknown => "?",
            }
        )
    }
}

/// A CURIE is a namespace + accession identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Curie {
    pub cv: ControlledVocabulary,
    pub accession: AccessionCode,
}

#[macro_export]
macro_rules! curie {
    ($ns:ident:$acc:tt) => {
        $crate::Curie {
            cv: $crate::ControlledVocabulary::$ns,
            accession: $crate::accession_code!($acc),
        }
    };
}

impl std::fmt::Display for Curie {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.cv, self.accession)
    }
}

#[derive(Debug)]
pub enum CURIEParsingError {
    UnknownControlledVocabulary,
    AccessionParsingError(AccessionCodeParseError),
    MissingNamespaceSeparator,
}

impl std::str::FromStr for Curie {
    type Err = CURIEParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((cv, accession)) = s.split_once(':').or(s.split_once('_')) {
            let cv = cv.parse::<ControlledVocabulary>().unwrap(); // Unknown CVs are handled with the CV::Unknown option
            let accession = accession
                .parse()
                .map_err(|e| CURIEParsingError::AccessionParsingError(e))?;
            Ok(Curie { cv, accession })
        } else {
            Err(CURIEParsingError::MissingNamespaceSeparator)
        }
    }
}

/// An accession code, Can either be a numeric code (u32 to 4 milion, so 9 fully utilised digits).
/// Or it can be an ASCII alphanumeric code (case-sensitive) of 1 to 8 characters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AccessionCode {
    Numeric(u32),
    Alphanumeric(NonZeroU8, [u8; 7]),
}

impl std::fmt::Display for AccessionCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Numeric(n) => write!(f, "{n:07}"),
            Self::Alphanumeric(first, bytes) => {
                write!(
                    f,
                    "{}{}",
                    char::from(first.get()),
                    std::str::from_utf8(bytes.trim_ascii_end()).unwrap()
                )
            }
        }
    }
}

/// Create an accession code, using this with a numeric code is match arm valid, using this with an alphanumeric code is const compatible.
/// ```rust no-run
/// match curie!(MS:000111) {
///     curie!(MS:0000111) => (),
///     x if x == curie!(MS:H000111G) => (),
///     _ => (),
/// }
/// ```
#[macro_export]
macro_rules! accession_code {
    ($acc:literal) => {
        // #[allow(clippy::zero_prefixed_literal)] cannot put allows on single expressions yet
        $crate::AccessionCode::Numeric($acc)
    };
    ($acc:tt) => {
        const {
            const BYTES: &[u8] = stringify!($acc).as_bytes();
            assert!(BYTES[0] != 0);
            assert!(BYTES.is_ascii());
            match BYTES {
                [a, b, c, d, e, f, g, h] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, *g, *h],
                ),
                [a, b, c, d, e, f, g] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, *g, 0x20],
                ),
                [a, b, c, d, e, f] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, *f, 0x20, 0x20],
                ),
                [a, b, c, d, e] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, *e, 0x20, 0x20, 0x20],
                ),
                [a, b, c, d] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, *d, 0x20, 0x20, 0x20, 0x20],
                ),
                [a, b, c] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, *c, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                [a, b] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [*b, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                [a] => $crate::AccessionCode::Alphanumeric(
                    std::num::NonZeroU8::new(*a).unwrap(),
                    [0x20, 0x20, 0x20, 0x20, 0x20, 0x20, 0x20],
                ),
                _ => panic!(concat!(
                    "Cannot convert ",
                    stringify!($acc),
                    " to accession code. Can only handle strings of 1 to 8 bytes"
                )),
            }
        }
    };
}

#[derive(Debug)]
pub enum AccessionCodeParseError {
    InvalidCharacters(String),
    TooLong(String),
    Empty,
}

impl std::str::FromStr for AccessionCode {
    type Err = AccessionCodeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.is_empty() {
            return Err(AccessionCodeParseError::Empty);
        }
        if s.chars().any(|c| !c.is_alphanumeric()) {
            return Err(AccessionCodeParseError::InvalidCharacters(s.to_string()));
        }
        s.parse::<u32>().map_or_else(
            |_| {
                if s.len() > 8 {
                    return Err(AccessionCodeParseError::TooLong(s.to_string()));
                }
                let mut bytes = [0x20; 7]; // Space
                bytes[..s.len() - 1].copy_from_slice(&s.as_bytes()[1..]);
                curie!(NCIT:P383);
                Ok(Self::Alphanumeric(
                    NonZeroU8::new(s.as_bytes()[0]).unwrap(),
                    bytes,
                ))
            },
            |u| Ok(Self::Numeric(u)),
        )
    }
}

fn test() {
    match curie!(MS:000111) {
        curie!(MS:0000111) => (),
        x if x == curie!(MS:H000111G) => (),
        _ => (),
    }
}

#[cfg(test)]
mod tests {
    use crate::Curie;

    #[test]
    fn parse_curies() {
        let options = [
            "NCIT:C25330",
            "NCIT:R100",
            "NCIT:P383",
            "MS:1000014",
            "UO:0000245",
            "BFO:0000015",
            "BAO_0000925", // Seen in one Obo file, but also quite common in URLs
            "GAZ:00002933",
            "AfPO_0000233",
            "BTO:0004947",
            "PRIDE:0000521",
            "MOD:01188",
            "XLMOD:07097",
            "GNO:G00001NT",
            "GNO:00000015",
        ];
        for option in options {
            let _curie: Curie = option.parse().unwrap();
        }
    }

    #[test]
    fn curie_macro() {
        assert_eq!("NCIT:C25330".parse::<Curie>().unwrap(), curie!(NCIT:C25330));
        assert_eq!("NCIT:R100".parse::<Curie>().unwrap(), curie!(NCIT:R100));
        assert_eq!("NCIT:P383".parse::<Curie>().unwrap(), curie!(NCIT:P383));
        assert_eq!("MS:1000014".parse::<Curie>().unwrap(), curie!(MS:1000014));
        assert_eq!("UO:0000245".parse::<Curie>().unwrap(), curie!(UO:0000245));
        assert_eq!("BFO:0000015".parse::<Curie>().unwrap(), curie!(BFO:0000015));
        assert_eq!("BTO:0004947".parse::<Curie>().unwrap(), curie!(BTO:0004947));
        assert_eq!(
            "PRIDE:0000521".parse::<Curie>().unwrap(),
            curie!(PRIDE:0000521)
        );
        assert_eq!("MOD:01188".parse::<Curie>().unwrap(), curie!(MOD:01188));
        assert_eq!("XLMOD:07097".parse::<Curie>().unwrap(), curie!(XLMOD:07097));
        assert_eq!(
            "GNO:G00001NT".parse::<Curie>().unwrap(),
            curie!(GNO:G00001NT)
        );
        assert_eq!(
            "GNO:00000015".parse::<Curie>().unwrap(),
            curie!(GNO:00000015)
        );
    }
}
