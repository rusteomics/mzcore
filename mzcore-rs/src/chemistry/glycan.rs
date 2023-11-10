#![allow(dead_code)]
/// This code originates from [rustyms](https://github.com/snijderlab/rustyms)
/// Copyright (c) 2023 Douwe Schulte and contributors
/// MIT License

use std::{fmt::Display, hash::Hash, ops::Range, sync::OnceLock};

use serde::{Deserialize, Serialize};

use crate::chemistry::{
    api::Chemical,
    composition::*,
    //helper_functions::*,
    element::{Element, ELEMENTS_SORTED_FOR_PARSING},
};
use crate::common::error::Context;
use crate::common::error::CustomError;

pub type GlycanComposition = Vec<(MonoSaccharide, i16)>;

/// A monosaccharide with all its complexity
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct MonoSaccharide {
    pub(super) base_sugar: BaseSugar,
    pub(super) substituents: Vec<GlycanSubstituent>,
    pub(super) furanose: bool,
    pub(super) proforma_name: Option<String>,
}

impl MonoSaccharide {
    /// Create a new monosaccharide
    pub fn new(sugar: BaseSugar, substituents: &[GlycanSubstituent]) -> Self {
        molecular_formula!();
        Self {
            base_sugar: sugar,
            substituents: substituents.to_owned(),
            furanose: false,
            proforma_name: None,
        }
    }

    /// Get this same monosaccharide but now with the given pro forma name
    #[must_use]
    #[allow(dead_code)]
    pub fn with_name(self, name: &str) -> Self {
        Self {
            proforma_name: Some(name.to_string()),
            ..self
        }
    }

    /// Set this saccharide up as to be a furanose
    #[must_use]
    #[allow(dead_code)]
    pub fn furanose(self) -> Self {
        Self {
            furanose: true,
            ..self
        }
    }

    /// Parse a short IUPAC name from this string starting at `start` and returning,
    /// if successful, a monosaccharide and the offset in the string where parsing ended.
    /// # Errors
    /// Fails if it finds a structure that does not fit the IUPAC glycan name.
    pub fn from_short_iupac(
        line: &str,
        start: usize,
        line_number: usize,
    ) -> Result<(Self, usize), CustomError> {
        let mut index = start;
        let bytes = line.as_bytes();
        let mut substituents = Vec::new();

        // ignore stuff
        index += line[index..].ignore(&["keto-"]);
        index += line[index..].ignore(&["D-", "L-", "?-"]);
        // Prefix mods
        let mut amount = 1;
        if bytes[index].is_ascii_digit() {
            match bytes[index + 1] {
                b',' if bytes[index + 3] == b':' => {
                    let start_index = index;
                    index += 7;
                    index += line[index..].ignore(&["-"]);
                    if !line[index..].starts_with("Anhydro") {
                        return Err(CustomError::error(
                            "Invalid iupac monosaccharide name",
                            "This internally linked glycan could not be parsed, expected Anhydro as modification",
                            Context::Line {
                                linenumber: line_number,
                                line: line.to_string(),
                                offset: start_index,
                                length: index-start_index+5,
                            },
                        ));
                    }
                    index += 7;
                    substituents.extend_from_slice(&[
                        GlycanSubstituent::Didehydro,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Deoxy,
                    ]);
                }
                b',' => {
                    let num = bytes[index + 1..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                        .count();
                    index += num + 1;
                    amount = num / 2;
                    // X,X{mod} (or 3/4/5/etc mods)
                }
                _ => index += 1, // X{mod}
            }
            index += line[index..].ignore(&["-"]);
        }
        // Detect & ignore epi state
        index += line[index..].ignore(&["e"]);
        // Get the prefix mods
        if let Some(o) = line[index..].take_any(PREFIX_SUBSTITUENTS, |e| {
            substituents.extend(std::iter::repeat(e.clone()).take(amount));
        }) {
            index += o;
        }
        index += line[index..].ignore(&["-"]);
        // Another optional isomeric state
        index += line[index..].ignore(&["D-", "L-", "?-"]);
        // Base sugar
        let mut sugar = None;
        for sug in BASE_SUGARS {
            if line[index..].starts_with(sug.0) {
                index += sug.0.len();
                sugar = Some((sug.1.clone(), sug.2));
                break;
            }
        }
        let mut sugar = sugar
            .map(|(b, s)| {
                let mut alo = Self {
                    base_sugar: b,
                    substituents,
                    furanose: false,
                    proforma_name: None,
                };
                alo.substituents.extend(s.iter().cloned());
                alo
            })
            .ok_or_else(|| {
                CustomError::error(
                    "Invalid iupac monosaccharide name",
                    "This name could not be recognised as a standard iupac glycan name",
                    Context::Line {
                        linenumber: line_number,
                        line: line.to_string(),
                        offset: index,
                        length: 3,
                    },
                )
            })?;
        // Furanose
        if index < bytes.len() && bytes[index] == b'f' {
            index += 1;
            sugar.furanose = true;
        }
        // Postfix mods
        while index < bytes.len() {
            index += line[index..].ignore(&["-"]);
            let mut single_amount = 0;
            let mut double_amount = 0;
            // Location
            let (offset, mut amount, mut double) = line[index..].parse_location();
            index += offset;
            if double {
                double_amount = amount;
            } else {
                single_amount = amount;
            }
            if bytes[index] == b':' {
                // additional place
                let (offset, amt, dbl) = line[index + 1..].parse_location();
                index += offset + 1;
                amount += amt;
                if double {
                    double_amount += amount;
                } else {
                    single_amount += amount;
                }
                double |= dbl; // if any is double
            }

            index += line[index..].ignore(&["-"]);
            index += line[index..].ignore(&["(X)", "(R)", "(S)"]);
            if double {
                if let Some(o) = line[index..].take_any(DOUBLE_LINKED_POSTFIX_SUBSTITUENTS, |e| {
                    sugar.substituents.extend(
                        e.iter()
                            .flat_map(|s| std::iter::repeat(s).take(double_amount))
                            .cloned(),
                    );
                    if single_amount > 0 {
                        sugar.substituents.extend(
                            e.iter()
                                .filter(|s| **s != GlycanSubstituent::Water)
                                .flat_map(|s| std::iter::repeat(s).take(single_amount))
                                .cloned(),
                        );
                    }
                }) {
                    index += o;
                } else {
                    return Err(CustomError::error(
                        "Invalid iupac monosaccharide name",
                        "No detected double linked glycan substituent was found, while the pattern for location is for a double linked substituent",
                        Context::Line {
                            linenumber: line_number,
                            line: line.to_string(),
                            offset: index,
                            length: 2,
                        },
                    ));
                }
            } else {
                // Mod or an element
                if let Some(o) = line[index..].take_any(POSTFIX_SUBSTITUENTS, |e| {
                    sugar
                        .substituents
                        .extend(std::iter::repeat(e.clone()).take(amount));
                }) {
                    index += o;
                } else if let Some(o) = line[index..].take_any(ELEMENTS_SORTED_FOR_PARSING, |e| {
                    sugar
                        .substituents
                        .extend(std::iter::repeat(GlycanSubstituent::Element(*e)).take(amount));
                }) {
                    index += o;
                } else {
                    break;
                }
            }
            // Amount
            if amount != 1 {
                // Ignore the amount number, already determined before
                index += 1;
            }
        }
        Ok((sugar, index))
    }

    // fn symbol(&self) -> char {
    //     // ⬠◇♢▭◮⬘
    //     // ⬟◆♦▬
    //     match self {
    //         Self::Hexose => '○',      //●
    //         Self::HexNAc => '□',      // ■
    //         Self::Deoxyhexose => '△', // ▲
    //         Self::Pentose => '☆',     // ★
    //         Self::HexN => '⬔',
    //         _ => '⬡', // ⬢
    //     }
    // }
}

trait ParseHelper {
    fn ignore(self, ignore: &[&str]) -> usize;
    fn take_any<T>(self, parse_list: &[(&str, T)], f: impl FnMut(&T)) -> Option<usize>;
    fn parse_location(self) -> (usize, usize, bool);
}

impl ParseHelper for &str {
    /// Ignore any of the given things, greedily ignores the first match
    fn ignore(self, ignore: &[&str]) -> usize {
        for i in ignore {
            if self.starts_with(i) {
                return i.len();
            }
        }
        0
    }

    fn take_any<T>(self, parse_list: &[(&str, T)], mut f: impl FnMut(&T)) -> Option<usize> {
        let mut found = None;
        for element in parse_list {
            if self.starts_with(element.0) {
                found = Some(element.0.len());
                f(&element.1);
                break;
            }
        }
        found
    }

    // Get a location, return the new index, the amount of the mod to place and if it is doubly linked or not
    fn parse_location(self) -> (usize, usize, bool) {
        let bytes = self.as_bytes();
        let mut index = 0;
        let mut amount = 1;
        let mut double = false;
        if bytes[0].is_ascii_digit() || bytes[0] == b'?' {
            match bytes[1] {
                b',' => {
                    let num = bytes[1..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                        .count();
                    index += num + 1;
                    amount = num / 2 + 1;
                    // X,X{mod} (or 3/4/5/etc mods)
                }
                b'-' if bytes[index] != b'?' => {
                    index += 7;
                    double = true;
                } // X-X,X-X (Py)
                b'/' => {
                    let num = bytes[2..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b'/')
                        .count();
                    index += num + 2;
                    // X/X/X...{mod} multiple possible locations
                }
                c if (c.is_ascii_digit() || c == b'?') && bytes[0] == b'?' => {
                    if bytes[2] == b',' {
                        let num = bytes[2..]
                            .iter()
                            .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                            .count();
                        index += num + 2;
                        amount = num / 2 + 1;
                        // ?X,X{mod} (or 3/4/5/etc mods)
                    } else if bytes[2] == b'/' {
                        let num = bytes[3..]
                            .iter()
                            .take_while(|c| c.is_ascii_digit() || **c == b'/')
                            .count();
                        index += num + 3;
                        // ?X/X{mod} multiple possible locations
                    } else {
                        index += 2; // ?X{mod}
                    }
                }
                _ => index += 1, // X{mod}
            }
        }
        (index, amount, double)
    }
}

impl Chemical for MonoSaccharide {
    fn composition(&self) -> ElementalComposition {
        self.base_sugar.composition() + self.substituents.as_slice().composition()
    }
}

impl Display for MonoSaccharide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.proforma_name.clone().unwrap_or_else(|| format!(
                "{}{}{}",
                self.base_sugar,
                if self.furanose { "f" } else { "" },
                self.substituents
                    .iter()
                    .map(ToString::to_string)
                    .collect::<String>()
            ))
        )
    }
}

/// The base sugar of a monosaccharide, optionally with the isomeric state saved as well.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum BaseSugar {
    /// Edge case, no sugar at all, because ProForma enforces that a separate phosphate and sulphate have to be handled.
    None,
    /// 2 carbon base sugar
    Sugar,
    /// 3 carbon base sugar
    Triose,
    /// 4 carbon base sugar
    Tetrose(Option<TetroseIsomer>),
    /// 5 carbon base sugar
    Pentose(Option<PentoseIsomer>),
    /// 6 carbon base sugar
    Hexose(Option<HexoseIsomer>),
    /// 7 carbon base sugar
    Heptose(Option<HeptoseIsomer>),
    /// 8 carbon base sugar
    Octose,
    /// 9 carbon base sugar
    Nonose,
    /// 10 carbon base sugar
    Decose,
}

impl Display for BaseSugar {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::None => "None",
                Self::Sugar => "Sug",
                Self::Triose => "Tri",
                Self::Tetrose(_) => "Tet",
                Self::Pentose(_) => "Pen",
                Self::Hexose(_) => "Hex",
                Self::Heptose(_) => "Hep",
                Self::Octose => "Oct",
                Self::Nonose => "Non",
                Self::Decose => "Dec",
            }
        )
    }
}

impl Chemical for BaseSugar {
    fn composition(&self) -> ElementalComposition {
        match self {
            Self::None => ElementalComposition::default(),
            Self::Sugar => molecular_formula!(H 2 C 2 O 1),
            Self::Triose => molecular_formula!(H 4 C 3 O 2),
            Self::Tetrose(_) => molecular_formula!(H 6 C 4 O 3),
            Self::Pentose(_) => molecular_formula!(H 8 C 5 O 4),
            Self::Hexose(_) => molecular_formula!(H 10 C 6 O 5),
            Self::Heptose(_) => molecular_formula!(H 12 C 7 O 6),
            Self::Octose => molecular_formula!(H 14 C 8 O 7),
            Self::Nonose => molecular_formula!(H 16 C 9 O 8),
            Self::Decose => molecular_formula!(H 18 C 10 O 9),
        }
    }
}

/// Any 4 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum TetroseIsomer {
    /// Ery
    Erythrose,
    /// Tho
    Threose,
}

/// Any 5 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PentoseIsomer {
    /// Rib
    Ribose,
    /// Ara
    Arabinose,
    /// Xyl
    Xylose,
    /// Lyx
    Lyxose,
    /// Xul
    Xylulose,
}

/// Any 6 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HexoseIsomer {
    /// glc
    Glucose,
    /// Gal
    Galactose,
    /// Man
    Mannose,
    /// All
    Allose,
    /// Alt
    Altrose,
    /// Gul
    Gulose,
    /// Ido
    Idose,
    /// Tal
    Talose,
    /// Psi
    Psicose,
    /// Fru
    Fructose,
    /// Sor
    Sorbose,
    /// Tag
    Tagatose,
}

/// Any 7 carbon glycan
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum HeptoseIsomer {
    /// gro-manHep
    GlyceroMannoHeptopyranose, // TODO: Does this indicate some mods?
    /// Sed
    Sedoheptulose,
}

const BASE_SUGARS: &[(&str, BaseSugar, &[GlycanSubstituent])] = &[
    ("Sug", BaseSugar::Sugar, &[]),
    ("Tri", BaseSugar::Triose, &[]),
    ("Tet", BaseSugar::Tetrose(None), &[]),
    (
        "Ery",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[],
    ),
    ("Tho", BaseSugar::Tetrose(Some(TetroseIsomer::Threose)), &[]),
    ("Pen", BaseSugar::Pentose(None), &[]),
    ("Rib", BaseSugar::Pentose(Some(PentoseIsomer::Ribose)), &[]),
    (
        "Ara",
        BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)),
        &[],
    ),
    ("Xyl", BaseSugar::Pentose(Some(PentoseIsomer::Xylose)), &[]),
    ("Lyx", BaseSugar::Pentose(Some(PentoseIsomer::Lyxose)), &[]),
    ("Hex", BaseSugar::Hexose(None), &[]),
    ("Glc", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[]),
    ("Gal", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[]),
    ("Man", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[]),
    ("All", BaseSugar::Hexose(Some(HexoseIsomer::Allose)), &[]),
    ("Alt", BaseSugar::Hexose(Some(HexoseIsomer::Altrose)), &[]),
    ("Gul", BaseSugar::Hexose(Some(HexoseIsomer::Gulose)), &[]),
    ("Ido", BaseSugar::Hexose(Some(HexoseIsomer::Idose)), &[]),
    ("Tal", BaseSugar::Hexose(Some(HexoseIsomer::Talose)), &[]),
    ("Hep", BaseSugar::Heptose(None), &[]),
    (
        "gro-manHep",
        BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)),
        &[],
    ),
    (
        "Neu",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Sia",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Kdn",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Kdo",
        BaseSugar::Octose,
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Acid],
    ),
    (
        "Fuc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Rha",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Qui",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "Oli",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Tyv",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Asc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Abe",
        BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Par",
        BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Dig",
        BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "Col",
        BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    ("Psi", BaseSugar::Hexose(Some(HexoseIsomer::Psicose)), &[]),
    ("Fru", BaseSugar::Hexose(Some(HexoseIsomer::Fructose)), &[]),
    ("Sor", BaseSugar::Hexose(Some(HexoseIsomer::Sorbose)), &[]),
    ("Tag", BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[]),
    (
        "Xul",
        BaseSugar::Pentose(Some(PentoseIsomer::Xylulose)),
        &[],
    ),
    (
        "Sed",
        BaseSugar::Heptose(Some(HeptoseIsomer::Sedoheptulose)),
        &[],
    ),
    (
        "MurNAc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::NAcetyl, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "MurNGc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::NGlycolyl,
            GlycanSubstituent::OCarboxyEthyl,
        ],
    ),
    (
        "Mur",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Amino, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "Api",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[GlycanSubstituent::HydroxyMethyl],
    ),
    (
        "Dha",
        BaseSugar::Heptose(None),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Bac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
        ],
    ),
    (
        "Pse",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "Leg",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "Aci",
        BaseSugar::Nonose,
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
];

/// Any substituent on a monosaccharide.
/// Source: <https://www.ncbi.nlm.nih.gov/glycans/snfg.html> table 3.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum GlycanSubstituent {
    ///Am N-acetimidoyl
    Acetimidoyl,
    ///Ac acetyl
    Acetyl,
    ///Ala2Ac N-acetyl-D-alanyl
    AcetylAlanyl,
    ///Gln2Ac N-acetyl-glutaminyl
    AcetylGlutaminyl,
    ///A acid
    Acid,
    ///Ala D-alanyl
    Alanyl,
    ///ol alcohol
    Alcohol,
    ///N amino
    Amino,
    ///aric ??
    Aric,
    ///Pyr 1-carboxyethylidene
    CargoxyEthylidene,
    ///d Deoxy
    Deoxy,
    ///3,4Hb 3,4-dihydroxybutyryl
    DiHydroxyButyryl,
    ///DiMe two methyl
    DiMethyl,
    ///AmMe2 N-(N,N-dimethyl-acetimidoyl)
    DiMethylAcetimidoyl,
    ///Gr2,3Me2 2,3-di-O-methyl-glyceryl
    DiMethylGlyceryl,
    ///en didehydro an addition of a double bond
    Didehydro,
    ///An element that replaces a side chain
    Element(Element),
    ///Etn Ethanolamine
    Ethanolamine,
    ///EtOH O linked ethanol
    EtOH,
    ///Fo formyl
    Formyl,
    ///Gr glyceryl
    Glyceryl,
    ///Gc glycolyl
    Glycolyl,
    ///Gly glycyl
    Glycyl,
    ///4Hb 4-hydroxybutyryl, 3RHb (R)-3-hydroxybutyryl, 3SHb (S)-3-hydroxybutyryl
    HydroxyButyryl,
    ///HydroxyMethyl
    HydroxyMethyl,
    ///Lac
    Lac,
    ///Lt lactyl
    Lactyl,
    ///Me methyl
    Methyl,
    ///AmMe N-(N-methyl-acetimidoyl)
    MethylAcetimidoyl,
    ///5Glu2Me N-methyl-5-glutamyl
    MethylGlutamyl,
    ///NAc N-acetyl
    NAcetyl,
    ///N2DiMe N linked double methyl
    NDiMe,
    ///NFo N linked formyl
    NFo,
    ///NGc N linked glycolyl
    NGlycolyl,
    ///carboxyethyl used in Mur
    OCarboxyEthyl,
    ///PCho phosphate linked choline
    PCholine,
    ///P phosphate
    Phosphate,
    ///Py pyruvyl
    Pyruvyl,
    ///Suc ??
    Suc,
    ///S sulfate
    Sulfate,
    ///Tau tauryl
    Tauryl,
    ///ulo ??
    Ulo,
    ///ulof ??
    Ulof,
    ///water H2O
    Water,
}

const POSTFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("Ac", GlycanSubstituent::Acetyl),
    ("Ala2Ac", GlycanSubstituent::AcetylAlanyl),
    ("Ala", GlycanSubstituent::Alanyl),
    ("AmMe2", GlycanSubstituent::DiMethylAcetimidoyl),
    ("AmMe", GlycanSubstituent::MethylAcetimidoyl),
    ("Am", GlycanSubstituent::Acetimidoyl),
    ("en", GlycanSubstituent::Didehydro),
    ("Fo", GlycanSubstituent::Formyl),
    ("Gc", GlycanSubstituent::Glycolyl),
    ("Gln2Ac", GlycanSubstituent::AcetylGlutaminyl),
    ("5Glu2Me", GlycanSubstituent::MethylGlutamyl), // TODO: Does this number has to be there? (Not found in the glycosmos data)
    ("Gly", GlycanSubstituent::Glycyl),
    ("Gr", GlycanSubstituent::Glyceryl),
    ("Gr2,3Me2", GlycanSubstituent::DiMethylGlyceryl),
    ("4Hb", GlycanSubstituent::HydroxyButyryl),
    ("3RHb", GlycanSubstituent::HydroxyButyryl),
    ("3SHb", GlycanSubstituent::HydroxyButyryl),
    ("3,4Hb", GlycanSubstituent::DiHydroxyButyryl),
    ("Lt", GlycanSubstituent::Lactyl),
    ("Lac", GlycanSubstituent::Lac),
    ("Me", GlycanSubstituent::Methyl),
    ("CMe", GlycanSubstituent::Methyl), // unsure about the difference with Me
    ("NAc", GlycanSubstituent::NAcetyl),
    ("Pyr", GlycanSubstituent::CargoxyEthylidene),
    ("Tau", GlycanSubstituent::Tauryl),
    ("onic", GlycanSubstituent::Acid),
    ("uronic", GlycanSubstituent::Acid),
    ("aric", GlycanSubstituent::Aric),
    ("ol", GlycanSubstituent::Alcohol),
    ("Etn", GlycanSubstituent::Ethanolamine),
    ("EtOH", GlycanSubstituent::Ethanolamine),
    ("ulof", GlycanSubstituent::Ulof),
    ("ulo", GlycanSubstituent::Ulo),
    ("N2DiMe", GlycanSubstituent::NDiMe),
    ("NDiMe", GlycanSubstituent::NDiMe),
    ("PCho", GlycanSubstituent::PCholine),
    ("CE", GlycanSubstituent::Glycyl), // Same molecular formula
    ("Suc", GlycanSubstituent::Suc),
    ("NFo", GlycanSubstituent::NFo),
    ("DiMe", GlycanSubstituent::DiMethyl),
    ("A", GlycanSubstituent::Acid),
    ("P", GlycanSubstituent::Phosphate),
    ("p", GlycanSubstituent::Phosphate),
    ("S", GlycanSubstituent::Sulfate),
    ("N", GlycanSubstituent::Amino),
];

const DOUBLE_LINKED_POSTFIX_SUBSTITUENTS: &[(&str, &[GlycanSubstituent])] = &[
    ("Py", &[GlycanSubstituent::Pyruvyl]),
    ("N", &[GlycanSubstituent::Water]),
    (
        "P",
        &[GlycanSubstituent::Water, GlycanSubstituent::Phosphate],
    ),
];

const PREFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("deoxy", GlycanSubstituent::Deoxy),
    ("Anhydro", GlycanSubstituent::Deoxy),
    ("d", GlycanSubstituent::Deoxy),
];

impl Display for GlycanSubstituent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Acetimidoyl => "Am".to_string(),
                Self::Acetyl => "Ac".to_string(),
                Self::AcetylAlanyl => "Ala2Ac".to_string(),
                Self::AcetylGlutaminyl => "Gln2Ac".to_string(),
                Self::Acid => "A".to_string(),
                Self::Alanyl => "Ala".to_string(),
                Self::Alcohol => "ol".to_string(),
                Self::Amino => "N".to_string(),
                Self::Aric => "aric".to_string(),
                Self::CargoxyEthylidene => "Pyr".to_string(),
                Self::Deoxy => "d".to_string(),
                Self::Didehydro => "en".to_string(),
                Self::DiHydroxyButyryl => "3,4Hb".to_string(),
                Self::DiMethyl => "DiMe".to_string(),
                Self::DiMethylAcetimidoyl => "AmMe2".to_string(),
                Self::DiMethylGlyceryl => "Gr2,3Me2".to_string(),
                Self::Ethanolamine => "Etn".to_string(),
                Self::Element(el) => el.to_string(),
                Self::EtOH => "EtOH".to_string(),
                Self::Formyl => "Fo".to_string(),
                Self::Glyceryl => "Gr".to_string(),
                Self::Glycolyl => "Gc".to_string(),
                Self::Glycyl => "Gly".to_string(),
                Self::HydroxyButyryl => "Hb".to_string(),
                Self::HydroxyMethyl => "HMe".to_string(),
                Self::Lac => "Lac".to_string(),
                Self::Lactyl => "Lt".to_string(),
                Self::Methyl => "Me".to_string(),
                Self::MethylAcetimidoyl => "AmMe".to_string(),
                Self::MethylGlutamyl => "5Glu2Me".to_string(),
                Self::NAcetyl => "NAc".to_string(),
                Self::NDiMe => "NDiMe".to_string(),
                Self::NFo => "NFo".to_string(),
                Self::NGlycolyl => "NGc".to_string(),
                Self::OCarboxyEthyl => "carboxyethyl".to_string(),
                Self::PCholine => "PCho".to_string(),
                Self::Phosphate => "P".to_string(),
                Self::Pyruvyl => "Py".to_string(),
                Self::Suc => "Suc".to_string(),
                Self::Sulfate => "S".to_string(),
                Self::Tauryl => "Tau".to_string(),
                Self::Ulo => "ulo".to_string(),
                Self::Ulof => "ulof".to_string(),
                Self::Water => "waterloss".to_string(),
            }
        )
    }
}

impl Chemical for GlycanSubstituent {
    fn composition(&self) -> ElementalComposition {
        let side = match self {
            Self::Acetimidoyl => molecular_formula!(H 5 C 2 N 1),
            Self::Acetyl => molecular_formula!(H 3 C 2 O 1),
            Self::AcetylAlanyl => molecular_formula!(H 8 C 5 N 1 O 2),
            Self::AcetylGlutaminyl => molecular_formula!(H 11 C 7 N 2 O 3),
            Self::Acid => molecular_formula!(H -1 O 2), // Together with the replacement below this is H-2 O+1
            Self::Alanyl => molecular_formula!(H 6 C 3 N 1 O 1),
            Self::Alcohol => molecular_formula!(H 3 O 1), // Together with the replacement below this is H+2
            Self::Amino => molecular_formula!(H 2 N 1),
            Self::Aric => molecular_formula!(H 3 O 3), // Together with replacement below this is H2O2
            Self::CargoxyEthylidene => molecular_formula!(H 3 C 3 O 3), // double substituent, calculated to work with the additional side chain deletion
            Self::Deoxy => molecular_formula!(H 1), // Together with the replacement below this is O-1
            Self::Didehydro => molecular_formula!(H -1 O 1), // Together with the replacement below this is H-2
            Self::DiHydroxyButyryl => molecular_formula!(H 7 C 4 O 3),
            Self::DiMethyl => molecular_formula!(H 5 C 2), // assumed to replace the both the OH and H on a single carbon
            Self::DiMethylAcetimidoyl => molecular_formula!(H 9 C 4 N 1),
            Self::DiMethylGlyceryl => molecular_formula!(H 9 C 5 O 3),
            Self::Ethanolamine => molecular_formula!(H 6 C 2 N 1 O 1),
            Self::EtOH => molecular_formula!(H 5 C 2 O 2),
            Self::Element(el) => ElementalComposition::from_tuples(&[(*el, 0, 1)]),
            Self::Formyl => molecular_formula!(H 1 C 1 O 1),
            Self::Glyceryl | Self::Lac => molecular_formula!(H 5 C 3 O 3),
            Self::Glycolyl => molecular_formula!(H 3 C 2 O 2),
            Self::Glycyl | Self::NAcetyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::HydroxyButyryl => molecular_formula!(H 7 C 4 O 2),
            Self::HydroxyMethyl | Self::Ulo => molecular_formula!(H 3 C 1 O 2), // Ulo: replaces H, together with replacement below this is H2C1O1
            Self::Lactyl => molecular_formula!(H 5 C 3 O 2),
            Self::Methyl => molecular_formula!(H 3 C 1),
            Self::MethylAcetimidoyl => molecular_formula!(H 7 C 3 N 1),
            Self::MethylGlutamyl => molecular_formula!(H 10 C 6 N 1 O 3),
            Self::NDiMe => molecular_formula!(H 6 C 2 N 1),
            Self::NFo => molecular_formula!(H 2 C 1 N 1 O 1),
            Self::NGlycolyl => molecular_formula!(H 4 C 2 N 1 O 2),
            Self::OCarboxyEthyl => molecular_formula!(H 6 C 3 O 3), // Replaces H, together with replacement below this is H5C3O2
            Self::PCholine => molecular_formula!(H 14 C 5 N 1 O 4 P 1),
            Self::Phosphate => molecular_formula!(H 2 O 4 P 1),
            Self::Pyruvyl => molecular_formula!(H 3 C 3 O 2),
            Self::Suc => molecular_formula!(H 6 C 4 N 1 O 3),
            Self::Sulfate => molecular_formula!(H 1 O 4 S 1),
            Self::Tauryl => molecular_formula!(H 6 C 2 N 1 O 3 S 1),
            Self::Ulof => molecular_formula!(H 4 C 1 O 2), // Replaces H, together with replacement below this is H3C1O1
            Self::Water => molecular_formula!(H - 1),
        };
        side - molecular_formula!(O 1 H 1) // substituent so replaces a standard oxygen side chain
    }
}

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
#[allow(dead_code)]
pub fn glycan_parse_list() -> &'static Vec<(String, MonoSaccharide)> {
    GLYCAN_PARSE_CELL.get_or_init(|| {
        vec![
            (
                "phosphate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("phosphate".to_string()),
                    furanose: false,
                },
            ),
            (
                "sulfate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("sulfate".to_string()),
                    furanose: false,
                },
            ),
            (
                "Sug".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Sugar,
                    substituents: vec![],
                    proforma_name: Some("Sug".to_string()),
                    furanose: false,
                },
            ),
            (
                "Tri".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Triose,
                    substituents: vec![],
                    proforma_name: Some("Tri".to_string()),
                    furanose: false,
                },
            ),
            (
                "Tet".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Tetrose(None),
                    substituents: vec![],
                    proforma_name: Some("Tet".to_string()),
                    furanose: false,
                },
            ),
            (
                "Pen".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                },
            ),
            (
                "a-Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Acid],
                    proforma_name: Some("a-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "en,a-Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Didehydro,
                        GlycanSubstituent::Deoxy,
                    ],
                    proforma_name: Some("en,a-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "d-Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("d-Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNAc(S)".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNAc(S)".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNAc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexNS".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNS".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexN".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino],
                    proforma_name: Some("HexN".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexS".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexS".to_string()),
                    furanose: false,
                },
            ),
            (
                "HexP".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("HexP".to_string()),
                    furanose: false,
                },
            ),
            (
                "Hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "Hep".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Heptose(None),
                    substituents: vec![],
                    proforma_name: Some("Hep".to_string()),
                    furanose: false,
                },
            ),
            (
                "Oct".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![],
                    proforma_name: Some("Oct".to_string()),
                    furanose: false,
                },
            ),
            (
                "Non".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                },
            ),
            (
                "Dec".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Decose,
                    substituents: vec![],
                    proforma_name: Some("Dec".to_string()),
                    furanose: false,
                },
            ),
            (
                "Neu5Ac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "Neu5Gc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Gc".to_string()),
                    furanose: false,
                },
            ),
            (
                "Neu".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu".to_string()),
                    furanose: false,
                },
            ),
            (
                "Fuc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
            (
                "P".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("Hexphosphate".to_string()), // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements,  both getting counts after them
                    furanose: false,
                },
            ),
            (
                "H".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                },
            ),
            (
                "N".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                },
            ),
            (
                "F".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                },
            ),
            (
                "S".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "A".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Ac".to_string()),
                    furanose: false,
                },
            ),
            (
                "G".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose,
                    substituents: vec![
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                    proforma_name: Some("Neu5Gc".to_string()),
                    furanose: false,
                },
            ),
        ]
    })
}
#[allow(dead_code)]
static GLYCAN_PARSE_CELL: OnceLock<Vec<(String, MonoSaccharide)>> = OnceLock::new();

/// Rose tree representation of glycan structure
#[allow(dead_code)]
#[derive(Eq, PartialEq, Clone, Hash, Debug, Serialize, Deserialize)]
pub struct GlycanStructure {
    pub(super) sugar: MonoSaccharide,
    pub(super) branches: Vec<GlycanStructure>,
}

impl GlycanStructure {
    /// Parse a short IUPAC glycan structure
    /// # Panics
    /// Panics if there is no single sugar found
    /// # Errors
    /// Errors when the format is not correct, could be unknown monosaccharide, or a open brace
    pub fn from_short_iupac(
        line: &str,
        range: Range<usize>,
        line_number: usize,
    ) -> Result<Self, CustomError> {
        // GlcNAc(?1-?)Man(?1-?)[Man(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc-ol
        let mut offset = range.start;
        let mut branch = Self {
            sugar: MonoSaccharide::new(BaseSugar::Decose, &[]),
            branches: Vec::new(),
        }; // Starting sugar, will be removed
        let mut last_branch: &mut Self = &mut branch;
        let bytes = line.as_bytes();
        while offset < range.end {
            while bytes[offset] == b'[' {
                let end = end_of_enclosure(bytes, offset + 1, b'[', b']').ok_or_else(|| {
                    CustomError::error(
                        "Invalid iupac short glycan",
                        "No closing brace for branch",
                        Context::line(line_number, line, offset, range.end - offset),
                    )
                })?;
                last_branch.branches.push(Self::from_short_iupac(
                    line,
                    offset + 1..end,
                    line_number,
                )?);

                offset = end + 1;
            }
            let (sugar, new_offset) = MonoSaccharide::from_short_iupac(line, offset, line_number)?;
            offset = new_offset;

            last_branch.branches.push(Self {
                sugar: sugar.clone(),
                branches: Vec::new(),
            });
            last_branch = last_branch.branches.last_mut().unwrap();
            if bytes[offset] == b'(' {
                if let Some(end) = index_of_char(bytes, offset + 1, b')') {
                    offset = end + 1; // just ignore all linking stuff I do not care
                } else {
                    // This only happens for incomplete branches where the last parts of the branch are unknown.
                    assert!(range.end - offset < 10); // make sure it is the last part
                    offset = range.end; // assume it is the last not closed brace
                }
            }
        }
        branch.branches.pop().map_or_else(
            || {
                Err(CustomError::error(
                    "Invalid iupac short glycan",
                    "No glycan found",
                    Context::line(line_number, line.to_string(), range.start, range.len()),
                ))
            },
            Ok,
        )
    }
}

impl Chemical for GlycanStructure {
    fn composition(&self) -> ElementalComposition {
        self.sugar.composition()
            + self
                .branches
                .iter()
                .map(Chemical::composition)
                .sum::<ElementalComposition>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pro_forma_compliance() {
        let cases = &[
            ("Hep", molecular_formula!(H 12 C 7 O 6)),
            ("phosphate", molecular_formula!(H 1 O 3 P 1)),
            ("a-Hex", molecular_formula!(H 8 C 6 O 6)),
            ("Sug", molecular_formula!(H 2 C 2 O 1)),
            ("HexN", molecular_formula!(H 11 C 6 N 1 O 4)),
            ("Pen", molecular_formula!(H 8 C 5 O 4)),
            ("Tet", molecular_formula!(H 6 C 4 O 3)),
            ("HexP", molecular_formula!(H 11 C 6 O 8 P 1)),
            ("Neu5Ac", molecular_formula!(H 17 C 11 N 1 O 8)),
            ("Non", molecular_formula!(H 16 C 9 O 8)),
            ("HexNAc(S)", molecular_formula!(H 13 C 8 N 1 O 8 S 1)),
            ("Dec", molecular_formula!(H 18 C 10 O 9)),
            ("en,a-Hex", molecular_formula!(H 6 C 6 O 5)),
            ("Neu5Gc", molecular_formula!(H 17 C 11 N 1 O 9)),
            ("Neu", molecular_formula!(H 15 C 9 N 1 O 7)),
            ("HexNAc", molecular_formula!(H 13 C 8 N 1 O 5)),
            ("Fuc", molecular_formula!(H 10 C 6 O 4)),
            ("HexNS", molecular_formula!(H 11 C 6 N 1 O 7 S 1)),
            ("Tri", molecular_formula!(H 4 C 3 O 2)),
            ("Oct", molecular_formula!(H 14 C 8 O 7)),
            ("sulfate", molecular_formula!(O 3 S 1)),
            ("d-Hex", molecular_formula!(H 10 C 6 O 5)),
            ("Hex", molecular_formula!(H 10 C 6 O 5)),
            ("HexS", molecular_formula!(H 10 C 6 O 8 S 1)),
        ];
        for (name, composition) in cases {
            assert_eq!(
                glycan_parse_list()
                    .iter()
                    .find(|p| p.0 == *name)
                    .unwrap_or_else(|| panic!("Assumed {name} would be defined"))
                    .1
                    .composition(),
                *composition,
                "{name}",
            );
        }
    }

    #[test]
    fn iupac_short_names() {
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal2,3Ac24-1,6-1Py", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Pyruvyl,
                    ]
                ),
                18
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcNAc", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[GlycanSubstituent::NAcetyl]
                ),
                6
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal6S", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[GlycanSubstituent::Sulfate]
                ),
                5
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcN2Gc", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[GlycanSubstituent::Amino, GlycanSubstituent::Glycolyl,]
                ),
                7
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GalNAc3S", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate]
                ),
                8
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("GlcN2,6S2", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Sulfate,
                        GlycanSubstituent::Sulfate
                    ]
                ),
                9
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Tagf1,6P2", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)),
                    &[GlycanSubstituent::Phosphate, GlycanSubstituent::Phosphate]
                )
                .furanose(),
                9
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal2,3Ac24-1,6-1Py", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Pyruvyl,
                    ]
                ),
                18
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("D-Araf", 0, 0),
            Ok((
                MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)), &[])
                    .furanose(),
                6
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Xyl-onic", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Pentose(Some(PentoseIsomer::Xylose)),
                    &[GlycanSubstituent::Acid]
                ),
                8
            ))
        );
        assert_eq!(
            MonoSaccharide::from_short_iupac("Glc2,3,4,6Ac4", 0, 0),
            Ok((
                MonoSaccharide::new(
                    BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    &[
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acetyl
                    ]
                ),
                13
            ))
        );
    }

    /*
    #[allow(clippy::float_cmp)] // Handled in a different way
    fn iupac_masses() {
        assert_eq!(
            MonoSaccharide::from_short_iupac("Gal3DiMe(b1-4)GlcNAc(b1-", 0, 0)
                .unwrap()
                .0
                .composition()
                .monoisotopic_mass()
                .unwrap()
                .value
                .round(),
            411.0
        );
    }*/
}

// TODO: put in common utils
/// Get the index of the next copy of the given char
fn index_of_char(chars: &[u8], start: usize, char: u8) -> Option<usize> {
    for (i, ch) in chars[start..].iter().enumerate() {
        if *ch == char {
            return Some(start + i);
        }
    }
    None
}

// TODO: put in common utils
/// Find the enclosed text by the given symbols, assumes a single open is already read just before the start
fn end_of_enclosure(chars: &[u8], start: usize, open: u8, close: u8) -> Option<usize> {
    let mut state = 1;
    for (i, ch) in chars[start..].iter().enumerate() {
        if *ch == open {
            state += 1;
        } else if *ch == close {
            state -= 1;
            if state == 0 {
                return Some(start + i);
            }
        }
    }
    None
}