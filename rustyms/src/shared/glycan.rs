use std::{fmt::Display, hash::Hash, sync::OnceLock};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    formula::{Chemical, MolecularFormula},
    Element, SequencePosition, ELEMENT_PARSE_LIST,
};

/// Glycan absolute configuration
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum Configuration {
    /// D configuration
    D,
    /// L configuration
    L,
    /// Double configuration D and D
    DD,
    /// Double configuration L and L
    LL,
    /// Double configuration D and L
    DL,
    /// Double configuration L and D
    LD,
}

/// A monosaccharide with all its complexity
#[derive(Clone, Ord, PartialOrd, Debug, Serialize, Deserialize)]
pub struct MonoSaccharide {
    pub(super) base_sugar: BaseSugar,
    pub(super) substituents: Vec<GlycanSubstituent>,
    pub(super) furanose: bool,
    pub(super) configuration: Option<Configuration>,
    pub(super) proforma_name: Option<String>,
}

impl MonoSaccharide {
    fn equivalent(&self, other: &Self, precise: bool) -> bool {
        self.base_sugar.equivalent(&other.base_sugar, precise)
            && self.substituents == other.substituents
            && (!precise
                || (self.furanose == other.furanose && self.configuration == other.configuration))
    }
}

impl std::cmp::PartialEq for MonoSaccharide {
    fn eq(&self, other: &Self) -> bool {
        self.equivalent(other, true)
    }
}

impl std::cmp::Eq for MonoSaccharide {}

impl std::hash::Hash for MonoSaccharide {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        self.base_sugar.hash(hasher);
        self.substituents.hash(hasher);
        self.furanose.hash(hasher);
        self.configuration.hash(hasher);
    }
}

impl MonoSaccharide {
    /// Check if this is a fucose
    pub fn is_fucose(&self) -> bool {
        self.base_sugar == BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
            && self.substituents.contains(&GlycanSubstituent::Deoxy)
    }

    /// Create a new monosaccharide
    pub fn new(sugar: BaseSugar, substituents: &[GlycanSubstituent]) -> Self {
        Self {
            base_sugar: sugar,
            substituents: substituents.to_owned(),
            furanose: false,
            configuration: None,
            proforma_name: None,
        }
    }

    /// Get this same monosaccharide but now with the given ProForma name
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

    /// Set this saccharide up to be a certain configuration
    #[must_use]
    #[allow(dead_code)]
    pub fn configuration(self, configuration: Configuration) -> Self {
        Self {
            configuration: Some(configuration),
            ..self
        }
    }

    /// Simplify a glycan composition to be sorted and deduplicated.
    /// Returns None if overflow occurred, meaning that there where more than `isize::MAX` or less then `isize::MIN` monosaccharides for one species.
    pub(crate) fn simplify_composition(
        mut composition: Vec<(Self, isize)>,
    ) -> Option<Vec<(Self, isize)>> {
        // Sort on monosaccharide
        composition.retain(|el| el.1 != 0);
        composition.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate
        let mut max = composition.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &composition[index];
            let next = &composition[index + 1];
            if this.0 == next.0 {
                composition[index].1 = composition[index].1.checked_add(next.1)?;
                composition.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        composition.retain(|el| el.1 != 0);
        Some(composition)
    }

    /// Parse the given text (will be changed to lowercase) as a glycan composition.
    /// # Errors
    /// When the composition could not be read. Or when any of the glycans occurs outside of the valid range
    pub fn from_composition(text: &str) -> Result<Vec<(Self, isize)>, CustomError> {
        let basic_error =
            CustomError::error("Invalid glycan composition", "..", Context::show(text));
        Self::simplify_composition(
            crate::helper_functions::parse_named_counter(
                &text.to_ascii_lowercase(),
                glycan_parse_list(),
                false,
            )
            .map_err(|e| {
                basic_error.with_long_description(format!(
                    "This modification cannot be read as a valid glycan: {e}"
                ))
            })?,
        )
        .ok_or_else(|| {
            basic_error.with_long_description(format!(
                "The occurrence of one monosaccharide species is outside of the range {} to {}",
                isize::MIN,
                isize::MAX
            ))
        })
    }

    /// Parse a short IUPAC name from this string starting at `start` and returning,
    /// if successful, a monosaccharide and the offset in the string where parsing ended.
    /// # Errors
    /// Fails if it finds a structure that does not fit the IUPAC glycan name.
    pub fn from_short_iupac(
        original_line: &str,
        start: usize,
        line_index: usize,
    ) -> Result<(Self, usize), CustomError> {
        let mut index = start;
        let line = original_line.to_ascii_lowercase();
        let bytes = line.as_bytes();
        let mut substituents = Vec::new();
        let mut configuration = None;
        let mut epi = false;

        // ignore stuff
        index += line[index..].ignore(&["keto-"]);
        if line[index..].starts_with("d-") {
            configuration = Some(Configuration::D);
            index += 2;
        } else if line[index..].starts_with("l-") {
            configuration = Some(Configuration::L);
            index += 2;
        } else if line[index..].starts_with("?-") {
            configuration = None;
            index += 2;
        }
        // Prefix mods
        let mut amount = 1;
        if bytes[index].is_ascii_digit() {
            match bytes.get(index + 1) {
                Some(b',') if bytes.get(index + 3).copied() == Some(b':') => {
                    let start_index = index;
                    index += 7;
                    index += line[index..].ignore(&["-"]);
                    if !line[index..].starts_with("anhydro") {
                        return Err(CustomError::error(
                            "Invalid iupac monosaccharide name",
                            "This internally linked glycan could not be parsed, expected Anhydro as modification",
                            Context::Line {
                                line_index: Some(line_index),
                                line: original_line.to_string(),
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
                Some(b',') => {
                    let num = bytes[index + 1..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b',' || **c == b'?')
                        .count();
                    index += num + 1;
                    amount = num / 2;
                    // X,X{mod} (or 3/4/5/etc mods)
                }
                Some(_) => index += 1, // X{mod}
                None => (),
            }
            index += line[index..].ignore(&["-"]);
        }
        // Detect epi state
        if line[index..].starts_with('e') {
            epi = true;
            index += 1;
        }
        // Get the prefix mods
        if !line[index..].starts_with("dig") && !line[index..].starts_with("dha") {
            if let Some(o) = line[index..].take_any(PREFIX_SUBSTITUENTS, |e| {
                substituents.extend(std::iter::repeat(*e).take(amount));
            }) {
                index += o;
            }
            index += line[index..].ignore(&["-"]);
        }
        // Another optional isomeric state
        if line[index..].starts_with("d-") {
            configuration = Some(Configuration::D);
            index += 2;
        } else if line[index..].starts_with("l-") {
            configuration = Some(Configuration::L);
            index += 2;
        } else if line[index..].starts_with("?-") {
            configuration = None;
            index += 2;
        }
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
                    base_sugar: match b {
                        BaseSugar::Nonose(Some(NonoseIsomer::Leg)) if epi => {
                            BaseSugar::Nonose(Some(NonoseIsomer::ELeg))
                        }
                        other => other,
                    },
                    substituents,
                    furanose: false,
                    configuration,
                    proforma_name: None,
                };
                alo.substituents.extend(s.iter().copied());
                alo
            })
            .ok_or_else(|| {
                CustomError::error(
                    "Invalid iupac monosaccharide name",
                    "This name could not be recognised as a standard iupac glycan name",
                    Context::Line {
                        line_index: Some(line_index),
                        line: original_line.to_string(),
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
            index += line[index..].ignore(&["(x)", "(r)", "(s)"]);
            if double {
                if let Some(o) = line[index..].take_any(DOUBLE_LINKED_POSTFIX_SUBSTITUENTS, |e| {
                    sugar.substituents.extend(
                        e.iter()
                            .flat_map(|s| std::iter::repeat(s).take(double_amount))
                            .copied(),
                    );
                    if single_amount > 0 {
                        sugar.substituents.extend(
                            e.iter()
                                .filter(|s| **s != GlycanSubstituent::Water)
                                .flat_map(|s| std::iter::repeat(s).take(single_amount))
                                .copied(),
                        );
                    }
                }) {
                    index += o;
                } else {
                    return Err(CustomError::error(
                        "Invalid iupac monosaccharide name",
                        "No detected double linked glycan substituent was found, while the pattern for location is for a double linked substituent",
                        Context::Line {
                            line_index: Some(line_index),
                            line: original_line.to_string(),
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
                        .extend(std::iter::repeat(*e).take(amount));
                }) {
                    index += o;
                } else if let Some(o) = line[index..].take_any(ELEMENT_PARSE_LIST, |e| {
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
        index += line[index..].ignore(&["?"]); // I guess to indicate partial structures
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
        let possibly_unknown_number = |n: u8| n.is_ascii_digit() || n == b'?';
        let number_or_slash = |n: &u8| n.is_ascii_digit() || *n == b'/';
        let possibly_unknown_number_or_comma = |n: &u8| possibly_unknown_number(*n) || *n == b',';

        if possibly_unknown_number(bytes[0]) && bytes.len() > 1 {
            match bytes[1] {
                b',' => {
                    let num = bytes[1..]
                        .iter()
                        .copied()
                        .take_while(possibly_unknown_number_or_comma)
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
                        .copied()
                        .take_while(number_or_slash)
                        .count();
                    index += num + 2;
                    // X/X/X...{mod} multiple possible locations
                }
                c if possibly_unknown_number(c) && bytes[0] == b'?' => {
                    if bytes[2] == b',' {
                        let num = bytes[2..]
                            .iter()
                            .copied()
                            .take_while(possibly_unknown_number_or_comma)
                            .count();
                        index += num + 2;
                        amount = num / 2 + 1;
                        // ?X,X{mod} (or 3/4/5/etc mods)
                    } else if bytes[2] == b'/' {
                        let num = bytes[3..]
                            .iter()
                            .copied()
                            .take_while(number_or_slash)
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
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.base_sugar
            .formula_inner(sequence_index, peptidoform_index)
            + self
                .substituents
                .as_slice()
                .formula_inner(sequence_index, peptidoform_index)
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
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
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
    Nonose(Option<NonoseIsomer>),
    /// 10 carbon base sugar
    Decose,
}

impl BaseSugar {
    fn equivalent(&self, other: &Self, precise: bool) -> bool {
        match (self, other) {
            (Self::None, Self::None)
            | (Self::Sugar, Self::Sugar)
            | (Self::Octose, Self::Octose)
            | (Self::Decose, Self::Decose)
            | (Self::Triose, Self::Triose) => true,
            (Self::Tetrose(a), Self::Tetrose(b)) => !precise || a == b,
            (Self::Pentose(a), Self::Pentose(b)) => !precise || a == b,
            (Self::Hexose(a), Self::Hexose(b)) => !precise || a == b,
            (Self::Heptose(a), Self::Heptose(b)) => !precise || a == b,
            (Self::Nonose(a), Self::Nonose(b)) => !precise || a == b,
            _ => false,
        }
    }
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
                Self::Nonose(_) => "Non",
                Self::Decose => "Dec",
            }
        )
    }
}

impl Chemical for BaseSugar {
    fn formula_inner(
        &self,
        _sequence_index: SequencePosition,
        _peptidoform_index: usize,
    ) -> MolecularFormula {
        match self {
            Self::None => MolecularFormula::default(),
            Self::Sugar => molecular_formula!(H 2 C 2 O 1),
            Self::Triose => molecular_formula!(H 4 C 3 O 2),
            Self::Tetrose(_) => molecular_formula!(H 6 C 4 O 3),
            Self::Pentose(_) => molecular_formula!(H 8 C 5 O 4),
            Self::Hexose(_) => molecular_formula!(H 10 C 6 O 5),
            Self::Heptose(_) => molecular_formula!(H 12 C 7 O 6),
            Self::Octose => molecular_formula!(H 14 C 8 O 7),
            Self::Nonose(_) => molecular_formula!(H 16 C 9 O 8),
            Self::Decose => molecular_formula!(H 18 C 10 O 9),
        }
    }
}

/// Any 4 carbon glycan
#[allow(dead_code)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum TetroseIsomer {
    /// Ery
    Erythrose,
    /// Tho
    Threose,
}

/// Any 5 carbon glycan
#[allow(dead_code)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
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
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
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
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum HeptoseIsomer {
    /// gro-manHep
    GlyceroMannoHeptopyranose, // TODO: Does this indicate some mods?
    /// Sed
    Sedoheptulose,
}

/// Any 9 carbon glycan, these isomers are modification specific (need the correct substituents
/// applied to be meaningful). These are to be used only to store isomeric state that was inferred
/// from other sources that cannot be tracked in other ways in the current structure. Any isomer
/// used that does not have the correct monosaccharide substituents applied is meaningless.
#[allow(dead_code)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum NonoseIsomer {
    /// 3-Deoxy-D-glycero-D-galacto-non-2-ulopyranosonic acid
    Kdn,
    /// 5,7-Diamino-3,5,7,9-tetradeoxy-L-glycero-L-manno-non-2-ulopyranosonic acid
    Pse,
    /// 5,7-Diamino-3,5,7,9-tetradeoxy-D-glycero-D-galacto-non-2-ulopyranosonic acid
    Leg,
    /// 4 or 8 eLeg
    ELeg,
    /// 5,7-Diamino-3,5,7,9-tetradeoxy-L-glycero-L-altro-non-2-ulopyranosonic acid
    Aci,
}

/// Any substituent on a monosaccharide.
/// Source: <https://www.ncbi.nlm.nih.gov/glycans/snfg.html> table 3.
#[allow(dead_code)]
#[derive(Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum GlycanSubstituent {
    ///`Am` N-acetimidoyl
    Acetimidoyl,
    ///`Ac` acetyl
    Acetyl,
    ///`A` acid
    Acid,
    ///`Ala` D-alanyl
    Alanyl,
    ///`ol` alcohol
    Alcohol,
    ///`N` amino
    Amino,
    ///`aric` ??
    Aric,
    ///`Pyr` 1-carboxyethylidene
    CargoxyEthylidene,
    ///`d` Deoxy
    Deoxy,
    ///`DiMe` two methyl
    DiMethyl,
    ///`en` didehydro an addition of a double bond
    Didehydro,
    ///`An` element that replaces a side chain
    Element(Element),
    ///`Etn` Ethanolamine
    Ethanolamine,
    ///`EtOH` O linked ethanol
    EtOH,
    ///`Fo` formyl
    Formyl,
    ///`Gr` glyceryl
    Glyceryl,
    ///`Gc` glycolyl
    Glycolyl,
    ///`Gly` glycyl
    Glycyl,
    ///`4Hb` 4-hydroxybutyryl, 3RHb (R)-3-hydroxybutyryl, 3SHb (S)-3-hydroxybutyryl
    HydroxyButyryl,
    ///`HydroxyMethyl`
    HydroxyMethyl,
    ///`Lac`
    Lac,
    ///`Lt` lactyl
    Lactyl,
    ///`Me` methyl
    Methyl,
    ///`NAc` N-acetyl
    NAcetyl,
    ///`N2DiMe` N linked double methyl
    NDiMe,
    ///`NFo` N linked formyl
    NFo,
    ///`NGc` N linked glycolyl
    NGlycolyl,
    ///`carboxyethyl` used in Mur
    OCarboxyEthyl,
    ///`PCho` phosphate linked choline
    PCholine,
    ///`P` phosphate
    Phosphate,
    ///`Py` pyruvyl
    Pyruvyl,
    ///`Suc` ??
    Suc,
    ///`S` sulfate
    Sulfate,
    ///`Tau` tauryl
    Tauryl,
    ///`ulo` ??
    Ulo,
    ///`ulof` ??
    Ulof,
    ///`water` H2O
    Water,
}

impl GlycanSubstituent {
    /// Get the symbol used to denote this substituent
    pub const fn notation(self) -> &'static str {
        match self {
            Self::Acetimidoyl => "Am",
            Self::Acetyl => "Ac",
            Self::Acid => "A",
            Self::Alanyl => "Ala",
            Self::Alcohol => "ol",
            Self::Amino => "N",
            Self::Aric => "aric",
            Self::CargoxyEthylidene => "Pyr",
            Self::Deoxy => "d",
            Self::Didehydro => "en",
            Self::DiMethyl => "Me2",
            Self::Ethanolamine => "Etn",
            Self::Element(el) => el.symbol(),
            Self::EtOH => "EtOH",
            Self::Formyl => "Fo",
            Self::Glyceryl => "Gr",
            Self::Glycolyl => "Gc",
            Self::Glycyl => "Gly",
            Self::HydroxyButyryl => "Hb",
            Self::HydroxyMethyl => "HMe",
            Self::Lac => "Lac",
            Self::Lactyl => "Lt",
            Self::Methyl => "Me",
            Self::NAcetyl => "NAc",
            Self::NDiMe => "NDiMe",
            Self::NFo => "NFo",
            Self::NGlycolyl => "NGc",
            Self::OCarboxyEthyl => "carboxyethyl",
            Self::PCholine => "PCho",
            Self::Phosphate => "P",
            Self::Pyruvyl => "Py",
            Self::Suc => "Suc",
            Self::Sulfate => "S",
            Self::Tauryl => "Tau",
            Self::Ulo => "ulo",
            Self::Ulof => "ulof",
            Self::Water => "water_loss",
        }
    }
}

impl Display for GlycanSubstituent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.notation())
    }
}

impl Chemical for GlycanSubstituent {
    fn formula_inner(
        &self,
        _sequence_index: SequencePosition,
        _peptidoform_index: usize,
    ) -> MolecularFormula {
        let side = match self {
            Self::Acetimidoyl => molecular_formula!(H 5 C 2 N 1),
            Self::Acetyl => molecular_formula!(H 3 C 2 O 1),
            Self::Acid => molecular_formula!(H -1 O 2), // Together with the replacement below this is H-2 O+1
            Self::Alanyl => molecular_formula!(H 6 C 3 N 1 O 1),
            Self::Alcohol => molecular_formula!(H 3 O 1), // Together with the replacement below this is H+2
            Self::Amino => molecular_formula!(H 2 N 1),
            Self::Aric => molecular_formula!(H 3 O 3), // Together with replacement below this is H2O2
            Self::CargoxyEthylidene => molecular_formula!(H 3 C 3 O 3), // double substituent, calculated to work with the additional side chain deletion
            Self::Deoxy => molecular_formula!(H 1), // Together with the replacement below this is O-1
            Self::Didehydro => molecular_formula!(H -1 O 1), // Together with the replacement below this is H-2
            Self::DiMethyl => molecular_formula!(H 5 C 2), // assumed to replace the both the OH and H on a single carbon
            Self::Ethanolamine => molecular_formula!(H 6 C 2 N 1 O 1),
            Self::EtOH => molecular_formula!(H 5 C 2 O 2),
            Self::Element(el) => MolecularFormula::new(&[(*el, None, 1)], &[]).unwrap(),
            Self::Formyl => molecular_formula!(H 1 C 1 O 1),
            Self::Glyceryl | Self::Lac => molecular_formula!(H 5 C 3 O 3),
            Self::Glycolyl => molecular_formula!(H 3 C 2 O 2),
            Self::Glycyl | Self::NAcetyl => molecular_formula!(H 4 C 2 N 1 O 1),
            Self::HydroxyButyryl => molecular_formula!(H 7 C 4 O 2),
            Self::HydroxyMethyl | Self::Ulo => molecular_formula!(H 3 C 1 O 2), // Ulo: replaces H, together with replacement below this is H2C1O1
            Self::Lactyl => molecular_formula!(H 5 C 3 O 2),
            Self::Methyl => molecular_formula!(H 3 C 1),
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
