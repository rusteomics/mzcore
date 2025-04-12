use crate::{
    error::{Context, CustomError},
    formula::Chemical,
    glycan::MonoSaccharide,
    helper_functions::{explain_number_error, RangeExtension},
    Element, MolecularFormula,
};
use std::{
    num::NonZeroU16,
    ops::{Range, RangeBounds},
};

enum Brick {
    Element(Element),
    Formula(MolecularFormula),
}

/// # Errors
/// Errors if the provided text is not a Unimod composition brick.
fn parse_unimod_composition_brick(text: &str, range: Range<usize>) -> Result<Brick, CustomError> {
    match text[range.clone()].to_lowercase().as_str() {
        "ac" => Ok(Brick::Formula(molecular_formula!(C 2 H 2 O 1))),
        "me" => Ok(Brick::Formula(molecular_formula!(C 1 H 2))),
        "kdn" => Ok(Brick::Formula(molecular_formula!(C 9 H 14 O 8))),
        "kdo" => Ok(Brick::Formula(molecular_formula!(C 8 H 12 O 7))),
        "sulf" => Ok(Brick::Formula(molecular_formula!(S 1 O 3))),
        "phos" => Ok(Brick::Formula(molecular_formula!(P 1 O 3))),
        "water" => Ok(Brick::Formula(molecular_formula!(H 2 O 1))),
        _ => {
            Element::try_from(text[range.clone()].to_lowercase().as_str()).map_or_else(|()| if let Ok((ms, _)) =
                MonoSaccharide::from_short_iupac(text, range.start_index(), range.len())
            {
                Ok(Brick::Formula(ms.formula_inner(crate::SequencePosition::default(),0)))
            } else {
                Err(CustomError::error(
                    "Invalid Unimod chemical formula", 
                    "Unknown Unimod composition brick, use an element or one of the unimod shorthands. Eg: 'H(13) C(12) N O(3)'.",
                     Context::line_range(None, text, range)))
            }, |el| Ok(Brick::Element(el)))
        }
    }
}

impl MolecularFormula {
    /// Parses Unimod compositions into molecular formulas. As Unimod compositions can have glycans in them these are reported as molecular formula.
    /// ```text
    /// H(25) C(8) 13C(7) N 15N(2) O(3)
    /// H(6) C(4) N(2) dHex
    /// ```
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It panics if the string contains not UTF8 symbols.
    pub fn from_unimod(value: &str, range: impl RangeBounds<usize>) -> Result<Self, CustomError> {
        let (mut index, end) = range.bounds(value.len());
        assert!(value.is_ascii());

        let mut formula = Self::default();

        let mut isotope = None;
        let mut last_name_index = -1_isize;
        let mut last_name = String::new();
        while index < end {
            match value.as_bytes()[index] {
            b'(' => {
                let length = value.chars().skip(index+1).take_while(|c| *c == '-' || *c == '+' || c.is_ascii_digit()).count();
                let num = value[index+1..index+1+length].parse::<i32>()
                    .map_err(|err|
                        CustomError::error(
                            "Invalid Unimod chemical formula", 
                            format!("The element amount {}", explain_number_error(&err)),
                            Context::line(None, value, index+1, length)))?;
                match parse_unimod_composition_brick(value, last_name_index as usize..last_name_index as usize+last_name.len())? {
                    Brick::Element(el) => {if !formula.add((el, isotope.take(), num)) {
                        return Err(CustomError::error("Invalid Unimod chemical formula", "An element or isotope without a defined mass was found", Context::line_range(None, value,last_name_index as usize..last_name_index as usize+last_name.len())));
                    }},
                    Brick::Formula(f) => formula += f*num,
                }
                last_name.clear();
                last_name_index = -1;
                index += length + 2;
                if value.as_bytes()[index-1] != b')' {
                    return Err(CustomError::error("Invalid Unimod chemical formula", "The amount of an element should be closed by ')'", Context::line(None, value, index-1, 1)));
                }
            }
            b' ' => {
                if !last_name.is_empty() {
                    match parse_unimod_composition_brick(value, last_name_index as usize..last_name_index as usize+last_name.len())? {
                        Brick::Element(el) => {if !formula.add((el, isotope.take(), 1)) {
                            return Err(CustomError::error("Invalid Unimod chemical formula", "An element or isotope without a defined mass was found", Context::line_range(None, value,last_name_index as usize..last_name_index as usize+last_name.len())));
                        }},
                        Brick::Formula(f) => formula += f,
                        }
                    last_name.clear();
                    last_name_index = -1;
                }
                index += 1;
            }
            n if n.is_ascii_digit() => {
                let length = value.chars().skip(index).take_while(char::is_ascii_digit).count();
                isotope = Some(value[index..index+length].parse::<NonZeroU16>()
                    .map_err(|err|
                        CustomError::error(
                            "Invalid Unimod chemical formula", 
                            format!("The isotope {}", explain_number_error(&err)),
                            Context::line(None, value, index, length)))?);
                index += length;
            },
            n if n.is_ascii_alphabetic() => {
                last_name.push(n as char);
                if last_name_index == -1 {
                    last_name_index = isize::try_from(index).unwrap();
                }
                index += 1;
            }
            _ => return Err(CustomError::error(
            "Invalid Unimod chemical formula", 
            "Unexpected character, use an element or one of the unimod shorthands. Eg: 'H(13) C(12) N O(3)'.",
                        Context::line(None, value, index, 1))),
        }
        }
        if !last_name.is_empty() {
            match parse_unimod_composition_brick(
                value,
                last_name_index as usize..last_name_index as usize + last_name.len(),
            )? {
                Brick::Element(el) => {
                    if !formula.add((el, isotope.take(), 1)) {
                        return Err(CustomError::error(
                            "Invalid Unimod chemical formula",
                            "An element or isotope without a defined mass was found",
                            Context::line_range(
                                None,
                                value,
                                last_name_index as usize
                                    ..last_name_index as usize + last_name.len(),
                            ),
                        ));
                    }
                }
                Brick::Formula(f) => formula += f,
            }
        }
        Ok(formula)
    }
}
