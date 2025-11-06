use std::{num::NonZeroU16, ops::Range};

use context_error::*;

use crate::{
    chemistry::{Chemical, Element, MolecularFormula},
    glycan::MonoSaccharide,
    helper_functions::{RangeExtension, explain_number_error},
    sequence::SequencePosition,
};

enum Brick {
    Element(Element),
    Formula(MolecularFormula),
}

/// # Errors
/// Errors if the provided text is not a Unimod composition brick.
fn parse_unimod_composition_brick(
    text: &str,
    range: Range<usize>,
) -> Result<Brick, BoxedError<'_, BasicKind>> {
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
                MonoSaccharide::from_short_iupac(text, range.start_index(), 0)
            {
                Ok(Brick::Formula(ms.formula_inner(SequencePosition::default(),0)))
            } else {
                Err(BoxedError::new(BasicKind::Error,
                    "Invalid Unimod chemical formula",
                    "Unknown Unimod composition brick, use an element or one of the unimod shorthands. Eg: 'H(13) C(12) N O(3)'.",
                     Context::line_range(None, text, range)))
            }, |el| Ok(Brick::Element(el)))
        }
    }
}

impl MolecularFormula {
    /// Parse a molecular formula from a Unimod formula.
    /// # Errors
    /// If the formula is not valid according to the Unimod molecular formula format, with some help on what is going wrong.
    /// ```rust
    /// use mzcore::prelude::*;
    /// assert!(MolecularFormula::unimod("H(25) C(8) 13C(7) N 15N(2) O(3)").is_ok());
    /// assert!(MolecularFormula::unimod("H(6) C(4) N(2) dHex").is_ok());
    /// assert_eq!(MolecularFormula::unimod("C(1) 13C(1) H(6)"), Ok(molecular_formula!(C 1 [13 C 1] H 6)));
    ///
    /// ```
    pub fn unimod(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        Self::unimod_inner(&Context::none().lines(0, value), value, 0..value.len())
    }

    /// This parses a substring of the given string as a Unimod molecular formula definition.
    /// Additionally, this allows passing a base context to allow to set the line index and source
    /// and other properties. Note that the base context is assumed to contain the full line at
    /// line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid Unimod molecular formula string.
    pub fn unimod_inner<'a>(
        base_context: &Context<'a>,
        value: &'a str,
        range: Range<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let (mut index, end) = range.bounds(value.len());

        let mut formula = Self::default();

        let mut isotope = None;
        let mut last_name: Option<(usize, String)> = None;
        while index < end {
            match (value.as_bytes()[index], last_name.as_ref()) {
                (b'(', Some((last_name_i, last_name_s))) => {
                    let length = value
                        .chars()
                        .skip(index + 1)
                        .take_while(|c| *c == '-' || *c == '+' || c.is_ascii_digit())
                        .count();
                    let num = value[index + 1..index + 1 + length]
                        .parse::<i32>()
                        .map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Unimod chemical formula",
                                format!("The element amount {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((0, index + 1, length)),
                            )
                        })?;
                    match parse_unimod_composition_brick(
                        value,
                        *last_name_i..last_name_i + last_name_s.len(),
                    )? {
                        Brick::Element(el) => {
                            if let Err(err) = formula.add((el, isotope.take(), num)) {
                                return Err(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid Unimod chemical formula",
                                    err.reason(),
                                    base_context.clone().add_highlight((
                                        0,
                                        *last_name_i..last_name_i + last_name_s.len(),
                                    )),
                                ));
                            }
                        }
                        Brick::Formula(f) => formula += f * num,
                    }
                    last_name = None;
                    index += length + 2;
                    if value.as_bytes()[index - 1] != b')' {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid Unimod chemical formula",
                            "The amount of an element should be closed by ')'",
                            base_context.clone().add_highlight((0, index - 1, 1)),
                        ));
                    }
                }
                (b' ', Some((last_name_i, last_name_s))) => {
                    match parse_unimod_composition_brick(
                        value,
                        *last_name_i..last_name_i + last_name_s.len(),
                    )? {
                        Brick::Element(el) => {
                            if let Err(err) = formula.add((el, isotope.take(), 1)) {
                                return Err(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid Unimod chemical formula",
                                    err.reason(),
                                    base_context.clone().add_highlight((
                                        0,
                                        *last_name_i..last_name_i + last_name_s.len(),
                                    )),
                                ));
                            }
                        }
                        Brick::Formula(f) => formula += f,
                    }
                    last_name = None;
                    index += 1;
                }
                (b' ', None) => {
                    index += 1;
                }
                (n, _) if n.is_ascii_digit() => {
                    let length = value
                        .chars()
                        .skip(index)
                        .take_while(char::is_ascii_digit)
                        .count();
                    isotope = Some(value[index..index + length].parse::<NonZeroU16>().map_err(
                        |err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Unimod chemical formula",
                                format!("The isotope {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((0, index, length)),
                            )
                        },
                    )?);
                    index += length;
                }
                (n, _) if n.is_ascii_alphabetic() => {
                    if let Some((_, name)) = last_name.as_mut() {
                        name.push(n as char);
                    } else {
                        last_name = Some((index, (n as char).to_string()));
                    }
                    index += 1;
                }
                _ => {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid Unimod chemical formula",
                        "Unexpected character, use an element or one of the unimod shorthands. Eg: 'H(13) C(12) N O(3)'.",
                        base_context.clone().add_highlight((0, index, 1)),
                    ));
                }
            }
        }
        if let Some((last_name_i, last_name_s)) = last_name {
            match parse_unimod_composition_brick(
                value,
                last_name_i..last_name_i + last_name_s.len(),
            )? {
                Brick::Element(el) => {
                    if let Err(err) = formula.add((el, isotope.take(), 1)) {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid Unimod chemical formula",
                            err.reason(),
                            base_context
                                .clone()
                                .add_highlight((0, last_name_i..last_name_i + last_name_s.len())),
                        ));
                    }
                }
                Brick::Formula(f) => formula += f,
            }
        }
        Ok(formula)
    }
}
