use std::num::NonZeroU16;

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularFormula},
    helper_functions::{self, RangeExtension, explain_number_error},
};

impl MolecularFormula {
    /// Parse a molecular formula from an XL-MOD formula.
    /// # Errors
    /// If the formula is not valid according to the XL-MOD molecular formula format, with some help on what is going wrong.
    /// ```rust
    /// use mzcore::prelude::*;
    /// assert_eq!(
    ///     MolecularFormula::xlmod("C7 D10 H2 N4").unwrap(),
    ///     molecular_formula!(C 7 [2 H 10] H 2 N 4)
    /// );
    /// assert_eq!(
    ///     MolecularFormula::xlmod("-C1 -H2 O1").unwrap(),
    ///     molecular_formula!(C -1 H -2 O 1)
    /// );
    /// assert_eq!(
    ///     MolecularFormula::xlmod("13C6 H6 O2").unwrap(),
    ///     molecular_formula!([13 C 6] H 6 O 2)
    /// );    ///
    /// ```
    pub fn xlmod(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        Self::xlmod_inner(&Context::none().lines(0, value), value, 0..value.len())
    }

    /// This parses a substring of the given string as an XL-MOD molecular formula definition.
    /// Additionally, this allows passing a base context to allow to set the line index and source
    /// and other properties. Note that the base context is assumed to contain the full line at
    /// line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid XL-MOD molecular formula string.
    #[expect(clippy::missing_panics_doc)] // 2 is not zero, but that cannot be proven at compile time
    pub fn xlmod_inner<'a>(
        base_context: &Context<'a>,
        value: &'a str,
        range: std::ops::Range<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let (start, end) = range.bounds(value.len());
        let mut formula = Self::default();
        for (offset, block) in helper_functions::split_ascii_whitespace(&value[start..=end]) {
            let negative = block.starts_with('-');
            let isotope_len = block
                .chars()
                .skip(usize::from(negative))
                .take_while(char::is_ascii_digit)
                .count();
            let number_length = block.chars().rev().take_while(char::is_ascii_digit).count();
            if number_length + isotope_len + usize::from(negative) >= block.len() {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid Xlmod molecular formula",
                    "No element is defined",
                    base_context.clone().add_highlight((0, offset, block.len())),
                ));
            }
            let element_len = block.len() - number_length - isotope_len - usize::from(negative);
            let element = &block[isotope_len + usize::from(negative)..block.len() - number_length];
            let (mut isotope, element) = if element == "D" {
                if isotope_len == 0 {
                    (Some(NonZeroU16::new(2).unwrap()), Element::H)
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid Xlmod molecular formula",
                        "A deuterium cannot have a defined isotope as deuterium is by definition always isotope 2 of hydrogen",
                        base_context.clone().add_highlight((
                            0,
                            offset + usize::from(negative),
                            isotope_len,
                        )),
                    ));
                }
            } else {
                let mut found = None;
                for possible in ELEMENT_PARSE_LIST {
                    if element == possible.0 {
                        found = Some(possible.1);
                        break;
                    }
                }
                if let Some(element) = found {
                    (None, element)
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid Xlmod molecular formula",
                        "Not a valid character in formula",
                        base_context.clone().add_highlight((
                            0,
                            offset + usize::from(negative) + isotope_len,
                            element_len,
                        )),
                    ));
                }
            };
            if isotope.is_none() && isotope_len > 0 {
                isotope = Some(
                    block[usize::from(negative)..usize::from(negative) + isotope_len]
                        .parse::<NonZeroU16>()
                        .map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Xlmod molecular formula",
                                format!("The isotope number {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((
                                    0,
                                    offset + usize::from(negative),
                                    isotope_len,
                                )),
                            )
                        })?,
                );
            }
            let number = if negative { -1 } else { 1 }
                * if number_length == 0 {
                    1
                } else {
                    block[block.len() - number_length..]
                        .parse::<i32>()
                        .map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Xlmod molecular formula",
                                format!("The element count {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((
                                    0,
                                    offset + usize::from(negative) + isotope_len + element_len,
                                    number_length,
                                )),
                            )
                        })?
                };
            if let Err(err) = Self::add(&mut formula, (element, isotope, number)) {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid Xlmod molecular formula",
                    err.reason(),
                    base_context.clone().add_highlight((
                        0,
                        offset + usize::from(negative),
                        isotope_len + element_len,
                    )),
                ));
            }
        }
        Ok(formula)
    }
}
