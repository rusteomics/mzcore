use std::{num::NonZeroU16, ops::RangeBounds};

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularFormula},
    helper_functions::{self, RangeExtension, explain_number_error},
};

impl MolecularFormula {
    /// XLMOD: `C7 D10 H2 N4`
    /// `-C1 -H2 O1` `C7 D10 H2 N4` `13C6 H6 O2`
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains not UTF8 symbols.
    pub fn from_xlmod(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Self, BoxedError<'_, BasicKind>> {
        let (start, end) = range.bounds(value.len());
        let mut formula = Self::default();
        for (offset, block) in helper_functions::split_ascii_whitespace(&value[start..end]) {
            let negative = block.starts_with('-');
            let isotope_len = block
                .chars()
                .skip(usize::from(negative))
                .take_while(char::is_ascii_digit)
                .count();
            let number = block.chars().rev().take_while(char::is_ascii_digit).count();
            if number + isotope_len + usize::from(negative) >= block.len() {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid Xlmod molecular formula",
                    "No element is defined",
                    Context::line(None, value, offset, block.len()),
                ));
            }
            let element_len = block.len() - number - isotope_len - usize::from(negative);
            let element = &block[isotope_len + usize::from(negative)..block.len() - number];
            let (mut isotope, element) = if element == "D" {
                if isotope_len == 0 {
                    (Some(NonZeroU16::new(2).unwrap()), Element::H)
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid Xlmod molecular formula",
                        "A deuterium cannot have a defined isotope as deuterium is by definition always isotope 2 of hydrogen",
                        Context::line(None, value, offset + usize::from(negative), isotope_len),
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
                        Context::line(
                            None,
                            value,
                            offset + usize::from(negative) + isotope_len,
                            element_len,
                        ),
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
                                Context::line(
                                    None,
                                    value,
                                    offset + usize::from(negative),
                                    isotope_len,
                                ),
                            )
                        })?,
                );
            }
            let number = if negative { -1 } else { 1 }
                * if number == 0 {
                    1
                } else {
                    block[block.len() - number..block.len()]
                        .parse::<i32>()
                        .map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Xlmod molecular formula",
                                format!("The element count {}", explain_number_error(&err)),
                                Context::line(
                                    None,
                                    value,
                                    offset + usize::from(negative) + isotope_len + element_len,
                                    number,
                                ),
                            )
                        })?
                };
            if let Err(err) = Self::add(&mut formula, (element, isotope, number)) {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid Xlmod molecular formula",
                    err.reason(),
                    Context::line(
                        None,
                        value,
                        offset + usize::from(negative),
                        isotope_len + element_len,
                    ),
                ));
            }
        }
        Ok(formula)
    }
}
