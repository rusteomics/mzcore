use std::num::NonZeroU16;

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, MolecularFormula},
    helper_functions::{explain_number_error, str_starts_with},
};

impl MolecularFormula {
    /// Parse a molecular formula from a PSI-MOD formula.
    /// # Errors
    /// If the formula is not valid according to the PSI-MOD molecular formula format, with some help on what is going wrong.
    /// ```rust
    /// use mzcore::prelude::*;
    /// assert!(MolecularFormula::psi_mod("(12)C -5 (13)C 5 H 1 N 3 O -1 S 9").is_ok());
    /// assert!(MolecularFormula::psi_mod("C 6 H 10 N 0 O 5").is_ok());
    ///
    /// ```
    pub fn psi_mod(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        Self::psi_mod_inner(&Context::none().lines(0, value), value, 0..value.len())
    }

    /// This parses a substring of the given string as a PSI-MOD molecular formula definition.
    /// Additionally, this allows passing a base context to allow to set the line index and source
    /// and other properties. Note that the base context is assumed to contain the full line at
    /// line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid PSI-MOD molecular formula string.
    pub fn psi_mod_inner<'a>(
        base_context: &Context<'a>,
        value: &'a str,
        range: std::ops::Range<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let mut index = range.start;
        let end = range.end.min(value.len());
        let mut isotope = None;
        let mut element = None;
        let bytes = value.as_bytes();
        let mut result = Self::default();
        while index < end {
            match (bytes[index], element) {
                (b'(', _) if isotope.is_none() => {
                    let len = bytes
                        .iter()
                        .skip(index)
                        .position(|c| *c == b')')
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid PSI-MOD molecular formula",
                                "No closing round bracket found",
                                base_context.clone().add_highlight((0, index, 1)),
                            )
                        })?;
                    isotope = Some(
                        value[index + 1..index + len]
                            .parse::<NonZeroU16>()
                            .map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid PSI-MOD molecular formula",
                                    format!("The isotope number {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((0, index + 1, len)),
                                )
                            })?,
                    );
                    index += len + 1;
                }
                (b'-' | b'0'..=b'9', Some(ele)) => {
                    let length = value[index..end]
                        .char_indices()
                        .take_while(|(_, c)| c.is_ascii_digit() || *c == '-')
                        .last()
                        .map_or(0, |(i, c)| i + c.len_utf8());
                    let num = value[index..index + length].parse::<i32>().map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            format!("The element number {}", explain_number_error(&err)),
                            base_context.clone().add_highlight((0, index, length)),
                        )
                    })?;

                    if num != 0
                        && let Err(err) = Self::add(&mut result, (ele, isotope, num))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            err.reason(),
                            base_context.clone().add_highlight((0, index - 1, 1)),
                        ));
                    }
                    element = None;
                    isotope = None;
                    index += length;
                }
                (b' ' | b'\t', _) => index += 1,
                _ => {
                    if let Some(element) = element
                        && let Err(err) = Self::add(&mut result, (element, None, 1))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            err.reason(),
                            base_context.clone().add_highlight((0, index - 1, 1)),
                        ));
                    }

                    let mut found = false;
                    for possible in ELEMENT_PARSE_LIST {
                        if str_starts_with::<true>(&value[index..], possible.0) {
                            element = Some(possible.1);
                            index += possible.0.len();
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            "Not a valid character in formula",
                            base_context.clone().add_highlight((0, index, 1)),
                        ));
                    }
                }
            }
        }
        if isotope.is_some() || element.is_some() {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid PSI-MOD molecular formula",
                "Last element missed a count",
                base_context.clone().add_highlight((0, index, 1)),
            ))
        } else {
            Ok(result)
        }
    }
}
