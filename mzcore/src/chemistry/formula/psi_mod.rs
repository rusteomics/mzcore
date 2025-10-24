use std::{num::NonZeroU16, ops::RangeBounds};

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, MolecularFormula},
    helper_functions::{RangeExtension, explain_number_error, str_starts_with},
};

impl MolecularFormula {
    /// PSI-MOD: `(12)C -5 (13)C 5 H 1 N 3 O -1 S 9`
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains not UTF8 symbols.
    pub fn from_psi_mod(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Self, BoxedError<'_, BasicKind>> {
        let (mut index, end) = range.bounds(value.len());
        let mut isotope = None;
        let mut element = None;
        let bytes = value.as_bytes();
        let mut result = Self::default();
        while index < end {
            match bytes[index] {
                b'(' if isotope.is_none() => {
                    let len = bytes
                        .iter()
                        .skip(index)
                        .position(|c| *c == b')')
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid PSI-MOD molecular formula",
                                "No closing round bracket found",
                                Context::line(None, value, index, 1),
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
                                    Context::line(None, value, index + 1, len),
                                )
                            })?,
                    );
                    index += len + 1;
                }
                b'-' | b'0'..=b'9' if element.is_some() => {
                    let (num, len) = std::str::from_utf8(
                        &bytes
                            .iter()
                            .skip(index)
                            .take_while(|c| c.is_ascii_digit() || **c == b'-')
                            .copied()
                            .collect::<Vec<_>>(),
                    )
                    .map_or_else(
                        |e| panic!("Non UTF8 in PSI-MOD molecular formula, error: {e}"),
                        |v| {
                            (
                                v.parse::<i32>().map_err(|err| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid PSI-MOD molecular formula",
                                        format!(
                                            "The isotope number {}",
                                            explain_number_error(&err)
                                        ),
                                        Context::line(None, value, index, v.len()),
                                    )
                                }),
                                v.len(),
                            )
                        },
                    );
                    let num = num?;
                    if num != 0
                        && let Err(err) = Self::add(&mut result, (element.unwrap(), isotope, num))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            err.reason(),
                            Context::line(None, value, index - 1, 1),
                        ));
                    }
                    element = None;
                    isotope = None;
                    index += len;
                }
                b' ' => index += 1,
                _ => {
                    if let Some(element) = element
                        && let Err(err) = Self::add(&mut result, (element, None, 1))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD molecular formula",
                            err.reason(),
                            Context::line(None, value, index - 1, 1),
                        ));
                    }

                    let mut found = false;
                    for possible in ELEMENT_PARSE_LIST {
                        if str_starts_with(&value[index..], possible.0, true) {
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
                            Context::line(None, value, index, 1),
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
                Context::line(None, value, index, 1),
            ))
        } else {
            Ok(result)
        }
    }
}
