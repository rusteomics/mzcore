use std::{num::NonZeroU16, ops::RangeBounds};

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularFormula},
    helper_functions::{RangeExtension, explain_number_error},
};

impl MolecularFormula {
    /// Parse ProForma formulas: `[13C2][12C-2]H2N`.
    /// # The specification (copied from ProForma v2)
    /// As no widely accepted specification exists for expressing elemental formulas, we have adapted a standard with the following rules (taken from <https://github.com/rfellers/chemForma>):
    /// ## Formula Rule 1
    /// A formula will be composed of pairs of atoms and their corresponding cardinality (two Carbon atoms: C2). Pairs SHOULD be separated by spaces but are not required to be.
    /// Atoms and cardinality SHOULD NOT be. Also, the Hill system for ordering (<https://en.wikipedia.org/wiki/Chemical_formula#Hill_system>) is preferred, but not required.
    /// ```text
    /// Example: C12H20O2 or C12 H20 O2
    /// ```
    /// ## Formula Rule 2
    /// Cardinalities must be positive or negative integer values. Zero is not supported. If a cardinality is not included with an atom, it is assumed to be +1.
    /// ```text
    /// Example: HN-1O2
    /// ```
    /// ## Formula Rule 3
    /// Isotopes will be handled by prefixing the atom with its isotopic number in square brackets. If no isotopes are specified, previous rules apply. If no isotope is specified, then it is
    /// assumed the natural isotopic distribution for a given element applies.
    /// ```text
    /// Example: [13C2][12C-2]H2N
    /// Example: [13C2]C-2H2N
    /// ```
    /// ## Allow charge
    /// Allows `:z{x}` to define the charge of a formula, eg `:z+1`, `:z-3`. As defined in ProForma 2.1.
    /// ## Allow empty
    /// Allows the string `(empty)` to be used to denote an empty formula
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains not UTF8 symbols.
    pub fn pro_forma<const ALLOW_CHARGE: bool, const ALLOW_EMPTY: bool>(
        value: &str,
    ) -> Result<Self, BoxedError<'_, BasicKind>> {
        Self::pro_forma_inner::<ALLOW_CHARGE, ALLOW_EMPTY>(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
        )
    }

    pub fn pro_forma_inner<'a, const ALLOW_CHARGE: bool, const ALLOW_EMPTY: bool>(
        base_context: &Context<'a>,
        value: &'a str,
        range: impl RangeBounds<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let (mut index, end) = range.bounds(value.len().saturating_sub(1));
        if index > end || end >= value.len() || &value[index..=end] == "(empty)" {
            return if ALLOW_EMPTY {
                Ok(Self::default())
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid ProForma molecular formula",
                    "The formula is empty",
                    base_context.clone().add_highlight((0, range)),
                ))
            };
        }

        let mut element = None;
        let bytes = value.as_bytes();
        let mut result = Self::default();
        'main_parse_loop: while index <= end {
            match bytes[index] {
                b'[' => {
                    // Skip the open square bracket and leading spaces
                    index += 1 + bytes[index + 1..]
                        .iter()
                        .take_while(|b| **b == b' ')
                        .count();
                    let len = bytes
                        .iter()
                        .skip(index)
                        .position(|c| *c == b']')
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ProForma molecular formula",
                                "No closing square bracket found",
                                base_context.clone().add_highlight((0, index, 1)),
                            )
                        })?;
                    let isotope = bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit())
                        .count();
                    let ws1 = bytes[index + isotope..]
                        .iter()
                        .take_while(|b| **b == b' ')
                        .count();
                    let ele = bytes
                        .iter()
                        .skip(index + isotope + ws1)
                        .take_while(|c| c.is_ascii_alphabetic())
                        .count();

                    for possible in ELEMENT_PARSE_LIST {
                        if &value[index + isotope + ws1..index + isotope + ws1 + ele] == possible.0
                        {
                            element = Some(possible.1);
                            break;
                        }
                    }
                    if let Some(parsed_element) = element {
                        let ws2 = bytes[index + isotope + ws1 + ele..]
                            .iter()
                            .take_while(|c| **c == b' ')
                            .count();
                        let num_len = bytes[index + isotope + ws1 + ele + ws2..]
                            .iter()
                            .take_while(|c| **c == b'-' || **c == b'+' || c.is_ascii_digit())
                            .count();
                        let num = value[index + isotope + ws1 + ele + ws2
                            ..index + isotope + ws1 + ele + ws2 + num_len]
                            .parse::<i32>()
                            .map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid ProForma molecular formula",
                                    format!("The element number {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((
                                        0,
                                        index + isotope + ws1 + ele + ws2,
                                        num_len,
                                    )),
                                )
                            })?;
                        let isotope = value[index..index + isotope]
                            .parse::<NonZeroU16>()
                            .map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid ProForma molecular formula",
                                    format!("The isotope number {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((0, index, isotope)),
                                )
                            })?;

                        if let Err(err) =
                            Self::add(&mut result, (parsed_element, Some(isotope), num))
                        {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid ProForma molecular formula",
                                err.reason(),
                                base_context.clone().add_highlight((0, index, len)),
                            ));
                        }
                        element = None;
                        index += len + 1;
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ProForma molecular formula",
                            "Invalid element",
                            base_context
                                .clone()
                                .add_highlight((0, index + isotope, ele)),
                        ));
                    }
                }
                b'-' | b'0'..=b'9' if element.is_some() => {
                    let (num, len) = std::str::from_utf8(
                        &bytes
                            .iter()
                            .skip(index)
                            .take(end - index + 1) // Bind the maximal length if this is used as part of the molecular charge parsing
                            .take_while(|c| c.is_ascii_digit() || **c == b'-')
                            .copied()
                            .collect::<Vec<_>>(),
                    )
                    .map_or_else(
                        |e| panic!("Non UTF8 in ProForma molecular formula, error: {e}"),
                        |v| {
                            (
                                v.parse::<i32>().map_err(|err| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid ProForma molecular formula",
                                        format!(
                                            "The element number {}",
                                            explain_number_error(&err)
                                        ),
                                        base_context.clone().add_highlight((0, index, v.len())),
                                    )
                                }),
                                v.len(),
                            )
                        },
                    );
                    let num = num?;
                    if num != 0
                        && let Err(err) = Self::add(&mut result, (element.unwrap(), None, num))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ProForma molecular formula",
                            err.reason(),
                            base_context.clone().add_highlight((
                                0,
                                index - element.unwrap().symbol().len(),
                                element.unwrap().symbol().len(),
                            )),
                        ));
                    }
                    element = None;
                    index += len;
                }
                b' ' => index += 1,
                b':' if ALLOW_CHARGE => {
                    if Some(&b'z') == bytes.get(index + 1) {
                        index += 2;
                        let num = value[index..=end].parse::<i32>().map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ProForma molecular formula",
                                format!("The charge number is {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((
                                    0,
                                    index,
                                    end.saturating_sub(index),
                                )),
                            )
                        })?;
                        let _ = result.add((Element::Electron, None, -num));
                        break 'main_parse_loop;
                    }
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid ProForma molecular formula",
                        "A charge tag was not set up properly, a charge tag should be formed as ':z<sign><number>'",
                        base_context.clone().add_highlight((
                            0,
                            index.saturating_sub(1),
                            if bytes.len() < index { 1 } else { 2 },
                        )),
                    ));
                }
                _ => {
                    if let Some(element) = element
                        && let Err(err) = Self::add(&mut result, (element, None, 1))
                    {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ProForma molecular formula",
                            err.reason(),
                            base_context.clone().add_highlight((
                                0,
                                index - element.symbol().len(),
                                element.symbol().len(),
                            )),
                        ));
                    }
                    let element_text: String = value[index..].chars().take(2).collect::<String>();
                    for possible in ELEMENT_PARSE_LIST {
                        if element_text.starts_with(possible.0) {
                            element = Some(possible.1);
                            index += possible.0.len();
                            continue 'main_parse_loop;
                        }
                    }
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid ProForma molecular formula",
                        "Not a valid character in formula",
                        base_context.clone().add_highlight((
                            0,
                            index,
                            value[index..]
                                .chars()
                                .next()
                                .map(char::len_utf8)
                                .unwrap_or_default(),
                        )),
                    ));
                }
            }
        }
        if let Some(element) = element
            && let Err(err) = Self::add(&mut result, (element, None, 1))
        {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid ProForma molecular formula",
                err.reason(),
                base_context.clone().add_highlight((
                    0,
                    index - element.symbol().len(),
                    element.symbol().len(),
                )),
            ));
        }
        // Simplify
        result.elements.retain(|el| el.2 != 0);
        if !ALLOW_EMPTY && result.is_empty() {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid ProForma molecular formula",
                "The formula is empty",
                base_context.clone().add_highlight((0, range)),
            ))
        } else {
            Ok(result)
        }
    }
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn fuzz() {
    let _a = MolecularFormula::pro_forma::<true, true>(":");
    let _a = MolecularFormula::pro_forma::<true, true>(":1002\\[d2C-2]H2N");
    let _a = MolecularFormula::pro_forma::<true, true>("+Wv:z-,33U");
    assert!(MolecularFormula::pro_forma::<true, false>("").is_err());
    assert!(MolecularFormula::pro_forma::<true, false>("f{}").is_err());
}
