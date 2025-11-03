use std::ops::RangeBounds;

use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularFormula},
    helper_functions::{RangeExtension, next_num, str_starts_with},
    quantities::Multi,
};

impl MolecularFormula {
    /// Parse RESID formulas: `C 2 H 3 N 1 O 1 +` or `C 4 H 5 N 1 O 3, C 4 H 6 N 2 O 2`. If the `range` (byte range in the given line) is not specified it defaults to the full line.
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains non UTF8 symbols.
    pub fn from_resid(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Multi<Self>, BoxedError<'_, BasicKind>> {
        let mut multi = Vec::new();
        let mut start = 0;
        for part in value[range.start_index()..range.end_index(value.len())].split(',') {
            multi.push(Self::from_resid_single(
                value,
                range.start_index() + start..range.start_index() + start + part.len(),
            )?);
            start += part.len() + 1;
        }
        Ok(multi.into())
    }

    /// Parse RESID formulas: `C 2 H 3 N 1 O 1 +` but does not allow multi formulas (split with commas). If the `range` (byte range in the given line) is not specified it defaults to the full line.
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains non UTF8 symbols.
    pub fn from_resid_single(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Self, BoxedError<'_, BasicKind>> {
        let (mut index, end) = range.bounds(value.len().saturating_sub(1));
        let mut result = Self::default();
        while index <= end {
            trim(&mut index, value);
            let mut element = None;
            let mut amount: i32 = 1;
            for possible in ELEMENT_PARSE_LIST {
                if str_starts_with::<true>(&value[index..], possible.0) {
                    element = Some(possible.1);
                    index += possible.0.len();
                    break;
                }
            }
            if element.is_none() {
                if value[index..].starts_with('+') {
                    element = Some(Element::Electron);
                    index += 1;
                    amount = -1;
                } else if value[index..].starts_with('-') {
                    element = Some(Element::Electron);
                    index += 1;
                }
            }
            if let Some(element) = element.take() {
                trim(&mut index, value);
                if let Some(number) = next_num(value.as_bytes(), index, false) {
                    index += number.0;
                    amount *= number.1 as i32;
                }
                if let Err(err) = result.add((element, None, amount)) {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid RESID molecular formula",
                        err.reason(),
                        Context::line(
                            None,
                            value,
                            index,
                            value[index..]
                                .chars()
                                .next()
                                .map(char::len_utf8)
                                .unwrap_or_default(),
                        ),
                    ));
                }
                trim(&mut index, value);
            } else {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid RESID molecular formula",
                    format!("Not a valid character in formula, now has: {result:?}"),
                    Context::line(
                        None,
                        value,
                        index,
                        value[index..]
                            .chars()
                            .next()
                            .map(char::len_utf8)
                            .unwrap_or_default(),
                    ),
                ));
            }
        }
        Ok(result)
    }
}

fn trim(index: &mut usize, text: &str) {
    *index = *index
        + text[*index..]
            .chars()
            .take_while(char::is_ascii_whitespace) // defined to be one byte each
            .count();
}
