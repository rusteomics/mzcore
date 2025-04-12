use crate::{
    error::{Context, CustomError},
    helper_functions::{next_num, RangeExtension},
    MolecularFormula, Multi, ELEMENT_PARSE_LIST,
};
use std::ops::RangeBounds;

impl MolecularFormula {
    /// Parse RESID formulas: `C 2 H 3 N 1 O 1 +` or `C 4 H 5 N 1 O 3, C 4 H 6 N 2 O 2`. If the `range` (byte range in the given line) is not specified it defaults to the full line.
    /// # Errors
    /// If the formula is not valid according to the above specification, with some help on what is going wrong.
    /// # Panics
    /// It can panic if the string contains non UTF8 symbols.
    pub fn from_resid(
        value: &str,
        range: impl RangeBounds<usize>,
    ) -> Result<Multi<Self>, CustomError> {
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
    ) -> Result<Self, CustomError> {
        let (mut index, end) = range.bounds(value.len().saturating_sub(1));
        let mut result = Self::default();
        while index <= end {
            trim(&mut index, value);
            let element_text: String = value[index..]
                .chars()
                .take(2)
                .collect::<String>()
                .to_ascii_lowercase();
            let mut element = None;
            let mut amount: i32 = 1;
            for possible in ELEMENT_PARSE_LIST {
                if element_text.starts_with(possible.0) {
                    element = Some(possible.1);
                    index += possible.0.len();
                    break;
                }
            }
            if element.is_none() {
                if element_text.starts_with('+') {
                    element = Some(crate::Element::Electron);
                    index += 1;
                    amount = -1;
                } else if element_text.starts_with('-') {
                    element = Some(crate::Element::Electron);
                    index += 1;
                }
            }
            if let Some(element) = element.take() {
                trim(&mut index, value);
                if let Some(number) = next_num(value.as_bytes(), index, false) {
                    index += number.0;
                    amount *= number.1 as i32;
                }
                if !result.add((element, None, amount)) {
                    return Err(CustomError::error(
                        "Invalid RESID molecular formula",
                        "An element with undefined mass was used",
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
                return Err(CustomError::error(
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
