use context_error::*;

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularFormula},
    helper_functions::{next_num, str_starts_with},
    quantities::Multi,
};

impl MolecularFormula {
    /// Parse a molecular formula from a RESID formula.
    /// # Errors
    /// If the formula is not valid according to the PSI-MOD molecular formula format, with some help on what is going wrong.
    /// ```rust
    /// use mzcore::prelude::*;
    /// assert_eq!(MolecularFormula::resid("C 2 H 3 N 1 O 1 +"), Ok(molecular_formula!(C 2 H 3 N 1 O 1 :z+1).into()));
    /// assert!(dbg!(MolecularFormula::resid("C 4 H 5 N 1 O 3, C 4 H 6 N 2 O 2")).is_ok());
    ///
    /// ```
    pub fn resid(value: &str) -> Result<Multi<Self>, BoxedError<'_, BasicKind>> {
        Self::resid_inner(&Context::none().lines(0, value), value, 0..value.len())
    }

    /// This parses a substring of the given string as a PSI-MOD molecular formula definition.
    /// Additionally, this allows passing a base context to allow to set the line index and source
    /// and other properties. Note that the base context is assumed to contain the full line at
    /// line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid PSI-MOD molecular formula string.
    pub fn resid_inner<'a>(
        base_context: &Context<'a>,
        value: &'a str,
        range: std::ops::Range<usize>,
    ) -> Result<Multi<Self>, BoxedError<'a, BasicKind>> {
        let mut multi = Vec::new();
        let mut start = 0;
        for part in value[range.clone()].split(',') {
            multi.push(Self::resid_single_inner(
                base_context,
                value,
                range.start + start..range.start + start + part.len(),
            )?);
            start += part.len() + 1;
        }
        Ok(multi.into())
    }

    /// Parse a molecular formula from a RESID formula.
    /// # Errors
    /// If the formula is not valid according to the PSI-MOD molecular formula format, with some help on what is going wrong.
    /// ```rust
    /// use mzcore::prelude::*;
    /// assert!(MolecularFormula::resid_single("C 2 H 3 N 1 O 1 +").is_ok());
    /// assert!(MolecularFormula::resid_single("C 4 H 5 N 1 O 3, C 4 H 6 N 2 O 2").is_err());
    ///
    /// ```
    pub fn resid_single(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        Self::resid_single_inner(&Context::none().lines(0, value), value, 0..value.len())
    }

    /// This parses a substring of the given string as a PSI-MOD molecular formula definition.
    /// Additionally, this allows passing a base context to allow to set the line index and source
    /// and other properties. Note that the base context is assumed to contain the full line at
    /// line index 0.
    ///
    /// # Errors
    /// It fails when the string is not a valid PSI-MOD molecular formula string.
    pub fn resid_single_inner<'a>(
        base_context: &Context<'a>,
        value: &'a str,
        range: std::ops::Range<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let mut index = range.start;
        let end = range.end.min(value.len());
        let mut result = Self::default();
        while index < end {
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
                        base_context
                            .clone()
                            .add_highlight((0, index, element.symbol().len())),
                    ));
                }
                trim(&mut index, value);
            } else {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid RESID molecular formula",
                    format!("Not a valid character in formula, now has: {result:?}"),
                    base_context.clone().add_highlight((0, index, 1)),
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
