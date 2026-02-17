use std::{
    num::{IntErrorKind, ParseIntError},
    ops::{Bound, Range, RangeBounds},
    path::Path,
    str::FromStr,
};

pub(crate) trait InvertResult<T, E> {
    /// # Errors
    /// If any of the errors contained within has an error.
    fn invert(self) -> Result<Option<T>, E>;
}

impl<T, E> InvertResult<T, E> for Option<Result<T, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.map_or_else(|| Ok(None), |o| o.map(|v| Some(v)))
    }
}
impl<T, E> InvertResult<T, E> for Option<Option<Result<T, E>>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.flatten()
            .map_or_else(|| Ok(None), |o| o.map(|v| Some(v)))
    }
}
impl<T, E> InvertResult<T, E> for Option<Result<Option<T>, E>> {
    fn invert(self) -> Result<Option<T>, E> {
        self.unwrap_or_else(|| Ok(None))
    }
}

pub(crate) trait RangeExtension
where
    Self: Sized,
{
    fn start_index(&self) -> usize;
    // Give the max of the end index (inclusive) or the upper bound
    fn end_index(&self, upper_bound: usize) -> usize;
}

impl<Ra: RangeBounds<usize>> RangeExtension for Ra {
    fn start_index(&self) -> usize {
        match self.start_bound() {
            Bound::Unbounded => 0,
            Bound::Included(s) => *s,
            Bound::Excluded(s) => s + 1,
        }
    }

    fn end_index(&self, upper_bound: usize) -> usize {
        match self.end_bound() {
            Bound::Unbounded => upper_bound,
            Bound::Included(s) => *s.min(&upper_bound),
            Bound::Excluded(s) => ((*s).saturating_sub(1)).min(upper_bound),
        }
    }
}

/// Helper function to check extensions in filenames
pub(crate) fn check_extension(filename: impl AsRef<Path>, extension: impl AsRef<Path>) -> bool {
    filename
        .as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case(extension.as_ref()))
}

/// Find the enclosed text by the given symbols, assumes a single open is already read just before the start, guarantees to only pick full characters
pub(crate) fn end_of_enclosure(text: &str, start: usize, open: u8, close: u8) -> Option<usize> {
    let mut state = 1;
    for (i, ch) in text.as_bytes()[start..].iter().enumerate() {
        // Check if this byte is a full character (is_char_boundary also works on index==len)
        if text.is_char_boundary(start + i) && text.is_char_boundary(start + i + 1) {
            if *ch == open {
                state += 1;
            } else if *ch == close {
                state -= 1;
                if state == 0 {
                    return Some(start + i);
                }
            }
        }
    }
    None
}

/// Split the given range based on the separator.
/// This also takes brackets into account and these take precedence over the separator searched for.
pub(crate) fn split_with_brackets(
    text: &str,
    range: Range<usize>,
    separator: u8,
    open: u8,
    close: u8,
) -> Vec<Range<usize>> {
    let mut state: usize = 0;
    let mut index = range.start;
    let mut last_field = range.start;
    let mut fields = Vec::new();
    while index < range.end {
        if !text.is_char_boundary(index) {
            index += 1;
            continue;
        }
        if index + 1 < text.len() && !text.is_char_boundary(index + 1) {
            index += 1;
            continue;
        }
        let ch = text.as_bytes()[index];
        if ch == open {
            state += 1;
        } else if ch == close {
            state = state.saturating_sub(1);
        } else if ch == separator && state == 0 {
            fields.push(last_field..index);
            last_field = index + 1;
        }
        index += 1;
    }
    fields.push(last_field..index);
    fields
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn test_split_with_brackets() {
    assert_eq!(
        split_with_brackets(
            "23-CHEMMOD:+15.995,23-[MS, MS:1001524, fragment neutral loss, 63.998285]",
            0..72,
            b',',
            b'[',
            b']'
        ),
        vec![0..18, 19..72]
    );
    assert_eq!(
        split_with_brackets(
            "0[MS,MS:1001876, modification probability, 0.1]|23[MS,MS:1001876, modification probability, 0.9]-UNIMOD:35",
            0..106,
            b',',
            b'[',
            b']'
        ),
        vec![0..106]
    );
    assert_eq!(
        split_with_brackets("0[,,,[,,]],,[,,l;]hj", 0..20, b',', b'[', b']'),
        vec![0..10, 11..11, 12..20]
    );
}

/// Get the next number starting at the byte range given, returns length in bytes, boolean indicating if the number is positive, and the number.
/// # Errors
/// Returns none if the number is too big to fit in a `Number`.
pub(crate) fn next_number<const ALLOW_SIGN: bool, const FLOATING_POINT: bool, Number: FromStr>(
    line: &str,
    range: impl RangeBounds<usize>,
) -> Option<(usize, bool, Result<Number, Number::Err>)> {
    let start = range.start_index();
    let end = range.end_index(line.len() - 1);
    let mut positive = true;
    let mut sign_set = false;
    let mut chars = line[start..=end].char_indices().peekable();
    if ALLOW_SIGN {
        match chars.peek() {
            Some((_, '-')) => {
                positive = false;
                sign_set = true;
            }
            Some((_, '+')) => {
                sign_set = true;
            }
            _ => (),
        }
        if sign_set {
            let _ = chars.next();
        }
    }

    let mut consumed = usize::from(sign_set);
    chars
        .take_while(|(_, c)| {
            if c.is_ascii_digit() || (FLOATING_POINT && ".eE+-".contains(*c)) {
                consumed += 1;
                true
            } else {
                false
            }
        })
        .last()
        .map(|(end_index, c)| {
            (
                consumed,
                positive,
                line[start..start + end_index + c.len_utf8()].parse::<Number>(),
            )
        })
}

/// To be used as `The xx number ` + the explanation from here (does not have a dot).
pub(crate) const fn explain_number_error(error: &ParseIntError) -> &'static str {
    match error.kind() {
        IntErrorKind::Empty => "is empty",
        IntErrorKind::InvalidDigit => "contains an invalid character",
        IntErrorKind::NegOverflow => "is too small to fit in the internal representation",
        IntErrorKind::PosOverflow => "is too big to fit in the internal representation",
        IntErrorKind::Zero => "is zero, which is not allowed here",
        _ => "is not a valid number",
    }
}

/// Count the number of digits, but only if it did not use exponential notation
pub(crate) fn float_digits(text: &str) -> Option<u8> {
    if text.is_ascii()
        && !text.eq_ignore_ascii_case("nan")
        && !text.eq_ignore_ascii_case("+nan")
        && !text.eq_ignore_ascii_case("-nan")
        && !text.eq_ignore_ascii_case("inf")
        && !text.eq_ignore_ascii_case("+inf")
        && !text.eq_ignore_ascii_case("-inf")
        && !text.eq_ignore_ascii_case("infinity")
        && !text.eq_ignore_ascii_case("+infinity")
        && !text.eq_ignore_ascii_case("-infinity")
    {
        let text = text
            .as_bytes()
            .iter()
            .position(|b| *b == b'e' || *b == b'E')
            .map_or(text, |pos| &text[..pos]);
        text.as_bytes()
            .iter()
            .position(|b| *b == b'.')
            .map_or(Some(0), |pos| u8::try_from(text[pos + 1..].len()).ok())
    } else {
        None
    }
}
