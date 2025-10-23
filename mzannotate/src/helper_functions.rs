use std::{
    collections::HashMap,
    hash::Hash,
    num::{IntErrorKind, ParseIntError},
    ops::{Bound, Range, RangeBounds},
    str::FromStr,
};

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

pub(crate) trait RangeMaths<Other>
where
    Self: Sized,
{
    fn add_start(&self, amount: Other) -> Self;
}

impl RangeMaths<isize> for Range<usize> {
    fn add_start(&self, amount: isize) -> Self {
        let new_start = self.start.saturating_add_signed(amount);
        Self {
            start: new_start,
            end: self.end.max(new_start),
        }
    }
}

impl RangeMaths<usize> for Range<usize> {
    fn add_start(&self, amount: usize) -> Self {
        let new_start = self.start.saturating_add(amount);
        Self {
            start: new_start,
            end: self.end.max(new_start),
        }
    }
}

/// Find the enclosed text by the given symbols, assumes a single open is already read just before
/// the start, guarantees to only pick full characters. Returns the byte offset (including start)
/// into the text where the closing bracket was found.
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

pub(crate) fn merge_hashmap<K, V>(one: HashMap<K, V>, two: HashMap<K, V>) -> HashMap<K, V>
where
    V: std::ops::MulAssign + Default,
    K: Eq + Hash,
{
    let mut new = one;
    for (key, value) in two {
        let v = new.entry(key).or_default();
        *v *= value;
    }
    new
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
