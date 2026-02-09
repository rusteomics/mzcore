use std::{
    collections::HashMap,
    hash::Hash,
    num::{IntErrorKind, ParseIntError},
    ops::{Bound, RangeBounds},
    path::Path,
    str::FromStr,
};

use crate::sequence::SequencePosition;

/// Handle a `ParserResult` by combining the errors into the `$errors` and returning from the enclosing scope if necessary.
macro_rules! handle {
    ($errors:ident, $call:expr) => {
        match $call {
            Ok((res, w)) => {
                combine_errors(&mut $errors, w);
                res
            }
            Err(errs) => {
                combine_errors(&mut $errors, errs);
                return Err($errors);
            }
        }
    };
    (single $errors:ident, $call:expr) => {
        match $call {
            Ok(res) => res,
            Err(err) => {
                combine_error(&mut $errors, err);
                return Err($errors);
            }
        }
    };
    (fail $errors:ident, $error:expr) => {
        combine_error(&mut $errors, $error);
        return Err($errors);
    };
}

pub(crate) fn peptide_range_contains(
    range: &impl RangeBounds<usize>,
    peptide_length: usize,
    position: SequencePosition,
) -> bool {
    match position {
        SequencePosition::NTerm => range.start_index() == 0,
        SequencePosition::Index(i) => range.contains(&i),
        SequencePosition::CTerm => range.end_index(peptide_length) == peptide_length,
    }
}

pub(crate) trait ResultExtensions<T, E> {
    /// # Errors
    /// If any of the errors contained within has an error.
    fn flat_err(self) -> Result<T, E>;
}

impl<T, E> ResultExtensions<T, E> for Result<T, Result<T, E>> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => Ok(o),
            Err(r) => r,
        }
    }
}

impl<T, E> ResultExtensions<T, E> for Result<Result<T, E>, E> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => o,
            Err(r) => Err(r),
        }
    }
}

pub(crate) trait RangeExtension
where
    Self: Sized,
{
    fn start_index(&self) -> usize;
    // Give the max of the end index (inclusive) or the upper bound
    fn end_index(&self, upper_bound: usize) -> usize;
    fn bounds(&self, upper_bound: usize) -> (usize, usize) {
        (self.start_index(), self.end_index(upper_bound))
    }
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

/// Split a string into chunks of text separated by whitespace with the offset before each chunk returned for nice error generation.
pub(crate) fn split_ascii_whitespace(input: &str) -> Vec<(usize, &str)> {
    let mut index = input.chars().take_while(char::is_ascii_whitespace).count();
    let mut chunks = Vec::new();
    while index < input.len() {
        let chunk_len = input[index..]
            .char_indices()
            .take_while(|(_, c)| !c.is_ascii_whitespace())
            .last()
            .map_or(0, |(i, c)| i + c.len_utf8());
        chunks.push((index, &input[index..index + chunk_len]));
        index += chunk_len;
        index += input[index..]
            .char_indices()
            .take_while(|(_, c)| c.is_ascii_whitespace())
            .last()
            .map_or(0, |(i, c)| i + c.len_utf8());
    }
    chunks
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn test_split_ascii_whitespace() {
    assert_eq!(
        split_ascii_whitespace("C7 D10 H2 N4"),
        vec![(0, "C7"), (3, "D10"), (7, "H2"), (10, "N4")]
    );
}

/// Get the index of the next copy of the given char (looking at the byte value, does not guarantee full character)
pub(crate) fn next_char(chars: &[u8], start: usize, char: u8) -> Option<usize> {
    for (i, ch) in chars[start..].iter().enumerate() {
        if *ch == char {
            return Some(start + i);
        }
    }
    None
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

/// Find the enclosed text by the given symbols, assumes a single open is already read just before the start.
/// This also takes brackets '[]' into account and these take precedence over the enclosure searched for.
pub(crate) fn end_of_enclosure_with_brackets(
    text: &str,
    start: usize,
    open: u8,
    close: u8,
) -> Option<usize> {
    let mut state = 1;
    let mut index = start;
    while index < text.len() {
        if !text.is_char_boundary(index) {
            index += 1;
            continue;
        }
        if index + 1 < text.len() && !text.is_char_boundary(index + 1) {
            index += 1;
            continue;
        }
        let ch = text.as_bytes()[index];
        if ch == b'[' {
            index = end_of_enclosure(text, index + 1, b'[', b']')?;
        }
        if ch == open {
            state += 1;
        } else if ch == close {
            state -= 1;
            if state == 0 {
                return Some(index);
            }
        }
        index += 1;
    }
    None
}

/// Get the next number, returns length in bytes and the number.
/// # Panics
/// If the text is not valid UTF-8.
/// # Errors
/// Returns none if the number is too big to fit in a `isize`.
pub(crate) fn next_num(
    chars: &[u8],
    mut start: usize,
    allow_only_sign: bool,
) -> Option<(usize, isize)> {
    let mut sign = 1;
    let mut sign_set = false;
    if chars.get(start) == Some(&b'-') {
        sign = -1;
        start += 1;
        sign_set = true;
    } else if chars.get(start) == Some(&b'+') {
        start += 1;
        sign_set = true;
    }
    let len = chars[start..]
        .iter()
        .take_while(|c| c.is_ascii_digit())
        .count();
    if len == 0 {
        if allow_only_sign && sign_set {
            Some((1, sign))
        } else {
            None
        }
    } else {
        let num: isize = std::str::from_utf8(&chars[start..start + len])
            .unwrap()
            .parse()
            .ok()?;
        Some((usize::from(sign_set) + len, sign * num))
    }
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

/// Get a canonicalised u64 for f64 to be able to hash f64, based on the `ordered_float` crate (MIT license)
pub(crate) fn f64_bits(value: f64) -> u64 {
    if value.is_nan() {
        0x7ff8_0000_0000_0000_u64 // CANONICAL_NAN_BITS
    } else {
        (value + 0.0).to_bits() // The +0.0 is to guarantee even handling of negative and positive zero
    }
}

pub(crate) fn merge_hashmap<K, V>(
    mut into: HashMap<K, V>,
    from: &HashMap<K, V>,
    default_into: &V,
    default_from: &V,
) -> HashMap<K, V>
where
    V: std::ops::MulAssign + Default + Clone,
    K: Eq + Hash + Clone,
{
    for (key, value) in from {
        let v = into
            .entry(key.clone())
            .or_insert_with(|| default_into.clone());
        *v *= value.clone();
    }
    for (key, value) in &mut into {
        if from.get(key).is_none() {
            *value *= default_from.clone();
        }
    }
    into
}

/// Implement a binary operator for all ref cases after the implementation for the ref-ref case (assumes deref operator works)
macro_rules! impl_binop_ref_cases {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty, $o:ty) => {
        impl $imp<$u> for &'_ $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(self, &other)
            }
        }

        impl<'a> $imp<&'a $u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: &'a $u) -> $o {
                $imp::$method(&self, other)
            }
        }

        impl $imp<$u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(&self, &other)
            }
        }
    };
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

/// Check if 'a' starts with 'b' with or without ignoring casing
pub(crate) fn str_starts_with<const IGNORE_CASING: bool>(a: &str, b: &str) -> bool {
    if a.len() >= b.len() {
        for (a, b) in a.chars().zip(b.chars()) {
            if IGNORE_CASING && !a.eq_ignore_ascii_case(&b) || !IGNORE_CASING && a != b {
                return false;
            }
        }
        true
    } else {
        false
    }
}

#[allow(clippy::missing_panics_doc)]
#[test]
fn starts_with() {
    assert!(str_starts_with::<false>("aaabbb", "a"));
    assert!(str_starts_with::<false>("aaabbb", "aa"));
    assert!(str_starts_with::<false>("aaabbb", "aaa"));
    assert!(!str_starts_with::<false>("aaabbb", "b"));
    assert!(!str_starts_with::<false>("aaabbb", "ab"));
    assert!(!str_starts_with::<false>("aaabbb", "aab"));
    assert!(str_starts_with::<true>("aaabbb", "a"));
    assert!(str_starts_with::<true>("aaabbb", "aa"));
    assert!(str_starts_with::<true>("aaabbb", "aaa"));
    assert!(str_starts_with::<true>("aaabbb", "A"));
    assert!(str_starts_with::<true>("aaabbb", "AA"));
    assert!(str_starts_with::<true>("aaabbb", "AAA"));
    assert!(str_starts_with::<true>("aaabbb", "aaA"));
    assert!(!str_starts_with::<false>("aaabbb", "A"));
    assert!(!str_starts_with::<false>("aaabbb", "AA"));
    assert!(!str_starts_with::<false>("aaabbb", "AAA"));
    assert!(!str_starts_with::<false>("aaabbb", "aaA"));
}

/// Helper function to check extensions in filenames
pub(crate) fn check_extension(filename: impl AsRef<Path>, extension: impl AsRef<Path>) -> bool {
    filename
        .as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case(extension.as_ref()))
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

#[allow(clippy::missing_panics_doc)]
#[test]
fn test_float_digits() {
    assert_eq!(float_digits("+3.14"), Some(2));
    assert_eq!(float_digits("+3.1415"), Some(4));
    assert_eq!(float_digits(".14"), Some(2));
    assert_eq!(float_digits("1"), Some(0));
    assert_eq!(float_digits("1."), Some(0));
    assert_eq!(float_digits("1E2"), Some(0));
    assert_eq!(float_digits("1.2E2"), Some(1));
    assert_eq!(float_digits("1.2E-2"), Some(1));
    assert_eq!(float_digits("1.21e2"), Some(2));
    assert_eq!(float_digits("1.2E2222"), Some(1));
    assert_eq!(float_digits("0006"), Some(0));
    assert_eq!(float_digits("inf"), None);
    assert_eq!(float_digits("NaN"), None);
}
