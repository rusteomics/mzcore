use std::{
    collections::HashMap,
    hash::Hash,
    num::{IntErrorKind, ParseIntError},
    ops::{Bound, Range, RangeBounds},
    path::Path,
    str::FromStr,
};

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
        self.map_or_else(|| Ok(None), |o| o)
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

pub(crate) trait RangeMaths<Other>
where
    Self: Sized,
{
    fn add_start(&self, amount: Other) -> Self;
    fn add_end(&self, amount: Other) -> Self;
    fn sub_start(&self, amount: Other) -> Self;
    fn sub_end(&self, amount: Other) -> Self;
}

impl RangeMaths<isize> for Range<usize> {
    fn add_start(&self, amount: isize) -> Self {
        let new_start = self.start.saturating_add_signed(amount);
        Self {
            start: new_start,
            end: self.end.max(new_start),
        }
    }
    fn add_end(&self, amount: isize) -> Self {
        let new_end = self.end.saturating_add_signed(amount);
        Self {
            start: self.start.min(new_end),
            end: new_end,
        }
    }
    fn sub_start(&self, amount: isize) -> Self {
        let new_start = self.start.saturating_add_signed(-amount);
        Self {
            start: new_start,
            end: self.end.max(new_start),
        }
    }
    fn sub_end(&self, amount: isize) -> Self {
        let new_end = self.end.saturating_add_signed(-amount);
        Self {
            start: self.start.min(new_end),
            end: new_end,
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
    fn add_end(&self, amount: usize) -> Self {
        let new_end = self.end.saturating_add(amount);
        Self {
            start: self.start.min(new_end),
            end: new_end,
        }
    }
    fn sub_start(&self, amount: usize) -> Self {
        let new_start = self.start.saturating_add(amount);
        Self {
            start: new_start,
            end: self.end.max(new_start),
        }
    }
    fn sub_end(&self, amount: usize) -> Self {
        let new_end = self.end.saturating_sub(amount);
        Self {
            start: self.start.min(new_end),
            end: new_end,
        }
    }
}

/// # Errors
/// If the name cannot be recognised or a number is not valid.
pub(crate) fn parse_named_counter<T: Clone>(
    value: &str,
    names: &[(String, T)],
    allow_negative: bool,
) -> Result<Vec<(T, isize)>, String> {
    let mut index = 0;
    let mut output = Vec::new();
    while index < value.len() {
        if value[index..].starts_with(' ') {
            index += 1;
        } else {
            let mut found = false;
            for name in names {
                if value[index..].starts_with(&name.0) {
                    index += name.0.len();
                    let num = &value[index..]
                        .chars()
                        .skip_while(char::is_ascii_whitespace)
                        .take_while(|c| c.is_ascii_digit() || (allow_negative && *c == '-'))
                        .collect::<String>()
                        .trim()
                        .to_string();
                    if num.is_empty() {
                        output.push((name.1.clone(), 1));
                    } else {
                        output.push((
                            name.1.clone(),
                            num.parse()
                                .map_err(|_| format!("Not a valid number '{num}'"))?,
                        ));
                        index += num.len()
                            + value[index..]
                                .chars()
                                .take_while(char::is_ascii_whitespace)
                                .count();
                    }
                    found = true;
                    break; // Names loop
                }
            }
            if !found {
                return Err(format!("Name not recognised {}", &value[index..]));
            }
        }
    }
    Ok(output)
}

/// Split a string into chunks of text separated by whitespace with the offset before each chunk returned for nice error generation.
pub(crate) fn split_ascii_whitespace(input: &str) -> Vec<(usize, &str)> {
    let mut index = input.chars().take_while(char::is_ascii_whitespace).count();
    let mut chunks = Vec::new();
    while index < input.len() {
        let chunk_len = input[index..]
            .chars()
            .take_while(|c| !c.is_ascii_whitespace())
            .count();
        chunks.push((index, &input[index..index + chunk_len]));
        index += chunk_len;
        index += input[index..]
            .chars()
            .take_while(char::is_ascii_whitespace)
            .count();
    }
    chunks
}

/// Helper function to check extensions in filenames
pub(crate) fn check_extension(filename: impl AsRef<Path>, extension: impl AsRef<Path>) -> bool {
    filename
        .as_ref()
        .extension()
        .is_some_and(|ext| ext.eq_ignore_ascii_case(extension.as_ref()))
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

/// A number of characters, used as length or index
pub(crate) type Characters = usize;

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

/// Check if two strings are equal with or without ignoring casing
pub(crate) fn str_eq(a: &str, b: &str, ignore_casing: bool) -> bool {
    if ignore_casing {
        a.eq_ignore_ascii_case(b)
    } else {
        a == b
    }
}

/// Check if 'a' starts with 'b' with or without ignoring casing
pub(crate) fn str_starts_with(a: &str, b: &str, ignore_casing: bool) -> bool {
    for (a, b) in a.chars().zip(b.chars()) {
        if ignore_casing && !a.eq_ignore_ascii_case(&b) || !ignore_casing && a != b {
            return false;
        }
    }
    a.len() >= b.len()
}

#[allow(clippy::missing_panics_doc)]
#[test]
fn starts_with() {
    assert!(str_starts_with("aaabbb", "a", false));
    assert!(str_starts_with("aaabbb", "aa", false));
    assert!(str_starts_with("aaabbb", "aaa", false));
    assert!(!str_starts_with("aaabbb", "b", false));
    assert!(!str_starts_with("aaabbb", "ab", false));
    assert!(!str_starts_with("aaabbb", "aab", false));
    assert!(str_starts_with("aaabbb", "a", true));
    assert!(str_starts_with("aaabbb", "aa", true));
    assert!(str_starts_with("aaabbb", "aaa", true));
    assert!(str_starts_with("aaabbb", "A", true));
    assert!(str_starts_with("aaabbb", "AA", true));
    assert!(str_starts_with("aaabbb", "AAA", true));
    assert!(str_starts_with("aaabbb", "aaA", true));
    assert!(!str_starts_with("aaabbb", "A", false));
    assert!(!str_starts_with("aaabbb", "AA", false));
    assert!(!str_starts_with("aaabbb", "AAA", false));
    assert!(!str_starts_with("aaabbb", "aaA", false));
}
