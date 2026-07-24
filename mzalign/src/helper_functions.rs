/// Get the next number, returns length in bytes and the number.
/// # Panics
/// If the text is not valid UTF-8.
/// # Errors
/// Returns none if the number is too big to fit in a `u16`.
pub(crate) fn next_num(chars: &[u8], start: usize) -> Option<(usize, u16)> {
    let len = chars[start..].iter().take_while(|c| c.is_ascii_digit()).count();
    if len == 0 {
        None
    } else {
        let num: u16 = std::str::from_utf8(&chars[start..start + len]).unwrap().parse().ok()?;
        Some((len, num))
    }
}
