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
