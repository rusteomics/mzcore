//! Text handling utilities

#[cfg(feature = "search-index")]
pub(crate) fn tags(text: &str) -> impl Iterator<Item = [u8; 3]> {
    text.char_indices().filter_map(move |(from, _)| {
        let mut chars = text[from..].chars();
        let a = simplify_char(chars.next()?);
        let b = simplify_char(chars.next()?);
        let c = simplify_char(chars.next()?);
        Some([a, b, c])
    })
}

#[cfg(feature = "search-index")]
const fn simplify_char(c: char) -> u8 {
    if c.is_ascii() { c as u32 as u8 } else { 0 }
}

pub(crate) fn levenshtein_distance(word1: &str, word2: &str) -> usize {
    let word2_len = word2.chars().count();
    let mut row1 = (0..).take(word2_len + 1).collect::<Vec<_>>();
    let mut row2 = vec![0; word2_len + 1];

    for (index1, char1) in word1.chars().enumerate() {
        row2[0] = index1 + 1;

        for (index2, char2) in word2.chars().enumerate() {
            let del = row1[index2 + 1] + 1;
            let ins = row2[index2] + 1;
            let sub = row1[index2] + usize::from(char1 != char2);
            row2[index2 + 1] = del.min(ins).min(sub);
        }

        row1.clone_from_slice(&row2);
    }

    row2[word2_len - 1]
}
