use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _ = rustyms::sequence::Peptidoform::sloppy_pro_forma(
                s,
                0..s.len(),
                None,
                &rustyms::sequence::SloppyParsingParameters::default(),
            );
        }
    });
}
