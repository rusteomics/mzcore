//! Fuzz target for parsing mzPAF sequences
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused =
                mzannotate::mzspeclib::MzSpecLibTextParser::open(s.as_bytes(), None, None);
        }
    });
}
