//! Fuzz target for parsing mzSpecLib files
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data)
            && let Ok(parser) =
                mzannotate::mzspeclib::MzSpecLibTextParser::open(s.as_bytes(), None, None)
        {
            let _unused: Vec<_> = parser.collect();
        }
    });
}
