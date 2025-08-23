//! Fuzz target for parsing mzPAF sequences
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused = rustyms::fragment::parse_mzPAF(s, None);
        }
    });
}
