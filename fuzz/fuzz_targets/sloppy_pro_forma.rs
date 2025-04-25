//! Fuzz target for parsing general peptide sequences aka 'sloppy ProForma'
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused = rustyms::sequence::Peptidoform::sloppy_pro_forma(
                s,
                0..s.len(),
                None,
                &rustyms::sequence::SloppyParsingParameters::default(),
            );
        }
    });
}
