//! Fuzz target for parsing ProForma sequences
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused = mzcore::sequence::CompoundPeptidoformIon::pro_forma(s, None);
        }
    });
}
