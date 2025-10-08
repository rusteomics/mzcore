//! Fuzz target for parsing mzPAF sequences
use afl::*;
use mzcore::prelude::CompoundPeptidoformIon;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused =
                mzannotate::fragment::parse_mz_paf(s, None, &CompoundPeptidoformIon::default());
        }
    });
}
