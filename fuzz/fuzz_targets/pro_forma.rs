//! Fuzz target for parsing ProForma sequences
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data)
            && let Ok(def) = mzcore::sequence::CompoundPeptidoformIon::pro_forma(s, None)
        {
            // Enforce that all displayed peptides are actually valid ProForma according to the parser
            mzcore::sequence::CompoundPeptidoformIon::pro_forma(&def.to_string(), None).unwrap();
        }
    });
}
