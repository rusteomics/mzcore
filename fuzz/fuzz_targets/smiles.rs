//! Fuzz target for parsing OpenSMILES sequences
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let _unused = mzcore::chemistry::StructuralFormula::from_smiles(s);
        }
    });
}
