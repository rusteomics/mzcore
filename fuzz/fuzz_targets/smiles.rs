//! Fuzz target for parsing OpenSMILES sequences
use afl::*;

fn main() {
    fuzz!(|data: &str| {
        let _unused = mzcore::chemistry::StructuralFormula::from_smiles(data);
    });
}
