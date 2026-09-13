//! Fuzz target for parsing general peptide sequences aka 'sloppy ProForma'
use afl::*;

fn main() {
    fuzz!(|data: &str| {
        let _unused = mzcore::sequence::Peptidoform::sloppy_pro_forma(
            data,
            &mzcore::ontology::STATIC_ONTOLOGIES,
            &mzcore::sequence::SloppyParsingParameters::default(),
        );
    });
}
