//! Fuzz target for parsing mzSpecLib files
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(parser) = mzannotate::mzspeclib::MzSpecLibTextParser::open(
            data,
            None,
            &mzcore::ontology::STATIC_ONTOLOGIES,
        ) {
            let _unused: Vec<_> = parser.collect();
        }
    });
}
