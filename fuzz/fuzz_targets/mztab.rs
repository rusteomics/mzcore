//! Fuzz target for parsing mzSpecLib files
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        let parser = mzident::MzTabPSM::parse_reader(
            data,
            &mzcore::ontology::STATIC_ONTOLOGIES,
            context_error::Context::default(),
        )
        .unwrap()
        .2;
        let _unused: Vec<_> = parser.collect();
    });
}
