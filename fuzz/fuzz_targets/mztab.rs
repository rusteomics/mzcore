//! Fuzz target for parsing mzSpecLib files
use afl::*;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(_) = std::str::from_utf8(data) {
            let parser =
                mzident::MZTabData::parse_reader(data, &mzcore::ontology::STATIC_ONTOLOGIES);
            let _unused: Vec<_> = parser.collect();
        }
    });
}
