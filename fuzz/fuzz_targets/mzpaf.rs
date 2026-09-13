//! Fuzz target for parsing mzPAF sequences
use afl::*;

fn main() {
    fuzz!(|data: &str| {
        let _unused = mzannotate::fragment::Fragment::mz_paf_strict(
            data,
            &mzcore::ontology::STATIC_ONTOLOGIES,
            &[],
        );
    });
}
