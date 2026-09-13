//! Fuzz target for parsing ProForma sequences
use afl::*;

fn main() {
    fuzz!(|data: &str| {
        if let Ok((def, _warnings)) = mzcore::sequence::PeptidoformIonSet::pro_forma(
            data,
            &mzcore::ontology::STATIC_ONTOLOGIES,
        ) {
            // Enforce that all displayed peptides are actually valid ProForma according to the
            // parser
            mzcore::sequence::PeptidoformIonSet::pro_forma(
                &def.to_string(),
                &mzcore::ontology::STATIC_ONTOLOGIES,
            )
            .unwrap();
        }
    });
}
