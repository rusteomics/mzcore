//! Fuzz target for Peaks CSV file parsing
use std::io::BufReader;

use afl::*;
use mzident::{IdentifiedPeptidoformSource, PeaksData};

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data)
            && let Ok(csv) =
                mzcore::csv::parse_csv_raw(BufReader::new(s.as_bytes()), b',', None, None)
        {
            let _unused: Vec<_> =
                PeaksData::parse_many(csv, &mzcore::ontology::STATIC_ONTOLOGIES, false, None)
                    .collect();
        }
    });
}
