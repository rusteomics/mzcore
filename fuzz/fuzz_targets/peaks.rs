//! Fuzz target for Peaks CSV file parsing
use std::io::BufReader;

use afl::*;
use mzident::{PSMSource, PeaksPSM};

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(csv) = mzcore::csv::parse_csv_raw(BufReader::new(data), b',', None, None) {
            let _unused: Vec<_> =
                PeaksPSM::parse_many(csv, &mzcore::ontology::STATIC_ONTOLOGIES, false, None)
                    .collect();
        }
    });
}
