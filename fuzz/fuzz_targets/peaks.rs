//! Fuzz target for Peaks CSV file parsing
use std::io::BufReader;

use afl::*;
use mzident::{IdentifiedPeptidoformSource, PeaksData, csv};

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data)
            && let Ok(csv) = csv::parse_csv_raw(BufReader::new(s.as_bytes()), b',', None)
        {
            let _unused: Vec<_> = PeaksData::parse_many(csv, None, false, None).collect();
        }
    });
}
