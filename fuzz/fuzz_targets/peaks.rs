use std::io::BufReader;

use afl::*;
use rustyms::identification::IdentifiedPeptideSource;

use std::io::Read;
use std::io::Result;
use std::slice::Iter;

// Thanks to crate 'stringreader'
pub struct StringReader<'a> {
    iter: Iter<'a, u8>,
}

impl<'a> StringReader<'a> {
    /// Wrap a string in a `StringReader`, which implements `std::io::Read`.
    pub fn new(data: &'a str) -> Self {
        Self {
            iter: data.as_bytes().iter(),
        }
    }
}

impl Read for StringReader<'_> {
    fn read(&mut self, buf: &mut [u8]) -> Result<usize> {
        for (i, item) in buf.iter_mut().enumerate() {
            if let Some(x) = self.iter.next() {
                *item = *x;
            } else {
                return Ok(i);
            }
        }
        Ok(buf.len())
    }
}

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(s) = std::str::from_utf8(data) {
            let csv = rustyms::csv::parse_csv_raw(BufReader::new(StringReader::new(s)), b',', None)
                .unwrap();
            let _: Vec<_> =
                rustyms::identification::PeaksData::parse_many(csv, None, false).collect();
        }
    });
}
