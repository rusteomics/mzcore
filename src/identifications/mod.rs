//! Read in the annotations from peptide identification sources

mod common_parser;
#[path = "../shared/csv.rs"]
mod csv;
mod fasta;
mod identified_peptide;
mod novor;
mod opair;
mod peaks;

use self::csv::CsvLine;
use crate::error::CustomError;
pub use fasta::*;
pub use identified_peptide::*;
pub use novor::*;
pub use opair::*;
pub use peaks::*;

#[cfg(test)]
mod novor_tests;
#[cfg(test)]
mod opair_tests;
#[cfg(test)]
mod peaks_tests;

use std::str::FromStr;

impl CsvLine {
    pub fn column_context(&self, column: usize) -> crate::error::Context {
        crate::error::Context::line(
            self.line_index,
            self.line.clone(),
            self.fields[column].start,
            self.fields[column].len(),
        )
    }

    pub fn range_context(&self, range: std::ops::Range<usize>) -> crate::error::Context {
        crate::error::Context::line(self.line_index, self.line.clone(), range.start, range.len())
    }

    pub fn full_context(&self) -> crate::error::Context {
        crate::error::Context::full_line(self.line_index, self.line.clone())
    }
    /// Parse a column into the given format, if erroneous extend the base error with the correct context and return that
    pub fn parse_column<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<F, CustomError> {
        self[column]
            .parse()
            .map_err(|_| base_error.with_context(self.column_context(column)))
    }
    /// Parse a column into the given format, if erroneous extend the base error with the correct context and return that
    pub fn parse_column_or_empty<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<Option<F>, CustomError> {
        let text = &self[column];
        if text.is_empty() || text == "-" {
            Ok(None)
        } else {
            Ok(Some(text.parse().map_err(|_| {
                base_error.with_context(self.column_context(column))
            })?))
        }
    }
}
