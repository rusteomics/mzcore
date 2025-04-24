//! Read in the annotations from peptide identification sources

#[macro_use]
mod common_parser;

pub mod csv;
mod file_format;
mod formats;
mod general;
mod identified_peptide;

use crate::*;
pub use file_format::*;
pub use formats::*;
pub use general::*;
pub use identified_peptide::*;
