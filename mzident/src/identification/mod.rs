//! Read in the annotations from peptide identification sources

#[macro_use]
mod common_parser;
mod file_format;
mod formats;
mod general;
mod identified_peptidoform;
mod metadata;
mod peaks_family_id;
mod peptidoform_availability;
mod source;
mod spectrum_id;
#[cfg(test)]
mod test;

pub mod csv;

#[cfg(test)]
use test::*;

pub use file_format::*;
pub use formats::*;
pub use general::*;
pub use identified_peptidoform::*;
pub use metadata::*;
pub use peaks_family_id::*;
pub use peptidoform_availability::*;
pub use source::*;
pub use spectrum_id::*;
