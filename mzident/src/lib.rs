#![doc = include_str!("../README.md")]

#[macro_use]
mod common_parser;
#[cfg(feature = "mzannotate")]
mod annotated_spectrum;
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

mod helper_functions;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::{IdentifiedPeptidoform, MetaData, open_identified_peptidoforms_file};
}
