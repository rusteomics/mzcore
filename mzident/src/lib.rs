#![doc = include_str!("../README.md")]

#[macro_use]
mod common_parser;
#[cfg(feature = "mzannotate")]
mod annotated_spectrum;
mod file_format;
mod formats;
mod general;
pub mod mztab_writer;
mod peaks_family_id;
mod peptidoform_availability;
mod protein_metadata;
mod psm;
mod psm_metadata;
mod source;
mod spectrum_id;
#[cfg(test)]
mod test;

#[cfg(test)]
use test::*;

pub use file_format::*;
pub use formats::*;
pub use general::*;
pub use peaks_family_id::*;
pub use peptidoform_availability::*;
pub use protein_metadata::*;
pub use psm::*;
pub use psm_metadata::*;
pub use source::*;
pub use spectrum_id::*;

mod helper_functions;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::{PSM, PSMMetaData, ProteinMetaData, open_psm_file};
}
