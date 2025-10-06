#![doc = include_str!("../README.md")]

mod helper_functions;
/// Only available with feature `identification`.
pub mod identification;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::identification::{IdentifiedPeptidoform, MetaData};
}
