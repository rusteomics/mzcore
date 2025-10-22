//! Handle mzSpecLib files. For now only the text encoding is supported.
mod analyte;
mod attribute;
mod header;
mod interpretation;
mod protein_description;
mod read;
mod spectrum_description;
mod write;

pub use analyte::*;
pub use attribute::*;
pub use header::*;
pub use interpretation::*;
pub use protein_description::*;
pub use read::*;
pub use spectrum_description::*;
pub use write::*;

/// An ID
pub type Id = u32;
