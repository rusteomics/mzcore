#![doc = include_str!("../README.md")]

mod curie;
mod cv_error;
mod cv_index;
mod cv_source;
mod hash_buf_reader;
mod load;
// Currently LZW decoding is only accessible when downloading files, if this is also needed when
// feeding local files directly the API will have to be changed.
#[cfg(feature = "http")]
mod lzw;
mod obo;
pub mod text;

pub use curie::*;
pub use cv_error::*;
pub use cv_index::*;
pub use cv_source::*;
pub use hash_buf_reader::*;
pub use obo::*;
