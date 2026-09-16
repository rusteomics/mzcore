#![doc = include_str!("../README.md")]

//! This crate handles parsing the [IMGT LIGM-DB database](https://www.imgt.org/) into structures compatible with mzcore.
//! It additionally stores all regions and annotations. There are two main ways of selecting
//! germline(s), by using the [`IMGT`] CV access or by building a
//! query over the data [`Selection`].
//!
//! <details><summary>Data present per species</summary>
#![doc = include_str!("germlines.md")]
//!
//! </details>

use flate2 as _;

mod combine;
mod cv;
mod fancy;
mod imgt_gene;
mod parse;
mod regions;
mod select;
mod species;
mod structs;

pub use cv::{IMGT, STATIC_IMGT};
pub use fancy::*;
pub use regions::*;
pub use select::*;
pub use species::*;
