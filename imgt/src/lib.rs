#![doc = include_str!("../README.md")]

//! This crate handles parsing the [IMGT LIGM-DB database](https://www.imgt.org/) into structures compatible with mzcore.
//! It additionally stores all regions and annotations. There are two main ways of selecting germline(s), specified by name
//! [`get_germline`](crate::imgt::get_germline) or by building a query over the data [`Selection`](crate::imgt::Selection).
//!
//! <details><summary>Data present per species</summary>
//!
#![doc = include_str!("germlines.md")]
//!
//! </details>

mod combine;
mod cv;
mod fancy;
mod imgt_gene;
mod parse;
mod regions;
mod select;
mod species;
mod structs;

pub use fancy::*;

pub use cv::IMGT;
pub use regions::*;
pub use select::*;
pub use species::*;
