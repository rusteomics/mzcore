//! # Handling CVs
//! This library helps to handle CVs (controlled vocabularies) easily. Start of by implementing
//! [`CVSource`] and [`CVData`] for the needed types.
//!
//! Handles:
//! * Retrieving data items ([`CVIndex::get_by_index`], [`CVIndex::get_by_name`], [`CVIndex::search`])
//! * Statically included CV data ([`CVSource::static_data`])
//! * Updating the CV ([`CVIndex::update`], [`CVIndex::update_from_path`], [`CVIndex::update_from_url`])
//!
//! Automatically stores any runtime given data in a binary index for quick caching for the next
//! time the tool is run. Because all tools that use that same CV definition cache to the same
//! location, this gives automatic syncing between all rusteomics derived tools on the same system.
//!
//! # Features
//! * `search-index` turns on the creation of a trigram index for faster fuzzy match searching.
//!   With this function turned off fuzzy match searching is still available, but slower. The
//!   creation of the index does take initial startup time and memory.
//! * `http` turns on downloading from the internet, as this pulls in `reqwest` this is made
//!   optional to prevent relying on too many dependencies.

mod curie;
mod cv_error;
mod cv_index;
mod cv_source;
mod hash_buf_reader;
mod load;
mod obo;
mod text;

pub use curie::*;
pub use cv_error::*;
pub use cv_index::*;
pub use cv_source::*;
pub use hash_buf_reader::*;
pub use obo::*;
