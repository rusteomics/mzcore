//! Ontologies or controlled vocabularies or databases of modifications
//!
//! Use [`Ontologies`] to interact with these in th most pleasent way. Or use [`STATIC_ONTOLOGIES`]
//! if only the static data is needed.
mod custom;
mod gnome;
mod ontologies;
mod ontology_modification;
mod psimod;
mod resid;
mod unimod;
mod xlmod;

pub use custom::*;
pub use gnome::*;
pub use ontologies::*;
pub use psimod::*;
pub use resid::*;
pub use unimod::*;
pub use xlmod::*;
