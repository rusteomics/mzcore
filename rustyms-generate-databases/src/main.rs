use std::path::Path;

#[macro_use]
extern crate uom;

#[path = "../../rustyms/src/system.rs"]
mod system;
#[macro_use]
#[path = "../../rustyms/src/helper_functions.rs"]
mod helper_functions;
#[path = "../../rustyms/src/shared/element.rs"]
mod element;
#[macro_use]
#[path = "../../rustyms/src/shared/formula/mod.rs"]
mod formula;

#[path = "../../rustyms/src/shared/multi.rs"]
mod multi;
#[path = "../../rustyms/src/shared/sequence_position.rs"]
mod sequence_position;

mod atomic_masses;
#[path = "../../rustyms/src/shared/csv.rs"]
mod csv;
#[path = "../../rustyms/src/error/mod.rs"]
pub mod error;
pub mod glycan;
mod gnome;
mod obo;
mod ontology_modification;
mod psi_mod;
mod resid;
mod unimod;
mod xlmod;

use atomic_masses::*;
use gnome::*;
use ontology_modification::*;
use psi_mod::*;
use resid::*;
use unimod::*;
use xlmod::*;

use serde::{Deserialize, Serialize};

use crate::formula::MultiChemical;
use crate::glycan::{GlycanStructure, MonoSaccharide};
use crate::system::OrderedMass;

use ordered_float::OrderedFloat;
use std::cmp::Ordering;

include!("../../rustyms/src/shared/neutral_loss.rs");
include!("../../rustyms/src/shared/modification.rs");
include!("../../rustyms/src/shared/aminoacid.rs");

mod fragment {
    use super::*;
    use serde::{Deserialize, Serialize};
    /// The index in the branches as stored in the structure
    pub type GlycanBranchIndex = usize;
    /// The index in the branches when the branches are sorted on mass, this is used to properly render the names of the branches for human consumption
    pub type GlycanBranchMassIndex = usize;
    include!("../../rustyms/src/shared/glycan_position.rs");
}

impl crate::Element {
    pub fn is_valid(self, _isotope: Option<std::num::NonZeroU16>) -> bool {
        true
    }
}

pub use crate::element::*;
pub use crate::formula::{AmbiguousLabel, MolecularFormula};
pub use crate::multi::Multi;
pub use crate::sequence_position::*;

fn main() {
    let out_dir = Path::new("rustyms/src/databases");
    build_atomic_masses(out_dir);
    build_gnome_ontology(out_dir);
    build_psi_mod_ontology(out_dir);
    build_resid_ontology(out_dir);
    build_unimod_ontology(out_dir);
    build_xlmod_ontology(out_dir);
}
