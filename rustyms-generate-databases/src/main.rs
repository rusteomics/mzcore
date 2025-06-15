//! Parse modification ontologies and generate the binary blobs for rustyms
use std::path::Path;

mod atomic_masses;
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

fn main() {
    let out_dir = Path::new("rustyms/src/databases");
    build_atomic_masses(out_dir);
    build_gnome_ontology(out_dir);
    build_psi_mod_ontology(out_dir);
    build_resid_ontology(out_dir);
    build_xlmod_ontology(out_dir);
    build_unimod_ontology(out_dir);
}
