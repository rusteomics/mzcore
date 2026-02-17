//! Parse modification ontologies and generate the binary blobs for rustyms
use std::path::Path;

mod atomic_masses;

use atomic_masses::*;

/// # Panics
/// If parsing the atomic masses table did not work, or if the output directory could not be made.
fn main() {
    let out_dir = Path::new("mzcore/src/databases");
    if !out_dir.exists() {
        std::fs::create_dir(out_dir).unwrap();
    }
    build_atomic_masses(out_dir);
}
