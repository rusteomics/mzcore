//! Parse modification ontologies and generate the binary blobs for rustyms
use std::path::Path;

mod atomic_masses;

use atomic_masses::*;
fn main() {
    let out_dir = Path::new("mzcore/src/databases");
    if !out_dir.exists() {
        std::fs::create_dir(out_dir).unwrap();
    }
    build_atomic_masses(out_dir);
}
