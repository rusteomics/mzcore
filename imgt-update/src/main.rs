//! Update the ontologies used in mzcore
use imgt::IMGT;
use itertools::Itertools;
use mzcv::CVIndex;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

fn main() {
    let mut index = CVIndex::<IMGT>::empty();
    index.update_from_url(&[]).unwrap();

    println!(
        "IMGT version: {}, last updated: {}, germlines: {}",
        index.version().version.as_deref().unwrap_or("-"),
        index.version().last_updated().as_deref().unwrap_or("-"),
        index.len()
    );
    index
        .save_to_cache_at(std::path::Path::new("imgt/src/IMGT.dat"))
        .unwrap();

    // Build docs
    let mut docs = BufWriter::new(File::create("imgt/src/germlines.md").unwrap());
    for (species, germlines) in index.data().iter().sorted_by_key(|(s, _)| **s) {
        writeln!(
            docs,
            "## {} / {}

| Kind | V | J | C |
|------|---|---|---|
|IGHV{}
|IGKV{}
|IGLV{}
|IGIV{}

_Number of genes / number of alleles_
",
            species.scientific_name(),
            species.common_name(),
            germlines.h.doc_row(),
            germlines.k.doc_row(),
            germlines.l.doc_row(),
            germlines.i.doc_row(),
        )
        .unwrap();
    }
}
