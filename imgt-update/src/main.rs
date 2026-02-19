//! Update the ontologies used in mzcore
use imgt::IMGT;
use itertools::Itertools;
use mzcv::CVIndex;
use std::{
    fs::File,
    io::{BufWriter, Write},
};

fn main() {
    let args: Vec<String> = std::env::args()
        .skip(1)
        .map(|v| v.to_ascii_lowercase())
        .collect();
    let mut index = CVIndex::<IMGT>::empty();

    let errs = if let Some(path) = args.first() {
        index
            .update_from_path([Some(std::path::Path::new(path))], false)
            .unwrap()
    } else {
        index.update_from_url(&[]).unwrap()
    };

    for err in errs {
        println!("{err}");
    }

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
