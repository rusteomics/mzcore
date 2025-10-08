🦀 Rust: [![Crates.io](https://img.shields.io/crates/v/rustyms.svg)](https://crates.io/crates/rustyms) [![rustyms documentation](https://docs.rs/rustyms/badge.svg)](https://docs.rs/rustyms)
🐍 Python: [![PyPI version](https://badge.fury.io/py/rustyms.svg)](https://badge.fury.io/py/rustyms) [![Python Docs](https://readthedocs.org/projects/rustyms/badge/?version=latest)](https://rustyms.readthedocs.io/)

# Match those fragments!

A set of libraries to handle peptide centric mass spectrometry calculations. Built to handle very complex peptidoforms in a sensible way. Centered around the following HUPO-PSI standards:
- [ProForma](https://www.psidev.info/proforma) A standard notation for proteo/peptidoforms allowing for highly complex definitions
- [mzSpecLib](https://www.psidev.info/mzspeclib) A standard notation for spectral libraries
- [mzPAF](https://www.psidev.info/mzpaf) A standard notation of peak fragment annotation
- [mzTab](https://www.psidev.info/mztab-specifications) A standard notation for matched peptidoforms from database and _de novo_ searches

For raw data centered HUPO-PSI standards support (eg mzML, USI) see [mzdata](https://crates.io/crates/mzdata). 

## Features

- mzcore
  - Read [ProForma](https://github.com/HUPO-PSI/ProForma) sequences (complete 2.0 specification supported: 'level 2-ProForma + top-down compliant + cross-linking compliant + glycans compliant + mass spectrum compliant')
  - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking
  - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
  - Extensive support for glycans, including generating bitmap and vector images
- mzannotate
  - Generate theoretical fragments with control over the fragmentation model from any ProForma peptidoform
    - Generate theoretical fragments for chimeric spectra
    - Generate theoretical fragments for cross-links (also disulfides)
    - Generate theoretical fragments for modifications of unknown position
    - Generate peptide backbone (a, b, c, x, y, and z) and satellite ion fragments (d, v, and w)
    - Generate glycan fragments (B, Y, and internal fragments)
  - Integrated with [mzdata](https://crates.io/crates/mzdata)
  - Read and write [mzSpecLib](https://www.psidev.info/mzspeclib) and [mzPAF](https://www.psidev.info/mzpaf)
  - Match spectra to the generated fragments
- mzalign
  - [Align peptides based on mass](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00188)
  - Consecutive alignment of one sequence on a stretch of multiple sequences
  - Indexed alignment for fast alignments for big datasets
- imgt
  - Fast access to the IMGT database of antibody germlines
- mzident
  - Reading of multiple identified peptide file formats (amongst others: [mzTab](https://www.psidev.info/mztab-specifications), Fasta, MaxQuant, MSFragger, Novor, OPair, Peaks, and Sage)
- rustyms-py
  - Python bindings are provided to several core components of the libraries. Go to the [Python documentation](https://rustyms.readthedocs.io/) for more information.

# Folder organisation

## mzcore/mzalign/mzannotate/mzcv/mzident/imgt

These are the main librares. This contains all source code, databases (Unimod etc) and example data.

## examples

Some examples on how to use the libraries provided here, see the readme file in the examples themselves for more details.

## fuzz

The harness to fuzz test the libraries for increased stability, see the readme for more details.

## rustyms-py

This Rust library provides python bindings (using pyO3) for rustyms.

## rustyms-generate-databases

Using the `rustyms-generate-databases` the definitions for the databases can be updated. See the readme on the download locations for all databases. Then run `cargo run -p rustyms-generate-databases` (from the root folder of this repository).

## rustyms-generate-imgt

Using the `rustyms-generate-imgt` the definitions for the germlines can be updated. Put the imgt.dat.Z file in the `rustyms-generate-imgt/data` directory and unpack it (this can be downloaded from https://www.imgt.org/download/LIGM-DB/imgt.dat.Z). Then run `cargo run -p rustyms-generate-imgt` (from the root folder of this repository).

# Contributing

Any contribution is welcome (especially adding/fixing documentation as that is very hard to do as main developer).