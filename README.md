# Oxidise your MS data analysis

A set of libraries to handle peptide centric mass spectrometry calculations. Built to handle very complex peptidoforms in a sensible way. Centered around the following HUPO-PSI standards:
- [ProForma](https://www.psidev.info/proforma) A standard notation for proteo/peptidoforms allowing for highly complex definitions
- [mzSpecLib](https://www.psidev.info/mzspeclib) A standard notation for spectral libraries
- [mzPAF](https://www.psidev.info/mzpaf) A standard notation of peak fragment annotation
- [mzTab](https://www.psidev.info/mztab-specifications) A standard notation for matched peptidoforms from database and _de novo_ searches

For raw data centered HUPO-PSI standards support (eg mzML, USI) see [mzdata](https://crates.io/crates/mzdata).

| Crate        | Crates.io                                                                                           | Docs                                                                                                              |
| ------------ | --------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------- |
| 🦀 mzcore     | [![Crates.io](https://img.shields.io/crates/v/mzcore.svg)](https://crates.io/crates/mzcore)         | [![mzcore documentation](https://docs.rs/mzcore/badge.svg)](https://docs.rs/mzcore)                               |
| 🦀 mzannotate | [![Crates.io](https://img.shields.io/crates/v/mzannotate.svg)](https://crates.io/crates/mzannotate) | [![mzannotate documentation](https://docs.rs/mzannotate/badge.svg)](https://docs.rs/mzannotate)                   |
| 🦀 mzalign    | [![Crates.io](https://img.shields.io/crates/v/mzalign.svg)](https://crates.io/crates/mzalign)       | [![mzalign documentation](https://docs.rs/mzalign/badge.svg)](https://docs.rs/mzalign)                            |
| 🦀 imgt       | [![Crates.io](https://img.shields.io/crates/v/imgt.svg)](https://crates.io/crates/imgt)             | [![imgt documentation](https://docs.rs/imgt/badge.svg)](https://docs.rs/imgt)                                     |
| 🦀 mzident    | [![Crates.io](https://img.shields.io/crates/v/mzident.svg)](https://crates.io/crates/mzident)       | [![mzident documentation](https://docs.rs/mzident/badge.svg)](https://docs.rs/mzident)                            |
| 🦀 mzcv       | [![Crates.io](https://img.shields.io/crates/v/mzcv.svg)](https://crates.io/crates/mzcv)             | [![mzcv documentation](https://docs.rs/mzcv/badge.svg)](https://docs.rs/mzcv)                                     |
| 🐍 rustyms-py | [![PyPI version](https://badge.fury.io/py/rustyms.svg)](https://badge.fury.io/py/rustyms)           | [![Python Docs](https://readthedocs.org/projects/rustyms/badge/?version=latest)](https://rustyms.readthedocs.io/) |

## Features

- mzcore
  - Read [ProForma](https://github.com/HUPO-PSI/ProForma) sequences (complete 2.0 and nearly complete 2.1)
  - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking
  - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
  - Extensive support for glycans, including generating bitmap and vector images
- mzannotate
  - Generate theoretical fragments with control over the fragmentation model from any ProForma peptidoform
    - Complex features supported: chimeric spectra, cross-links (also disulfides), modifications of unknown position
    - Generate peptide backbone (a, b, c, x, y, and z) and satellite ion fragments (d, v, and w)
    - Generate glycan fragments (B, Y, and internal fragments)
  - Integrated with [mzdata](https://crates.io/crates/mzdata)
  - Read and write [mzSpecLib](https://www.psidev.info/mzspeclib) and [mzPAF](https://www.psidev.info/mzpaf)
  - Match spectra to the generated fragments
- mzalign
  - [Align peptides based on mass](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00188)
  - Consecutive alignment of one sequence on a stretch of multiple sequences
  - Indexed alignment for fast alignments for big datasets
  - Multiple sequence alignment based on the same mass-based alignment 
- imgt
  - Fast access to the IMGT database of antibody germlines
- mzident
  - Reading of multiple PSM file formats (amongst others: [mzTab](https://www.psidev.info/mztab-specifications), Fasta, MaxQuant, MSFragger, Novor, OPair, Peaks, and Sage)
  - Writing of mzTab files
- mzcv
  - Handle ontologies both statically included and runtime updating
- rustyms-py
  - Python bindings are provided to several core components of the libraries. Go to the [Python documentation](https://rustyms.readthedocs.io/) for more information.

## Supported formats 

The final goal would be to support all open standards (or at least the ones that are (widely) used) for both reading and writing. Below is the list of formats that are currently supported.

| Format | Version |  crate | Reading | Writing | Comment |
| --- | --- | ---| --- | --- | --- |
| [ProForma](https://github.com/HUPO-PSI/ProForma) | 2.0 & 2.1 | mzcore | ✅ | ✅ | Nearly full 2.1 support (full support is planned) |
| [mzPAF](https://www.psidev.info/mzpaf) | 1.0 |mzannotate | ✅ | ✅ | |
| [mzSpecLib](https://www.psidev.info/mzspeclib) | 1.0 | mzannotate | ✅ | ✅ | Not all metadata is used |
| FASTA | - | mzident | ✅ | ❌ | |
| [mzTab](https://www.psidev.info/mztab-specifications) | 1.0 | mzident | ✅ | ✅ | Not all metadata is accessible, peptides and small molecules are ignored |
| [Spectrum Sequence List (SSL)](https://skyline.ms/home/software/BiblioSpec/wiki-page.view?name=BiblioSpec%20input%20and%20output%20file%20formats) | - | mzident | ✅ | ❌ | Small molecules are ignored |

For raw data related formats (MGF/mzML/USI) see [mzdata](https://crates.io/crates/mzdata).

These formats are envisioned to have support for. Open an issue if you have a need for these or if you have some thoughts on the implementation. PRs to add support for these are also very welcome. Regardless of inclusion in this list any (open) standard can be suggested for inclusion.

| Format | Version |  crate | Comment |
| --- | --- | ---| --- | 
| [mzIdentML](https://www.psidev.info/peff) | 1.3 | mzident | Including support for cross-linked identifications and mzSpecLib like annotated spectra |
| [PEFF](https://www.psidev.info/mzidentml) | 1.0 | mzident |  |

# Folder organisation

## mzcore/mzalign/mzannotate/mzcv/mzident/imgt

These are the main libraries. This contains all source code, databases (Unimod etc) and example data.

## examples

Some examples on how to use the libraries provided here, see the readme file in the examples themselves for more details.

## fuzz

The harness to fuzz test the libraries for increased stability, see the readme for more details.

## rustyms-py

This Rust library provides python bindings (using pyO3) for rustyms.

## rustyms-generate-databases

Using the `rustyms-generate-databases` the definitions for the elemental data can be updated. See the readme on the download locations for all databases. Then run `cargo run -p rustyms-generate-databases` (from the root folder of this repository).

## mzcore-update

Using the `mzcore-update` the modification databases can be updated. All databases expect RESID will be downloaded and updated. For RESID download the file `ftp://ftp.proteininformationresource.org/pir_databases/other_databases/resid/RESIDUES.XML` and place this at `mzcore-update/data/RESID.xml`, a version is already provided there are as RESID is not developed any more you likely do not need to download the file. Then run `cargo run --releases -p mzcore-update` (from the root folder of this repository).

## imgt-update

Using the `imgt-update` the definitions for the germlines can be updated. Put the imgt.dat.Z file in the `imgt-update/data` directory and unpack it (this can be downloaded from https://www.imgt.org/download/LIGM-DB/imgt.dat.Z). Then run `cargo run --release -p imgt-update` (from the root folder of this repository).

# Contributing

Any contribution is welcome (especially adding/fixing documentation as that is very hard to do as main developer).