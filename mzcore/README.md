# (Proteomics) core for mass spectrometry in Rust

Handle mass spectrometry calculations in Rust. This crate is mostly used for proteomics data in the broadest sense: bottom-up, top-down, cross-linking, glycopeptidomics, and much more. Support for other mass spectrometry fields is not missing by design but by lack of expertise of the authors, so if you are doing something mass spectrometry related and need some extensions to this crate please reach out. 

## Library features

 - Read [ProForma](https://github.com/HUPO-PSI/ProForma) sequences (complete specification supported: 'level 2-ProForma + top-down compliant + cross-linking compliant + glycans compliant + mass spectrum compliant version 2.1')
 - Handle glycans
  - Has GNOme database built-in with all structures
  - Can render glycan structures as well as fragments of structures (use feature `glycan-render`)
 - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
 - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking
 - Handle molecular formula, mass + average weight + most abundant mass (for the latter use feature `isotopes`)
 - Generate isotopic patterns (use feature `isotopes`)
 - Access to Unimod, PSI-MOD, RESID, GNOme, and XL-MOD databases statically and possible to update at runtime

## Compilation features

* `isotopes` - gives access to generation of an averagine model for isotopes, also enables two additional dependencies.
* `glycan-render` - enables the rendering to SVGs for glycans and glycan fragments
* `glycan-render-bitmap` - enables the rendering to bitmaps for glycans, by enabling the optional dependencies zeno and swash
