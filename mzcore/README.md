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

## Changelog
### 0.2.0 

- Made mass calculations generic, they can now either calculate as MolecularFormula or as Mass directly, the latter is faster (2.5 times so in some benchmarks) but makes it impossible to generate isotopic distributions or to calculate most abundant isotope. To do this the trait Chemical has been updated with a lot of additional functions, but still needs only one function to be imnplemented. No changes should be needed if the previous behaviour is to be retained, except that the trait Chemical has been renamed to Molecule, and MultiChemical to AmbiguousMolecule. Peptidoforms, ions, and sets now implement the Molecule trait as well which should mean that it is easier to calculate the masses there as well.
- CompoundPeptidoformIon has been renamed to PeptidoformIonSet as decided for ProForma 2.1.
- ProForma 2.1 is now released and fully supported in mzcore.
- PlacementRule Anywhere and Terminal have been combined into PlacementRule::Position.
- OpenSMILES v1.0 is now supported (only reading and converting to MolecularFormula).
- Molecular formulas can now be formatted in PSI-MOD and XLMOD notation alongside the already in place ProForma notation.
- Tracked positions of modifications for monosaccharides (not all positions are parsed yet from short IUPAC, open an issue if this is needed).
- Removed all styling options from glycan rendering functions in favour of a struct with these options. New options are control over the modifications, if they should be rendered, and if the positions should be rendered.
- The CrossId has been broken into the different options leading to more precise parsing, and the ability to fix formatting issues in Obo files.
- Added shrink_to_fit functions to trim any unnecessary empty space from Vecs, Strings, and the like.
- Used wider dependency definitions to be less likely to lead to duplicated dependencies.
- Added charges when defined in PSI-MOD
- Added stubs when defined in XLMOD
- Updated static modification databases 
