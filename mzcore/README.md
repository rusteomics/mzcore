# Match those fragments!

Handle mass spectrometry data in Rust. This crate is set up to handle very complex peptides with
loads of ambiguity and complexity. It pivots around the [`CompoundPeptidoformIon`](crate::sequence::CompoundPeptidoformIon),
[`PeptidoformIon`](crate::sequence::PeptidoformIon) and [`Peptidoform`](crate::sequence::Peptidoform) which encode the
[ProForma](https://github.com/HUPO-PSI/ProForma) specification. Additionally, this crate enables the
reading of [mgf](spectrum::mgf), doing [spectrum annotation](crate::annotation::AnnotatableSpectrum::annotate)
(BU/MD/TD), finding [isobaric sequences](crate::prelude::find_isobaric_sets), doing [alignments of peptides](crate::align::align)
, accessing the [IMGT germline database](crate::imgt), and [reading identified peptide files](crate::identification).

## Library features

 - Read [ProForma](https://github.com/HUPO-PSI/ProForma) sequences (complete specification supported: 'level 2-ProForma + top-down compliant + cross-linking compliant + glycans compliant + mass spectrum compliant')
 - Generate theoretical fragments with control over the fragmentation model from any ProForma peptidoform/proteoform
   - Generate theoretical fragments for chimeric spectra
   - Generate theoretical fragments for cross-links (also disulfides)
   - Generate theoretical fragments for modifications of unknown position
   - Generate peptide backbone (a, b, c, x, y, and z) and satellite ion fragments (w, d, and v)
   - Generate glycan fragments (B, Y, and internal fragments)
 - Integrated with [mzdata](https://crates.io/crates/mzdata) for reading raw data files
 - Match spectra to the generated fragments
 - [Align peptides based on mass](https://pubs.acs.org/doi/10.1021/acs.jproteome.4c00188)
 - Fast access to the IMGT database of antibody germlines
 - Reading of multiple identified peptide file formats (Fasta, MaxQuant, MSFragger, Novor, OPair, Peaks, Sage, and many more)
 - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
 - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking

## Compilation features

Rustyms ties together multiple smaller modules into one cohesive structure.
It has multiple features which allow you to slim it down if needed (all are enabled by default).
* `align` - gives access to mass based alignment of peptides.
* `identification` - gives access to methods reading many different identified peptide formats.
* `imgt` - enables access to the IMGT database of antibodies germline sequences, with annotations.
* `isotopes` - gives access to generation of an averagine model for isotopes, also enables two additional dependencies.
* `rand` - allows the generation of random peptides.
* `rayon` - enables parallel iterators using rayon, mostly for `imgt` but also in consecutive align.
* `mzdata` - enables integration with [mzdata](https://github.com/mobiusklein/mzdata) which has more advanced raw file support.
* `glycan-render` - enables the rendering to SVGs for glycans and glycan fragments
* `glycan-render-bitmap` - enables the rendering to bitmaps for glycans, by enabling the optional dependencies zeno and swash
