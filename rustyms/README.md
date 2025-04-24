# Match those fragments!

Handle mass spectrometry data in Rust. This crate is set up to handle very complex peptides with
loads of ambiguity and complexity. It pivots around the [`CompoundPeptidoformIon`], [`PeptidoformIon`] and [`Peptidoform`]
which encode the [ProForma](https://github.com/HUPO-PSI/ProForma) specification. Additionally
this crate enables the reading of [mgf](rawfile::mgf), doing [spectrum annotation](RawSpectrum::annotate)
(BU/MD/TD), finding [isobaric sequences](find_isobaric_sets), doing [alignments of peptides](align::align)
, accessing the [IMGT germline database](imgt), and [reading identified peptide files](identification).

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
 - Reading of multiple identified peptide file formats (Fasta, MaxQuant, MSFragger, Novor, OPair, Peaks, and Sage)
 - Exhaustively fuzz tested for reliability (using [cargo-afl](https://crates.io/crates/cargo-afl))
 - Extensive use of [uom](https://docs.rs/uom/latest/uom/) for compile time unit checking

## Example usage

```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
# let raw_file_path = "data/annotated_example.mgf";
use rustyms::{prelude::*, system::{usize::Charge, e}};
// Open example raw data (this is the built in mgf reader, look into mzdata for more advanced raw file readers)
let spectrum = rustyms::spectrum::mgf::open(raw_file_path)?;
// Parse the given ProForma definition
let peptide = CompoundPeptidoformIon::pro_forma("[Gln->pyro-Glu]-QVQEVSERTHGGNFD", None)?;
// Generate theoretical fragments for this peptide given EThcD fragmentation
let model = FragmentationModel::ethcd();
let fragments = peptide.generate_theoretical_fragments(Charge::new::<e>(2), model);
let parameters = MatchingParameters::default();
// Annotate the raw data with the theoretical fragments
let annotated = spectrum[0].annotate(peptide, &fragments, &parameters, MassMode::Monoisotopic);
// Calculate a peak false discovery rate for this annotation 
let (fdr, _) = annotated.fdr(&fragments, &parameters, MassMode::Monoisotopic);
// This is the incorrect sequence for this spectrum so the peak FDR will indicate this
# dbg!(&fdr, fdr.peaks_sigma(), fdr.peaks_fdr(), fdr.peaks_score());
assert!(fdr.peaks_sigma() > 2.0);
# Ok(()) }
```

```rust
# fn main() -> Result<(), rustyms::error::CustomError> {
use rustyms::{prelude::*, sequence::SimpleLinear, align::*};
// Check how this peptide compares to a similar peptide (using the feature `align`)
let first_peptide = Peptidoform::pro_forma("IVQEVT", None)?.into_simple_linear().unwrap();
let second_peptide = Peptidoform::pro_forma("LVQVET", None)?.into_simple_linear().unwrap();
// Align the two peptides using mass based alignment
// IVQEVT A
// LVQVET B
// ─  ╶╴
let alignment = align::<4, SimpleLinear, SimpleLinear>(
                  &first_peptide, 
                  &second_peptide,
                  AlignScoring::default(), 
                  AlignType::GLOBAL);
# dbg!(&alignment);
// Calculate some more statistics on this alignment
let stats = alignment.stats();
assert_eq!(stats.mass_similar, 6); // 6 out of the 6 positions are mass similar
# Ok(()) }
```

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
