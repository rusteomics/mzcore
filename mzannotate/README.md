# Match those fragments!

Handle fragment spectra annotation. Has support for [generating theoretical fragments](crate::prelude::PeptidoformFragmentation::generate_theoretical_fragments), [matching theoretical fragments to spectra](crate::annotation::AnnotatableSpectrum::annotate), and calculating many types of scores on [annotated spectra](crate::spectrum::AnnotatedSpectrum).

## Library features

 - [Read](crate::mzspeclib::MzSpecLibTextParser) and [write](crate::mzspeclib::MzSpecLibTextWriter) [mzSpecLib](https://www.psidev.info/mzspeclib) text files
 - [Read](crate::fragment::Fragment::mz_paf) and [write](crate::fragment::ToMzPAF::to_mz_paf) [mzPAF](https://www.psidev.info/mzpaf) peak annotations
 - Generate theoretical fragments with [control over the fragmentation model](crate::annotation::model::FragmentationModel) from any [ProForma](https://www.psidev.info/proforma) peptidoform/proteoform
   - Generate theoretical fragments for chimeric spectra
   - Generate theoretical fragments for cross-links (also disulfides)
   - Generate theoretical fragments for modifications of unknown position
   - Generate peptide backbone (a, b, c, x, y, and z) and satellite (w, d, and v) ion fragments
   - Generate glycan fragments (B, Y, and internal fragments)
 - Integrated with [mzdata](https://crates.io/crates/mzdata)
 - Match spectra to the generated fragments

## Compilation features

* `coloured-errors` - writes out error messages in with colours.

## Changelog
### 0.2.0

- Support calculating the fragments as MolecularFormula or directly as mass. See mzcore 0.3.0 for more details.
- Uses the CV structures from mzcv now.
- Added isotope pattern matching with additional filtering based on cosine similarity
- Added internal fragment generation and matching.