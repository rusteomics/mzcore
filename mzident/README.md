# Handle PSM files

Handling many different formats of PSM files. Supports these formats:

- mzTab
- Fasta
- Spectrum Sequence List (SSL)
- mzSpecLib (only with feature `mzannotate`)

And output from the following programs:

- BiatNovo
- DeepNovo
- InstaNovo
- MaxQuant
- MetaMorpheus
- MSFragger
- NovoB
- Novor
- OPair
- Peaks
- PepNet
- PGPointNovo
- PLGS
- pLink
- PointNovo
- PowerNovo
- Proteoscape
- pUniFind
- Sage
- π-HelixNovo
- π-PrimeNovo

## Compilation features

* `mzannotate` - Adds mzannotate as a dependency and allow mzSpecLib spectra to be used as PSM and allow other formats to parse annotated spectra

## Changelog
### 0.2.0

- Better InstaNovo parsing and version detection thanks to @BioGeek
- Full mzTab metadata support in reading and writing
- Renamed identified peptidoform to PSM to align better with general terms
- Better protein handling with the trait ProteinMetaData
- All software now have a PSI-MS term