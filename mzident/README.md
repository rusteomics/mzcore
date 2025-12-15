# Handle PSM files

Handling many different formats of PSM files. Supports these formats:

- mzTab
- Fasta
- Spectrum Sequence List (SSL)
- mzSpecLib (only with feature `mzannotate`)

And output from the following programs:

- DeepNovo
- PointNovo
- BiatNovo
- PGPointNovo
- InstaNovo
- MaxQuant
- MetaMorpheus
- MSFragger
- NovoB
- Novor
- OPair
- Peaks
- PepNet
- π-HelixNovo
- π-PrimeNovo
- pLink
- PLGS
- PowerNovo
- Proteoscape
- pUniFind
- Sage

## Compilation features

* `mzannotate` - Adds mzannotate as a dependency and allow mzSpecLib spectra to be used as IdentifiedPeptiform and allow other formats to parse annotated spectra
