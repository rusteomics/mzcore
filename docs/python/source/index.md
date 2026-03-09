# rustyms Python bindings

Python bindings to the [rustyms](https://docs.rs/rustyms/) library for proteomics
and mass spectrometry.

## Quickstart

Python bindings are provided to several core components of the rustyms library, including:

- {py:class}`~rustyms.Element` Chemical elements
- {py:class}`~rustyms.MolecularFormula` Molecular formulas
- {py:class}`~rustyms.AminoAcid` Amino acids
- {py:class}`~rustyms.Modification` Amino acid modifications, with support for mass shifts,
  chemical formulas, Unimod, and PSI-MOD labels
- {py:class}`~rustyms.AmbiguousModification` {py:class}`~rustyms.Modification` with ambiguous
  localization.
- {py:class}`~rustyms.Fragment` Theoretical fragment ion
- {py:class}`~rustyms.SequenceElement` One position in a peptide sequence with amino acid and
  modifications
- {py:class}`~rustyms.PeptidoformIonSet` Peptide sequence, modifications, and charge, using
  [ProForma 2.0](https://proforma.readthedocs.io) (see {ref}`ProForma support` for more
  information)
- {py:class}`~rustyms.RawPeak` A single peak in a mass spectrum
- {py:class}`~rustyms.RawSpectrum` A mass spectrum without any annotations
- {py:class}`~rustyms.AnnotatedPeak` A single peak in a mass spectrum with annotations
- {py:class}`~rustyms.AnnotatedSpectrum` A mass spectrum with annotations

{py:mod}`rustyms` can, for example, annotate a mass spectrum from a
[ProForma 2.0](https://proforma.readthedocs.io) peptidoforms string:

```python
import rustyms

# Create a new spectrum
raw_spectrum = rustyms.RawSpectrum(
    title="spectrum_1",
    num_scans=6,
    rt=10,
    precursor_mass=436.12634,
    precursor_charge=2,
    mz_array=[72.04444, 148.06048, 175.05362, 263.08742, 290.08056, 366.09661],
    intensity_array=[100, 600, 300, 400, 500, 200],
)

# Create a new peptide from a ProForma 2.0 string
peptide = rustyms.PeptidoformIonSet("ACDE/2")

# Annotate the spectrum with the peptide
annotated_spectrum = raw_spectrum.annotate(peptide, "cid_hcd")
```

## Installation

Install with pip:

```bash
pip install rustyms
```

## Contributing and development

See {ref}`Contributing` for information on the development setup and how to contribute to the
rustyms Python bindings.

```{toctree}
:hidden:
:includehidden:
:maxdepth: 2

About <self>
api
contributing
```
