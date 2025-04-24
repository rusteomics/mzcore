#![doc = include_str!("../README.md")]
#![warn(clippy::all, clippy::pedantic, clippy::nursery, missing_docs)]
#![allow(
    clippy::must_use_candidate,
    clippy::cast_precision_loss,
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::wildcard_imports,
    clippy::module_name_repetitions,
    clippy::too_many_lines,
    clippy::too_long_first_doc_paragraph
)]

#[cfg(feature = "align")]
/// Only available with feature `align`.
pub mod align;

#[cfg(feature = "identification")]
/// Only available with feature `identification`.
pub mod identification;

#[cfg(feature = "imgt")]
/// Only available with feature `imgt`.
pub mod imgt;

#[cfg(test)]
mod fragmentation_tests;
#[macro_use]
mod helper_functions;

pub mod annotation;
pub mod chemistry;
pub mod error;
pub mod fragment;
pub mod glycan;
mod isobaric_sets;
pub mod ontology;
pub mod quantities;
#[cfg(feature = "rand")]
/// Only available with features `rand`.
mod rand;
pub mod sequence;
pub mod spectrum;
pub mod system;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::annotation::{
        model::{FragmentationModel, MatchingParameters},
        AnnotatableSpectrum,
    };
    pub use crate::chemistry::{
        Chemical, Element, MassMode, MolecularCharge, MolecularFormula, MultiChemical,
    };
    pub use crate::fragment::Fragment;
    pub use crate::isobaric_sets::{building_blocks, find_isobaric_sets};
    pub use crate::sequence::{
        AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, IsAminoAcid, Peptidoform,
        PeptidoformIon, Protease, SequenceElement, SequencePosition,
    };
    pub use crate::spectrum::RawSpectrum;
}

#[macro_use]
extern crate uom;

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use crate::{
        annotation::{
            model::{FragmentationModel, MatchingParameters},
            AnnotatableSpectrum,
        },
        prelude::*,
    };

    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = Peptidoform::pro_forma("WFWF", None)
            .unwrap()
            .into_linear()
            .unwrap();
        let fragments = peptide.generate_theoretical_fragments(
            system::usize::Charge::new::<system::e>(1),
            FragmentationModel::all(),
        );
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }

    #[test]
    fn simple_matching() {
        let model = FragmentationModel::all();
        let parameters = MatchingParameters::default();
        let spectrum = spectrum::mgf::open("data/example.mgf").unwrap();
        let peptide = CompoundPeptidoformIon::pro_forma("WFWF", None).unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(system::usize::Charge::new::<system::e>(1), model);
        let annotated =
            spectrum[0].annotate(peptide, &fragments, &parameters, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
