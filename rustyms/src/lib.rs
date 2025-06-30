#![doc = include_str!("../README.md")]

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

/// Contains all things related to annotations (MS2 spectrum annotations that is).
pub mod annotation;
/// Contains all things related to the underlying chemistry.
pub mod chemistry;
pub mod error;
/// Contains all things related to fragments and fragmentation.
pub mod fragment;
pub mod glycan;
mod isobaric_sets;
pub mod ontology;
/// Contains all things related to tolerances and structures to handle multiple mass/formula options.
pub mod quantities;
#[cfg(feature = "rand")]
/// Only available with features `rand`.
mod rand;
/// Contains all things related to sequences, amongst others amino acids and peptidoforms.
pub mod sequence;
pub mod spectrum;
pub mod system;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::annotation::{
        AnnotatableSpectrum,
        model::{FragmentationModel, MatchingParameters},
    };
    pub use crate::chemistry::{
        Chemical, Element, MassMode, MolecularCharge, MolecularFormula, MultiChemical,
    };
    pub use crate::fragment::Fragment;
    pub use crate::isobaric_sets::{
        BuildingBlocks, TerminalBuildingBlocks, building_blocks, find_isobaric_sets,
    };
    pub use crate::sequence::{
        AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, HasCompoundPeptidoformIon,
        HasPeptidoformImpl, HasPeptidoformIon, IsAminoAcid, Peptidoform, PeptidoformIon, Protease,
        SequenceElement, SequencePosition,
    };
    pub use crate::spectrum::RawSpectrum;
}

#[macro_use]
extern crate uom;

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use crate::prelude::*;

    use super::*;

    #[test]
    fn simple_fragments() {
        let peptide = Peptidoform::pro_forma("WFWF", None)
            .unwrap()
            .into_linear()
            .unwrap();
        let fragments = peptide.generate_theoretical_fragments(
            system::isize::Charge::new::<system::e>(1),
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
            .generate_theoretical_fragments(system::isize::Charge::new::<system::e>(1), model);
        let annotated =
            spectrum[0].annotate(peptide, &fragments, &parameters, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
