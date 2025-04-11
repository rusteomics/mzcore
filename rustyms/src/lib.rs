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
#[macro_use]
mod formula;

#[doc(hidden)]
#[path = "shared/csv.rs"]
pub mod csv;

/// Contains logic surrounding amino acids, see [`AminoAcid`] for the main structure.
pub mod aminoacid;
mod checked_aminoacid;
mod element;
pub mod error;
mod formula_search;
pub mod fragment;
pub mod glycan;
mod isobaric_sets;
#[cfg(feature = "isotopes")]
/// Only available with feature `isotopes`.
mod isotopes;
mod mass_mode;
pub mod model;
pub mod modification;
mod molecular_charge;
#[path = "shared/multi.rs"]
mod multi;
mod mzpaf;
mod neutral_loss;
pub mod ontologies;
pub mod peptidoform;
pub mod placement_rule;
mod protease;
#[cfg(feature = "rand")]
/// Only available with features `rand`.
mod rand;
pub mod rawfile;
mod sequence_element;
#[path = "shared/sequence_position.rs"]
mod sequence_position;
pub mod spectrum;
pub mod system;
mod tolerance;

pub use aminoacid::{AminoAcid, IsAminoAcid};
pub use checked_aminoacid::CheckedAminoAcid;
pub use element::*;
pub use formula::*;
pub use formula_search::find_formulas;
pub use fragment::Fragment;
pub use isobaric_sets::{building_blocks, find_isobaric_sets};
pub use mass_mode::MassMode;
pub use model::FragmentationModel;
pub use modification::{CrossLinkName, Modification};
pub use molecular_charge::MolecularCharge;
pub use multi::*;
pub use neutral_loss::*;
pub use peptidoform::*;
pub use peptidoform::{CompoundPeptidoformIon, Peptidoform, PeptidoformIon};
pub use protease::*;
pub use sequence_element::SequenceElement;
pub use sequence_position::*;
pub use spectrum::{AnnotatableSpectrum, AnnotatedSpectrum, RawSpectrum};
pub use tolerance::*;

#[macro_use]
extern crate uom;

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use crate::model::MatchingParameters;

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
        let spectrum = rawfile::mgf::open("data/example.mgf").unwrap();
        let peptide = CompoundPeptidoformIon::pro_forma("WFWF", None).unwrap();
        let fragments = peptide
            .generate_theoretical_fragments(system::usize::Charge::new::<system::e>(1), model);
        let annotated =
            spectrum[0].annotate(peptide, &fragments, &parameters, MassMode::Monoisotopic);
        println!("{annotated:?}");
    }
}
