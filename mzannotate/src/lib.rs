#![doc = include_str!("../README.md")]

#[cfg(test)]
mod fragmentation_tests;

mod aminoacid;
/// Contains all things related to annotations (fragment spectrum annotations that is).
pub mod annotation;
/// Contains all things related to fragments and fragmentation.
pub mod fragment;
pub mod glycan;
mod helper_functions;
mod modification;
mod monosaccharide;
pub mod mzspeclib;
mod peptidoform;
/// Defines annotated spectra
pub mod spectrum;

/// Reexport mzdata to make it easier to get the exact same version.
pub use mzdata;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::annotation::{
        AnnotatableSpectrum,
        model::{FragmentationModel, MatchingParameters},
    };
    pub use crate::fragment::{Fragment, ToMzPAF};
    pub use crate::glycan::GlycanFragmention;
    pub use crate::mzspeclib::{MzSpecLibTextParser, MzSpecLibTextWriter};
    pub use crate::peptidoform::PeptidoformFragmentation;
    pub use crate::spectrum::{AnnotatedPeak, AnnotatedSpectrum};
    pub use crate::term;
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use mzcore::prelude::*;

    use crate::prelude::*;

    #[test]
    fn simple_fragments() {
        let peptide = Peptidoform::pro_forma("WFWF", &mzcore::ontology::STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_linear()
            .unwrap();
        let fragments = peptide.generate_theoretical_fragments(
            mzcore::system::isize::Charge::new::<mzcore::system::e>(1),
            FragmentationModel::all(),
        );
        println!("{}", fragments.len());
        println!("{fragments:?}");
    }
}
