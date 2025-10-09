#![doc = include_str!("../README.md")]

#[cfg(test)]
mod fragmentation_tests;

mod aminoacid;
/// Contains all things related to annotations (MS2 spectrum annotations that is).
pub mod annotation;
/// Contains all things related to fragments and fragmentation.
pub mod fragment;
pub mod glycan;
mod helper_functions;
mod modification;
mod monosaccharide;
pub mod mzspeclib;
mod peptidoform;
pub mod spectrum;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::annotation::{
        AnnotatableSpectrum,
        model::{FragmentationModel, MatchingParameters},
    };
    pub use crate::fragment::{Fragment, ToMzPAF};
    pub use crate::glycan::GlycanFragmention;
    pub use crate::peptidoform::PeptidoformFragmentation;
    pub use crate::spectrum::{AnnotatedPeak, AnnotatedSpectrum};
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use mzcore::prelude::*;

    use crate::prelude::*;

    #[test]
    fn simple_fragments() {
        let peptide = Peptidoform::pro_forma("WFWF", None)
            .unwrap()
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
