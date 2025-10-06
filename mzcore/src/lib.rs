#![doc = include_str!("../README.md")]
#![expect(macro_use_extern_crate)] // Could not get uom to work without
use bincode as _;

#[macro_use]
mod helper_functions;

/// Contains all things related to the underlying chemistry.
pub mod chemistry;
pub mod glycan;
mod isobaric_sets;
pub mod ontology;
mod parse_json;
/// Contains all things related to tolerances and structures to handle multiple mass/formula options.
pub mod quantities;
#[cfg(feature = "rand")]
/// Only available with features `rand`.
mod rand;
/// Contains all things related to sequences, amongst others amino acids and peptidoforms.
pub mod sequence;
pub mod system;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::chemistry::{
        Chemical, Element, MassMode, MolecularCharge, MolecularFormula, MultiChemical,
    };
    pub use crate::isobaric_sets::{
        BuildingBlocks, TerminalBuildingBlocks, building_blocks, find_isobaric_sets,
    };
    pub use crate::sequence::{
        AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, HasCompoundPeptidoformIon,
        HasPeptidoformImpl, HasPeptidoformIon, IsAminoAcid, Peptidoform, PeptidoformIon, Protease,
        SequenceElement, SequencePosition,
    };
}

#[macro_use]
extern crate uom;
