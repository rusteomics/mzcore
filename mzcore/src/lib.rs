#![doc = include_str!("../README.md")]
#![expect(macro_use_extern_crate)] // Could not get uom to work without
use bincode as _;

#[macro_use]
mod helper_functions;

/// Contains all things related to the underlying chemistry.
pub mod chemistry;
/// Parse CSV files while keeping track of all necessary info to generate great error messages
pub mod csv;
pub mod glycan;
mod isobaric_sets;
pub mod ontology;
/// Contains logic to parse the custom modifications and custom models databases with backwards compatibility
pub mod parse_json;
/// Contains all things related to tolerances and structures to handle multiple mass/formula options.
pub mod quantities;
/// Contains all things related to sequences, amongst others amino acids and peptidoforms.
pub mod sequence;
mod space;
pub mod system;

/// A subset of the types and traits that are envisioned to be used the most, importing this is a good starting point for working with the crate
pub mod prelude {
    pub use crate::chemistry::{
        Chemical, Element, MassMode, MolecularCharge, MolecularFormula, MultiChemical,
    };
    pub use crate::isobaric_sets::{
        BuildingBlocks, TerminalBuildingBlocks, building_blocks, find_isobaric_sets,
    };
    pub use crate::molecular_formula;
    pub use crate::sequence::{
        AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, HasCompoundPeptidoformIon,
        HasPeptidoformImpl, HasPeptidoformIon, IsAminoAcid, Peptidoform, PeptidoformIon, Protease,
        SequenceElement, SequencePosition,
    };
}

#[macro_use]
extern crate uom;

/// The result of a parser, contains the result and a list of warnings if it succeeded and only a list of errors if it failed.
pub type ParserResult<'a, T, Kind> =
    Result<(T, Vec<context_error::BoxedError<'a, Kind>>), Vec<context_error::BoxedError<'a, Kind>>>;
