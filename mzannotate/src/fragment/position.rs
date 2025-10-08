use std::fmt::Debug;

use serde::{Deserialize, Serialize};

use mzcore::{
    glycan::{GlycanPosition, MonoSaccharide},
    prelude::*,
    sequence::{Modification, PeptidePosition},
};

/// Any position on a glycan or a peptide
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum DiagnosticPosition {
    /// A position on a glycan
    Glycan(GlycanPosition, MonoSaccharide),
    /// A position on a compositional glycan (attachment AA + sequence index + the sugar)
    GlycanCompositional(MonoSaccharide, Option<(AminoAcid, SequencePosition)>),
    /// A position on a peptide
    Peptide(PeptidePosition, AminoAcid),
    /// Labile modification
    Labile(Modification),
    /// Reporter ion
    Reporter,
}
