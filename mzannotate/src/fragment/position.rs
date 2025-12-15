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

impl mzcore::space::Space for DiagnosticPosition {
    fn space(&self) -> mzcore::space::UsedSpace {
        match self {
            Self::Glycan(p, s) => p.space() + s.space(),
            Self::GlycanCompositional(s, a) => s.space() + a.space(),
            Self::Peptide(p, a) => p.space() + a.space(),
            Self::Labile(m) => m.space(),
            Self::Reporter => mzcore::space::UsedSpace::default(),
        }
        .set_total::<Self>()
    }
}
