//! Defines the different levels of complexity a peptide can be.
//! Used for compile time checking for incorrect use of peptides.
use serde::{Deserialize, Serialize};

/// A [`crate::sequence::Peptidoform`] that (potentially) is linked, either with cross-links or branches
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct Linked;

/// A [`crate::sequence::Peptidoform`] that is not cross-linked or branched, but can use the whole breath of the complexity otherwise
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct Linear;

/// A [`crate::sequence::Peptidoform`] that does not have any of the following:
/// * Labile modifications
/// * Global isotope modifications
/// * Charge carriers, use of charged ions apart from protons
/// * Cyclic structures: inter/intra cross-links or branches
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct SimpleLinear;

/// A [`crate::sequence::Peptidoform`] that does not have any of the following:
/// * Ambiguous modifications
/// * Ambiguous amino acid sequence `(?AA)`
///
/// On top of the outlawed features in [`SimpleLinear`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct SemiAmbiguous;

/// A [`crate::sequence::Peptidoform`] that does not have any of the following:
/// * Ambiguous amino acids (B/Z)
///
/// On top of the outlawed features in [`SemiAmbiguous`].
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct UnAmbiguous;

/// Indicate that a peptide has at max this level of complexity, or lower.
/// `Self` will always be the highest (or identical) of the two complexities.
pub trait AtMax<T> {}
impl<T> AtMax<T> for T {}
impl AtMax<Linked> for Linear {}
impl AtMax<Linked> for SimpleLinear {}
impl AtMax<Linked> for SemiAmbiguous {}
impl AtMax<Linked> for UnAmbiguous {}
impl AtMax<Linear> for SimpleLinear {}
impl AtMax<Linear> for SemiAmbiguous {}
impl AtMax<Linear> for UnAmbiguous {}
impl AtMax<SimpleLinear> for SemiAmbiguous {}
impl AtMax<SimpleLinear> for UnAmbiguous {}
impl AtMax<SemiAmbiguous> for UnAmbiguous {}

/// Indicate that a peptide has at least this level of complexity, or higher.
/// The type parameter will always be the highest (or identical) of the two complexities.
pub trait AtLeast<T> {}
impl<T> AtLeast<T> for T {}
impl AtLeast<Linear> for Linked {}
impl AtLeast<SimpleLinear> for Linked {}
impl AtLeast<SemiAmbiguous> for Linked {}
impl AtLeast<UnAmbiguous> for Linked {}
impl AtLeast<SimpleLinear> for Linear {}
impl AtLeast<SemiAmbiguous> for Linear {}
impl AtLeast<UnAmbiguous> for Linear {}
impl AtLeast<SemiAmbiguous> for SimpleLinear {}
impl AtLeast<UnAmbiguous> for SimpleLinear {}
impl AtLeast<UnAmbiguous> for SemiAmbiguous {}

/// Type level max between two complexity levels. The highest of the two levels is stored in the
/// associated type `HighestLevel`.
pub trait HighestOf<T> {
    /// This is the highest complexity level out of Self and the type parameter
    type HighestLevel;
}
impl<T> HighestOf<T> for T {
    type HighestLevel = T;
}
impl HighestOf<Linked> for Linear {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for SimpleLinear {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for SemiAmbiguous {
    type HighestLevel = Linked;
}
impl HighestOf<Linked> for UnAmbiguous {
    type HighestLevel = Linked;
}
impl HighestOf<Linear> for SimpleLinear {
    type HighestLevel = Linear;
}
impl HighestOf<Linear> for SemiAmbiguous {
    type HighestLevel = Linear;
}
impl HighestOf<Linear> for UnAmbiguous {
    type HighestLevel = Linear;
}
impl HighestOf<SimpleLinear> for SemiAmbiguous {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<SimpleLinear> for UnAmbiguous {
    type HighestLevel = SimpleLinear;
}
impl HighestOf<SemiAmbiguous> for UnAmbiguous {
    type HighestLevel = SemiAmbiguous;
}
impl HighestOf<Linear> for Linked {
    type HighestLevel = Self;
}
impl HighestOf<SimpleLinear> for Linked {
    type HighestLevel = Self;
}
impl HighestOf<SemiAmbiguous> for Linked {
    type HighestLevel = Self;
}
impl HighestOf<UnAmbiguous> for Linked {
    type HighestLevel = Self;
}
impl HighestOf<SimpleLinear> for Linear {
    type HighestLevel = Self;
}
impl HighestOf<SemiAmbiguous> for Linear {
    type HighestLevel = Self;
}
impl HighestOf<UnAmbiguous> for Linear {
    type HighestLevel = Self;
}
impl HighestOf<SemiAmbiguous> for SimpleLinear {
    type HighestLevel = Self;
}
impl HighestOf<UnAmbiguous> for SimpleLinear {
    type HighestLevel = Self;
}
impl HighestOf<UnAmbiguous> for SemiAmbiguous {
    type HighestLevel = Self;
}
