use std::borrow::Cow;
use std::slice::Iter;
use anyhow::*;

use crate::chemistry::composition::ElementalComposition;

/// Any item that can provide its elemental composition
pub trait Chemical {
    /// Returns the elemental composition (molecular formula or percent composition)
    fn composition(&self) -> ElementalComposition;
}

impl<T: Chemical> Chemical for &[T] {
    fn composition(&self) -> ElementalComposition {
        self.iter().map(Chemical::composition).sum()
    }
}

impl<T: Chemical> Chemical for &Vec<T> {
    fn composition(&self) -> ElementalComposition {
        self.iter().map(Chemical::composition).sum()
    }
}

/// Any item that has both a name and a symbol
// FIXME: shall we return Result<&Cow<str>>? this would avoid unwrapping in some implems
pub trait HasNameAndSymbol: Sized {
    fn name(&self) -> Cow<str>;
    fn symbol(&self) -> Cow<str>;
}

/// Any item having a measurable mass (monoisotopic and average)
///
// FIXME: shall we return Result<f64>? this would avoid unwrapping in some implems
pub trait HasMass: Sized {
    fn mono_mass(&self) -> f64;
    fn average_mass(&self) -> Option<f64>;
}

pub trait MolecularEntity: HasNameAndSymbol + HasMass {}

pub trait IsAminoAcid: HasNameAndSymbol + HasMass {
    fn is_valid(&self) -> bool;
    fn single_letter_code(&self) -> u8;
    fn three_letter_code(&self) -> Cow<str>;
}

pub trait IsAminoAcidSeq {
    fn amino_acids_as_bytes(&self) -> Iter<u8>;
    fn length(&self) -> usize;
}

pub trait AminoAcidFactory<T: IsAminoAcid>: Sized {
    fn aa_from_byte(&self, aa: &u8) -> Result<&T>;

    fn aa_iter_from_bytes<'a>(
        &'a self,
        sequence: &'a Cow<[u8]>,
    ) -> Box<dyn Iterator<Item = Result<&'a T>> + 'a> {
        Box::new( sequence.iter().map( |aa| self.aa_from_byte(aa)))
    }
}
