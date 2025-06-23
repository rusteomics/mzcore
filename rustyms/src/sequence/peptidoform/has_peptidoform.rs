use std::sync::Arc;

use crate::{
    prelude::{CompoundPeptidoformIon, Peptidoform, PeptidoformIon},
    sequence::{AtLeast, AtMax},
};

/// A structure that has a reference to a peptidoform.
pub trait HasPeptidoformImpl {
    /// The complexity level for this peptidoform containing structure
    type Complexity;
    /// The peptidoform
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity>;
}

/// A structure that has a reference to a peptidoform. Implement [HasPeptidoformImpl] instead.
#[diagnostic::on_unimplemented(
    message = "Implement HasPeptidoformImpl on this type instead of using this type directly"
)]
pub trait HasPeptidoform<Complexity> {
    /// Get a reference to a peptidoform.
    fn peptidoform(&self) -> &Peptidoform<Complexity>;
}

impl<Complexity> HasPeptidoformImpl for Peptidoform<Complexity> {
    type Complexity = Complexity;
    fn peptidoform(&self) -> &Self {
        self
    }
}

impl<Complexity> HasPeptidoformImpl for &Peptidoform<Complexity> {
    type Complexity = Complexity;
    fn peptidoform(&self) -> &Peptidoform<Complexity> {
        self
    }
}

impl<T: HasPeptidoformImpl> HasPeptidoformImpl for Arc<T> {
    type Complexity = T::Complexity;
    fn peptidoform(&self) -> &Peptidoform<T::Complexity> {
        self.as_ref().peptidoform()
    }
}

impl<T: HasPeptidoformImpl, Complexity: AtLeast<T::Complexity>> HasPeptidoform<Complexity> for T {
    fn peptidoform(&self) -> &Peptidoform<Complexity> {
        self.peptidoform().as_ref()
    }
}

/// A structure that has a reference to a peptidoform ion.
pub trait HasPeptidoformIon {
    /// Get a reference to a peptidoform ion.
    fn peptidoform_ion(&self) -> &PeptidoformIon;
}

/// A structure that has a reference to a compound peptidoform ion.
pub trait HasCompoundPeptidoformIon {
    /// Get a reference to a compound peptidoform ion.
    fn compound_peptidoform_ion(&self) -> &CompoundPeptidoformIon;
}
