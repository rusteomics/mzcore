use std::{rc::Rc, sync::Arc};

use crate::{
    prelude::{CompoundPeptidoformIon, Peptidoform, PeptidoformIon},
    sequence::AtLeast,
};

/// A structure that has a reference to a peptidoform.
pub trait HasPeptidoformImpl {
    /// The complexity level for this peptidoform containing structure
    type Complexity;
    /// The peptidoform
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity>;
}

/// A structure that has a reference to a peptidoform. Which can be retrieved at a specific
/// [`Complexity`](crate::sequence::Complexity) level. Do not implement this trait directly
/// implement [`HasPeptidoformImpl`] instead.
#[diagnostic::on_unimplemented(
    message = "Implement HasPeptidoformImpl instead of HasPeptidoform<Complexity>"
)]
pub trait HasPeptidoform<Complexity> {
    /// Get a reference to a peptidoform.
    fn cast_peptidoform(&self) -> &Peptidoform<Complexity>;
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

impl<T: HasPeptidoformImpl> HasPeptidoformImpl for Rc<T> {
    type Complexity = T::Complexity;
    fn peptidoform(&self) -> &Peptidoform<T::Complexity> {
        self.as_ref().peptidoform()
    }
}

impl<T: HasPeptidoformImpl> HasPeptidoformImpl for Box<T> {
    type Complexity = T::Complexity;
    fn peptidoform(&self) -> &Peptidoform<T::Complexity> {
        self.as_ref().peptidoform()
    }
}

impl<T: HasPeptidoformImpl, Complexity: AtLeast<T::Complexity>> HasPeptidoform<Complexity> for T {
    fn cast_peptidoform(&self) -> &Peptidoform<Complexity> {
        self.peptidoform().as_ref()
    }
}

/// A structure that has a reference to a peptidoform ion.
pub trait HasPeptidoformIon {
    /// Get a reference to a peptidoform ion.
    fn peptidoform_ion(&self) -> &PeptidoformIon;
}

impl HasPeptidoformIon for PeptidoformIon {
    fn peptidoform_ion(&self) -> &Self {
        self
    }
}

impl HasPeptidoformIon for &PeptidoformIon {
    fn peptidoform_ion(&self) -> &PeptidoformIon {
        self
    }
}

impl<T: HasPeptidoformIon> HasPeptidoformIon for Arc<T> {
    fn peptidoform_ion(&self) -> &PeptidoformIon {
        self.as_ref().peptidoform_ion()
    }
}

/// A structure that has a reference to a compound peptidoform ion.
pub trait HasCompoundPeptidoformIon {
    /// Get a reference to a compound peptidoform ion.
    fn compound_peptidoform_ion(&self) -> &CompoundPeptidoformIon;
}

impl HasCompoundPeptidoformIon for CompoundPeptidoformIon {
    fn compound_peptidoform_ion(&self) -> &Self {
        self
    }
}

impl HasCompoundPeptidoformIon for &CompoundPeptidoformIon {
    fn compound_peptidoform_ion(&self) -> &CompoundPeptidoformIon {
        self
    }
}

impl<T: HasCompoundPeptidoformIon> HasCompoundPeptidoformIon for Arc<T> {
    fn compound_peptidoform_ion(&self) -> &CompoundPeptidoformIon {
        self.as_ref().compound_peptidoform_ion()
    }
}
