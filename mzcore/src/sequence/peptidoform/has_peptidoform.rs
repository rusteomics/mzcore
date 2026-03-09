use std::{rc::Rc, sync::Arc};

use crate::sequence::{AtLeast, Peptidoform, PeptidoformIon, PeptidoformIonSet};

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
    message = "This structure does not have a `Peptidoform<{Complexity}>`",
    note = "If this is a structure you control implement `HasPeptidoformImpl` instead of `HasPeptidoform<{Complexity}>`"
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

impl<T: HasPeptidoformImpl> HasPeptidoformImpl for &T {
    type Complexity = T::Complexity;
    fn peptidoform(&self) -> &Peptidoform<T::Complexity> {
        (*self).peptidoform()
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
    #[inline]
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

/// A structure that has a reference to a peptidoform ion set.
pub trait HasPeptidoformIonSet {
    /// Get a reference to a peptidoform ion set.
    fn peptidoform_ion_set(&self) -> &PeptidoformIonSet;
}

impl HasPeptidoformIonSet for PeptidoformIonSet {
    fn peptidoform_ion_set(&self) -> &Self {
        self
    }
}

impl HasPeptidoformIonSet for &PeptidoformIonSet {
    fn peptidoform_ion_set(&self) -> &PeptidoformIonSet {
        self
    }
}

impl<T: HasPeptidoformIonSet> HasPeptidoformIonSet for Arc<T> {
    fn peptidoform_ion_set(&self) -> &PeptidoformIonSet {
        self.as_ref().peptidoform_ion_set()
    }
}
