//! Defines the different levels of complexity a peptide can be.
//! Used for compile time checking for incorrect use of peptides.
use serde::{Deserialize, Serialize};

/// A trait to mark all options for availability of peptidoforms
pub trait PeptidoformAvailability {}
impl PeptidoformAvailability for MaybePeptidoform {}
impl PeptidoformAvailability for PeptidoformPresent {}

/// A structure where a peptidoform might be present
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct MaybePeptidoform;

/// A structure where a peptidoform definitely is present
#[derive(
    Debug, Default, Copy, Clone, PartialEq, PartialOrd, Ord, Eq, Hash, Serialize, Deserialize,
)]
pub struct PeptidoformPresent;

impl From<PeptidoformPresent> for MaybePeptidoform {
    fn from(_value: PeptidoformPresent) -> Self {
        Self
    }
}
