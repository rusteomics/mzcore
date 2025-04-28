use serde::{Deserialize, Serialize};

/// The mode of mass to use
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Default, Debug, Serialize, Deserialize,
)]
#[non_exhaustive]
pub enum MassMode {
    /// Monoisotopic mass, use the base isotope to calculate the mass (eg always 12C)
    #[default]
    Monoisotopic,
    /// The average weight, the average between all occurring isotopes (eg something in between 12C and 13C depending on the number of C)
    Average,
    #[cfg(feature = "isotopes")]
    /// The most abundant mass, the most abundant single isotopic species (eg 12C or 13C depending on the number of C).
    ///
    /// Only available with crate feature 'isotopes'.
    MostAbundant,
}

impl std::fmt::Display for MassMode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Average => write!(f, "average weight"),
            Self::Monoisotopic => write!(f, "monoisotopic mass"),
            #[cfg(feature = "isotopes")]
            Self::MostAbundant => write!(f, "most abundant mass"),
        }
    }
}
