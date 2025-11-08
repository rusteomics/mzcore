//! The available ontologies

use std::sync::LazyLock;

use context_error::*;
use mzcv::{CVIndex, CVVersion};
use serde::{Deserialize, Serialize};

use crate::{
    ontology::{Custom, Gnome, PsiMod, Resid, Unimod, XlMod},
    sequence::SimpleModification,
};

/// A single shared static access to the static data in the ontologies for cases where no runtime resolution is needed (like tests).
pub static STATIC_ONTOLOGIES: LazyLock<Ontologies> = LazyLock::new(Ontologies::init_static);

/// Handle all ProForma needed ontologies.
///
/// Get a copy via [`Self::init()`], [`Self::init_static()`] (or [`STATIC_ONTOLOGIES`]), or even
/// [`Self::empty`]. Then find modifications using either [`Self::get_by_name`],
/// , [`Self::search`], or the same methods on one particular ontology
/// as `Self::unimod().get_by_index()`.
///
/// ```rust
/// use mzcore::{molecular_formula, ontology::STATIC_ONTOLOGIES, prelude::*};
/// // Using the static version to not do IO in test, but more logical in library users:
/// // let ontologies = Self::init();
/// let ontologies = &STATIC_ONTOLOGIES;
/// // Get modifications by name
/// let modification = ontologies.get_by_name(&[], "Oxidation").unwrap();
/// assert_eq!(modification.formula(), molecular_formula!(O 1));
/// // or by index from a particular ontology
/// let modification2 = ontologies.unimod().get_by_index(&35).unwrap();
/// assert_eq!(modification, modification2);
/// // or search all (or a subset) for fuzzy matches
/// let search = ontologies.search(&[], "Oxidated");
/// assert!(search.contains(&modification));
/// ```
pub struct Ontologies {
    custom: CVIndex<Custom>,
    gnome: CVIndex<Gnome>,
    psimod: CVIndex<PsiMod>,
    resid: CVIndex<Resid>,
    unimod: CVIndex<Unimod>,
    xlmod: CVIndex<XlMod>,
}

impl std::fmt::Debug for Ontologies {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Ontologies")
            .field("Unimod", &self.unimod.len())
            .field("PSI-MOD", &self.psimod.len())
            .field("XL-MOD", &self.xlmod.len())
            .field("GNOme", &self.gnome.len())
            .field("RESID", &self.resid.len())
            .field("Custom", &self.custom.len())
            .finish()
    }
}

impl Ontologies {
    /// Initialize all ontologies, also returns all warnings detected when initialising the ontologies
    pub fn init() -> (Self, Vec<BoxedError<'static, mzcv::CVError>>) {
        let mut errors = Vec::new();
        let (unimod, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);
        let (psimod, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);
        let (xlmod, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);
        let (gnome, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);
        let (resid, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);
        let (custom, mut warnings) = CVIndex::init();
        errors.append(&mut warnings);

        (
            Self {
                custom,
                gnome,
                psimod,
                resid,
                unimod,
                xlmod,
            },
            errors,
        )
    }

    /// Initialize all ontologies with their static data (or empty in the case of Custom)
    pub fn init_static() -> Self {
        Self {
            custom: CVIndex::init_static(),
            gnome: CVIndex::init_static(),
            psimod: CVIndex::init_static(),
            resid: CVIndex::init_static(),
            unimod: CVIndex::init_static(),
            xlmod: CVIndex::init_static(),
        }
    }

    /// Initialize all ontologies with their empty indices
    pub fn empty() -> Self {
        Self {
            custom: CVIndex::empty(),
            gnome: CVIndex::empty(),
            psimod: CVIndex::empty(),
            resid: CVIndex::empty(),
            unimod: CVIndex::empty(),
            xlmod: CVIndex::empty(),
        }
    }

    /// Update the custom database with the given data, mostly there to quickly create some test data.
    /// This does not store the new data on the disk
    #[must_use]
    pub fn with_custom(mut self, data: impl IntoIterator<Item = SimpleModification>) -> Self {
        self.custom
            .update_do_not_save_to_disk(CVVersion::default(), data);
        self
    }

    /// Get Unimod
    pub const fn unimod(&self) -> &CVIndex<Unimod> {
        &self.unimod
    }

    /// Get Unimod
    pub const fn unimod_mut(&mut self) -> &mut CVIndex<Unimod> {
        &mut self.unimod
    }

    /// Get PSI-MOD
    pub const fn psimod(&self) -> &CVIndex<PsiMod> {
        &self.psimod
    }

    /// Get PSI-MOD
    pub const fn psimod_mut(&mut self) -> &mut CVIndex<PsiMod> {
        &mut self.psimod
    }

    /// Get XL-MOD
    pub const fn xlmod(&self) -> &CVIndex<XlMod> {
        &self.xlmod
    }

    /// Get XL-MOD
    pub const fn xlmod_mut(&mut self) -> &mut CVIndex<XlMod> {
        &mut self.xlmod
    }

    /// Get GNOme
    pub const fn gnome(&self) -> &CVIndex<Gnome> {
        &self.gnome
    }

    /// Get GNOme
    pub const fn gnome_mut(&mut self) -> &mut CVIndex<Gnome> {
        &mut self.gnome
    }

    /// Get RESID
    pub const fn resid(&self) -> &CVIndex<Resid> {
        &self.resid
    }

    /// Get RESID
    pub const fn resid_mut(&mut self) -> &mut CVIndex<Resid> {
        &mut self.resid
    }

    /// Get Custom
    pub const fn custom(&self) -> &CVIndex<Custom> {
        &self.custom
    }

    /// Get Custom
    pub const fn custom_mut(&mut self) -> &mut CVIndex<Custom> {
        &mut self.custom
    }

    /// Find the closest names in the given ontologies, or if empty in all ontologies
    pub fn search(&self, ontologies: &[Ontology], term: &str) -> Vec<(SimpleModification, bool)> {
        let ontologies = if ontologies.is_empty() {
            &[
                Ontology::Unimod,
                Ontology::Psimod,
                Ontology::Xlmod,
                Ontology::Gnome,
                Ontology::Resid,
                Ontology::Custom,
            ]
        } else {
            ontologies
        };

        let mut options = Vec::new();
        for ontology in ontologies {
            options.append(&mut match ontology {
                Ontology::Unimod => self.unimod.search(term, 5, 6),
                Ontology::Psimod => self.psimod.search(term, 5, 6),
                Ontology::Xlmod => self.xlmod.search(term, 5, 6),
                Ontology::Gnome => self.gnome.search(term, 5, 6),
                Ontology::Resid => self.resid.search(term, 5, 6),
                Ontology::Custom => self.custom.search(term, 5, 6),
            });
        }

        options.sort_unstable_by(|a, b| a.2.cmp(&b.2).then(a.1.cmp(&b.1)).then(a.0.cmp(&b.0)));

        options
            .into_iter()
            .map(|(m, n, _)| (m, n))
            .take(10)
            .collect()
    }

    /// Find the given name in this ontology.
    pub fn get_by_name(&self, ontologies: &[Ontology], term: &str) -> Option<SimpleModification> {
        let ontologies = if ontologies.is_empty() {
            &[
                Ontology::Unimod,
                Ontology::Psimod,
                Ontology::Xlmod,
                Ontology::Gnome,
                Ontology::Resid,
                Ontology::Custom,
            ]
        } else {
            ontologies
        };

        for ontology in ontologies {
            match ontology {
                Ontology::Unimod => {
                    if let Some(m) = self.unimod.get_by_name(term) {
                        return Some(m);
                    }
                }
                Ontology::Psimod => {
                    if let Some(m) = self.psimod.get_by_name(term) {
                        return Some(m);
                    }
                }
                Ontology::Xlmod => {
                    if let Some(m) = self.xlmod.get_by_name(term) {
                        return Some(m);
                    }
                }
                Ontology::Gnome => {
                    if let Some(m) = self.gnome.get_by_name(term) {
                        return Some(m);
                    }
                }
                Ontology::Resid => {
                    if let Some(m) = self.resid.get_by_name(term) {
                        return Some(m);
                    }
                }
                Ontology::Custom => {
                    if let Some(m) = self.custom.get_by_name(term) {
                        return Some(m);
                    }
                }
            }
        }

        None
    }

    /// Get all modifications in the selected ontologies (or in all if the list is empty).
    pub fn data(&self, ontologies: &[Ontology]) -> impl Iterator<Item = SimpleModification> {
        let ontologies = if ontologies.is_empty() {
            &[
                Ontology::Unimod,
                Ontology::Psimod,
                Ontology::Xlmod,
                Ontology::Gnome,
                Ontology::Resid,
                Ontology::Custom,
            ]
        } else {
            ontologies
        };

        ontologies.iter().flat_map(|ontology| match ontology {
            Ontology::Unimod => self.unimod().data(),
            Ontology::Psimod => self.psimod().data(),
            Ontology::Xlmod => self.xlmod().data(),
            Ontology::Gnome => self.gnome().data(),
            Ontology::Resid => self.resid().data(),
            Ontology::Custom => self.custom().data(),
        })
    }
}

/// All allowed ontologies for modification names
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum Ontology {
    #[default]
    /// Unimod
    Unimod,
    /// PSI-MOD
    Psimod,
    /// GNOme
    Gnome,
    /// XLMOD
    Xlmod,
    /// Resid
    Resid,
    /// Custom
    Custom,
}

impl Ontology {
    /// Get the prefix character for the ontology
    pub const fn char(self) -> char {
        match self {
            Self::Unimod => 'U',
            Self::Psimod => 'M',
            Self::Gnome => 'G',
            Self::Xlmod => 'X',
            Self::Resid => 'R',
            Self::Custom => 'C',
        }
    }

    /// Get the accession number name for the ontology
    pub const fn name(self) -> &'static str {
        match self {
            Self::Unimod => "UNIMOD",
            Self::Psimod => "MOD",
            Self::Gnome => "GNO",
            Self::Xlmod => "XLMOD",
            Self::Resid => "RESID",
            Self::Custom => "CUSTOM",
        }
    }
}

impl std::fmt::Display for Ontology {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Unimod => "Unimod",
                Self::Psimod => "PSI-MOD",
                Self::Gnome => "GNOme",
                Self::Xlmod => "XLMOD",
                Self::Resid => "Resid",
                Self::Custom => "Custom",
            },
        )
    }
}
