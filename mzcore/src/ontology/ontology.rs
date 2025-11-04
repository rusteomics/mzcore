//! The available ontologies

use context_error::*;
use itertools::Itertools;
use mzcv::CVIndex;
use serde::{Deserialize, Serialize};

use crate::{
    ontology::{Custom, CustomDatabase, Gnome, PsiMod, Resid, Unimod, XlMod},
    sequence::SimpleModification,
};

pub struct Ontologies {
    custom: CVIndex<Custom>,
    gnome: CVIndex<Gnome>,
    psimod: CVIndex<PsiMod>,
    resid: CVIndex<Resid>,
    unimod: CVIndex<Unimod>,
    xlmod: CVIndex<XlMod>,
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

    /// Find the closest names in the given ontologies, or if empty in all ontologies
    pub fn find_closest(&self, ontologies: &[Ontology], term: &str) -> Vec<SimpleModification> {
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

        options
    }

    /// Find the given name in this ontology.
    pub fn find_name(&self, ontologies: &[Ontology], term: &str) -> Option<SimpleModification> {
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

    /// Find the given name in this ontology.
    pub fn find_id(&self, ontologies: &[Ontology], id: usize) -> Option<SimpleModification> {
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
                    if let Some(m) = self.unimod.get_by_index(&id) {
                        return Some(m);
                    }
                }
                Ontology::Psimod => {
                    if let Some(m) = self.psimod.get_by_index(&id) {
                        return Some(m);
                    }
                }
                Ontology::Xlmod => {
                    if let Some(m) = self.xlmod.get_by_index(&id) {
                        return Some(m);
                    }
                }
                Ontology::Gnome => {
                    if let Some(m) = self.gnome.get_by_index(&id) {
                        return Some(m);
                    }
                }
                Ontology::Resid => {
                    if let Some(m) = self.resid.get_by_index(&id) {
                        return Some(m);
                    }
                }
                Ontology::Custom => {
                    if let Some(m) = self.custom.get_by_index(&id) {
                        return Some(m);
                    }
                }
            }
        }

        None
    }
}

/// An empty list of modifications (needed for lifetime reasons)
static EMPTY_LIST: OntologyModificationList = Vec::new();

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

/// The shared type for contact between the build and compile steps
pub type OntologyModificationList = Vec<(Option<usize>, String, SimpleModification)>;

impl Ontology {
    /// Get the modifications lookup list for this ontology
    #[allow(unused_variables, clippy::unused_self)]
    pub fn lookup(self, custom_database: Option<&CustomDatabase>) -> &OntologyModificationList {
        #[cfg(not(feature = "internal-no-data"))]
        {
            match self {
                Self::Gnome => &databases::GNOME,
                Self::Psimod => &databases::PSIMOD,
                Self::Unimod => &databases::UNIMOD,
                Self::Resid => &databases::RESID,
                Self::Xlmod => &databases::XLMOD,
                Self::Custom => custom_database.map_or(&EMPTY_LIST, |c| c),
            }
        }
        #[cfg(feature = "internal-no-data")]
        {
            &EMPTY_LIST
        }
    }

    /// Find the closest names in this ontology
    pub fn find_closest<'a>(
        self,
        code: &'a str,
        custom_database: Option<&CustomDatabase>,
    ) -> BoxedError<'a, BasicKind> {
        BoxedError::new(
            BasicKind::Error,
            "Invalid modification",
            format!("The provided name does not exists in {}", self.name()),
            Context::show(code),
        )
        .suggestions(Self::similar_names(&[self], code, custom_database))
    }

    /// # Panics
    /// Asserts that the ontology list is not empty.
    pub fn find_closest_many(
        ontologies: &[Self],
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> BoxedError<'static, BasicKind> {
        assert!(!ontologies.is_empty());
        let names = if ontologies.len() > 1 {
            let mut names = ontologies[..ontologies.len() - 1]
                .iter()
                .map(|o| o.name().to_string())
                .collect_vec();
            let last = names.len() - 1;
            names[last] = format!("{} or {}", names[last], ontologies.last().unwrap().name());
            names.join(", ")
        } else {
            ontologies[0].name().to_string()
        };
        BoxedError::new(
            BasicKind::Error,
            "Invalid modification",
            format!("The provided name does not exists in {names}"),
            Context::show(code.to_string()),
        )
        .suggestions(Self::similar_names(ontologies, code, custom_database))
    }

    /// Get the closest similar names in the given ontologies. Finds both modifications and linkers
    pub fn similar_names(
        ontologies: &[Self],
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Vec<String> {
        let mut resulting = Vec::new();
        for ontology in ontologies {
            let options: Vec<&str> = ontology
                .lookup(custom_database)
                .iter()
                .map(|option| option.1.as_str())
                .collect();
            resulting.extend(
                similar::get_close_matches(code, &options, 3, 0.7)
                    .iter()
                    .map(|o| format!("{}:{}", ontology.char(), o)),
            );
        }
        resulting
    }

    /// Find the given name in this ontology.
    pub fn find_name(
        self,
        code: &str,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        for option in self.lookup(custom_database) {
            if option.1.eq_ignore_ascii_case(code) {
                return Some(option.2.clone());
            }
        }
        None
    }

    /// Find the given id in this ontology
    pub fn find_id(
        self,
        id: usize,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        for option in self.lookup(custom_database) {
            if option.0.is_some_and(|i| i == id) {
                return Some(option.2.clone());
            }
        }
        None
    }
}
