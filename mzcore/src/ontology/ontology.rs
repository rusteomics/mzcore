//! The available ontologies

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::sequence::SimpleModification;

/// A database of custom modifications
pub type CustomDatabase = OntologyModificationList;

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
#[cfg(not(feature = "internal-no-data"))]
mod databases {
    use crate::ontology::OntologyModificationList;
    use bincode::config::Configuration;
    use std::sync::LazyLock;

    /// Get the unimod ontology
    /// # Panics
    /// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
    pub(super) static UNIMOD: LazyLock<OntologyModificationList> = LazyLock::new(|| {
        bincode::serde::decode_from_slice::<OntologyModificationList, Configuration>(
            include_bytes!("../databases/unimod.dat"),
            Configuration::default(),
        )
        .unwrap()
        .0
    });
    /// Get the PSI-MOD ontology
    /// # Panics
    /// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
    pub(super) static PSIMOD: LazyLock<OntologyModificationList> = LazyLock::new(|| {
        bincode::serde::decode_from_slice::<OntologyModificationList, Configuration>(
            include_bytes!("../databases/psimod.dat"),
            Configuration::default(),
        )
        .unwrap()
        .0
    });
    /// Get the Gnome ontology
    /// # Panics
    /// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
    pub(super) static GNOME: LazyLock<OntologyModificationList> = LazyLock::new(|| {
        bincode::serde::decode_from_slice::<OntologyModificationList, Configuration>(
            include_bytes!("../databases/gnome.dat"),
            Configuration::default(),
        )
        .unwrap()
        .0
    });
    /// Get the Resid ontology
    /// # Panics
    /// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
    pub(super) static RESID: LazyLock<OntologyModificationList> = LazyLock::new(|| {
        bincode::serde::decode_from_slice::<OntologyModificationList, Configuration>(
            include_bytes!("../databases/resid.dat"),
            Configuration::default(),
        )
        .unwrap()
        .0
    });
    /// Get the Xlmod ontology
    /// # Panics
    /// Panics when the modifications are not correctly provided at compile time, always report a panic if it occurs here.
    pub(super) static XLMOD: LazyLock<OntologyModificationList> = LazyLock::new(|| {
        bincode::serde::decode_from_slice::<OntologyModificationList, Configuration>(
            include_bytes!("../databases/xlmod.dat"),
            Configuration::default(),
        )
        .unwrap()
        .0
    });
}
