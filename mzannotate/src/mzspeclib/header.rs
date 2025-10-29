//! Library headers and the attribute sets they can contain
use std::{collections::HashMap, fmt::Display};

use crate::mzspeclib::{Attribute, Attributes};

/// The header for a spectral library
#[derive(Clone, Debug)]
pub struct LibraryHeader {
    /// The version of the format
    pub format_version: String,
    /// The attributes for this library
    pub attributes: Attributes,
    /// The attribute classes for this library
    pub attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
}

impl Default for LibraryHeader {
    fn default() -> Self {
        Self {
            format_version: "1.0".into(),
            attributes: vec![Vec::new(); 1],
            attribute_classes: HashMap::new(),
        }
    }
}

impl LibraryHeader {
    /// Create a new library header
    pub const fn new(
        format_version: String,
        attributes: Attributes,
        attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
    ) -> Self {
        Self {
            format_version,
            attributes,
            attribute_classes,
        }
    }

    /// Check if this attribute is already set.
    pub(crate) fn is_already_defined(
        &self,
        attribute: &Attribute,
        entry: EntryType,
        names: &[&str],
    ) -> bool {
        self.attribute_classes
            .get(&entry)
            .iter()
            .flat_map(|s| s.iter())
            .filter(|s| names.contains(&s.id.as_str()))
            .any(|s| {
                s.attributes
                    .iter()
                    .flatten()
                    .any(|a| a.name == attribute.name && a.value.equivalent(&attribute.value))
            })
    }
}

/// A set of attributes in the library header that contains global settings for all or a subset of spectra
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct AttributeSet {
    /// The id
    pub id: String,
    /// The type of attributes
    pub namespace: EntryType,
    /// The attributes themselves
    pub attributes: Attributes,
}

impl PartialOrd for AttributeSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.id.partial_cmp(&other.id) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.namespace.partial_cmp(&other.namespace) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.attributes.len().partial_cmp(&other.attributes.len())
    }
}

/// All types of entries for global attribute sets
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum EntryType {
    /// Spectrum attributes
    Spectrum,
    /// Analyte attributes
    Analyte,
    /// Interpretation attributes
    Interpretation,
    /// A cluster
    Cluster,
}

impl Display for EntryType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Spectrum => "Spectrum",
                Self::Analyte => "Analyte",
                Self::Interpretation => "Interpretation",
                Self::Cluster => "Cluster",
            }
        )
    }
}
