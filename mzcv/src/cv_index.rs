//! The [`CVIndex`] itself with the core functionality.

use std::{collections::HashMap, sync::Arc};

use crate::{CVData, CVSource, CVVersion, text::*};

/// An index into a CV which contains the main ways of handling CVs.
///
/// Data can be accessed using the following main ways:
/// * Unique indexing (constant time) [`Self::get_by_index`] and [`Self::get_by_name`].
/// * Searching for fuzzy matches in names and synonyms using [`Self::search`].
/// * Iterate through all data elements using [`Self::data`].
///
/// Data can be updated in the following ways:
/// * From a direct memory representation [`Self::update`].
/// * From a file [`Self::update_from_path`].
/// * From a URL [`Self::update_from_url`].
///
/// Each of these stores the full data file at the default location
/// ([`CVSource::default_stem`].[`CVSource::cv_extension`].gz). And updates the binary cache that
/// is used as the fastest way of loading the data in [`Self::init`].
#[derive(Debug)]
pub struct CVIndex<CV: CVSource> {
    /// All data elements
    data: Vec<Arc<CV::Data>>,
    /// Index number index
    index: HashMap<<CV::Data as CVData>::Index, Arc<CV::Data>>,
    /// Lowercased name index
    name: HashMap<Box<str>, Arc<CV::Data>>,
    /// Lowercased synonym index, multiple synonyms can link to the same data, but every synonym can only link to one data
    synonyms: HashMap<Box<str>, Arc<CV::Data>>,
    #[cfg(feature = "search-index")]
    /// All trigrams for the names and synonyms
    trigram_index: HashMap<[u8; 3], Vec<Box<str>>>,
    /// The version
    version: CVVersion,
}

impl<CV: CVSource> CVIndex<CV> {
    /// See if this CV does not contain any elements
    pub const fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// Get the number of data elements in this CV
    pub const fn len(&self) -> usize {
        self.data.len()
    }

    /// Load a data item from index
    pub fn get_by_index(&self, index: &<CV::Data as CVData>::Index) -> Option<Arc<CV::Data>> {
        self.index.get(index).cloned()
    }

    /// Load a data item by name, names are matched in a case insensitive manner.
    pub fn get_by_name(&self, name: &str) -> Option<Arc<CV::Data>> {
        let name = name.to_ascii_lowercase().into_boxed_str();
        self.name.get(&name).cloned()
    }

    /// Load a data item by name or if that fails by synonym, names are matched in a case insensitive manner. Returns a boolean indicating if it matches a name `true` or a synonym `false`
    pub fn get_by_name_or_synonym(&self, name: &str) -> Option<(bool, Arc<CV::Data>)> {
        let name = name.to_ascii_lowercase().into_boxed_str();
        self.name
            .get(&name)
            .map(|m| (true, m.clone()))
            .or_else(|| self.synonyms.get(&name).map(|m| (false, m.clone())))
    }

    /// Search trough the name and synonyms lists to find the closest terms.
    /// This can be limited to a maximum amount of entries sent back and to a maximum edit distance with the search term.
    /// The results are sorted unstably on edit distance so multiple calls could return different lists.
    ///
    /// Returns the found mod, if this is based on name `true` or synonym `false` and the edit distance
    pub fn search(
        &self,
        term: &str,
        limit: usize,
        max_distance: usize,
    ) -> Vec<(Arc<CV::Data>, Option<String>, usize)> {
        // Convert to lowercase, see if any name or synonym exactly matches before going over the trigram index and doing distance calculations
        let lowercase = term.to_ascii_lowercase().into_boxed_str();
        self.name
            .get(&lowercase)
            .map(|v| (v.clone(), None, 0))
            .or_else(|| {
                self.synonyms
                    .get(&lowercase)
                    .map(|v| (v.clone(), Some(term.to_string()), 0))
            })
            .map_or_else(
                || {
                    let mut results: Vec<(&str, usize)> = Vec::with_capacity(limit);
                    #[cfg(feature = "search-index")]
                    let mut set = std::collections::HashSet::new();
                    for (distance, t) in {
                        #[cfg(feature = "search-index")]
                        {
                            tags(&lowercase)
                                .filter_map(|tag| self.trigram_index.get(&tag))
                                .flatten()
                                .filter(|term| set.insert(lowercase.clone()))
                        }
                        #[cfg(not(feature = "search-index"))]
                        {
                            self.name.keys().chain(self.synonyms.keys())
                        }
                    }
                    .map(|t| (levenshtein_distance(&lowercase, t), t))
                    .filter(|(distance, _)| *distance <= max_distance)
                    {
                        let index = results
                            .binary_search_by(|item| item.1.cmp(&distance))
                            .unwrap_or_else(|i| i);
                        if index < limit {
                            if results.len() >= limit {
                                results.remove(limit - 1);
                            }
                            results.insert(index, (t, distance));
                        }
                    }

                    results
                        .into_iter()
                        .filter_map(|(name, distance)| {
                            self.name
                                .get(name)
                                .map(|m| (m.clone(), None, distance))
                                .or_else(|| {
                                    self.synonyms
                                        .get(name)
                                        .map(|m| (m.clone(), Some(name.to_string()), distance))
                                })
                        })
                        .collect()
                },
                |v| vec![v],
            )
    }

    /// Get the underlying data in insertion order. Note that the iterator is spelled out fully to
    /// allow merging multiple CVs with the same data type but this should not be relied upon.
    pub fn data(&self) -> std::iter::Cloned<std::slice::Iter<'_, Arc<<CV as CVSource>::Data>>> {
        self.data.iter().cloned()
    }

    /// Get the version
    pub const fn version(&self) -> &CVVersion {
        &self.version
    }

    /// Update without overwriting the cache (used when the data is loaded from the cache anyways)
    pub(crate) fn update_skip_rebuilding_cache(
        &mut self,
        data: impl IntoIterator<Item = Arc<CV::Data>>,
        version: CVVersion,
    ) {
        self.data.clear();
        self.index.clear();
        self.name.clear();
        self.synonyms.clear();
        #[cfg(feature = "search-index")]
        self.trigram_index.clear();
        self.version = version;

        for element in data {
            self.add(element);
        }
    }

    /// Create an empty CV
    pub fn empty() -> Self {
        Self {
            data: Vec::new(),
            index: HashMap::new(),
            name: HashMap::new(),
            synonyms: HashMap::new(),
            #[cfg(feature = "search-index")]
            trigram_index: HashMap::new(),
            version: CVVersion::default(),
        }
    }

    /// Remove a modification. Returns true if this element was sucessfully deletd and false if this element could not be found.
    pub fn remove(&mut self, index: &<CV::Data as CVData>::Index) -> bool {
        if let Some(pos) = self
            .data
            .iter()
            .position(|m| m.index().is_some_and(|id| id == *index))
        {
            let m = self.data[pos].clone();
            if let Some(name) = m.name() {
                self.name.remove(name);
                #[cfg(feature = "search-index")]
                for tag in tags(&name) {
                    self.trigram_index
                        .entry(tag)
                        .and_modify(|v| v.retain(|i| **i != *name));
                }
            }
            for synonym in m.synonyms() {
                self.synonyms.remove(synonym);
                #[cfg(feature = "search-index")]
                for tag in tags(&synonym) {
                    self.trigram_index
                        .entry(tag)
                        .and_modify(|v| v.retain(|i| **i != *synonym));
                }
            }
            self.index.remove(index);
            self.data.remove(pos);
            true
        } else {
            false
        }
    }

    /// Add a single element to the ontology. This updates all indices but does not save it to disk.
    // TODO: what to do on duplicate insertions?
    #[allow(clippy::needless_pass_by_value)] // This fits the use case of update_skip_rebuilding_cache
    pub fn add(&mut self, element: Arc<CV::Data>) {
        self.data.push(element.clone());
        if let Some(index) = element.index() {
            self.index.insert(index, element.clone());
        }
        if let Some(name) = element.name() {
            let name = name.trim_ascii().to_ascii_lowercase().into_boxed_str();
            #[cfg(feature = "search-index")]
            for tag in tags(&name) {
                self.trigram_index
                    .entry(tag)
                    .or_default()
                    .push(name.clone());
            }
            self.name.insert(name, element.clone());
        }
        for keyword in element.synonyms() {
            let keyword = keyword.trim_ascii().to_ascii_lowercase().into_boxed_str();
            #[cfg(feature = "search-index")]
            for tag in tags(&keyword) {
                self.trigram_index
                    .entry(tag)
                    .or_default()
                    .push(keyword.clone());
            }
            self.synonyms.insert(keyword, element.clone());
        }
    }
}
