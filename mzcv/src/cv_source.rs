//! The main traits any user of the CV logic needs to implement.

use bincode::{Decode, Encode};
use context_error::BoxedError;
use directories::{BaseDirs, ProjectDirs};
use sha2::Digest;

use crate::{CVError, Curie, hash_buf_reader::HashBufReader};

/// Implement this trait to create a new CV. The best way of using this is with a ZST (zero sized type).
/// ```rust ignore
/// struct Unimod {}
///
/// impl CVSource for Unimod {
///     ...
/// }
///
/// fn main() {
///     Unimod::get_by_name("deamidated");
/// }
/// ```
pub trait CVSource {
    /// Set this constant to true to enable automatic writing of the CV when the cahe is updated
    const AUTOMATICALLY_WRITE_UNCOMPRESSED: bool = false;
    /// The data item that is stored in the CV
    type Data: CVData + 'static;
    /// The type of the main datastructure to keep all data items (used to build any kind of hierarchy necessary)
    type Structure: CVStructure<Self::Data> + Encode + Decode<()>;
    /// The name of the CV, used to create the paths to store intermediate files and caches so has to be valid in that context
    fn cv_name() -> &'static str;
    /// The source files for the
    fn files() -> &'static [CVFile];
    /// The static data of this CV
    fn static_data() -> Option<(CVVersion, Self::Structure)>;
    /// The default file stem (no extension).
    /// # Panics
    /// If both [`ProjectDirs::from`] and [`BaseDirs::new`] failed to retrieve a suitable base folder.
    fn default_stem() -> std::path::PathBuf {
        let proj_dirs = ProjectDirs::from("org", "rusteomics", "mzcore");
        let base = BaseDirs::new();
        let folder = proj_dirs
            .map(|p| p.data_dir().to_owned())
            .or_else(|| base.map(|b| b.home_dir().to_owned()))
            .expect("No suitable base directory could be found");
        if !folder.exists() {
            drop(std::fs::create_dir_all(&folder));
        }
        folder.join(Self::cv_name())
    }
    /// Parse the textual representation of this CV. Return the version, the data, and possibly a
    /// list of warnings or non critical errors encountered while parsing.
    /// # Errors
    /// If the parsing failed.
    fn parse(
        reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl Digest>>,
    ) -> Result<
        (
            CVVersion,
            Self::Structure,
            Vec<BoxedError<'static, CVError>>,
        ),
        Vec<BoxedError<'static, CVError>>,
    >;

    /// Write out the data to the standard file, only need to implement this if
    /// `AUTOMATICALLY_WRITE_UNCOMPRESSED` is set to true. This is used to write out data when the
    /// main way of interacting with this structure is via the [`crate::CVIndex::update`] method.
    /// So for custom modification CVs for example.
    /// # Errors
    /// When the underlying writer failed.
    fn write_uncompressed<W: std::io::Write>(
        _writer: W,
        _version: &CVVersion,
        _data: impl Iterator<Item = std::sync::Arc<Self::Data>>,
    ) -> Result<(), BoxedError<'static, CVError>> {
        Ok(())
    }
}

/// The description of a file that is used to built a controlled vocabulary.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct CVFile {
    /// The name of the CV, used to create the paths to store intermediate files and caches so has to be valid in that context
    pub name: &'static str,
    /// The file extension of the CV
    pub extension: &'static str,
    /// The url where this CV can be found (note that only http/https downloading is supported)
    pub url: Option<&'static str>,
    /// The source compression of this file (note that only uncompressed is currently supported)
    pub compression: CVCompression,
}

/// Version information for a CV
#[derive(Clone, Debug, Decode, Default, Encode)]
pub struct CVVersion {
    /// The last updated date as reported by the CV (year, month, day, hour, minute)
    pub last_updated: Option<(u16, u8, u8, u8, u8)>,
    /// The version of the CV
    pub version: Option<String>,
    /// The hash of the uncompressed file
    pub hash: Vec<u8>,
}

impl CVVersion {
    /// A nice representation of the date
    pub fn last_updated(&self) -> Option<String> {
        self.last_updated.map(|(year, month, day, hour, minute)| {
            format!("{year:04}-{month:02}-{day:02} {hour:02}:{minute:02}")
        })
    }

    /// A nice string representation of the hash
    pub fn hash_hex(&self) -> String {
        use std::fmt::Write;
        let mut res = String::with_capacity(self.hash.len() * 2);
        for n in &self.hash {
            write!(&mut res, "{n:02x}").unwrap();
        }
        res
    }
}

/// A data element from a CV. Note that technically an implementation could be made that does not
/// have any index, name, or keywords for a data element. Such an element would be kept and stored
/// but would not be accessible via anything else then [`crate::CVIndex::data`].
pub trait CVData: Clone {
    /// The type used for the index
    type Index: std::fmt::Debug + Clone + std::hash::Hash + Eq;
    /// The numerical index or id of the item.
    fn index(&self) -> Option<Self::Index>;

    /// The CURIE for this entry. This behavior is strictly *opt-in*.
    fn curie(&self) -> Option<Curie> {
        None
    }

    /// The name of the item, this will be stored in a case insensitive manner in the
    /// [`crate::CVIndex`] but should be reported in the correct casing here.
    fn name(&self) -> Option<std::borrow::Cow<'_, str>>;

    /// Any synonyms that can be uniquely attributed to this data element. So EXACT synonyms from Obo files.
    fn synonyms(&self) -> impl Iterator<Item = &str>;

    /// Enumerate all the immediate parents of this term, e.g. those with an `is_a` relationship with it. This behavior is strictly *opt-in*.
    fn parents(&self) -> impl Iterator<Item = &Self::Index> {
        let value: Option<&Self::Index> = None;
        value.into_iter()
    }
}

/// The used compression of the source CV.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub enum CVCompression {
    /// No compression
    #[default]
    None,
    /// For LZW (.Z) compressed files like IMGT
    Lzw,
}

/// A structure to contain [`CVData`] elements but leave the implementation up to the needs of the specific CV.
pub trait CVStructure<Data>: Default {
    /// See if no data items are present
    fn is_empty(&self) -> bool;
    /// The total number of data items present
    fn len(&self) -> usize;
    /// Clear this entire structure
    fn clear(&mut self);
    /// The iterator to use when the index and data items are requested
    type IterIndexed<'a>: Iterator<Item = (Self::Index, std::sync::Arc<Data>)>
    where
        Self: 'a;
    /// Iterate over all data items and return both the index and the data item
    fn iter_indexed(&self) -> Self::IterIndexed<'_>;
    /// The iterator to use when only the data items are requested
    type IterData<'a>: Iterator<Item = std::sync::Arc<Data>>
    where
        Self: 'a;
    /// Iterate over all data items
    fn iter_data(&self) -> Self::IterData<'_>;
    /// Add a sinlge data item to the structure
    fn add(&mut self, data: std::sync::Arc<Data>);
    /// The indexing type
    type Index: Clone;
    /// Get an element by the index
    fn index(&self, index: Self::Index) -> Option<std::sync::Arc<Data>>;
    /// Remove an element based on the index
    fn remove(&mut self, index: Self::Index);
}

impl<T> CVStructure<T> for Vec<std::sync::Arc<T>> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }
    fn len(&self) -> usize {
        self.len()
    }
    fn clear(&mut self) {
        self.clear();
    }
    type IterIndexed<'a>
        = std::iter::Enumerate<std::iter::Cloned<std::slice::Iter<'a, std::sync::Arc<T>>>>
    where
        T: 'a;
    fn iter_indexed(&self) -> Self::IterIndexed<'_> {
        self.iter().cloned().enumerate()
    }
    type IterData<'a>
        = std::iter::Cloned<std::slice::Iter<'a, std::sync::Arc<T>>>
    where
        T: 'a;
    fn iter_data(&self) -> Self::IterData<'_> {
        self.iter().cloned()
    }
    fn add(&mut self, data: std::sync::Arc<T>) {
        self.push(data);
    }
    type Index = usize;
    fn index(&self, index: Self::Index) -> Option<std::sync::Arc<T>> {
        self.get(index).cloned()
    }
    fn remove(&mut self, index: Self::Index) {
        self.remove(index);
    }
}
