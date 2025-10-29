//! The main traits any user of the CV logic needs to implement.

use std::sync::Arc;

use bincode::{Decode, Encode};
use context_error::BoxedError;
use directories::{BaseDirs, ProjectDirs};
use sha2::Digest;

use crate::{CVError, hash_buf_reader::HashBufReader};

/// Implement this trait to create a new CV. The best way of using this might be with a ZST (zero sized type).
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
    /// The data item that is stored in the CV
    type Data: CVData + 'static;
    /// The name of the CV, used to create the paths to store intermediate files and caches so has to be valid in that context
    fn cv_name() -> &'static str;
    /// The source files for the
    fn files() -> &'static [CVFile];
    /// The static data of this CV
    fn static_data() -> Option<(CVVersion, &'static [Arc<Self::Data>])>;
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
    /// Parse the textual representation of this CV
    /// # Errors
    /// If the parsing failed.
    fn parse(
        reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>;
}

pub struct CVFile {
    /// The name of the CV, used to create the paths to store intermediate files and caches so has to be valid in that context
    pub name: &'static str,
    /// The file extension of the CV
    pub extension: &'static str,
    /// The url where this CV can be found
    pub url: Option<&'static str>,
    /// The source compression of this file
    pub compression: CVCompression,
}

/// Version information for a CV
#[derive(Debug, Decode, Default, Encode)]
pub struct CVVersion {
    /// The last updated date as reported by the CV (year, month, day, hour, minute)
    pub last_updated: Option<(u16, u8, u8, u8, u8)>,
    /// The version of the CV
    pub version: Option<String>,
    /// The hash of the uncompressed file
    pub hash: Vec<u8>,
}

impl CVVersion {
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
pub trait CVData: Encode + Decode<()> + Clone {
    /// The type used for the index
    // Example IDs:
    // NCIT:C25330
    // NCIT:R100
    // NCIT:P383
    // MS:1000014
    // UO:0000245
    // BAO_0000925
    // BFO:0000015
    // GAZ:00002933
    // AfPO_0000233 // In URL
    // BTO:0004947
    // PRIDE:0000521
    // MOD:01188
    // XLMOD:07097
    // GNO:G00001NT // Name and ID are identical
    // GNO:00000015
    type Index: std::fmt::Debug + Clone + std::hash::Hash + Eq;
    /// The numerical index or id of the item.
    fn index(&self) -> Option<Self::Index>;
    /// The name of the item, this will be stored in a case insensitive manner in the
    /// [`crate::CVIndex`] but should be reported in the correct casing here.
    fn name(&self) -> Option<&str>;
    /// Any synonyms that can be uniquely attributed to this data element.
    fn synonyms(&self) -> impl Iterator<Item = &str>;
}

/// The used compression of the source CV.
#[derive(Clone, Copy, Debug, Default)]
pub enum CVCompression {
    /// No compression
    #[default]
    None,
    /// For LZW (.Z) compressed files like IMGT TODO: does not work yet
    LZW,
}
