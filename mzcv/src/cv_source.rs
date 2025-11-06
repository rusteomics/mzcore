//! The main traits any user of the CV logic needs to implement.

use bincode::{Decode, Encode};
use context_error::BoxedError;
use directories::{BaseDirs, ProjectDirs};
use sha2::Digest;

use crate::{CVError, hash_buf_reader::HashBufReader};

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
    /// The name of the CV, used to create the paths to store intermediate files and caches so has to be valid in that context
    fn cv_name() -> &'static str;
    /// The source files for the
    fn files() -> &'static [CVFile];
    /// The static data of this CV
    fn static_data() -> Option<(CVVersion, Vec<Self::Data>)>;
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
    /// The name of the item, this will be stored in a case insensitive manner in the
    /// [`crate::CVIndex`] but should be reported in the correct casing here.
    fn name(&self) -> Option<&str>;
    /// Any synonyms that can be uniquely attributed to this data element.
    fn synonyms(&self) -> impl Iterator<Item = &str>;

    /// The cache type, needed to be generic to handle serde (de)serialisation nicely
    type Cache: CVCache<Self>;
}

/// A trait to help setting the encoding/decoding for the [`CVData`]
pub trait CVCache<Data>: Encode + Decode<()> {
    /// Construct a cache
    fn construct(version: CVVersion, data: Vec<Data>) -> Self;
    /// Deconstruct a cache
    fn deconstruct(self) -> (CVVersion, Vec<Data>);
}

/// A cache using [`bincode`]
#[derive(Debug, Decode, Encode)]
pub struct CVCacheBincode<D: Decode<()> + Encode> {
    version: CVVersion,
    data: Vec<D>,
}

impl<T: Encode + Decode<()>> CVCache<T> for CVCacheBincode<T> {
    fn construct(version: CVVersion, data: Vec<T>) -> Self {
        Self { version, data }
    }
    fn deconstruct(self) -> (CVVersion, Vec<T>) {
        (self.version, self.data)
    }
}

/// A cache using [`serde`]
#[cfg(feature = "serde")]
#[derive(Debug)]
pub struct CVCacheSerde<D: serde::Serialize + for<'de> serde::Deserialize<'de>> {
    version: CVVersion,
    data: Vec<D>,
}

impl<D: serde::Serialize + for<'de> serde::Deserialize<'de>> Encode for CVCacheSerde<D> {
    fn encode<E: bincode::enc::Encoder>(
        &self,
        encoder: &mut E,
    ) -> Result<(), bincode::error::EncodeError> {
        Encode::encode(&self.version, encoder)?;
        Encode::encode(
            &bincode::serde::encode_to_vec(&self.data, *encoder.config())?,
            encoder,
        )?;
        Ok(())
    }
}

impl<Context, T: serde::Serialize + for<'de> serde::Deserialize<'de>> Decode<Context>
    for CVCacheSerde<T>
{
    fn decode<D: bincode::de::Decoder<Context = Context>>(
        decoder: &mut D,
    ) -> Result<Self, bincode::error::DecodeError> {
        let version = Decode::decode(decoder)?;
        let data: Vec<u8> = Decode::decode(decoder)?;
        let (data, _) = bincode::serde::decode_from_slice(&data, *decoder.config())?;
        Ok(Self { version, data })
    }
}

#[cfg(feature = "serde")]
impl<T: serde::Serialize + for<'de> serde::Deserialize<'de>> CVCache<T> for CVCacheSerde<T> {
    fn construct(version: CVVersion, data: Vec<T>) -> Self {
        Self { version, data }
    }
    fn deconstruct(self) -> (CVVersion, Vec<T>) {
        (self.version, self.data)
    }
}

/// The used compression of the source CV.
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub enum CVCompression {
    /// No compression
    #[default]
    None,
    /// For LZW (.Z) compressed files like IMGT TODO: does not work yet
    LZW,
}
