//! The [`CVError`] which makes it easy for downstream users of the error type to match on the exact error.

use context_error::ErrorKind;

/// An error encountered while parsing a CV
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum CVError {
    /// The binary cache does not exist
    CacheDoesNotExist,
    /// The binary cache could not be opened
    CacheCouldNotBeOpenend,
    /// The binary cache could not be made
    CacheCouldNotBeMade,
    /// The binary cache could not be parsed (likely the result of different versions of the same CV with different struct definitions)
    CacheCouldNotBeParsed,
    /// The CV file does not exist
    #[default]
    FileDoesNotExist,
    /// The CV file could not be made
    FileCouldNotBeMade,
    /// The CV file could not be opened
    FileCouldNotBeOpenend,
    /// The CV file could not be moved
    FileCouldNotBeMoved,
    /// The CV file could not be parsed
    FileCouldNotBeParsed,
    /// The CV does not have a URL set and neither was a URL supplied
    CVUrlNotSet,
    /// The CV could not be read from the given URL, the URL could be wrong, or not up, or fail to download the file
    CVUrlCouldNotBeRead,
}

impl ErrorKind for CVError {
    type Settings = ();
    fn descriptor(&self) -> &'static str {
        "error"
    }
    fn ignored(&self, _settings: Self::Settings) -> bool {
        false
    }
    fn is_error(&self, _settings: Self::Settings) -> bool {
        true
    }
}
