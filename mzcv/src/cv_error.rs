//! The [`CVError`] which makes it easy for downstream users of the error type to match on the exact error.

use context_error::ErrorKind;

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum CVError {
    CacheDoesNotExist,
    CacheCouldNotBeOpenend,
    CacheCouldNotBeMade,
    CacheCouldNotBeParsed,
    #[default]
    FileDoesNotExist,
    FileCouldNotBeMade,
    FileCouldNotBeOpenend,
    FileCouldNotBeMoved,
    FileCouldNotBeParsed,
    CVUrlNotSet,
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
