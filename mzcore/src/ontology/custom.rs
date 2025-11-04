//! Code to handle the Custom ontology
use context_error::{BoxedError, Context, CreateError, FullErrorContent};

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader};

use crate::{
    parse_json::ParseJson,
    sequence::{SimpleModification, SimpleModificationInner},
};

/// The shared type for contact between the build and compile steps
pub type CustomDatabase = Vec<(Option<usize>, String, SimpleModification)>;

pub struct Custom {}

impl CVSource for Custom {
    type Data = SimpleModificationInner;
    fn cv_name() -> &'static str {
        "CUSTOM"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "custom_modifications",
            extension: "json",
            url: None,
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Vec<Self::Data>)> {
        None
    }

    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
        let reader = reader.next().unwrap();

        let json: serde_json::Value = serde_json::de::from_reader(reader).map_err(|err| {
            vec![BoxedError::new(
                CVError::FileCouldNotBeParsed,
                "Could not parse custom modifications file",
                err.to_string(),
                Context::default(),
            )]
        })?;

        let db = CustomDatabase::from_json_value(json)
            .map_err(|err| vec![err.convert(|_| CVError::FileCouldNotBeParsed)])?;

        Ok((
            CVVersion {
                hash: Vec::new(),
                last_updated: None,
                version: None,
            },
            db.into_iter()
                .map(|(_, _, s)| std::sync::Arc::unwrap_or_clone(s)),
        ))
    }
}
