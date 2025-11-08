//! Code to handle the Custom ontology
use context_error::{BoxedError, Context, CreateError, FullErrorContent};

use mzcv::{CVError, CVFile, CVSource, CVVersion, HashBufReader};

use crate::{
    parse_json::ParseJson,
    sequence::{SimpleModification, SimpleModificationInner},
};

/// The shared type for contact between the build and compile steps
pub type CustomDatabase = Vec<(Option<usize>, String, SimpleModification)>;

/// Custom modifications
///
/// These can be used to handle `C:modification` and `Custom:0001` modifications in ProForma
/// definitions. As these are custom made by the user there is no way of automatically updating
/// from the internet. To parse these files the `custom_modifications.json` file is parsed with
/// [`CustomDatabase::from_json_value`] which generates the modifications.
///
/// To update the modifications use [`mzcv::CVIndex::update`] or [`mzcv::CVIndex::update_from_path`].
/// Both automatically write the cache and the original JSON file.
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct Custom {}

impl CVSource for Custom {
    const AUTOMATICALLY_WRITE_UNCOMPRESSED: bool = true;

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
        let mut reader = reader.next().unwrap();

        let json: serde_json::Value = serde_json::de::from_reader(&mut reader).map_err(|err| {
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
                hash: reader.hash().to_vec(),
                last_updated: None,
                version: None,
            },
            db.into_iter()
                .map(|(_, _, s)| std::sync::Arc::unwrap_or_clone(s)),
        ))
    }

    fn write_uncompressed<W: std::io::Write>(
        writer: W,
        _version: &CVVersion,
        data: impl Iterator<Item = std::sync::Arc<Self::Data>>,
    ) -> Result<(), BoxedError<'static, CVError>> {
        let value = data
            .map(|m| {
                (
                    m.description().and_then(|i| i.id()),
                    m.description()
                        .map_or_else(|| "UNNAMED".to_string(), |i| i.name.to_string()),
                    m,
                )
            })
            .collect::<Vec<_>>();
        serde_json::to_writer_pretty(writer, &value).map_err(|e| {
            BoxedError::new(
                CVError::FileCouldNotBeMade,
                "Could not save custom modifications file",
                e.to_string(),
                Context::default(),
            )
        })
    }
}
