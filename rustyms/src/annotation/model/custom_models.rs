use std::path::Path;

use crate::{
    error::{Context, CustomError},
    parse_json::ParseJson,
    prelude::FragmentationModel,
};

/// Parse a custom models JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom models in previous version of
/// the library.
/// # Errors
/// If the file could not be opened or parsed.
pub fn parse_custom_models(path: &Path) -> Result<Vec<(String, FragmentationModel)>, CustomError> {
    let string = std::fs::read_to_string(path).map_err(|err| {
        CustomError::error(
            "Could not parse custom models file",
            err,
            Context::show(path.display()),
        )
    })?;
    Vec::from_json(&string)
}

/// Parse a custom models JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom models in previous version of
/// the library.
/// # Errors
/// If the string could not be parsed.
pub fn parse_custom_models_str(
    value: &str,
) -> Result<Vec<(String, FragmentationModel)>, CustomError> {
    Vec::from_json(value)
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_reading_custom_models_json_2024() {
    let data = include_str!("custom_model_2024.json");
    let mods = parse_custom_models_str(data).unwrap();
    assert!(mods.len() > 1);
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_reading_custom_models_json_2025() {
    let data = include_str!("custom_model_20250528.json");
    let mods = parse_custom_models_str(data).unwrap();
    assert!(mods.len() > 1);
}
