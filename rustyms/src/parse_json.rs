use std::{any::type_name, ops::RangeInclusive, sync::Arc};

use itertools::Itertools;
use serde::de::DeserializeOwned;
use serde_json::Value;

use crate::error::{Context, CustomError};

/// Custom JSON parser, needed to allow backwards compatibility to older version of JSON formats.
pub(crate) trait ParseJson: Sized {
    /// Parse a JSON value element into this structure
    /// # Errors
    /// If the JSON is not valid to the format
    fn from_json_value(value: Value) -> Result<Self, CustomError>;

    /// Parse a string containing JSON into this structure
    /// # Errors
    /// If the JSON is not valid to the format
    fn from_json(value: &str) -> Result<Self, CustomError> {
        let value = serde_json::from_str::<Value>(value).map_err(|err| {
            CustomError::error(
                format!("Invalid JSON (for {})", type_name::<Self>()),
                err,
                Context::show(value.to_string()),
            )
        })?;
        Self::from_json_value(value)
    }
}

/// Parse a JSON value element into this structure using the serde JSON parser
/// # Errors
/// If the JSON is not valid to the format
#[expect(clippy::needless_pass_by_value)]
pub(crate) fn use_serde<T: DeserializeOwned>(value: Value) -> Result<T, CustomError> {
    serde_json::from_value(value.clone()).map_err(|err| {
        CustomError::error(
            format!("Could not parse JSON into {}", type_name::<T>()),
            err,
            Context::show(value.to_string()),
        )
    })
}

impl ParseJson for bool {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for String {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for isize {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for usize {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for u16 {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for u8 {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for i8 {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for f64 {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl<T: DeserializeOwned> ParseJson for RangeInclusive<T> {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl<T: ParseJson> ParseJson for Option<T> {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if value == Value::Null {
            Ok(None)
        } else {
            T::from_json_value(value).map(Some)
        }
    }
}

impl<T: ParseJson> ParseJson for Arc<T> {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        T::from_json_value(value).map(Self::new)
    }
}

impl<T: ParseJson> ParseJson for thin_vec::ThinVec<T> {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Array(arr) = value {
            arr.into_iter()
                .map(|element| T::from_json_value(element))
                .collect()
        } else {
            Err(CustomError::error(
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a list of {}",
                    type_name::<T>()
                ),
                Context::show(value),
            ))
        }
    }
}

impl<T: ParseJson> ParseJson for Vec<T> {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Array(arr) = value {
            arr.into_iter()
                .map(|element| T::from_json_value(element))
                .collect()
        } else {
            Err(CustomError::error(
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a list of {}",
                    type_name::<T>()
                ),
                Context::show(value),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson> ParseJson for (A, B) {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 2 {
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b))
            } else {
                Err(CustomError::error(
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 2 children",
                    Context::show(arr.iter().join(",")),
                ))
            }
        } else {
            Err(CustomError::error(
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                ),
                Context::show(value),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson, C: ParseJson> ParseJson for (A, B, C) {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 3 {
                let c = C::from_json_value(arr.pop().unwrap())?;
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b, c))
            } else {
                Err(CustomError::error(
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 3 children",
                    Context::show(arr.iter().join(",")),
                ))
            }
        } else {
            Err(CustomError::error(
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                    type_name::<C>(),
                ),
                Context::show(value),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson, C: ParseJson, D: ParseJson> ParseJson for (A, B, C, D) {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 4 {
                let d = D::from_json_value(arr.pop().unwrap())?;
                let c = C::from_json_value(arr.pop().unwrap())?;
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b, c, d))
            } else {
                Err(CustomError::error(
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 4 children",
                    Context::show(arr.iter().join(",")),
                ))
            }
        } else {
            Err(CustomError::error(
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {}, {}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                    type_name::<C>(),
                    type_name::<D>(),
                ),
                Context::show(value),
            ))
        }
    }
}
