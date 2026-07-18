use std::{any::type_name, ops::RangeInclusive, str::FromStr, sync::Arc};

use context_error::*;
use itertools::Itertools;
use mzcv::{AccessionCode, SynonymScope};
use serde::de::DeserializeOwned;
use serde_json::Value;

use crate::{
    ontology::Ontology,
    sequence::{GnoComposition, GnoSubsumption},
};

/// Custom JSON parser, needed to allow backwards compatibility to older version of JSON formats.
pub trait ParseJson: Sized {
    /// Parse a JSON value element into this structure
    /// # Errors
    /// If the JSON is not valid to the format
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>>;

    /// Parse a string containing JSON into this structure
    /// # Errors
    /// If the JSON is not valid to the format
    fn from_json(value: &str) -> Result<Self, BoxedError<'static, BasicKind>> {
        let value = serde_json::from_str::<Value>(value).map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                format!("Invalid JSON (for {})", type_name::<Self>()),
                err.to_string(),
                Context::default().lines(0, value.to_string()),
            )
        })?;
        Self::from_json_value(value)
    }
}

/// Parse a JSON value element into this structure using the serde JSON parser
/// # Errors
/// If the JSON is not valid to the format
#[expect(clippy::needless_pass_by_value)]
pub fn use_serde<T: DeserializeOwned>(value: Value) -> Result<T, BoxedError<'static, BasicKind>> {
    serde_json::from_value(value.clone()).map_err(|err| {
        BoxedError::new(
            BasicKind::Error,
            format!("Could not parse JSON into {}", type_name::<T>()),
            err.to_string(),
            Context::default().lines(0, value.to_string()),
        )
    })
}

macro_rules! use_serde {
    ($ty:ty) => {
        impl ParseJson for $ty {
            fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
                use_serde(value)
            }
        }
    };
}

use_serde!(bool);
use_serde!(String);
use_serde!(Box<str>);
use_serde!(usize);
use_serde!(u128);
use_serde!(u64);
use_serde!(u32);
use_serde!(u16);
use_serde!(u8);
use_serde!(isize);
use_serde!(i128);
use_serde!(i64);
use_serde!(i32);
use_serde!(i16);
use_serde!(i8);
use_serde!(f64);
use_serde!(f32);
use_serde!(Ontology);
use_serde!(SynonymScope);
use_serde!(GnoComposition);
use_serde!(GnoSubsumption);

impl<T: DeserializeOwned> ParseJson for RangeInclusive<T> {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

impl<T: ParseJson> ParseJson for Option<T> {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if value == Value::Null {
            Ok(None)
        } else {
            T::from_json_value(value).map(Some)
        }
    }
}

impl<T: ParseJson> ParseJson for Arc<T> {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        T::from_json_value(value).map(Self::new)
    }
}

impl<T: ParseJson> ParseJson for thin_vec::ThinVec<T> {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(arr) = value {
            arr.into_iter().map(|element| T::from_json_value(element)).collect()
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a list of {}",
                    type_name::<T>()
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl<T: ParseJson> ParseJson for Vec<T> {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(arr) = value {
            arr.into_iter().map(|element| T::from_json_value(element)).collect()
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a list of {}",
                    type_name::<T>()
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson> ParseJson for (A, B) {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 2 {
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 2 children",
                    Context::default().lines(0, arr.iter().join(",")),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson, C: ParseJson> ParseJson for (A, B, C) {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 3 {
                let c = C::from_json_value(arr.pop().unwrap())?;
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b, c))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 3 children",
                    Context::default().lines(0, arr.iter().join(",")),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                    type_name::<C>(),
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson, C: ParseJson, D: ParseJson> ParseJson for (A, B, C, D) {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 4 {
                let d = D::from_json_value(arr.pop().unwrap())?;
                let c = C::from_json_value(arr.pop().unwrap())?;
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b, c, d))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 4 children",
                    Context::default().lines(0, arr.iter().join(",")),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {}, {}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                    type_name::<C>(),
                    type_name::<D>(),
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl<A: ParseJson, B: ParseJson, C: ParseJson, D: ParseJson, E: ParseJson> ParseJson
    for (A, B, C, D, E)
{
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Array(mut arr) = value {
            if arr.len() == 5 {
                let e = E::from_json_value(arr.pop().unwrap())?;
                let d = D::from_json_value(arr.pop().unwrap())?;
                let c = C::from_json_value(arr.pop().unwrap())?;
                let b = B::from_json_value(arr.pop().unwrap())?;
                let a = A::from_json_value(arr.pop().unwrap())?;
                Ok((a, b, c, d, e))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid JSON",
                    "The JSON is a sequence but does not have 5 children",
                    Context::default().lines(0, arr.iter().join(",")),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid JSON",
                format!(
                    "The JSON has to be a sequence to parse a ({}, {}, {}, {}, {})",
                    type_name::<A>(),
                    type_name::<B>(),
                    type_name::<C>(),
                    type_name::<D>(),
                    type_name::<E>(),
                ),
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

impl ParseJson for AccessionCode {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        match value {
            Value::Number(n)
                if let Some(n) = n.as_u64()
                    && let Ok(n) = u32::try_from(n) =>
            {
                Ok(Self::Numeric(n))
            }
            Value::String(n) => Self::from_str(&n).map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid alphanumeric ModificationID",
                    err.to_string(),
                    Context::default().lines(0, n),
                )
            }),
            Value::Object(map) => map.get("Numeric").map_or_else(
                || {
                    map.get("Alphanumeric").map_or_else(
                        || {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid ModificationID",
                                "The 'id' can only be 'Numeric' or 'Alphanumeric'",
                                Context::default().lines(
                                    0,
                                    map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                                ),
                            ))
                        },
                        |n| {
                            Self::from_str(&n.to_string()).map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid alphanumeric ModificationID",
                                    err.to_string(),
                                    Context::default().lines(0, n.to_string()),
                                )
                            })
                        },
                    )
                },
                |n| {
                    u32::from_json_value(n.clone()).map_or_else(
                        |_| {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid numeric ModificationID",
                                format!(
                                    "The id has to be a positive integer within 0 to {}",
                                    u32::MAX
                                ),
                                Context::default().lines(0, n.to_string()),
                            ))
                        },
                        |n| Ok(Self::Numeric(n)),
                    )
                },
            ),
            value => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid ModificationID",
                "The 'id' is in an invalid format",
                Context::default().lines(0, value.to_string()),
            )),
        }
    }
}
