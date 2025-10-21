use std::{borrow::Cow, collections::HashMap, fmt::Display, str::FromStr};

use context_error::Context;
use mzdata::{
    Param,
    params::{CURIE, CURIEParsingError, ParamValue, ParamValueParseError, Unit, Value},
};

use crate::mzspeclib::Id;

#[derive(Debug)]
pub enum TermParserError {
    CURIEError(CURIEParsingError),
    MissingPipe(String),
}

/// A term, a CURIE plus its name
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Term {
    /// The CURIE (eg MS:0000000)
    pub accession: CURIE,
    /// The name
    pub name: Cow<'static, str>, // Static Cow to allow to have term in match arms
}

impl PartialOrd for Term {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Term {
    // TODO: only does sorting on accession, might be better to also sort on CV
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.accession.accession.cmp(&other.accession.accession)
    }
}

/// Create a new term `term!(MS:1002357|PSM-level probability)`. The accession/name combination
/// is not validated.
#[macro_export]
macro_rules! term {
    ($ns:ident:$accession:literal|$($name:tt)+) =>  {
        $crate::mzspeclib::Term {
            accession: mzdata::curie!($ns:$accession),
            name: std::borrow::Cow::Borrowed(stringify!($($name)+))
        }
    };
}

impl FromStr for Term {
    type Err = TermParserError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((curie, name)) = s.split_once('|') {
            let curie: CURIE = curie.parse().map_err(TermParserError::CURIEError)?;
            Ok(Self::new(curie, name.to_string()))
        } else {
            Err(TermParserError::MissingPipe(s.to_string()))
        }
    }
}

impl TryFrom<Param> for Term {
    type Error = ();
    fn try_from(value: Param) -> Result<Self, ()> {
        Ok(Self {
            accession: value.curie().ok_or(())?,
            name: value.name.into(),
        })
    }
}

impl TryFrom<&Param> for Term {
    type Error = ();
    fn try_from(value: &Param) -> Result<Self, ()> {
        Ok(Self {
            accession: value.curie().ok_or(())?,
            name: value.name.clone().into(),
        })
    }
}

impl Term {
    /// Create a new term, the term is not validated
    pub const fn new(accession: CURIE, name: String) -> Self {
        Self {
            accession,
            name: Cow::Owned(name),
        }
    }

    /// Create a [`Param`] from this term for use in mzdata.
    pub fn into_param(self, value: Value, unit: Unit) -> Param {
        Param {
            name: self.name.to_string(),
            value,
            accession: Some(self.accession.accession),
            controlled_vocabulary: Some(self.accession.controlled_vocabulary),
            unit,
        }
    }
}

impl Display for Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}|{}", self.accession, self.name)
    }
}

/// A value for an attribute
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttributeValue {
    /// A single value
    Scalar(Value),
    /// A list of values
    List(Vec<Value>),
    /// A term
    Term(Term),
}

impl AttributeValue {
    /// Coerce this value into a scalar value, making a string of the list and term options.
    pub fn scalar(&self) -> Cow<'_, Value> {
        match self {
            Self::Scalar(value) => Cow::Borrowed(value),
            Self::List(_) => Cow::Owned(Value::String(self.to_string())),
            Self::Term(term) => Cow::Owned(Value::String(term.to_string())),
        }
    }

    /// Create a scalar attribute
    pub fn from_scalar(value: impl Into<Value>) -> Self {
        Self::Scalar(value.into())
    }
}

impl FromStr for AttributeValue {
    type Err = ParamValueParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(term) = s.parse() {
            Ok(Self::Term(term))
        } else if s.contains(',') {
            let vals: Vec<Value> = s.split(',').flat_map(Value::from_str).collect();
            if vals.iter().all(|v| v.is_numeric() || v.is_boolean()) {
                Ok(Self::List(vals))
            } else {
                let v = Value::from_str(s)?;
                Ok(Self::Scalar(v))
            }
        } else {
            let v = if s.starts_with('\"') && s.ends_with('\"') && s.len() > 1 {
                Value::from_str(&s[1..s.len() - 1])
            } else {
                Value::from_str(s)
            }?;
            Ok(Self::Scalar(v))
        }
    }
}

impl From<Value> for AttributeValue {
    fn from(value: Value) -> Self {
        Self::Scalar(value)
    }
}

impl From<Term> for AttributeValue {
    fn from(value: Term) -> Self {
        Self::Term(value)
    }
}

impl From<AttributeValue> for Value {
    fn from(value: AttributeValue) -> Self {
        value.scalar().into_owned()
    }
}

impl Display for AttributeValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Scalar(value) => write!(f, "{value}"),
            Self::List(values) => {
                for (i, value) in values.iter().enumerate() {
                    if i != 0 {
                        write!(f, ",")?;
                    }
                    write!(f, "{value}")?;
                }
                Ok(())
            }
            Self::Term(term) => write!(f, "{term}"),
        }
    }
}

impl From<Vec<Value>> for AttributeValue {
    fn from(value: Vec<Value>) -> Self {
        Self::List(value)
    }
}

#[derive(Debug)]
pub enum AttributeParseError {
    TermParserError(TermParserError, String),
    ValueParseError(ParamValueParseError, String),
    GroupIDParseError(std::num::ParseIntError, String),
    MissingValueSeparator(String),
    Malformed(String),
}

/// An attribute as present in a mzSpecLib file
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attribute {
    /// The term
    pub name: Term,
    /// The value
    pub value: AttributeValue,
    /// The group id (or non if not part of a group)
    pub group_id: Option<u32>,
}

impl FromStr for Attribute {
    type Err = AttributeParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some(s) = s.strip_prefix('[') {
            if let Some((group_id, rest)) = s.split_once(']') {
                let group_id = match group_id.parse::<u32>() {
                    Ok(group_id) => group_id,
                    Err(e) => return Err(AttributeParseError::GroupIDParseError(e, s.to_string())),
                };
                if let Some((name, val)) = rest.split_once('=') {
                    let term: Term = name
                        .parse()
                        .map_err(|e| AttributeParseError::TermParserError(e, s.to_string()))?;
                    let value: AttributeValue = val
                        .parse()
                        .map_err(|e| AttributeParseError::ValueParseError(e, s.to_string()))?;
                    Ok(Self::new(term, value, Some(group_id)))
                } else {
                    Err(AttributeParseError::MissingValueSeparator(s.to_string()))
                }
            } else {
                Err(Self::Err::Malformed(s.to_string()))
            }
        } else if let Some((name, val)) = s.split_once('=') {
            let term: Term = name
                .parse()
                .map_err(|e| AttributeParseError::TermParserError(e, s.to_string()))?;
            let value: AttributeValue = val
                .parse()
                .map_err(|e| AttributeParseError::ValueParseError(e, s.to_string()))?;
            Ok(Self::new(term, value, None))
        } else {
            Err(AttributeParseError::MissingValueSeparator(s.to_string()))
        }
    }
}

impl Attribute {
    /// Create a new attribute
    pub fn new(name: Term, value: impl Into<AttributeValue>, group_id: Option<u32>) -> Self {
        Self {
            name,
            value: value.into(),
            group_id,
        }
    }

    /// Create a new attribute describing the given unit. This returns `None` if the unit is [`Unit::Unknown`].
    pub fn unit(unit: Unit, group_id: Option<u32>) -> Option<Self> {
        let (_, name) = unit.for_param();
        let accession = unit.to_curie()?;
        Some(Self::new(
            term!(UO:0000000|unit),
            Term::new(accession, name.to_string()),
            group_id,
        ))
    }
}

impl From<Attribute> for Param {
    fn from(value: Attribute) -> Self {
        Self {
            name: value.name.name.to_string(),
            value: value.value.into(),
            accession: Some(value.name.accession.accession),
            controlled_vocabulary: Some(value.name.accession.controlled_vocabulary),
            unit: Unit::Unknown,
        }
    }
}

impl Display for Attribute {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.group_id {
            Some(i) => {
                write!(f, "[{i}]{}={}", self.name, self.value)
            }
            None => {
                write!(f, "{}={}", self.name, self.value)
            }
        }
    }
}

/// A set of attributes in the library header that contains global settings for all or a subset of spectra
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AttributeSet {
    /// The id
    pub id: String,
    /// The type of attributes
    pub namespace: EntryType,
    /// The attributes themselves grouped by group id and with the original context to provide good error messages
    pub attributes: HashMap<Option<Id>, Vec<(Attribute, Context<'static>)>>,
}

impl PartialOrd for AttributeSet {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.id.partial_cmp(&other.id) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        match self.namespace.partial_cmp(&other.namespace) {
            Some(core::cmp::Ordering::Equal) => {}
            ord => return ord,
        }
        self.attributes.len().partial_cmp(&other.attributes.len())
    }
}

/// All types of entries for global attribute sets
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntryType {
    /// Spectrum attributes
    Spectrum,
    /// Analyte attributes
    Analyte,
    /// Interpretation attributes
    Interpretation,
    /// A cluster
    Cluster,
}

impl Display for EntryType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Spectrum => "Spectrum",
                Self::Analyte => "Analyte",
                Self::Interpretation => "Interpretation",
                Self::Cluster => "Cluster",
            }
        )
    }
}
