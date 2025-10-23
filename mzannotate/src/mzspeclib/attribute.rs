//! All code to make [`Attribute`]s and [`Term`]s
use std::{borrow::Cow, fmt::Display, str::FromStr};

use mzdata::{
    Param,
    params::{
        CURIE, CURIEParsingError, ControlledVocabulary, ParamValue, ParamValueParseError, Unit,
        Value,
    },
};

/// An error when parsing a term.
#[derive(Debug)]
pub enum TermParserError {
    /// The CURIE is invalid
    CURIEError(CURIEParsingError),
    /// The pipe symbol is missing
    MissingPipe,
    /// An unknown CV was used
    UnknownControlledVocabulary,
}

/// A term, a CURIE plus its name
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Term {
    /// The CURIE (e.g. MS:0000000)
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

/// Create a new term `term!(MS:1002357|PSM-level probability)`.
/// The accession + name combination is not validated.
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
            if curie.controlled_vocabulary == ControlledVocabulary::Unknown {
                Err(TermParserError::UnknownControlledVocabulary)
            } else {
                Ok(Self::new(curie, name.to_string()))
            }
        } else {
            Err(TermParserError::MissingPipe)
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
        // Displaying a CURIE with an unknown CV crashes
        if self.accession.controlled_vocabulary == ControlledVocabulary::Unknown {
            write!(f, "??:{}|{}", self.accession.accession, self.name)
        } else {
            write!(f, "{}|{}", self.accession, self.name)
        }
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

    /// Check if these are equal, plus if one if a term and the other a string check if they contain the same info
    pub(crate) fn equivalent(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Term(t), Self::Scalar(Value::String(s)))
            | (Self::Scalar(Value::String(s)), Self::Term(t)) => t.to_string() == *s,
            (Self::Scalar(Value::Int(0)), Self::Scalar(Value::Float(0.0)))
            | (Self::Scalar(Value::Float(0.0)), Self::Scalar(Value::Int(0))) => true,
            (a, b) => a == b,
        }
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

/// An error when parsing an attribute
#[derive(Debug)]
pub enum AttributeParseError {
    /// The term could not be parsed
    TermParserError(TermParserError),
    /// The value could not be parsed
    ValueParseError(ParamValueParseError),
    /// The group id was not numerical
    GroupIDParseError(std::num::ParseIntError),
    /// The equals '=' was missing
    MissingValueSeparator,
    /// The group id closing brace ']' was missing
    MissingClosingGroupBracket,
}

impl Display for AttributeParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::TermParserError(TermParserError::CURIEError(
                    CURIEParsingError::AccessionParsingError(_),
                )) => "Accession not numeric",
                Self::TermParserError(TermParserError::CURIEError(
                    CURIEParsingError::MissingNamespaceSeparator,
                )) => "There is no namespace",
                Self::TermParserError(TermParserError::CURIEError(
                    CURIEParsingError::UnknownControlledVocabulary(_),
                )) => "The controlled vocabulary is unknown",
                Self::TermParserError(TermParserError::MissingPipe) =>
                    "The term pipe symbol is missing, a term should be written like 'MS:1000896|normalized retention time'",
                Self::TermParserError(TermParserError::UnknownControlledVocabulary) =>
                    "An unknown controlled vocabulary was used",
                Self::ValueParseError(ParamValueParseError::FailedToExtractBuffer) =>
                    "Invalid value for the attribute, expected a buffer",
                Self::ValueParseError(ParamValueParseError::FailedToExtractFloat(_)) =>
                    "Invalid value for the attribute, expected a floating point",
                Self::ValueParseError(ParamValueParseError::FailedToExtractInt(_)) =>
                    "Invalid value for the attribute, expected an integer",
                Self::ValueParseError(ParamValueParseError::FailedToExtractString) =>
                    "Invalid value for the attribute, expected a string",
                Self::GroupIDParseError(_) => "The group id is not a valid number",
                Self::MissingValueSeparator => "The value separator '=' is missing",
                Self::MissingClosingGroupBracket => "The group id closing bracket ']' is missing",
            }
        )
    }
}

/// An attribute as present in a mzSpecLib file
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Attribute {
    /// The term
    pub name: Term,
    /// The value
    pub value: AttributeValue,
}

impl Attribute {
    /// Create a new attribute
    pub fn new(name: Term, value: impl Into<AttributeValue>) -> Self {
        Self {
            name,
            value: value.into(),
        }
    }

    /// Create a new attribute describing the given unit. This returns `None` if the unit is [`Unit::Unknown`].
    pub fn unit(unit: Unit) -> Option<Self> {
        let (_, name) = unit.for_param();
        let accession = unit.to_curie()?;
        Some(Self::new(
            term!(UO:0000000|unit),
            Term::new(accession, name.to_string()),
        ))
    }

    /// Parse an attribute from a string. Returns the group id (or None) and the attribute.
    /// # Errors
    /// If there is no valid attribute at this string.
    // TODO: this would be nice to have as an error with context.
    pub fn parse(s: &str) -> Result<(Option<u32>, Self), AttributeParseError> {
        if let Some(s) = s.strip_prefix('[') {
            if let Some((group_id, rest)) = s.split_once(']') {
                let group_id = match group_id.parse::<u32>() {
                    Ok(group_id) => group_id,
                    Err(e) => return Err(AttributeParseError::GroupIDParseError(e)),
                };
                if let Some((name, val)) = rest.split_once('=') {
                    let term: Term = name.parse().map_err(AttributeParseError::TermParserError)?;
                    let value: AttributeValue =
                        val.parse().map_err(AttributeParseError::ValueParseError)?;
                    Ok((Some(group_id), Self::new(term, value)))
                } else {
                    Err(AttributeParseError::MissingValueSeparator)
                }
            } else {
                Err(AttributeParseError::MissingClosingGroupBracket)
            }
        } else if let Some((name, val)) = s.split_once('=') {
            let term: Term = name.parse().map_err(AttributeParseError::TermParserError)?;
            let value: AttributeValue =
                val.parse().map_err(AttributeParseError::ValueParseError)?;
            Ok((None, Self::new(term, value)))
        } else {
            Err(AttributeParseError::MissingValueSeparator)
        }
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
        write!(f, "{}={}", self.name, self.value)
    }
}

/// A collection of attributes. They are grouped by group ID, the first is no group ID, after that it is 0 based index.
pub type Attributes = Vec<Vec<Attribute>>;

pub(crate) fn merge_attributes(a: &mut Attributes, b: &Attributes) {
    if a.is_empty() {
        a.push(Vec::new());
    }
    if !b.is_empty() {
        a[0].extend_from_slice(&b[0]);
        a.extend_from_slice(&b[1..]);
    }
}
