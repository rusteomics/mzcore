use std::{borrow::Cow, fmt::Display, str::FromStr};

use mzdata::{
    curie,
    params::{CURIE, CURIEParsingError, ParamValue, ParamValueParseError, Unit, Value},
};

#[derive(Debug)]
pub enum TermParserError {
    CURIEError(CURIEParsingError),
    MissingPipe(String),
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct Term {
    pub accession: CURIE,
    pub name: Cow<'static, str>, // Static Cow to allow to have term in match arms
}

impl PartialOrd for Term {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Term {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        if self.accession == other.accession {
            return std::cmp::Ordering::Equal;
        }
        self.name.cmp(&other.name)
    }
}

#[macro_export]
macro_rules! term {
    ($ns:ident:$accession:literal|$name:literal) =>  {
        $crate::mzspeclib::Term {
            accession: curie!($ns:$accession),
            name: std::borrow::Cow::Borrowed($name)
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

impl TryFrom<mzdata::Param> for Term {
    type Error = ();
    fn try_from(value: mzdata::Param) -> Result<Self, ()> {
        Ok(Self {
            accession: value.curie().ok_or(())?,
            name: value.name.into(),
        })
    }
}

impl Term {
    pub const fn new(accession: CURIE, name: String) -> Self {
        Self {
            accession,
            name: Cow::Owned(name),
        }
    }
}

impl Display for Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}|{}", self.accession, self.name)
    }
}

#[derive(Debug)]
pub enum AttributeValueParseError {
    ParamValueParseError(ParamValueParseError),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum AttributeValue {
    Scalar(Value),
    List(Vec<Value>),
    Term(Term),
}

impl AttributeValue {
    pub fn scalar(&self) -> Cow<'_, Value> {
        match self {
            Self::Scalar(value) => Cow::Borrowed(value),
            Self::List(_) => Cow::Owned(Value::String(self.to_string())),
            Self::Term(term) => Cow::Owned(Value::String(term.to_string())),
        }
    }

    pub fn from_scalar(value: impl Into<Value>) -> Self {
        Self::Scalar(value.into())
    }
}

impl FromStr for AttributeValue {
    type Err = ParamValueParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(term) = s.parse() {
            Ok(Self::Term(term))
        } else if s.contains(",") {
            let vals: Vec<Value> = s.split(",").flat_map(Value::from_str).collect();
            if vals.iter().all(|v| v.is_numeric() || v.is_boolean()) {
                Ok(Self::List(vals))
            } else {
                let v = Value::from_str(s)?;
                Ok(Self::Scalar(v))
            }
        } else {
            let v = Value::from_str(s)?;
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
            Self::Term(attribute_name) => f.write_str(attribute_name.to_string().as_str()),
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

#[derive(Debug, Clone, PartialEq)]
pub struct Attribute {
    pub name: Term,
    pub value: AttributeValue,
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
    pub fn new(name: Term, value: impl Into<AttributeValue>, group_id: Option<u32>) -> Self {
        Self {
            name,
            value: value.into(),
            group_id,
        }
    }

    pub fn unit(unit: Unit, group_id: Option<u32>) -> Option<Self> {
        let (_, name) = unit.for_param();
        let accession = unit.to_curie()?;
        Some(Self::new(
            term!(UO:0_000_000|"unit"),
            Term::new(accession, name.to_string()),
            group_id,
        ))
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

pub trait Attributed {
    fn attributes(&self) -> &[Attribute];

    fn find(&self, term: &Term) -> Option<(usize, &Attribute)> {
        self.attributes()
            .iter()
            .position(|v| {
                if v.name == *term {
                    true
                } else if v.name.accession == curie!(MS:1003275) {
                    if let AttributeValue::Term(v) = &v.value {
                        *v == *term
                    } else {
                        false
                    }
                } else {
                    false
                }
            })
            .map(|i| (i, &self.attributes()[i]))
    }

    fn find_by_id(&self, id: CURIE) -> Option<(usize, &Attribute)> {
        self.attributes()
            .iter()
            .position(|v| {
                if v.name.accession == id {
                    true
                } else if v.name.accession == curie!(MS:1003275) {
                    if let AttributeValue::Term(v) = &v.value {
                        v.accession == id
                    } else {
                        false
                    }
                } else {
                    false
                }
            })
            .map(|i| (i, &self.attributes()[i]))
    }

    fn find_group(&self, group_id: u32) -> Vec<&Attribute> {
        self.attributes()
            .iter()
            .filter(|a| a.group_id == Some(group_id))
            .collect()
    }

    fn find_all(&self, term: &Term) -> impl Iterator<Item = &Attribute> {
        self.attributes().iter().filter(|v| {
            if v.name == *term {
                true
            } else if v.name.accession == curie!(MS:1003275) {
                if let AttributeValue::Term(v) = &v.value {
                    *v == *term
                } else {
                    false
                }
            } else {
                false
            }
        })
    }

    fn iter_attribute_groups(&self) -> impl Iterator<Item = Vec<&Attribute>> {
        let last_group_id = self.find_last_group_id().unwrap_or_default();
        (0..=last_group_id).map(|v| self.find_group(v))
    }

    fn find_last_group_id(&self) -> Option<u32> {
        self.attributes()
            .iter()
            .map(|v| v.group_id)
            .reduce(|prev, new| {
                new.map(|new| prev.map_or(new, |prev| prev.max(new)))
                    .or(prev)
            })
            .unwrap_or_default()
    }
}

pub trait AttributedMut: Attributed {
    fn attributes_mut(&mut self) -> &mut [Attribute];

    fn find_group_mut(&mut self, group_id: u32) -> impl Iterator<Item = &mut Attribute> {
        self.attributes_mut()
            .iter_mut()
            .filter(move |a| a.group_id == Some(group_id))
    }

    fn add_attribute(&mut self, attr: Attribute);

    fn remove_attribute(&mut self, index: usize) -> Attribute;

    fn extend_attributes(&mut self, attrs: impl IntoIterator<Item = Attribute>) {
        for attr in attrs {
            self.add_attribute(attr);
        }
    }

    fn add_attribute_with_unit(
        &mut self,
        mut attr: Attribute,
        unit: Unit,
        group_id: Option<u32>,
    ) -> Option<u32> {
        let group_id =
            Some(group_id.unwrap_or_else(|| self.find_last_group_id().unwrap_or_default() + 1));
        attr.group_id = group_id;
        self.add_attribute(attr);
        if let Some(attr) = Attribute::unit(unit, group_id) {
            self.add_attribute(attr);
        }
        group_id
    }
}

macro_rules! impl_attributed {
    ($t:ty) => {
        impl Attributed for $t {
            fn attributes(&self) -> &[Attribute] {
                &self.attributes
            }
        }
    };

    (mut $t:ty) => {
        impl_attributed!($t);

        impl AttributedMut for $t {
            fn attributes_mut(&mut self) -> &mut [Attribute] {
                &mut self.attributes
            }

            fn add_attribute(&mut self, attr: Attribute) {
                self.attributes.push(attr);
            }

            fn remove_attribute(&mut self, index: usize) -> Attribute {
                self.attributes.remove(index)
            }
        }
    };
}

pub(crate) use impl_attributed;

impl Attributed for Vec<Attribute> {
    fn attributes(&self) -> &[Attribute] {
        self
    }
}

impl AttributedMut for Vec<Attribute> {
    fn attributes_mut(&mut self) -> &mut [Attribute] {
        self
    }

    fn add_attribute(&mut self, attr: Attribute) {
        self.push(attr);
    }

    fn remove_attribute(&mut self, index: usize) -> Attribute {
        self.remove(index)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct AttributeSet {
    pub id: String,
    pub namespace: EntryType,
    pub attributes: Vec<Attribute>,
}

impl_attributed!(mut AttributeSet);

impl AttributeSet {
    pub fn new(id: String, namespace: EntryType, attributes: Vec<Attribute>) -> Self {
        Self {
            id,
            namespace,
            attributes,
        }
    }
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EntryType {
    Spectrum,
    Analyte,
    Interpretation,
    Cluster,
}
