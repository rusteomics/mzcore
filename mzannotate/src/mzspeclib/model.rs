use std::{collections::HashMap, fmt::Display};

use mzcore::{prelude::PeptidoformIon, system::isize::Charge};
use mzdata::{curie, params::Value};

use crate::{
    mzspeclib::{
        Attribute, AttributeSet, AttributeValue, Attributed, AttributedMut, EntryType, Term,
        impl_attributed,
    },
    term,
};

/// The header for a spectral library
#[derive(Debug, Clone)]
pub struct LibraryHeader {
    /// The version of the format
    pub format_version: String,
    /// The attributes for this library
    pub attributes: Vec<Attribute>,
    /// The attribute classes for this library
    pub attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
}

impl Default for LibraryHeader {
    fn default() -> Self {
        Self {
            format_version: "1.0".into(),
            attributes: Vec::new(),
            attribute_classes: HashMap::new(),
        }
    }
}

impl_attributed!(mut LibraryHeader);

impl Display for LibraryHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let version = Attribute::new(
            Term::new(curie!(MS:1_003_186), "library format version".into()),
            AttributeValue::Scalar(Value::String(self.format_version.clone())),
            None,
        );
        writeln!(f, "<MzSpecLib>\n{version}")?;
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl LibraryHeader {
    /// Create a new library header
    pub const fn new(
        format_version: String,
        attributes: Vec<Attribute>,
        attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
    ) -> Self {
        Self {
            format_version,
            attributes,
            attribute_classes,
        }
    }
}

/// An ID
pub type IdType = u32;

/// An analyte that resulted in the spectrum.
#[derive(Default, Debug, Clone)]
pub struct Analyte {
    /// The numeric id of this analyte
    pub id: IdType,
    /// The peptidoform, this can be missing if the analyte is not a peptidoform (uses term MS:1003270)
    pub peptidoform_ion: Option<PeptidoformIon>,
    /// The charge (MS:1000041)
    pub charge: Option<Charge>,
    /// Other attributes for this analyte
    pub attributes: Vec<Attribute>,
}

impl_attributed!(mut Analyte);

impl Display for Analyte {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Analyte={}>", self.id)?;
        if let Some(sequence) = &self.peptidoform_ion {
            writeln!(
                f,
                "{}={sequence}",
                term!(MS:1_003_270|"proforma peptidoform ion notation")
            )?;
        }
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl Analyte {
    /// Create a new analyte
    pub const fn new(
        id: IdType,
        peptidoform_ion: Option<PeptidoformIon>,
        charge: Option<Charge>,
        attributes: Vec<Attribute>,
    ) -> Self {
        Self {
            id,
            peptidoform_ion,
            charge,
            attributes,
        }
    }
}

#[derive(Default, Debug, Clone)]
pub struct Interpretation {
    pub id: IdType,
    pub attributes: Vec<Attribute>,
    pub analyte_refs: Vec<IdType>,
    pub members: Vec<InterpretationMember>,
}

impl Interpretation {
    pub const fn new(
        id: IdType,
        attributes: Vec<Attribute>,
        analyte_refs: Vec<IdType>,
        members: Vec<InterpretationMember>,
    ) -> Self {
        Self {
            id,
            attributes,
            analyte_refs,
            members,
        }
    }
}

impl_attributed!(mut Interpretation);

impl Display for Interpretation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Interpretation={}>", self.id)?;
        if !self.analyte_refs.is_empty() {
            let val = AttributeValue::List(
                self.analyte_refs
                    .iter()
                    .map(|v| Value::Int(i64::from(*v)))
                    .collect(),
            );
            let mixture_ids = Attribute::new(
                Term::new(curie!(MS:1_003_163), "analyte mixture members".into()),
                val,
                None,
            );
            writeln!(f, "{mixture_ids}")?;
        }
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        for member in &self.members {
            write!(f, "{member}")?;
        }
        Ok(())
    }
}

#[derive(Default, Debug, Clone)]
pub struct InterpretationMember {
    pub id: IdType,
    pub attributes: Vec<Attribute>,
    pub analyte_ref: IdType,
}

impl_attributed!(mut InterpretationMember);

impl Display for InterpretationMember {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<InterpretationMember={}>", self.id)?;
        for attr in &self.attributes {
            f.write_str(attr.to_string().as_str())?;
        }
        Ok(())
    }
}
