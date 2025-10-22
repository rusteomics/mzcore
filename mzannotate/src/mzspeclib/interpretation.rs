use std::collections::HashMap;

use mzdata::params::Value;

use crate::{
    mzspeclib::{Attribute, AttributeValue, Id},
    term,
};

/// An interpretation, this gives details on how well the analytes match to the spectrum.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Interpretation {
    /// The ID
    pub id: Id,
    /// PSM level probability (MS:1002357) in range 0..=1
    pub probability: Option<f64>,
    /// Other attributes
    pub attributes: Vec<Attribute>,
    /// Analyte mixture members (MS:1003163)
    pub analyte_refs: Vec<Id>,
    /// The interpretation members, this links to a specific analyte and gives additional attributes
    /// about the interpretation of that analyte.
    pub members: HashMap<Id, Vec<Attribute>>,
}

impl Interpretation {
    pub fn attributes(&self) -> Vec<Attribute> {
        // TODO: add easily calculated attributes
        let mut attributes = Vec::new();
        if !self.analyte_refs.is_empty() {
            attributes.push(Attribute::new(
                term!(MS:1003163|analyte mixture members),
                AttributeValue::List(
                    self.analyte_refs
                        .iter()
                        .map(|v| Value::Int(i64::from(*v)))
                        .collect(),
                ),
                None,
            ));
        }
        if let Some(probability) = self.probability {
            attributes.push(Attribute::new(
                term!(MS:1002357|PSM-level probability),
                Value::Float(probability),
                None,
            ));
        }
        attributes.extend_from_slice(&self.attributes);

        attributes
    }
}
