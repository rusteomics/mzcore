use std::collections::HashMap;

use mzdata::params::Value;

use crate::{
    mzspeclib::{Attribute, AttributeValue, Attributes, Id, merge_attributes},
    term,
};

/// An interpretation, this gives details on how well the analytes match to the spectrum.
#[derive(Debug, Clone, PartialEq)]
pub struct Interpretation {
    /// The ID
    pub id: Id,
    /// PSM level probability (MS:1002357) in range 0..=1
    pub probability: Option<f64>,
    /// Other attributes
    pub attributes: Attributes,
    /// Analyte mixture members (MS:1003163)
    pub analyte_refs: Vec<Id>,
    /// The interpretation members, this links to a specific analyte and gives additional attributes
    /// about the interpretation of that analyte.
    pub members: HashMap<Id, Attributes>,
}

impl Default for Interpretation {
    fn default() -> Self {
        Self {
            id: 0,
            probability: None,
            attributes: vec![Vec::new(); 1],
            analyte_refs: Vec::new(),
            members: HashMap::new(),
        }
    }
}

impl Interpretation {
    pub fn attributes(&self) -> Attributes {
        // TODO: add easily calculated attributes
        let mut attributes = vec![Vec::new(); 1];
        if !self.analyte_refs.is_empty() {
            attributes[0].push(Attribute::new(
                term!(MS:1003163|analyte mixture members),
                AttributeValue::List(
                    self.analyte_refs
                        .iter()
                        .map(|v| Value::Int(i64::from(*v)))
                        .collect(),
                ),
            ));
        }
        if let Some(probability) = self.probability {
            attributes[0].push(Attribute::new(
                term!(MS:1002357|PSM-level probability),
                Value::Float(probability),
            ));
        }
        merge_attributes(&mut attributes, &self.attributes);

        attributes
    }
}
