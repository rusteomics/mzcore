use crate::{
    mzspeclib::{Attribute, Term},
    term,
};

use mzdata::{params::Value, spectrum::SpectrumDescription};

/// Create mzSpecLib attributes from a spectrum description
pub fn into_attributes(description: &SpectrumDescription) -> Vec<Attribute> {
    let mut attributes = Vec::new();
    let mut group_id = 0;

    // TODO: add all already parsed data elements to the list

    for param in &description.params {
        if let Ok(term) = param.try_into() {
            if let Some(curie) = param.unit.to_curie() {
                attributes.push(Attribute {
                    name: term,
                    value: param.value.clone().into(),
                    group_id: Some(group_id),
                });
                attributes.push(Attribute {
                    name: term!(UO:0000000|unit),
                    value: crate::mzspeclib::AttributeValue::Term(Term {
                        accession: curie,
                        name: std::borrow::Cow::Borrowed(param.unit.for_param().1),
                    }),
                    group_id: Some(group_id),
                });
                group_id += 1;
            } else {
                attributes.push(Attribute {
                    name: term,
                    value: param.value.clone().into(),
                    group_id: None,
                });
            }
        } else {
            attributes.push(Attribute {
                name: term!(MS:1003275|other attribute name),
                value: crate::mzspeclib::AttributeValue::Scalar(Value::String(param.name.clone())),
                group_id: Some(group_id),
            });
            attributes.push(Attribute {
                name: term!(MS:1003276|other attribute value ),
                value: param.value.clone().into(),
                group_id: Some(group_id),
            });
            if let Some(curie) = param.unit.to_curie() {
                attributes.push(Attribute {
                    name: term!(UO:0000000|unit),
                    value: crate::mzspeclib::AttributeValue::Term(Term {
                        accession: curie,
                        name: std::borrow::Cow::Borrowed(param.unit.for_param().1),
                    }),
                    group_id: Some(group_id),
                });
            }
            group_id += 1;
        }
    }
    attributes
}
