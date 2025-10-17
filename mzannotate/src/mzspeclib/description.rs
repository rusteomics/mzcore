use crate::{
    mzspeclib::{Attribute, Term},
    term,
};

use mzdata::{
    params::Value,
    spectrum::{IsolationWindowState, SpectrumDescription},
};

/// Create mzSpecLib attributes from a spectrum description
pub fn into_attributes(description: &SpectrumDescription) -> Vec<Attribute> {
    let mut attributes = Vec::new();
    let mut group_id = 0;

    // TODO: extend to all parsed attributes + all specific parameters at other levels
    if let Some(scan) = description.acquisition.scans.first() {
        if scan.start_time != 0.0 {
            attributes.push(Attribute::new(
                term!(MS:1000894|retention time),
                Value::Float(scan.start_time),
                None,
            ));
        }
        if let Some(window) = scan.scan_windows.first() {
            attributes.push(Attribute::new(
                term!(MS:1000501|scan window lower limit),
                Value::Float(window.lower_bound.into()),
                None,
            ));
            attributes.push(Attribute::new(
                term!(MS:1000500|scan window upper limit),
                Value::Float(window.upper_bound.into()),
                None,
            ));
        }
    }
    if let Some(precursor) = description.precursor.first() {
        if precursor.activation.energy != 0.0 {
            attributes.push(Attribute::new(
                term!(MS:1000509|activation energy),
                Value::Float(precursor.activation.energy.into()),
                None,
            ));
        }
        match precursor.isolation_window.flags {
            IsolationWindowState::Complete => {
                attributes.push(Attribute::new(
                    term!(MS:1003208|experimental precursor monoisotopic m/z),
                    Value::Float(precursor.isolation_window.target.into()),
                    None,
                ));
                attributes.push(Attribute::new(
                    term!(MS:1000828|isolation window lower offset),
                    Value::Float(precursor.isolation_window.lower_bound.into()),
                    None,
                ));
                attributes.push(Attribute::new(
                    term!(MS:1000829|isolation window upper offset),
                    Value::Float(precursor.isolation_window.upper_bound.into()),
                    None,
                ));
            }
            IsolationWindowState::Offset => {
                attributes.push(Attribute::new(
                    term!(MS:1000828|isolation window lower offset),
                    Value::Float(precursor.isolation_window.lower_bound.into()),
                    None,
                ));
                attributes.push(Attribute::new(
                    term!(MS:1000829|isolation window upper offset),
                    Value::Float(precursor.isolation_window.upper_bound.into()),
                    None,
                ));
            }
            _ => (),
        }
        if let Some(ion) = precursor.ions.first() {
            attributes.push(Attribute::new(
                term!(MS:1000744|selected ion m/z),
                Value::Float(ion.mz),
                None,
            ));
            if let Some(c) = ion.charge {
                attributes.push(Attribute::new(
                    term!(MS:1000041|charge state),
                    Value::Int(c.into()),
                    None,
                ));
            }
        }
    }

    if !description.id.is_empty() {
        attributes.push(Attribute::new(
            term!(MS:1003061|library spectrum name),
            Value::String(description.id.clone()),
            None,
        ));
    }
    attributes.push(Attribute::new(
        term!(MS:1000511|ms level),
        Value::Int(description.ms_level.into()),
        None,
    ));
    // TODO: add all already parsed data elements to the list

    for param in description
        .params
        .iter()
        .chain(description.acquisition.params.iter().flat_map(|p| p.iter()))
    {
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
