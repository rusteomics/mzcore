use std::borrow::Cow;

use crate::{
    mzspeclib::{Attribute, AttributeValue, Term},
    term,
};

use mzdata::{
    params::{CURIE, Value},
    spectrum::{IsolationWindowState, ScanPolarity, SpectrumDescription, SpectrumSummary},
};

/// Create mzSpecLib attributes from a spectrum description
pub fn get_attributes(
    description: &SpectrumDescription,
    summary: Option<&SpectrumSummary>,
) -> Vec<Attribute> {
    let mut attributes = Vec::new();
    let mut group_id = 0;

    if let Some(scan) = description.acquisition.scans.first() {
        if scan.start_time != 0.0 {
            attributes.push(Attribute::new(
                term!(MS:1000894|retention time),
                Value::Float(scan.start_time),
                None,
            ));
        }
        if let Some(filter) = scan.filter_string() {
            attributes.push(Attribute::new(
                term!(MS:1000512|filter string),
                Value::String(filter.to_string()),
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
        for dissociation in precursor.activation.methods() {
            attributes.push(Attribute::new(
                term!(MS:1000044|dissociation method),
                AttributeValue::Term(Term {
                    accession: CURIE {
                        controlled_vocabulary: dissociation.controlled_vocabulary(),
                        accession: dissociation.accession(),
                    },
                    name: Cow::Borrowed(dissociation.name()),
                }),
                None,
            ));
        }

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
                    Value::Float(
                        f64::from(precursor.isolation_window.target)
                            - f64::from(precursor.isolation_window.lower_bound),
                    ),
                    None,
                ));
                attributes.push(Attribute::new(
                    term!(MS:1000829|isolation window upper offset),
                    Value::Float(
                        f64::from(precursor.isolation_window.upper_bound)
                            - f64::from(precursor.isolation_window.target),
                    ),
                    None,
                ));
            }
            IsolationWindowState::Offset => {
                if precursor.isolation_window.lower_bound != 0.0 {
                    attributes.push(Attribute::new(
                        term!(MS:1000828|isolation window lower offset),
                        Value::Float(
                            f64::from(precursor.isolation_window.target)
                                - f64::from(precursor.isolation_window.lower_bound),
                        ),
                        None,
                    ));
                }
                if precursor.isolation_window.upper_bound != 0.0 {
                    attributes.push(Attribute::new(
                        term!(MS:1000829|isolation window upper offset),
                        Value::Float(
                            f64::from(precursor.isolation_window.upper_bound)
                                - f64::from(precursor.isolation_window.target),
                        ),
                        None,
                    ));
                }
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
            if ion.intensity != 0.0 {
                attributes.push(Attribute::new(
                    term!(MS:1003085|previous MSn-1 scan precursor intensity),
                    Value::Float(f64::from(ion.intensity)),
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
        term!(MS:1003062|library spectrum index),
        Value::Int(description.index as i64),
        None,
    ));
    attributes.push(Attribute::new(
        term!(MS:1000511|ms level),
        Value::Int(description.ms_level.into()),
        None,
    ));
    match description.polarity {
        ScanPolarity::Positive => attributes.push(Attribute::new(
            term!(MS:1000465|scan polarity),
            AttributeValue::Term(term!(MS:1000130|positive scan)),
            None,
        )),
        ScanPolarity::Negative => attributes.push(Attribute::new(
            term!(MS:1000465|scan polarity),
            AttributeValue::Term(term!(MS:1000129|negative scan)),
            None,
        )),
        ScanPolarity::Unknown => (),
    }

    if let Some(summary) = summary {
        attributes.push(Attribute::new(
            term!(MS:1000285|total ion current),
            Value::Float(f64::from(summary.tic)),
            None,
        ));
        attributes.push(Attribute::new(
            term!(MS:1000505|base peak intensity),
            Value::Float(f64::from(summary.base_peak.intensity)),
            None,
        ));
        attributes.push(Attribute::new(
            term!(MS:1000504|base peak m/z),
            Value::Float(summary.base_peak.mz),
            None,
        ));
        attributes.push(Attribute::new(
            term!(MS:1003059|number of peaks),
            Value::Int(summary.count as i64),
            None,
        ));
        attributes.push(Attribute::new(
            term!(MS:1000528|lowest observed m/z),
            Value::Float(summary.mz_range.0),
            None,
        ));
        attributes.push(Attribute::new(
            term!(MS:1000527|highest observed m/z),
            Value::Float(summary.mz_range.1),
            None,
        ));
    }

    for param in description
        .params
        .iter()
        .chain(description.acquisition.params.iter().flat_map(|p| p.iter()))
        .chain(
            description
                .precursor
                .iter()
                .flat_map(|p| p.activation.params.iter()),
        )
        .chain(description.precursor.iter().flat_map(|p| {
            p.ions
                .iter()
                .flat_map(|i| i.params.iter().flat_map(|p| p.iter()))
        }))
        .chain(
            description
                .acquisition
                .scans
                .iter()
                .flat_map(|p| p.params.iter().flat_map(|p| p.iter())),
        )
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
                    value: AttributeValue::Term(Term {
                        accession: curie,
                        name: Cow::Borrowed(param.unit.for_param().1),
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
                value: AttributeValue::Scalar(Value::String(param.name.clone())),
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
                    value: AttributeValue::Term(Term {
                        accession: curie,
                        name: Cow::Borrowed(param.unit.for_param().1),
                    }),
                    group_id: Some(group_id),
                });
            }
            group_id += 1;
        }
    }
    attributes
}
