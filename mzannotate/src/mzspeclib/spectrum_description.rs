//! Handle converting [`Attribute`]s to and from [`SpectrumDescription`]
use std::borrow::Cow;

use crate::{
    mzspeclib::{Attribute, AttributeValue, MzSpecLibErrorKind, Term},
    term,
};

use context_error::{BoxedError, Context, CreateError};
use mzdata::{
    curie,
    params::{CURIE, ParamDescribed, ParamValue, Unit, Value},
    spectrum::{IsolationWindowState, ScanPolarity, SpectrumDescription, SpectrumSummary},
};

/// Create mzSpecLib attributes from a spectrum description
pub fn get_attributes_from_spectrum_description(
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

/// Parse a set of attributes and store them in the correct location of the spectrum description.
/// Any unrecognised attributes will be stored in the provided list of attributes.
/// # Errors
/// If an attributes contains an invalid value.
pub(crate) fn populate_spectrum_description_from_attributes<'a>(
    attributes: impl Iterator<Item = (&'a Option<u32>, &'a Vec<(Attribute, Context<'static>)>)>,
    description: &mut SpectrumDescription,
    spectrum_attributes: &mut Vec<Attribute>,
) -> Result<(), BoxedError<'static, MzSpecLibErrorKind>> {
    let mut last_group_id = spectrum_attributes
        .iter()
        .filter_map(|a| a.group_id)
        .max()
        .unwrap_or_default();
    for (id, group) in attributes {
        if id.is_some() {
            let unit = if let Some((
                Attribute {
                    value: AttributeValue::Term(term),
                    ..
                },
                _,
            )) = group
                .iter()
                .find(|a| a.0.name.accession == curie!(UO:0000000))
            {
                Some(Unit::from_curie(&term.accession))
            } else {
                None
            };

            if let Some((Attribute { value, .. }, context)) = group
                .iter()
                .find(|a| a.0.name.accession == curie!(MS:1000894))
                && group.len() == 2
            {
                description.acquisition.scans[0].start_time =
                    f64::from(value.scalar().to_f32().map_err(|v| {
                        BoxedError::new(
                            MzSpecLibErrorKind::Attribute,
                            "Invalid attribute",
                            v.to_string(),
                            context.clone(),
                        )
                    })?) / match unit.ok_or_else(|| {
                        BoxedError::new(
                            MzSpecLibErrorKind::MissingUnit,
                            "Invalid attribute",
                            "MS:1000894|retention time needs a unit",
                            context.clone(),
                        )
                    })? {
                        Unit::Minute => 60.0,
                        _ => 1.0, // Assume seconds for anything else
                    };
            } else if let Some((Attribute { value, .. }, context)) = group
                .iter()
                .find(|a| a.0.name.accession == curie!(MS:1000896))
                && group.len() == 2
            {
                let rt = f64::from(value.scalar().to_f32().map_err(|v| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid attribute",
                        v.to_string(),
                        context.clone(),
                    )
                })?) / match unit.ok_or_else(|| {
                    BoxedError::new(
                        MzSpecLibErrorKind::MissingUnit,
                        "Invalid attribute",
                        "MS:1000896|normalized retention time needs a unit",
                        context.clone(),
                    )
                })? {
                    Unit::Minute => 60.0,
                    _ => 1.0, // Assume seconds for anything else
                };
                // This is normalised time so if normal time is already set ignore this param
                if description.acquisition.scans[0].start_time == 0.0 {
                    description.acquisition.scans[0].start_time = rt;
                } else {
                    description.acquisition.scans[0].add_param(
                        term!(MS:1000896|normalized retention time)
                            .into_param(Value::Float(rt), Unit::Second),
                    );
                }
            } else if let Some((attr, context)) = group.iter().find(|a| {
                a.0.name.accession == curie!(MS:1000045) || a.0.name.accession == curie!(MS:1000509)
            }) && group.len() == 2
            {
                description.precursor[0].activation.energy =
                    attr.value.scalar().to_f32().map_err(|v| {
                        BoxedError::new(
                            MzSpecLibErrorKind::Attribute,
                            "Invalid attribute",
                            v.to_string(),
                            context.clone(),
                        )
                    })?;
            } else if let Some(unit) = unit
                && group.len() == 2
            {
                if let Some((other, _)) = group
                    .iter()
                    .find(|a| a.0.name.accession != curie!(UO:0000000))
                {
                    description.add_param(
                        other
                            .name
                            .clone()
                            .into_param(other.value.clone().into(), unit),
                    );
                } else {
                    // Two units as a single group is not a sensible thing to have so just ignore
                }
            } else if let Some((value, _)) = group
                .iter()
                .find(|a| a.0.name.accession == curie!(MS:1003276))
                && let Some((name, _)) = group
                    .iter()
                    .find(|a| a.0.name.accession == curie!(MS:1003275))
                && group.len() == 2
            {
                description.add_param(match &name.value {
                    AttributeValue::Term(term) => mzdata::params::Param {
                        accession: Some(term.accession.accession),
                        name: term.name.to_string(),
                        value: value.value.scalar().into_owned(),
                        controlled_vocabulary: Some(term.accession.controlled_vocabulary),
                        unit: Unit::Unknown,
                    },
                    e => mzdata::params::Param {
                        accession: None,
                        name: e.scalar().to_string(),
                        value: value.value.scalar().into_owned(),
                        controlled_vocabulary: None,
                        unit: Unit::Unknown,
                    },
                });
            } else {
                last_group_id += 1;
                for attr in group {
                    let mut attr = attr.0.clone();
                    attr.group_id = Some(last_group_id);
                    spectrum_attributes.push(attr);
                }
            }
        } else {
            for (attr, context) in group {
                #[allow(clippy::unnested_or_patterns)]
                match attr.name.accession {
                    curie!(MS:1003208) => {
                        let window = &mut description.precursor[0].isolation_window;
                        window.target = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                        if matches!(window.flags, IsolationWindowState::Offset) {
                            window.lower_bound = window.target - window.lower_bound;
                            window.upper_bound += window.target;
                        }
                        window.flags = IsolationWindowState::Complete;
                    }
                    curie!(MS:1000828) => {
                        let offset = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                        let window = &mut description.precursor[0].isolation_window;
                        match window.flags {
                            IsolationWindowState::Unknown => {
                                window.flags = IsolationWindowState::Offset;
                                window.lower_bound = offset;
                            }
                            IsolationWindowState::Complete => {
                                window.lower_bound = window.target - offset;
                            }
                            _ => {}
                        }
                    }
                    curie!(MS:1000829) => {
                        let offset = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                        let window = &mut description.precursor[0].isolation_window;
                        match window.flags {
                            IsolationWindowState::Unknown => {
                                window.flags = IsolationWindowState::Offset;
                                window.upper_bound = offset;
                            }
                            IsolationWindowState::Complete => {
                                window.upper_bound = window.target - offset;
                            }
                            _ => {}
                        }
                    }
                    curie!(MS:1000794) => { // TODO: Also handle these when they are in a group with a unit
                        let limit = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                        let window = &mut description.precursor[0].isolation_window;
                        if matches!(
                            window.flags,
                            IsolationWindowState::Unknown | IsolationWindowState::Explicit
                        ) {
                            window.flags = IsolationWindowState::Explicit;
                            window.lower_bound = limit;
                        }
                    }
                    curie!(MS:1000793) => {
                        let limit = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                        let window = &mut description.precursor[0].isolation_window;
                        if matches!(
                            window.flags,
                            IsolationWindowState::Unknown | IsolationWindowState::Explicit
                        ) {
                            window.flags = IsolationWindowState::Explicit;
                            window.upper_bound = limit;
                        }
                    }
                    curie!(MS:1000501) => {
                        description.acquisition.scans[0].scan_windows[0].lower_bound =
                            attr.value.scalar().to_f32().map_err(|v| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                )
                            })?;
                    }
                    curie!(MS:1000500) => {
                        description.acquisition.scans[0].scan_windows[0].upper_bound =
                            attr.value.scalar().to_f32().map_err(|v| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                )
                            })?;
                    }
                    curie!(MS:1000894) => {
                        description.acquisition.scans[0].start_time =
                        attr.value.scalar().to_f64().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })? / 60.0;
                    }
                    curie!(MS:1000744) => {
                        description.precursor[0].ions[0].mz =
                            attr.value.scalar().to_f64().map_err(|v| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                )
                            })?;
                    }
                    curie!(MS:1000041) => {
                        description.precursor[0].ions[0].charge =
                            Some(attr.value.scalar().to_i32().map_err(|v| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                )
                            })?);
                    }
                    curie!(MS:1002476) => {
                        description.acquisition.scans[0].add_param(term!(MS:1002476|ion mobility drift time).into_param(attr.value.scalar().into_owned(), Unit::Unknown));
                    }
                    curie!(MS:1002815) => {
                        description.acquisition.scans[0].add_param(term!(MS:1002815|inverse reduced ion mobility drift time).into_param(attr.value.scalar().into_owned(), Unit::Unknown));
                    }
                    curie!(MS:1001581) => {
                        description.acquisition.scans[0].add_param(term!(MS:1001581|FAIMS compensation voltage).into_param(attr.value.scalar().into_owned(), Unit::Unknown));
                    }
                    curie!(MS:1003371) => {
                        description.acquisition.scans[0].add_param(term!(MS:1003371|SELEXION compensation voltage).into_param(attr.value.scalar().into_owned(), Unit::Unknown));
                    }
                    curie!(MS:1003063) => {
                        description.add_param(
                            attr.name
                                .clone()
                                .into_param(attr.value.clone().into(), Unit::Unknown),
                        );
                    }
                    curie!(MS:1003061) => description.id = attr.value.to_string(),
                    curie!(MS:1000511) => {
                        description.ms_level =
                            u8::try_from(attr.value.scalar().to_u64().map_err(|v| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                )
                            })?)
                            .map_err(|_| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid attribute",
                                    "The MS level is too large of a number",
                                    context.clone(),
                                )
                            })?;
                    }
                    curie!(MS:1000512) => description.acquisition.scans[0].add_param(
                        term!(MS:1000512|filter string)
                            .into_param(attr.value.clone().into(), Unit::Unknown),
                    ),
                    curie!(MS:1000008) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            description
                                .add_param(term.into_param(Value::Empty, Unit::Unknown));
                        } else {
                            description.add_param(
                                term!(MS:1000008|ionization type)
                                    .into_param(attr.value.clone().into(), Unit::Unknown),
                            );
                        }
                    }
                    curie!(MS:1000465) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            match term.accession {
                                curie!(MS:1000130) => {
                                    description.polarity =
                                        ScanPolarity::Positive;
                                }
                                curie!(MS:1000129) => {
                                    description.polarity =
                                        ScanPolarity::Negative;
                                }
                                _ => {
                                    description.polarity =
                                        ScanPolarity::Unknown;
                                }
                            }
                        } else {
                            description.polarity = ScanPolarity::Unknown;
                        }
                    }
                    curie!(MS:1000044) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            description.precursor[0]
                                .activation
                                .add_param(term.into_param(Value::Empty, Unit::Dimensionless));
                        } else {
                            description.precursor[0].activation.add_param(
                                term!(MS:1000044|dissociation method)
                                    .into_param(attr.value.clone().into(), Unit::Dimensionless),
                            );
                        }
                    }
                    curie!(MS:1003062) => {
                        description.index = attr.value.scalar().to_u64().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })? as usize;
                    }
                    curie!(MS:1003085) => {
                        description.precursor[0].ions[0].intensity = attr.value.scalar().to_f32().map_err(|v| {
                            BoxedError::new(
                                MzSpecLibErrorKind::Attribute,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            )
                        })?;
                    }
                    curie!(MS:1003086) => {
                        description.precursor[0].ions[0].add_param(attr.name.clone().into_param(attr.value.clone().into(), Unit::Unknown));
                    }
                    curie!(MS:1000028) => { // detector resolution
                        description.acquisition.scans[0].add_param(attr.name.clone().into_param(attr.value.clone().into(), Unit::Unknown));
                    }
                    | curie!(MS:1000285) // TIC 
                    | curie!(MS:1000505) // Base peak intensity
                    | curie!(MS:1000504) // Base peak m/z
                    | curie!(MS:1003059) // Number of peaks
                    | curie!(MS:1000528) // Lowest observed m/z
                    | curie!(MS:1000527) // Highest observed m/z
                    => (), // Ignore, can easily be calculated on the fly
                    | curie!(MS:1000002) // sample name
                    | curie!(MS:1000031) // instrument model
                    | curie!(MS:1000043) // intensity unit
                    | curie!(MS:1000443) // mass analyzer type
                    | curie!(MS:1003057) // scan number
                    | curie!(MS:1003299) // contributing replicate spectrum USI
                    | curie!(MS:1003065) // spectrum aggregation type
                    | curie!(MS:1003069) // number of replicate spectra available
                    | curie!(MS:1003070) // number of replicate spectra used
                    | curie!(MS:1003072) // spectrum origin type
                    | curie!(MS:1003203) // constituent spectrum file
                    | curie!(MS:1003213) // mass spectrometry acquisition method
                    => {
                        // These are metadata that could best live in the description
                        description.params.push(attr.clone().into());
                    }
                    _ => {
                        spectrum_attributes.push(attr.clone());
                    }
                }
            }
        }
    }
    Ok(())
}
