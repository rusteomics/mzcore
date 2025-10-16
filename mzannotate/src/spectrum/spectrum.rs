use mzcore::system::MassOverCharge;
use mzdata::{
    mzpeaks::{MZPeakSetType, peak_set::PeakSetIter, prelude::*},
    params::{ParamDescribed, ParamLike, Unit, Value, ValueRef},
    prelude::{IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::{SignalContinuity, SpectrumDescription},
};

use crate::{
    fragment::Fragment,
    mzspeclib::{Analyte, Attribute, Id, Interpretation, into_attributes},
    spectrum::AnnotatedPeak,
    term,
};

/// An annotated spectrum.
#[derive(Default, Debug, Clone)]
pub struct AnnotatedSpectrum {
    /// The Id for a spectrum
    pub key: Id,
    /// The spectrum description
    pub description: SpectrumDescription,
    /// The attributes that could not be interpretated as elements on the spectrum description
    pub attributes: Vec<Attribute>,
    /// The analytes
    pub analytes: Vec<Analyte>,
    /// The interpretations
    pub interpretations: Vec<Interpretation>,
    /// The spectrum itself
    pub peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
}

impl AnnotatedSpectrum {
    /// Create a new annotated spectrum
    pub const fn new(
        key: Id,
        description: SpectrumDescription,
        attributes: Vec<Attribute>,
        analytes: Vec<Analyte>,
        interpretations: Vec<Interpretation>,
        peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
    ) -> Self {
        Self {
            key,
            description,
            attributes,
            analytes,
            interpretations,
            peaks,
        }
    }
}

impl SpectrumLike<AnnotatedPeak<Fragment>, mzdata::mzpeaks::DeconvolutedPeak>
    for AnnotatedSpectrum
{
    fn description(&self) -> &SpectrumDescription {
        &self.description
    }

    fn description_mut(&mut self) -> &mut SpectrumDescription {
        &mut self.description
    }

    fn peaks(
        &'_ self,
    ) -> mzdata::spectrum::RefPeakDataLevel<
        '_,
        AnnotatedPeak<Fragment>,
        mzdata::mzpeaks::DeconvolutedPeak,
    > {
        mzdata::spectrum::RefPeakDataLevel::Centroid(&self.peaks)
    }

    fn raw_arrays(&'_ self) -> Option<&'_ mzdata::spectrum::BinaryArrayMap> {
        None
    }

    fn into_peaks_and_description(
        self,
    ) -> (
        mzdata::spectrum::PeakDataLevel<AnnotatedPeak<Fragment>, mzdata::mzpeaks::DeconvolutedPeak>,
        SpectrumDescription,
    ) {
        (
            mzdata::spectrum::PeakDataLevel::Centroid(self.peaks),
            self.description,
        )
    }
}

// TODO: this now misses most of the params in the other storage facilities (analytes, attributes, etc)
impl ParamDescribed for AnnotatedSpectrum {
    fn params(&self) -> &[mzdata::params::Param] {
        <SpectrumDescription as ParamDescribed>::params(&self.description)
    }

    fn params_mut(&mut self) -> &mut mzdata::params::ParamList {
        <SpectrumDescription as ParamDescribed>::params_mut(&mut self.description)
    }
}

#[allow(clippy::fallible_impl_from)] // Cannot fail, but cannot be proven to the compiler. 
// Basically the Attribute::unit has to be fallible to handle Unit::Unknown, but the units are
// known at compile time so this cannot happen.
impl<S: SpectrumLike> From<S> for AnnotatedSpectrum {
    fn from(value: S) -> Self {
        let mut this = Self {
            key: (value.index() + 1) as Id,
            description: value.description().clone(),
            ..Default::default()
        };

        let mut group_id = this
            .attributes
            .iter()
            .map(|v| v.group_id)
            .reduce(|prev, new| {
                new.map(|new| prev.map_or(new, |prev| prev.max(new)))
                    .or(prev)
            })
            .unwrap_or_default()
            .unwrap_or_default()
            + 1;

        let handle_param = |param: &mzdata::Param, this: &mut Self, group_id: &mut u32| {
            if matches!(param.value(), ValueRef::Empty) {
                // Ignore empty params
            } else if let Ok(term) = param.clone().try_into() {
                if let Some(unit) = Attribute::unit(param.unit, Some(*group_id)) {
                    this.attributes.push(Attribute::new(
                        term,
                        param.value.clone(),
                        Some(*group_id),
                    ));
                    this.attributes.push(unit);
                    *group_id += 1;
                } else {
                    this.attributes
                        .push(Attribute::new(term, param.value.clone(), None));
                }
            } else {
                this.attributes.push(Attribute::new(
                    term!(MS:1003275|other attribute name),
                    Value::String(param.name().to_string()),
                    Some(*group_id),
                ));
                this.attributes.push(Attribute::new(
                    term!(MS:1003276|other attribute value),
                    param.value.clone(),
                    Some(*group_id),
                ));
                *group_id += 1;
            }
        };

        for scan in value.acquisition().iter() {
            if let Some(filter_string) = scan.filter_string() {
                let attr = Attribute::new(
                    term!(MS:1000512|filter string),
                    Value::String(filter_string.to_string()),
                    None,
                );
                this.attributes.push(attr);
            }

            let attr = Attribute::new(
                term!(MS:1000016|scan start time),
                Value::Float(scan.start_time),
                Some(group_id),
            );
            let unit = Attribute::unit(Unit::Minute, Some(group_id)).unwrap();
            this.attributes.push(attr);
            this.attributes.push(unit);
            group_id += 1;

            if scan.injection_time != 0.0 {
                this.attributes.push(Attribute::new(
                    term!(MS:1000927|ion injection time),
                    Value::Float(f64::from(scan.injection_time)),
                    Some(group_id),
                ));
                this.attributes
                    .push(Attribute::unit(Unit::Millisecond, Some(group_id)).unwrap());
                group_id += 1;
            }

            for param in scan.params() {
                handle_param(param, &mut this, &mut group_id);
            }
        }
        if let Some(precursor) = value.precursor() {
            let ion = precursor.ion();
            this.attributes.push(Attribute::new(
                term!(MS:1000744|selected ion m/z),
                Value::Float(ion.mz()),
                Some(group_id),
            ));
            this.attributes
                .push(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
            group_id += 1;
            if let Some(z) = ion.charge() {
                this.attributes.push(Attribute::new(
                    term!(MS:1000041|charge state),
                    Value::Int(i64::from(z)),
                    None,
                ));
            }
            this.attributes.push(Attribute::new(
                term!(MS:1000042|peak intensity),
                Value::Float(f64::from(ion.intensity)),
                Some(group_id),
            ));
            this.attributes
                .push(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
            group_id += 1;

            for param in ion.params() {
                handle_param(param, &mut this, &mut group_id);
            }
        }

        if value.signal_continuity() == SignalContinuity::Centroid {
            this.peaks = value
                .peaks()
                .iter()
                .map(|v| {
                    AnnotatedPeak::new(
                        MassOverCharge::new::<mzcore::system::thomson>(v.mz()),
                        v.intensity(),
                        0,
                        Vec::new(),
                        Vec::new(),
                    )
                })
                .collect();
        }

        this
    }
}

// TODO: actually output the parameters that are part of the spectrum description
impl crate::mzspeclib::MzSpecLibEncode for AnnotatedSpectrum {
    /// The peak type
    type P = AnnotatedPeak<Fragment>;

    /// The key for this spectrum
    fn key(&self) -> Id {
        self.key
    }

    /// The attributes for this spectrum
    fn spectrum(&self) -> impl IntoIterator<Item = Attribute> {
        self.attributes
            .iter()
            .cloned()
            .chain(into_attributes(&self.description))
    }

    type AnalyteIter = Vec<Attribute>;
    /// The attributes for the analytes
    fn analytes(&self) -> impl Iterator<Item = (Id, Self::AnalyteIter)> {
        self.analytes.iter().map(|a| (a.id, a.attributes()))
    }

    /// The attributes for the interpretations
    fn interpretations(&self) -> impl Iterator<Item = (Id, &[Attribute])> {
        self.interpretations
            .iter()
            .map(|i| (i.id, i.attributes.as_slice()))
    }

    /// The peaks
    fn peaks(&self) -> PeakSetIter<'_, Self::P> {
        self.peaks.iter()
    }
}
