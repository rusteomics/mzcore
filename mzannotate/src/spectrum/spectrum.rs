use std::fmt::Display;

use mzcore::system::MassOverCharge;
use mzdata::{
    curie,
    mzpeaks::{MZPeakSetType, prelude::*},
    params::{ParamDescribed, ParamLike, ParamValue, Unit, Value, ValueRef},
    prelude::{IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::{SignalContinuity, SpectrumDescription},
};

use crate::{
    fragment::Fragment,
    mzspeclib::{
        Analyte, Attribute, Attributed, AttributedMut, Id, Interpretation, impl_attributed,
    },
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
    pub attributes: Vec<Attribute>,
    pub analytes: Vec<Analyte>,
    pub interpretations: Vec<Interpretation>,
    pub peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
}

impl_attributed!(mut AnnotatedSpectrum);

impl AnnotatedSpectrum {
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

    pub fn scan_number(&self) -> Option<usize> {
        let (_, v) = self.find_by_id(curie!(MS:1_003_057))?;
        v.value.scalar().to_i64().ok().map(|v| v as usize)
    }

    pub fn analyte(&self, id: Id) -> Option<&Analyte> {
        self.analytes.iter().find(|v| v.id == id)
    }

    pub fn add_analyte(&mut self, mut analyte: Analyte) {
        let k = (self.analytes.len() + 1) as u32;
        analyte.id = k;
        self.analytes.push(analyte);
    }

    pub fn remove_analyte(&mut self, id: Id) -> Option<Analyte> {
        let i = self.analytes.iter().position(|v| v.id == id)?;
        Some(self.analytes.remove(i))
    }

    pub fn interpretation(&self, id: Id) -> Option<&Interpretation> {
        self.interpretations.iter().find(|v| v.id == id)
    }

    pub fn add_interpretation(&mut self, mut interpretation: Interpretation) {
        let k = (self.interpretations.len() + 1) as u32;
        interpretation.id = k;
        self.interpretations.push(interpretation);
    }

    pub fn remove_interpretation(&mut self, id: Id) -> Option<Interpretation> {
        let i = self.interpretations.iter().position(|v| v.id == id)?;
        Some(self.interpretations.remove(i))
    }
}

impl Display for AnnotatedSpectrum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Spectrum={}>", self.key)?;
        writeln!(
            f,
            "MS:1003061|library spectrum name={}",
            self.description
                .title()
                .map_or_else(|| self.description.id.clone(), |v| v.to_string())
        )?;

        for attr in self.attributes() {
            writeln!(f, "{attr}")?;
        }
        for analyte in &self.analytes {
            write!(f, "{analyte}")?;
        }
        for interp in &self.interpretations {
            write!(f, "{interp}")?;
        }
        writeln!(f, "<Peaks>")?;
        for p in &self.peaks {
            writeln!(f, "{p}")?;
        }
        Ok(())
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

        let mut group_id = this.find_last_group_id().unwrap_or_default() + 1;

        let handle_param = |param: &mzdata::Param, this: &mut Self, group_id: &mut u32| {
            if matches!(param.value(), ValueRef::Empty) {
                // Ignore empty params
            } else if let Ok(term) = param.clone().try_into() {
                if let Some(unit) = Attribute::unit(param.unit, Some(*group_id)) {
                    this.add_attribute(Attribute::new(term, param.value.clone(), Some(*group_id)));
                    this.add_attribute(unit);
                    *group_id += 1;
                } else {
                    this.add_attribute(Attribute::new(term, param.value.clone(), None));
                }
            } else {
                this.add_attribute(Attribute::new(
                    term!(MS:1003275|"other attribute name"),
                    Value::String(param.name().to_string()),
                    Some(*group_id),
                ));
                this.add_attribute(Attribute::new(
                    term!(MS:1003276|"other attribute value"),
                    param.value.clone(),
                    Some(*group_id),
                ));
                *group_id += 1;
            }
        };

        for scan in value.acquisition().iter() {
            if let Some(filter_string) = scan.filter_string() {
                let attr = Attribute::new(
                    term!(MS:1_000_512|"filter string"),
                    Value::String(filter_string.to_string()),
                    None,
                );
                this.add_attribute(attr);
            }

            let attr = Attribute::new(
                term!(MS:1_000_016|"scan start time"),
                Value::Float(scan.start_time),
                Some(group_id),
            );
            let unit = Attribute::unit(Unit::Minute, Some(group_id)).unwrap();
            this.add_attribute(attr);
            this.add_attribute(unit);
            group_id += 1;

            if scan.injection_time != 0.0 {
                this.add_attribute_with_unit(
                    Attribute::new(
                        term!(MS:1_000_927|"ion injection time"),
                        Value::Float(f64::from(scan.injection_time)),
                        Some(group_id),
                    ),
                    Unit::Millisecond,
                    Some(group_id),
                );
                group_id += 1;
            }

            for param in scan.params() {
                handle_param(param, &mut this, &mut group_id);
            }
        }
        if let Some(precursor) = value.precursor() {
            let ion = precursor.ion();
            this.add_attribute(Attribute::new(
                term!(MS:1_000_744|"selected ion m/z"),
                Value::Float(ion.mz()),
                Some(group_id),
            ));
            this.add_attribute(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
            group_id += 1;
            if let Some(z) = ion.charge() {
                this.add_attribute(Attribute::new(
                    term!(MS:1_000_041|"charge state"),
                    Value::Int(i64::from(z)),
                    None,
                ));
            }
            this.add_attribute(Attribute::new(
                term!(MS:1_000_042|"peak intensity"),
                Value::Float(f64::from(ion.intensity)),
                Some(group_id),
            ));
            this.add_attribute(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
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
