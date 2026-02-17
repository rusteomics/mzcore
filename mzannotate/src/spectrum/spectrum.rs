use mzcore::system::MassOverCharge;
use mzdata::{
    mzpeaks::{
        MZPeakSetType,
        peak_set::{PeakSetIter, PeakSetVec},
        prelude::*,
    },
    params::{ParamDescribed, ParamLike, Unit, Value, ValueRef},
    prelude::{IonProperties, SpectrumLike},
    spectrum::{SignalContinuity, SpectrumDescription},
};

use crate::{
    fragment::Fragment,
    mzspeclib::{
        Analyte, Attribute, Attributes, Id, Interpretation,
        get_attributes_from_spectrum_description, merge_attributes,
    },
    spectrum::AnnotatedPeak,
    term,
};

/// An annotated spectrum.
#[derive(Clone, Debug)]
pub struct AnnotatedSpectrum {
    /// The ID for a spectrum
    pub key: Id,
    /// The spectrum description
    pub description: SpectrumDescription,
    /// The attributes that could not be interpreted as elements on the spectrum description
    pub attributes: Attributes,
    /// The analytes
    pub analytes: Vec<Analyte>,
    /// The interpretations
    pub interpretations: Vec<Interpretation>,
    /// The spectrum itself
    pub peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
}

impl Default for AnnotatedSpectrum {
    fn default() -> Self {
        Self {
            key: 0,
            description: SpectrumDescription::default(),
            attributes: vec![Vec::new(); 1],
            analytes: Vec::new(),
            interpretations: Vec::new(),
            peaks: PeakSetVec::default(),
        }
    }
}

impl AnnotatedSpectrum {
    /// Create a new annotated spectrum
    pub const fn new(
        key: Id,
        description: SpectrumDescription,
        attributes: Attributes,
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

#[expect(clippy::fallible_impl_from)] // Cannot fail, but cannot be proven to the compiler.
// Basically the Attribute::unit has to be fallible to handle Unit::Unknown, but the units are
// known at compile time so this cannot happen.
impl<Spectrum: SpectrumLike> From<Spectrum> for AnnotatedSpectrum {
    fn from(value: Spectrum) -> Self {
        let mut this = Self {
            key: (value.index() + 1) as Id,
            description: value.description().clone(),
            ..Default::default()
        };

        let handle_param = |param: &mzdata::Param, this: &mut Self| {
            if matches!(param.value(), ValueRef::Empty) {
                // Ignore empty params
            } else if let Ok(term) = param.clone().try_into() {
                if let Some(unit) = Attribute::unit(param.unit) {
                    this.attributes
                        .push(vec![Attribute::new(term, param.value.clone()), unit]);
                } else {
                    this.attributes[0].push(Attribute::new(term, param.value.clone()));
                }
            } else {
                this.attributes.push(vec![
                    Attribute::new(
                        term!(MS:1003275|other attribute name),
                        Value::String(param.name().to_string()),
                    ),
                    Attribute::new(term!(MS:1003276|other attribute value), param.value.clone()),
                ]);
            }
        };

        for scan in value.acquisition().iter() {
            if let Some(filter_string) = scan.filter_string() {
                this.attributes[0].push(Attribute::new(
                    term!(MS:1000512|filter string),
                    Value::String(filter_string.to_string()),
                ));
            }

            let attr = Attribute::new(
                term!(MS:1000016|scan start time),
                Value::Float(scan.start_time),
            );
            let unit = Attribute::unit(Unit::Minute).unwrap();
            this.attributes.push(vec![attr, unit]);

            if scan.injection_time != 0.0 {
                this.attributes.push(vec![
                    Attribute::new(
                        term!(MS:1000927|ion injection time),
                        Value::Float(f64::from(scan.injection_time)),
                    ),
                    Attribute::unit(Unit::Millisecond).unwrap(),
                ]);
            }

            for param in scan.params() {
                handle_param(param, &mut this);
            }
        }
        if let Some(precursor) = value.precursor()
            && let Some(ion) = precursor.ions.first()
        {
            this.attributes.push(vec![
                Attribute::new(term!(MS:1000744|selected ion m/z), Value::Float(ion.mz())),
                Attribute::unit(Unit::MZ).unwrap(),
            ]);
            if let Some(z) = ion.charge() {
                this.attributes[0].push(Attribute::new(
                    term!(MS:1000041|charge state),
                    Value::Int(i64::from(z)),
                ));
            }
            this.attributes.push(vec![
                Attribute::new(
                    term!(MS:1000042|peak intensity),
                    Value::Float(f64::from(ion.intensity)),
                ),
                Attribute::unit(Unit::MZ).unwrap(),
            ]);

            for param in ion.params() {
                handle_param(param, &mut this);
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

impl crate::mzspeclib::MzSpecLibEncode for AnnotatedSpectrum {
    /// The peak type
    type Peak = AnnotatedPeak<Fragment>;

    /// The key for this spectrum
    fn key(&self) -> Id {
        self.key
    }

    /// The attributes for this spectrum
    fn spectrum(&self) -> Attributes {
        let mut own = self.attributes.clone();
        let other = get_attributes_from_spectrum_description(
            &self.description,
            Some(&{
                // The `SummaryOps` trait from mzdata is private, so this is the exact implementation copy pasted into here
                let state = self.peaks.iter().fold(
                    (
                        0.0,
                        mzdata::mzpeaks::CentroidPeak::default(),
                        (f64::INFINITY, f64::NEG_INFINITY),
                    ),
                    |(tic, bp, (mz_min, mz_max)), p| {
                        (
                            tic + p.intensity(),
                            if p.intensity() > bp.intensity {
                                p.as_centroid()
                            } else {
                                bp
                            },
                            (mz_min.min(p.mz()), mz_max.max(p.mz())),
                        )
                    },
                );

                mzdata::spectrum::SpectrumSummary::new(state.0, state.1, state.2, self.peaks.len())
            }),
        );
        merge_attributes(&mut own, &other);
        own
    }

    /// The attributes for the analytes
    fn analytes(&self) -> impl Iterator<Item = (Id, Attributes)> {
        self.analytes.iter().map(|a| (a.id, a.attributes()))
    }

    type InterpretationMemberIter = Vec<(Id, Attributes)>;

    /// The attributes for the interpretations
    fn interpretations(
        &self,
    ) -> impl Iterator<Item = (Id, Attributes, Self::InterpretationMemberIter)> {
        self.interpretations.iter().map(|i| {
            (
                i.id,
                i.attributes(),
                i.members
                    .iter()
                    .map(|(id, attributes)| (*id, attributes.clone()))
                    .collect(),
            )
        })
    }

    /// The peaks
    fn peaks(&self) -> PeakSetIter<'_, Self::Peak> {
        self.peaks.iter()
    }
}

impl mzcore::space::Space for AnnotatedSpectrum {
    fn space(&self) -> mzcore::space::UsedSpace {
        mzcore::space::UsedSpace::default() // TODO: this does not take any data into account because the Space trait cannot be applied to mzdata because of orphan rules
            .set_total::<Self>()
    }
}
