use std::{collections::HashMap, fmt::Display};

use mzdata::{
    curie,
    mzpeaks::{MZPeakSetType, prelude::*},
    params::{ParamDescribed, ParamLike, ParamValue, Unit, Value, ValueRef},
    prelude::{IonProperties, PrecursorSelection, SpectrumLike},
    spectrum::SignalContinuity,
};

use crate::{
    fragment::Fragment,
    mzspeclib::{
        Attribute, AttributeSet, AttributeValue, Attributed, AttributedMut, EntryType, Term,
        impl_attributed,
    },
    term,
};

#[derive(Debug, Clone)]
pub struct LibraryHeader {
    pub format_version: String,
    pub attributes: Vec<Attribute>,
    pub attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
}

impl Default for LibraryHeader {
    fn default() -> Self {
        Self {
            format_version: "1.0".into(),
            attributes: Default::default(),
            attribute_classes: Default::default(),
        }
    }
}

impl_attributed!(mut LibraryHeader);

impl Display for LibraryHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let version = Attribute::new(
            Term::new(curie!(MS:1003186), "library format version".into()),
            AttributeValue::Scalar(Value::String(self.format_version.clone())),
            None,
        );
        writeln!(f, "<MzSpecLib>\n{version}")?;
        for attr in self.attributes.iter() {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl LibraryHeader {
    pub fn new(
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

#[derive(Debug, Clone)]
pub struct AnnotatedPeak<Annotation> {
    pub mz: f64,
    pub intensity: f32,
    pub index: mzdata::mzpeaks::IndexType,
    pub annotations: Vec<Annotation>,
    pub aggregations: Vec<String>,
}

impl<A> Default for AnnotatedPeak<A> {
    fn default() -> Self {
        Self {
            mz: f64::default(),
            intensity: f32::default(),
            index: mzdata::mzpeaks::IndexType::default(),
            annotations: Vec::new(),
            aggregations: Vec::new(),
        }
    }
}

impl<Annotation> AnnotatedPeak<Annotation> {
    pub fn new(
        mz: f64,
        intensity: f32,
        index: mzdata::mzpeaks::IndexType,
        annotations: Vec<Annotation>,
        aggregations: Vec<String>,
    ) -> Self {
        Self {
            mz,
            intensity,
            index,
            annotations,
            aggregations,
        }
    }
}

impl<B> AnnotatedPeak<B> {
    pub fn from<A>(value: AnnotatedPeak<A>) -> Self
    where
        B: From<A>,
    {
        Self {
            mz: value.mz,
            intensity: value.intensity,
            index: value.index,
            annotations: value.annotations.into_iter().map(From::from).collect(),
            aggregations: value.aggregations,
        }
    }
}

impl<B> AnnotatedPeak<B> {
    pub fn try_from<A>(value: AnnotatedPeak<A>) -> Result<Self, B::Error>
    where
        B: TryFrom<A>,
    {
        Ok(Self {
            mz: value.mz,
            intensity: value.intensity,
            index: value.index,
            annotations: value
                .annotations
                .into_iter()
                .map(B::try_from)
                .collect::<Result<Vec<_>, _>>()?,
            aggregations: value.aggregations,
        })
    }
}

mzdata::mzpeaks::implement_centroidlike!(AnnotatedPeak<Fragment>, true);

impl<A: Display> Display for AnnotatedPeak<A> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}\t{}", self.mz, self.intensity)?;

        if !self.annotations.is_empty() {
            f.write_str("\t")?;
            for (i, a) in self.annotations.iter().enumerate() {
                if i == 0 {
                    write!(f, "{a}")?;
                } else {
                    write!(f, ",{a}")?;
                }
            }
        }
        if !self.aggregations.is_empty() {
            f.write_str("\t")?;
            write!(f, "{}", self.aggregations.join(","))?;
        }
        Ok(())
    }
}

pub type IdType = u32;

#[derive(Default, Debug, Clone)]
pub struct Analyte {
    pub id: IdType,
    pub attributes: Vec<Attribute>,
}

impl_attributed!(mut Analyte);

impl Display for Analyte {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Analyte={}>", self.id)?;
        for attr in self.attributes.iter() {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl Analyte {
    pub fn new(id: IdType, attributes: Vec<Attribute>) -> Self {
        Self { id, attributes }
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
    pub fn new(
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
                    .map(|v| Value::Int((*v) as i64))
                    .collect(),
            );
            let mixture_ids = Attribute::new(
                Term::new(curie!(MS:1003163), "analyte mixture members".into()),
                val,
                None,
            );
            writeln!(f, "{mixture_ids}")?;
        }
        for attr in self.attributes.iter() {
            writeln!(f, "{attr}")?;
        }
        for member in self.members.iter() {
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
        for attr in self.attributes.iter() {
            f.write_str(attr.to_string().as_str())?;
        }
        Ok(())
    }
}

#[derive(Default, Debug, Clone)]
pub struct LibrarySpectrum {
    pub key: IdType,
    pub index: usize,
    pub name: Box<str>,
    pub attributes: Vec<Attribute>,
    pub analytes: Vec<Analyte>,
    pub interpretations: Vec<Interpretation>,
    pub peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
}

impl_attributed!(mut LibrarySpectrum);

impl LibrarySpectrum {
    pub fn new(
        key: IdType,
        index: usize,
        name: Box<str>,
        attributes: Vec<Attribute>,
        analytes: Vec<Analyte>,
        interpretations: Vec<Interpretation>,
        peaks: MZPeakSetType<AnnotatedPeak<Fragment>>,
    ) -> Self {
        Self {
            key,
            index,
            name,
            attributes,
            analytes,
            interpretations,
            peaks,
        }
    }

    pub fn scan_number(&self) -> Option<usize> {
        let (_, v) = self.find_by_id(curie!(MS:1003057))?;
        v.value.scalar().to_i64().ok().map(|v| v as usize)
    }

    pub fn analyte(&self, id: IdType) -> Option<&Analyte> {
        self.analytes.iter().find(|v| v.id == id)
    }

    pub fn add_analyte(&mut self, mut analyte: Analyte) {
        let k = (self.analytes.len() + 1) as u32;
        analyte.id = k;
        self.analytes.push(analyte);
    }

    pub fn remove_analyte(&mut self, id: IdType) -> Option<Analyte> {
        let i = self.analytes.iter().position(|v| v.id == id)?;
        Some(self.analytes.remove(i))
    }

    pub fn interpretation(&self, id: IdType) -> Option<&Interpretation> {
        self.interpretations.iter().find(|v| v.id == id)
    }

    pub fn add_interpretation(&mut self, mut interpretation: Interpretation) {
        let k = (self.interpretations.len() + 1) as u32;
        interpretation.id = k;
        self.interpretations.push(interpretation);
    }

    pub fn remove_interpretation(&mut self, id: IdType) -> Option<Interpretation> {
        let i = self.interpretations.iter().position(|v| v.id == id)?;
        Some(self.interpretations.remove(i))
    }
}

impl Display for LibrarySpectrum {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Spectrum={}>", self.key)?;
        writeln!(f, "MS:1003061|library spectrum name={}", self.name)?;
        for attr in self.attributes() {
            writeln!(f, "{attr}")?;
        }
        for analyte in self.analytes.iter() {
            write!(f, "{analyte}")?;
        }
        for interp in self.interpretations.iter() {
            write!(f, "{interp}")?;
        }
        writeln!(f, "<Peaks>")?;
        for p in self.peaks.iter() {
            writeln!(f, "{p}")?;
        }
        Ok(())
    }
}

impl From<mzdata::Spectrum> for LibrarySpectrum {
    fn from(value: mzdata::Spectrum) -> Self {
        let mut this = Self::default();
        this.key = (value.index() + 1) as IdType;
        this.index = value.index();
        this.name = value
            .description()
            .title()
            .map(|v| v.to_string())
            .unwrap_or_else(|| value.id().to_string())
            .into_boxed_str();

        let mut group_id = this.find_last_group_id().unwrap_or_default() + 1;

        macro_rules! handle_param {
            ($param:expr) => {
                if matches!($param.value(), ValueRef::Empty) {
                    todo!(
                        "Don't know how to convert value-less params yet: {:?}",
                        $param
                    )
                } else {
                    if $param.is_controlled() {
                        if matches!($param.unit, Unit::Unknown) {
                            this.add_attribute(Attribute::new(
                                Term::new(
                                    $param.curie().unwrap(),
                                    $param.name.clone().into_boxed_str(),
                                ),
                                $param.value.clone(),
                                None,
                            ))
                        } else {
                            this.add_attribute(Attribute::new(
                                Term::new(
                                    $param.curie().unwrap(),
                                    $param.name.clone().into_boxed_str(),
                                ),
                                $param.value.clone(),
                                Some(group_id),
                            ));
                            this.add_attribute(
                                Attribute::unit($param.unit, Some(group_id)).unwrap(),
                            );
                            group_id += 1;
                        }
                    } else {
                        this.add_attribute(
                            Attribute::new(
                                term!(MS:1003275|"other attribute name"),
                                Value::String($param.name().to_string()),
                                Some(group_id)
                            )
                        );
                        this.add_attribute(
                            Attribute::new(
                                term!(MS:1003276|"other attribute value"),
                                $param.value.clone(),
                                Some(group_id)
                            )
                        );
                        group_id += 1;
                    }
                }
            };
        }

        for scan in value.acquisition().iter() {
            if let Some(filter_string) = scan.filter_string() {
                let attr = Attribute::new(
                    term!(MS:1000512|"filter string"),
                    Value::String(filter_string.to_string()),
                    None,
                );
                this.add_attribute(attr);
            }

            let attr = Attribute::new(
                term!(MS:1000016|"scan start time"),
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
                        term!(MS:1000927|"ion injection time"),
                        Value::Float(scan.injection_time as f64),
                        Some(group_id),
                    ),
                    Unit::Millisecond,
                    Some(group_id),
                );
                group_id += 1;
            }

            for param in scan.params() {
                handle_param!(param);
            }
        }
        if let Some(precursor) = value.precursor() {
            let ion = precursor.ion();
            this.add_attribute(Attribute::new(
                term!(MS:1000744|"selected ion m/z"),
                Value::Float(ion.mz()),
                Some(group_id),
            ));
            this.add_attribute(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
            group_id += 1;
            if let Some(z) = ion.charge() {
                this.add_attribute(Attribute::new(
                    term!(MS:1000041|"charge state"),
                    Value::Int(z as i64),
                    None,
                ));
            }
            this.add_attribute(Attribute::new(
                term!(MS:1000042|"peak intensity"),
                Value::Float(ion.intensity as f64),
                Some(group_id),
            ));
            this.add_attribute(Attribute::unit(Unit::MZ, Some(group_id)).unwrap());
            group_id += 1;

            for param in ion.params() {
                handle_param!(param);
            }
        }

        if value.signal_continuity() == SignalContinuity::Centroid {
            this.peaks = value
                .peaks()
                .iter()
                .map(|v| {
                    AnnotatedPeak::new(
                        v.mz(),
                        v.intensity(),
                        0,
                        Default::default(),
                        Default::default(),
                    )
                })
                .collect();
        }

        this
    }
}
