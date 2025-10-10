use std::{
    collections::{HashMap, VecDeque},
    fmt::Display,
    io::{self, prelude::*},
};

use context_error::{BasicKind, BoxedError, Context, CreateError};
use indexmap::IndexMap;
use mzcore::{
    prelude::*,
    system::{MassOverCharge, isize::Charge},
};
use mzdata::{
    curie,
    mzpeaks::{peak_set::PeakSetVec, prelude::PeakCollectionMut},
    params::{ParamDescribed, ParamValue, Unit, Value},
    spectrum::{
        IsolationWindowState, Precursor, ScanEvent, ScanWindow, SelectedIon, SpectrumDescription,
    },
};

use crate::{
    fragment::{Fragment, parse_mz_paf},
    helper_functions::explain_number_error,
    mzspeclib::{
        Analyte, Attribute, AttributeParseError, AttributeSet, AttributeValue, Attributed,
        AttributedMut, EntryType, Id, Interpretation, LibraryHeader,
    },
    spectrum::{AnnotatedPeak, AnnotatedSpectrum},
    term,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ParserState {
    Initial,
    Header,
    AttributeSet,
    Spectrum,
    Analyte,
    Interpretation,
    InterpretationMember,
    Peaks,
    Between,
    EOF,
}

impl Display for ParserState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

#[derive(Debug)]
pub enum MzSpecLibTextParseError {
    IOError(io::Error),
    InvalidContentAtState(String, ParserState, u64),
    AttributeParseError(AttributeParseError, ParserState, u64),
    EOF,
    InconsistentCharge,
    MissingUnit,
    InvalidProforma(BoxedError<'static, BasicKind>),
    InvalidMzPAF(BoxedError<'static, BasicKind>),
    RichError(BoxedError<'static, BasicKind>),
}

// impl From<MzSpecLibTextParseError> for io::Error {
//     fn from(value: MzSpecLibTextParseError) -> Self {
//         match value {
//             MzSpecLibTextParseError::IOError(error) => error,
//             e => io::Error::new(io::ErrorKind::Other, e),
//         }
//     }
// }

#[derive(Debug)]
pub struct MzSpecLibParser<R: Read> {
    inner: io::BufReader<R>,
    header: LibraryHeader,
    state: ParserState,
    line_cache: VecDeque<String>,
    line_number: u64,
    offsets: LibraryIndex,
    last_error: bool,
}

impl<R: Read> MzSpecLibParser<R> {
    pub fn new(reader: io::BufReader<R>) -> Result<Self, MzSpecLibTextParseError> {
        let mut this = Self {
            inner: reader,
            header: LibraryHeader::default(),
            state: ParserState::Initial,
            line_cache: VecDeque::new(),
            line_number: 0,
            offsets: LibraryIndex::default(),
            last_error: false,
        };
        this.read_header()?;
        Ok(this)
    }

    pub fn header(&self) -> &LibraryHeader {
        &self.header
    }

    pub fn header_mut(&mut self) -> &mut LibraryHeader {
        &mut self.header
    }

    pub fn line_number(&self) -> u64 {
        self.line_number
    }

    fn push_back_line(&mut self, line: String) {
        self.line_cache.push_front(line);
        self.line_number -= 1;
    }

    fn read_next_line(&mut self, buf: &mut String) -> io::Result<usize> {
        buf.clear();
        if self.line_cache.is_empty() {
            match self.inner.read_line(buf) {
                Ok(mut z) => {
                    if z != 0 {
                        self.line_number += 1;
                    } else {
                        self.state = ParserState::EOF;
                        return Ok(0);
                    }

                    while buf.trim().is_empty() || buf.starts_with('#') {
                        let z1 = self.read_next_line(buf)?;
                        if z1 == 0 {
                            return Ok(z1);
                        }
                        z += z1;
                    }
                    Ok(z)
                }
                Err(e) => Err(e),
            }
        } else {
            *buf = self.line_cache.pop_front().unwrap();
            self.line_number += 1;
            Ok(buf.len())
        }
    }

    fn read_attribute(&mut self, buf: &mut String) -> Result<Attribute, MzSpecLibTextParseError> {
        let z = self
            .read_next_line(buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        buf.trim_ascii_end().parse().map_err(|e| {
            MzSpecLibTextParseError::AttributeParseError(e, self.state, self.line_number)
        })
    }

    fn read_attribute_sets(&mut self) -> Result<(), MzSpecLibTextParseError> {
        let mut buf = String::new();
        let mut current_attribute_set: Option<AttributeSet> = None;

        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => {
                    current_attribute_set
                        .as_mut()
                        .unwrap()
                        .attributes
                        .entry(attr.group_id)
                        .or_default()
                        .push((
                            attr,
                            Context::none()
                                .line_index(self.line_number as u32)
                                .lines(0, buf.clone()),
                        ));
                }
                Err(e) => {
                    if buf.starts_with('<') {
                        if buf.starts_with("<Spectrum=") {
                            if current_attribute_set.is_some() {
                                let set = current_attribute_set.take().unwrap();
                                self.header
                                    .attribute_classes
                                    .entry(set.namespace)
                                    .or_default()
                                    .push(set);
                            }
                            self.state = ParserState::Spectrum;
                            self.push_back_line(buf);
                            break;
                        } else if buf.starts_with("<AttributeSet") {
                            self.state = ParserState::AttributeSet;
                            if let Some((_, rest)) = buf.trim().split_once(' ') {
                                if !rest.ends_with('>') {
                                    return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                        buf,
                                        self.state,
                                        self.line_number,
                                    ));
                                }
                                if let Some((entry_tp, id)) = rest[..rest.len() - 1].split_once('=')
                                {
                                    if current_attribute_set.is_some() {
                                        let set = current_attribute_set.take().unwrap();
                                        self.header
                                            .attribute_classes
                                            .entry(set.namespace)
                                            .or_default()
                                            .push(set);
                                    }
                                    let set_entry_tp = match entry_tp {
                                        "Spectrum" => EntryType::Spectrum,
                                        "Analyte" => EntryType::Analyte,
                                        "Interpretation" => EntryType::Interpretation,
                                        "Cluster" => EntryType::Cluster,
                                        _ => {
                                            return Err(
                                                MzSpecLibTextParseError::InvalidContentAtState(
                                                    buf,
                                                    self.state,
                                                    self.line_number,
                                                ),
                                            );
                                        }
                                    };
                                    let set_id = id.to_string();
                                    current_attribute_set = Some(AttributeSet::new(
                                        set_id,
                                        set_entry_tp,
                                        HashMap::new(),
                                    ));
                                } else {
                                    return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                        buf,
                                        self.state,
                                        self.line_number,
                                    ));
                                }
                            } else {
                                return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                    buf,
                                    self.state,
                                    self.line_number,
                                ));
                            }
                        } else {
                            return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                buf,
                                self.state,
                                self.line_number,
                            ));
                        }
                    } else if buf.is_empty() {
                        // continue
                    } else {
                        return Err(e);
                    }
                }
            }
        }
        Ok(())
    }

    fn read_header(&mut self) -> Result<(), MzSpecLibTextParseError> {
        let mut buf = String::new();
        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }

        if !buf.starts_with("<mzSpecLib>") {
            return Err(MzSpecLibTextParseError::InvalidContentAtState(
                buf,
                ParserState::Initial,
                self.line_number,
            ));
        }
        self.state = ParserState::Header;

        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => {
                    // TODO: switch to assigning these basic attributes to header fields
                    match attr.name.accession {
                        curie!(MS:1003186) => {
                            self.header.format_version = attr.value.to_string();
                        }
                        _ => {
                            self.header.add_attribute(attr);
                        }
                    }
                }
                Err(e) => {
                    if buf.starts_with('<') {
                        if buf.starts_with("<Spectrum=") {
                            self.state = ParserState::Spectrum;
                        } else if buf.starts_with("<AttributeSet") {
                            self.state = ParserState::AttributeSet;
                        } else {
                            return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                buf,
                                self.state,
                                self.line_number,
                            ));
                        }

                        self.push_back_line(buf);
                        break;
                    } else if buf.is_empty() {
                        // continue
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        if !matches!(self.state, ParserState::AttributeSet) {
            return Ok(());
        }

        self.read_attribute_sets()?;

        Ok(())
    }

    /// Parse an open declaration tag. Example: `<Spectrum=1>`
    /// # Errors
    /// * If the open declaration tag is not the first thing in the buffer.
    /// * If the declaration does not contain a '='.
    /// * If the declaration is not followed by a valid integer.
    fn parse_open_declaration(
        &mut self,
        buf: &str,
        declaration: &str,
        to_state: ParserState,
    ) -> Result<u32, MzSpecLibTextParseError> {
        let id = if buf.starts_with(declaration) {
            self.state = to_state;
            if let Some((before, val)) = buf[..buf.trim_ascii_end().len() - 1].split_once('=') {
                val.trim().parse::<Id>().map_err(|e| {
                    MzSpecLibTextParseError::RichError(BoxedError::new(
                        BasicKind::Error,
                        "Invalid declaration",
                        format!("The ID was {}", explain_number_error(&e)),
                        Context::line(
                            Some(self.line_number as u32),
                            buf,
                            before.len() + 1,
                            val.len(),
                        )
                        .to_owned(),
                    ))
                })?
            } else {
                return Err(MzSpecLibTextParseError::RichError(BoxedError::new(
                    BasicKind::Error,
                    "Invalid declaration",
                    "The ID needs to contain an equals symbol `=` and this was missing",
                    Context::full_line(self.line_number as u32, buf).to_owned(),
                )));
            }
        } else {
            return Err(MzSpecLibTextParseError::RichError(BoxedError::new(
                BasicKind::Error,
                "Invalid declaration",
                format!("The declaration `{declaration}` was expected but not found"),
                Context::full_line(self.line_number as u32, buf).to_owned(),
            )));
        };
        Ok(id)
    }

    /// Parse an analyte from the stream.
    /// # Errors
    /// If the analyte
    fn read_analyte(&mut self) -> Result<Analyte, MzSpecLibTextParseError> {
        let mut buf = String::new();
        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        let id = self.parse_open_declaration(&buf, "<Analyte=", ParserState::Analyte)?;

        let mut analyte = Analyte::new(id, None, None, Vec::new());
        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => {
                    if attr.name.accession == curie!(MS:1_003_270)
                        && let Value::String(value) = attr.value.scalar().as_ref()
                    {
                        let peptidoform_ion = PeptidoformIon::pro_forma(value, None)
                            .map_err(|e| MzSpecLibTextParseError::InvalidProforma(e.to_owned()))?;
                        if let Some(charge) = peptidoform_ion.get_charge_carriers() {
                            let charge = charge.charge();
                            if analyte.charge.is_some_and(|c| c != charge) {
                                return Err(MzSpecLibTextParseError::InconsistentCharge);
                            } else if analyte.charge.is_none() {
                                analyte.charge = Some(charge);
                            }
                        }
                        analyte.peptidoform_ion = Some(peptidoform_ion);
                    } else if attr.name.accession == curie!(MS:1_000_041)
                        && let Value::Int(value) = attr.value.scalar().as_ref()
                    {
                        let charge = Charge::new::<mzcore::system::e>(*value as isize);
                        if analyte.charge.is_some_and(|c| c != charge) {
                            return Err(MzSpecLibTextParseError::InconsistentCharge);
                        }
                        if let Some(p) = analyte.peptidoform_ion.as_mut()
                            && p.get_charge_carriers().is_none()
                        {
                            p.set_charge_carriers(Some(MolecularCharge::proton(charge.value)));
                        }
                        analyte.charge = Some(charge);
                    } else {
                        analyte.add_attribute(attr);
                    }
                }
                Err(e) => {
                    if buf.starts_with('<') {
                        self.push_back_line(buf);
                        break;
                    } else if buf.trim().is_empty() {
                        // continue;
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let mut next_group_id = analyte.find_last_group_id().unwrap_or_default();
        let mut last_group_id = 0;

        let attr_set_name = term!(MS:1_003_212|"library attribute set name");
        let attr_sets: Vec<_> = analyte
            .find_all(&attr_set_name)
            .map(|v| v.value.to_string())
            .collect();
        for name in attr_sets {
            for attr_set in self
                .header
                .attribute_classes
                .get(&EntryType::Analyte)
                .into_iter()
                .flatten()
            {
                if attr_set.id == name || attr_set.id == "all" {
                    for (group_id, group) in &attr_set.attributes {
                        if group_id.is_some() {
                            next_group_id += 1;

                            for (attr, _) in group {
                                let mut attr = attr.clone();
                                attr.group_id = Some(next_group_id);
                                analyte.add_attribute(attr);
                            }
                        } else {
                            analyte
                                .attributes
                                .extend(group.iter().map(|(a, _)| a).cloned());
                        }
                    }
                }
            }
        }
        Ok(analyte)
    }

    fn read_interpretation(&mut self) -> Result<Interpretation, MzSpecLibTextParseError> {
        let mut buf = String::new();
        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        let id =
            self.parse_open_declaration(&buf, "<Interpretation=", ParserState::Interpretation)?;
        let mut interp = Interpretation::new(id, Vec::new(), Vec::new(), Vec::new());
        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => {
                    interp.add_attribute(attr);
                }
                Err(e) => {
                    if buf.starts_with("<InterpretationMember=") {
                        self.push_back_line(buf);
                        todo!();
                    } else if buf.starts_with('<') {
                        self.push_back_line(buf);
                        break;
                    } else if buf.trim().is_empty() {
                        // continue;
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let mut next_group_id = interp.find_last_group_id().unwrap_or_default();
        let mut last_group_id = 0;

        let attr_set_name = term!(MS:1003212|"library attribute set name");
        let attr_sets: Vec<_> = interp
            .find_all(&attr_set_name)
            .map(|v| v.value.to_string())
            .collect();
        for name in attr_sets {
            for attr_set in self
                .header
                .attribute_classes
                .get(&EntryType::Interpretation)
                .into_iter()
                .flatten()
            {
                if attr_set.id == name || attr_set.id == "all" {
                    for (group_id, group) in &attr_set.attributes {
                        if group_id.is_some() {
                            next_group_id += 1;

                            for (attr, _) in group {
                                let mut attr = attr.clone();
                                attr.group_id = Some(next_group_id);
                                interp.add_attribute(attr);
                            }
                        } else {
                            interp
                                .attributes
                                .extend(group.iter().map(|(a, _)| a).cloned());
                        }
                    }
                }
            }
        }
        Ok(interp)
    }

    fn parse_peak_line(
        &self,
        buf: &str,
    ) -> Result<AnnotatedPeak<Fragment>, MzSpecLibTextParseError> {
        let buf = buf.trim();
        let mut it = buf.split('\t');

        let mz = it
            .next()
            .ok_or_else(|| {
                MzSpecLibTextParseError::RichError(BoxedError::new(
                    BasicKind::Error,
                    "Peak m/z is missing",
                    "At least two columns are neccessary in the peaks data lines",
                    Context::full_line(self.line_number as u32, buf).to_owned(),
                ))
            })
            .and_then(|v| {
                v.parse::<f64>().map_err(|e| {
                    MzSpecLibTextParseError::RichError(BoxedError::new(
                        BasicKind::Error,
                        "Peak m/z is not a number",
                        e.to_string(),
                        Context::full_line(self.line_number as u32, buf).to_owned(),
                    ))
                })
            })?;

        let intensity = it
            .next()
            .ok_or_else(|| {
                MzSpecLibTextParseError::RichError(BoxedError::new(
                    BasicKind::Error,
                    "Peak intensity is missing",
                    "At least two columns are neccessary in the peaks data lines",
                    Context::full_line(self.line_number as u32, buf).to_owned(),
                ))
            })
            .and_then(|v| {
                v.parse::<f32>().map_err(|e| {
                    MzSpecLibTextParseError::RichError(BoxedError::new(
                        BasicKind::Error,
                        "Peak intensity is not a number",
                        e.to_string(),
                        Context::full_line(self.line_number as u32, buf).to_owned(),
                    ))
                })
            })?;

        let mut peak = AnnotatedPeak::new(
            MassOverCharge::new::<mzcore::system::thomson>(mz),
            intensity,
            0,
            Vec::new(),
            Vec::new(),
        );

        match it.next() {
            Some(v) => {
                let v = v.trim();
                if !v.is_empty() && v != "?" {
                    let annots = parse_mz_paf(v, None, &CompoundPeptidoformIon::default())
                        .map_err(|e| MzSpecLibTextParseError::InvalidMzPAF(e.to_owned()))?;
                    peak.annotations = annots;
                }
            }
            None => return Ok(peak),
        }

        match it.next() {
            Some(v) => {
                peak.aggregations = v.split(',').map(ToString::to_string).collect();
            }
            None => return Ok(peak),
        }

        Ok(peak)
    }

    fn read_spectrum(&mut self) -> Result<AnnotatedSpectrum, MzSpecLibTextParseError> {
        let mut buf = String::new();

        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        let id = self.parse_open_declaration(&buf, "<Spectrum=", ParserState::Spectrum)?;
        let mut spec = AnnotatedSpectrum::new(
            id,
            SpectrumDescription::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            PeakSetVec::default(),
        );

        // TODO: do better but should work for now
        spec.description
            .acquisition
            .scans
            .push(ScanEvent::default());
        spec.description.acquisition.scans[0]
            .scan_windows
            .push(ScanWindow::default());
        spec.description.precursor.push(Precursor::default());
        spec.description.precursor[0]
            .ions
            .push(SelectedIon::default());
        spec.description.signal_continuity = mzdata::spectrum::SignalContinuity::Centroid;

        let mut term_collection: HashMap<Option<u32>, Vec<(Attribute, Context<'static>)>> =
            HashMap::new();

        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => term_collection.entry(attr.group_id).or_default().push((
                    attr,
                    Context::none()
                        .line_index(self.line_number as u32)
                        .lines(0, buf.clone()),
                )),
                Err(e) => {
                    if buf.starts_with("<Analyte=") {
                        self.push_back_line(buf.clone());
                        let analyte = self.read_analyte()?;
                        spec.analytes.push(analyte);
                    } else if buf.starts_with("<Interpretation=") {
                        self.push_back_line(buf.clone());
                        let interp = self.read_interpretation()?;
                        spec.interpretations.push(interp);
                    } else if buf.starts_with("<Peaks>") {
                        self.state = ParserState::Peaks;
                        break;
                    } else if buf.trim().is_empty() {
                        // continue
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let attr_set_name = term!(MS:1003212|"library attribute set name");

        let set_names: Vec<_> = spec
            .find_all(&attr_set_name)
            .map(|v| v.value.to_string())
            .collect();

        parse_spectrum_attributes(term_collection.iter(), &mut spec)?;
        if let Some(sets) = self.header.attribute_classes.get(&EntryType::Spectrum) {
            parse_spectrum_attributes(
                sets.iter()
                    .filter_map(|attr_set| {
                        (set_names.contains(&attr_set.id) || attr_set.id == "all")
                            .then_some(attr_set.attributes.iter())
                    })
                    .flatten(),
                &mut spec,
            )?;
        }

        spec.description.precursor[0]
            .activation
            ._extract_methods_from_params();

        if matches!(self.state, ParserState::Peaks) {
            loop {
                let z = self
                    .read_next_line(&mut buf)
                    .map_err(MzSpecLibTextParseError::IOError)?;
                if z == 0 {
                    break;
                }
                let buf_trimmed = buf.trim();
                if buf_trimmed.starts_with('<') {
                    self.push_back_line(buf);
                    self.state = ParserState::Between;
                    break;
                } else if buf_trimmed.is_empty() {
                    self.state = ParserState::Between;
                    break;
                }

                let peak = self.parse_peak_line(buf_trimmed)?;
                spec.peaks.push(peak);
            }
        }
        Ok(spec)
    }

    pub fn read_next(&mut self) -> Result<AnnotatedSpectrum, MzSpecLibTextParseError> {
        if self.last_error {
            // If the last element went awry scan the buffer until the start of the next spectrum
            // is found, otherwise it generates an error for each skipped line.

            let mut buf = String::new();

            loop {
                let z = self
                    .read_next_line(&mut buf)
                    .map_err(MzSpecLibTextParseError::IOError)?;
                if z == 0 {
                    return Err(MzSpecLibTextParseError::EOF);
                }
                if buf.starts_with("<Spectrum=") {
                    self.push_back_line(buf);
                    break;
                }
            }
        }
        let res = self.read_spectrum();
        self.last_error = res.is_err();
        res
    }

    pub fn inner(&self) -> &io::BufReader<R> {
        &self.inner
    }

    pub fn as_mut(&mut self) -> &mut io::BufReader<R> {
        &mut self.inner
    }

    pub fn into_inner(self) -> io::BufReader<R> {
        self.inner
    }
}

fn parse_spectrum_attributes<'a>(
    attributes: impl Iterator<Item = (&'a Option<u32>, &'a Vec<(Attribute, Context<'static>)>)>,
    spectrum: &mut AnnotatedSpectrum,
) -> Result<(), MzSpecLibTextParseError> {
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
            {
                spectrum.description.acquisition.scans[0].start_time =
                    value.scalar().to_f32().map_err(|v| {
                        MzSpecLibTextParseError::RichError(BoxedError::new(
                            BasicKind::Error,
                            "Invalid attribute",
                            v.to_string(),
                            context.clone(),
                        ))
                    })? as f64
                        / match unit.ok_or(MzSpecLibTextParseError::MissingUnit)? {
                            Unit::Minute => 60.0,
                            _ => 1.0, // Assume seconds for anything else
                        };
            } else if let Some((attr, context)) = group.iter().find(|a| {
                a.0.name.accession == curie!(MS:1000045) || a.0.name.accession == curie!(MS:1000509)
            }) {
                spectrum.description.precursor[0].activation.energy =
                    attr.value.scalar().to_f32().map_err(|v| {
                        MzSpecLibTextParseError::RichError(BoxedError::new(
                            BasicKind::Error,
                            "Invalid attribute",
                            v.to_string(),
                            context.clone(),
                        ))
                    })?;
            } else if let Some(unit) = unit
                && group.len() == 2
            {
                let (other, _) = group
                    .iter()
                    .find(|a| a.0.name.accession != curie!(UO:0000000))
                    .unwrap(); // This assumes that the other attribute is NOT a unit
                spectrum.description.add_param(
                    other
                        .name
                        .clone()
                        .to_param(other.value.clone().into(), unit),
                );
            } else {
                // TODO: duplicate ids could happen due to merging of multiple attribute sets
                for attr in group {
                    spectrum.add_attribute(attr.0.clone());
                }
            }
        } else {
            for (attr, context) in group {
                match attr.name.accession {
                    curie!(MS:1003208) => {
                        let window = &mut spectrum.description.precursor[0].isolation_window;
                        window.target = attr.value.scalar().to_f32().map_err(|v| {
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })?;
                        if matches!(window.flags, IsolationWindowState::Offset) {
                            window.lower_bound = window.target - window.lower_bound;
                            window.upper_bound += window.target;
                        }
                        window.flags = IsolationWindowState::Complete;
                    }
                    curie!(MS:1000828) => {
                        let offset = attr.value.scalar().to_f32().map_err(|v| {
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })?;
                        let window = &mut spectrum.description.precursor[0].isolation_window;
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
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })?;
                        let window = &mut spectrum.description.precursor[0].isolation_window;
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
                    curie!(MS:1000794) => {
                        let limit = attr.value.scalar().to_f32().map_err(|v| {
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })?;
                        let window = &mut spectrum.description.precursor[0].isolation_window;
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
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })?;
                        let window = &mut spectrum.description.precursor[0].isolation_window;
                        if matches!(
                            window.flags,
                            IsolationWindowState::Unknown | IsolationWindowState::Explicit
                        ) {
                            window.flags = IsolationWindowState::Explicit;
                            window.upper_bound = limit;
                        }
                    }
                    curie!(MS:1000501) => {
                        spectrum.description.acquisition.scans[0].scan_windows[0].lower_bound =
                            attr.value.scalar().to_f32().map_err(|v| {
                                MzSpecLibTextParseError::RichError(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                ))
                            })?;
                    }
                    curie!(MS:1000500) => {
                        spectrum.description.acquisition.scans[0].scan_windows[0].upper_bound =
                            attr.value.scalar().to_f32().map_err(|v| {
                                MzSpecLibTextParseError::RichError(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                ))
                            })?;
                    }
                    curie!(MS:1003063) => {
                        spectrum.description.add_param(
                            attr.name
                                .clone()
                                .to_param(attr.value.clone().into(), Unit::Dimensionless),
                        );
                    }
                    curie!(MS:1003061) => spectrum.description.id = attr.value.to_string(),
                    curie!(MS:1000511) => {
                        spectrum.description.ms_level =
                            u8::try_from(attr.value.scalar().to_u64().map_err(|v| {
                                MzSpecLibTextParseError::RichError(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid attribute",
                                    v.to_string(),
                                    context.clone(),
                                ))
                            })?)
                            .map_err(|_| {
                                MzSpecLibTextParseError::RichError(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid attribute",
                                    "The MS level is too large of a number",
                                    context.clone(),
                                ))
                            })?;
                    }
                    curie!(MS:1000512) => spectrum.description.acquisition.scans[0].add_param(
                        term!(MS:1000512|"filter string")
                            .to_param(attr.value.clone().into(), Unit::Dimensionless),
                    ),
                    curie!(MS:1000008) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            spectrum
                                .description
                                .add_param(term.to_param(Value::Empty, Unit::Dimensionless));
                        } else {
                            spectrum.description.add_param(
                                term!(MS:1000008|"ionization type")
                                    .to_param(attr.value.clone().into(), Unit::Dimensionless),
                            );
                        }
                    }
                    curie!(MS:1000465) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            match term.accession {
                                curie!(MS:1000130) => {
                                    spectrum.description.polarity =
                                        mzdata::spectrum::ScanPolarity::Positive;
                                }
                                curie!(MS:1000129) => {
                                    spectrum.description.polarity =
                                        mzdata::spectrum::ScanPolarity::Negative;
                                }
                                _ => {
                                    spectrum.description.polarity =
                                        mzdata::spectrum::ScanPolarity::Unknown;
                                }
                            }
                        } else {
                            spectrum.description.polarity = mzdata::spectrum::ScanPolarity::Unknown;
                        }
                    }
                    curie!(MS:1000044) => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            spectrum.description.precursor[0]
                                .activation
                                .add_param(term.to_param(Value::Empty, Unit::Dimensionless));
                        } else {
                            spectrum.description.precursor[0].activation.add_param(
                                term!(MS:1000044|"dissociation method")
                                    .to_param(attr.value.clone().into(), Unit::Dimensionless),
                            );
                        }
                    }
                    curie!(MS:1003062) => {
                        spectrum.description.index = attr.value.scalar().to_u64().map_err(|v| {
                            MzSpecLibTextParseError::RichError(BoxedError::new(
                                BasicKind::Error,
                                "Invalid attribute",
                                v.to_string(),
                                context.clone(),
                            ))
                        })? as usize;
                    }
                    _ => {
                        if let AttributeValue::Term(term) = attr.value.clone() {
                            spectrum
                                .description
                                .add_param(term.to_param(Value::Empty, Unit::Dimensionless));
                        } else {
                            spectrum.description.add_param(
                                attr.name
                                    .clone()
                                    .to_param(attr.value.clone().into(), Unit::Dimensionless),
                            );
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

impl<R: Read> Iterator for MzSpecLibParser<R> {
    type Item = Result<AnnotatedSpectrum, MzSpecLibTextParseError>;

    fn next(&mut self) -> Option<Self::Item> {
        if matches!(self.state, ParserState::EOF) {
            return None;
        }
        Some(self.read_next())
    }
}

#[derive(Debug, Default, Clone, PartialEq, Eq, PartialOrd)]
pub struct LibraryIndexRecord {
    pub offset: u64,
    pub key: Id,
    pub index: usize,
    pub line_number: u64,
    pub name: Box<str>,
}

#[derive(Debug, Default, Clone)]
pub struct LibraryIndex {
    pub init: bool,
    primary_index: IndexMap<usize, LibraryIndexRecord>,
    key_index: HashMap<Id, usize>,
    name_index: HashMap<Box<str>, usize>,
}

impl LibraryIndex {
    pub fn insert(&mut self, index: usize, record: LibraryIndexRecord) {
        self.key_index.insert(record.key, index);
        self.name_index.insert(record.name.clone(), index);
        self.primary_index.insert(index, record);
    }

    pub fn get_by_key(&self, key: Id) -> Option<&LibraryIndexRecord> {
        self.primary_index.get(self.key_index.get(&key)?)
    }

    pub fn get_by_index(&self, index: usize) -> Option<&LibraryIndexRecord> {
        self.primary_index.get(&index)
    }

    pub fn get_by_name(&self, name: &str) -> Option<&LibraryIndexRecord> {
        self.primary_index.get(self.name_index.get(name)?)
    }

    #[allow(unused)]
    pub fn iter(&self) -> indexmap::map::Iter<'_, usize, LibraryIndexRecord> {
        self.primary_index.iter()
    }

    pub fn len(&self) -> usize {
        self.primary_index.len()
    }

    pub fn is_empty(&self) -> bool {
        self.primary_index.is_empty()
    }
}

impl<R: Read + Seek> MzSpecLibParser<R> {
    pub fn build_index(&mut self) -> io::Result<()> {
        self.inner.rewind()?;
        self.line_cache.clear();
        let mut buf = String::new();
        let mut line_count = 0u64;
        let mut offset_index = LibraryIndex::default();
        let mut offset = 0;
        let mut current_record: Option<LibraryIndexRecord> = None;
        let mut index = 0;
        while let Ok(z) = self.read_next_line(&mut buf) {
            if z == 0 {
                break;
            }

            if buf.starts_with("<Spectrum=") {
                let (_, rest) = buf.trim().split_once('=').unwrap();
                if let Some(key) = rest.strip_suffix('>') {
                    if let Some(current_record) = current_record {
                        offset_index.insert(current_record.index, current_record);
                    }
                    current_record = Some(LibraryIndexRecord::default());
                    if let Ok(key) = key.parse::<Id>() {
                        let r = current_record.as_mut().unwrap();
                        r.index = index;
                        r.key = key;
                        r.offset = offset;
                        r.line_number = line_count;
                        index += 1;
                    }
                }
            }
            line_count += 1;
            offset += z as u64;
        }
        if let Some(current_record) = current_record {
            offset_index.insert(current_record.index, current_record);
        }
        offset_index.init = true;
        self.offsets = offset_index;
        self.inner.rewind()?;
        Ok(())
    }

    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    pub fn iter(
        &mut self,
    ) -> impl Iterator<Item = Result<AnnotatedSpectrum, MzSpecLibTextParseError>> {
        (0..self.len()).map(|i| self.get_spectrum_by_index(i).unwrap())
    }

    pub fn get_spectrum_by_key(
        &mut self,
        key: Id,
    ) -> Option<Result<AnnotatedSpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init
            && let Err(e) = self.build_index()
        {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        let rec = self.offsets.get_by_key(key)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        self.line_number = rec.line_number;
        self.line_cache.clear();
        Some(self.read_next())
    }

    pub fn get_spectrum_by_index(
        &mut self,
        index: usize,
    ) -> Option<Result<AnnotatedSpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init
            && let Err(e) = self.build_index()
        {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        let rec = self.offsets.get_by_index(index)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        self.line_number = rec.line_number;
        self.line_cache.clear();
        Some(self.read_next())
    }

    pub fn get_spectrum_by_name(
        &mut self,
        name: &str,
    ) -> Option<Result<AnnotatedSpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init
            && let Err(e) = self.build_index()
        {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        let rec = self.offsets.get_by_name(name)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(MzSpecLibTextParseError::IOError(e)));
        }
        self.line_number = rec.line_number;
        self.line_cache.clear();
        Some(self.read_next())
    }
}

#[allow(clippy::missing_panics_doc, clippy::missing_errors_doc)]
#[cfg(test)]
mod test {
    use std::fs;

    use super::*;

    #[test]
    fn test_indexed_reader() {
        let buf = io::BufReader::new(
            fs::File::open("../data/chinese_hamster_hcd_selected_head.mzspeclib.txt").unwrap(),
        );

        let mut this = MzSpecLibParser::new(buf).unwrap();

        this.build_index().unwrap();

        eprintln!("{:#?}", this.offsets);

        assert_eq!(this.len(), 7);

        let spec = this.get_spectrum_by_index(3).unwrap().unwrap();
        assert_eq!(spec.key, 4);
    }

    #[test]
    fn test_header() -> Result<(), MzSpecLibTextParseError> {
        let buf = io::BufReader::new(
            fs::File::open("../data/chinese_hamster_hcd_selected_head.mzspeclib.txt").unwrap(),
        );

        let mut this = MzSpecLibParser::new(buf)?;
        let header = this.header();
        eprintln!("{header:?}");
        // TODO: switch to assigning these basic attributes to header fields
        assert_eq!(header.attributes().len(), 1);

        let spec = this.read_next()?;
        eprintln!("{spec}");

        Ok(())
    }
}
