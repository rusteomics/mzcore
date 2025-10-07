use std::{
    collections::{HashMap, VecDeque},
    fmt::Display,
    io::{self, prelude::*},
};

use indexmap::IndexMap;
use mzcore::prelude::CompoundPeptidoformIon;
use mzdata::{
    mzpeaks::prelude::PeakCollectionMut,
    params::{CURIE, ControlledVocabulary, ParamValue},
};

use crate::{
    fragment::{Fragment, parse_mz_paf},
    mzspeclib::{
        Analyte, AnnotatedPeak, Attribute, AttributeParseError, AttributeSet,
        AttributeValueParseError, Attributed, AttributedMut, EntryType, IdType, Interpretation,
        LibraryHeader, LibrarySpectrum, Term,
    },
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
                        if z1 == 0 { return Ok(z1) } else { z += z1 }
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
                    current_attribute_set.as_mut().unwrap().add_attribute(attr);
                }
                Err(e) => {
                    if buf.starts_with("<") {
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
                            if let Some((_, rest)) = buf.trim().split_once(" ") {
                                if !rest.ends_with(">") {
                                    return Err(MzSpecLibTextParseError::InvalidContentAtState(
                                        buf,
                                        self.state,
                                        self.line_number,
                                    ));
                                }
                                if let Some((entry_tp, id)) = rest[..rest.len() - 1].split_once("=")
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
                                    current_attribute_set =
                                        Some(AttributeSet::new(set_id, set_entry_tp, Vec::new()));
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
                        continue;
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
                        CURIE {
                            controlled_vocabulary: ControlledVocabulary::MS,
                            accession: 1003186,
                        } => {
                            self.header.format_version = attr.value.to_string();
                        }
                        _ => {
                            self.header.add_attribute(attr);
                        }
                    }
                }
                Err(e) => {
                    if buf.starts_with("<") {
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
                        continue;
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

    fn parse_decl(
        &mut self,
        buf: &str,
        decl: &str,
        to_state: ParserState,
    ) -> Result<u32, MzSpecLibTextParseError> {
        let id = if buf.starts_with(decl) {
            self.state = to_state;
            let buf_ = buf.trim();
            if let Some((_, val)) = buf_[..buf_.len() - 1].split_once("=") {
                let id = val.parse::<IdType>().map_err(|_| {
                    MzSpecLibTextParseError::InvalidContentAtState(
                        buf.to_string(),
                        self.state,
                        self.line_number,
                    )
                })?;
                id
            } else {
                return Err(MzSpecLibTextParseError::InvalidContentAtState(
                    buf.to_string(),
                    self.state,
                    self.line_number,
                ));
            }
        } else {
            return Err(MzSpecLibTextParseError::InvalidContentAtState(
                buf.to_string(),
                self.state,
                self.line_number,
            ));
        };
        Ok(id)
    }

    fn read_analyte(&mut self) -> Result<Analyte, MzSpecLibTextParseError> {
        let mut buf = String::new();
        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        let id = self.parse_decl(&buf, "<Analyte=", ParserState::Analyte)?;

        let mut analyte = Analyte::new(id, Vec::new());
        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => {
                    analyte.add_attribute(attr);
                }
                Err(e) => {
                    if buf.starts_with("<") {
                        self.push_back_line(buf);
                        break;
                    } else if buf.trim().is_empty() {
                        continue;
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let mut next_group_id = analyte.find_last_group_id().unwrap_or_default();
        let mut last_group_id = 0;

        let attr_set_name = Term::new(
            mzdata::curie!(MS:1003212),
            "library attribute set name".into(),
        );
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
                    for mut attr in attr_set.attributes.iter().cloned() {
                        if let Some(gi) = attr.group_id {
                            if gi != last_group_id {
                                next_group_id += 1;
                                last_group_id = gi;
                            }
                            attr.group_id = Some(next_group_id);
                        }
                        analyte.add_attribute(attr);
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
        let id = self.parse_decl(&buf, "<Interpretation=", ParserState::Interpretation)?;
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
                    } else if buf.starts_with("<") {
                        self.push_back_line(buf);
                        break;
                    } else if buf.trim().is_empty() {
                        continue;
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let mut next_group_id = interp.find_last_group_id().unwrap_or_default();
        let mut last_group_id = 0;

        let attr_set_name = Term::new(
            mzdata::curie!(MS:1003212),
            "library attribute set name".into(),
        );
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
                    for mut attr in attr_set.attributes.iter().cloned() {
                        if let Some(gi) = attr.group_id {
                            if gi != last_group_id {
                                next_group_id += 1;
                                last_group_id = gi;
                            }
                            attr.group_id = Some(next_group_id);
                        }
                        interp.add_attribute(attr);
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
        let mut it = buf.split("\t");

        let mz = match it.next().and_then(|v| v.parse::<f64>().ok()) {
            Some(v) => v,
            None => {
                return Err(MzSpecLibTextParseError::InvalidContentAtState(
                    buf.to_string(),
                    ParserState::Peaks,
                    self.line_number,
                ));
            }
        };

        let intensity = match it.next().and_then(|v| v.parse::<f32>().ok()) {
            Some(v) => v,
            None => {
                return Err(MzSpecLibTextParseError::InvalidContentAtState(
                    buf.to_string(),
                    ParserState::Peaks,
                    self.line_number,
                ));
            }
        };

        let mut peak = AnnotatedPeak::new(mz, intensity, 0, Vec::new(), Vec::new());

        match it.next() {
            Some(v) => {
                let annots = parse_mz_paf(v, None, &CompoundPeptidoformIon::default()).unwrap();
                peak.annotations = annots;
            }
            None => return Ok(peak),
        }

        match it.next() {
            Some(v) => {
                peak.aggregations = v.split(",").map(|i| i.to_string()).collect();
            }
            None => return Ok(peak),
        }

        Ok(peak)
    }

    fn read_spectrum(&mut self) -> Result<LibrarySpectrum, MzSpecLibTextParseError> {
        let mut buf = String::new();

        let z = self
            .read_next_line(&mut buf)
            .map_err(MzSpecLibTextParseError::IOError)?;
        if z == 0 {
            return Err(MzSpecLibTextParseError::EOF);
        }
        let id = self.parse_decl(&buf, "<Spectrum=", ParserState::Spectrum)?;
        let mut spec = LibrarySpectrum::new(
            id,
            0,
            Box::default(),
            Vec::new(),
            Vec::new(),
            Vec::new(),
            Default::default(),
        );

        loop {
            match self.read_attribute(&mut buf) {
                Ok(attr) => match attr.name.accession {
                    CURIE {
                        controlled_vocabulary: ControlledVocabulary::MS,
                        accession: 1003061,
                    } => spec.name = attr.value.to_string().into_boxed_str(),
                    CURIE {
                        controlled_vocabulary: ControlledVocabulary::MS,
                        accession: 1003062,
                    } => {
                        spec.index = attr.value.scalar().to_u64().map_err(|v| {
                            MzSpecLibTextParseError::AttributeParseError(
                                AttributeParseError::ValueParseError(v, buf.clone()),
                                self.state,
                                self.line_number,
                            )
                        })? as usize
                    }
                    _ => spec.add_attribute(attr),
                },
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
                        continue;
                    } else {
                        return Err(e);
                    }
                }
            }
        }

        let mut next_group_id = spec.find_last_group_id().unwrap_or_default();
        let mut last_group_id = 0;

        let attr_set_name = Term::new(
            mzdata::curie!(MS:1003212),
            "library attribute set name".into(),
        );
        let attr_sets: Vec<_> = spec
            .find_all(&attr_set_name)
            .map(|v| v.value.to_string())
            .collect();
        for name in attr_sets {
            for attr_set in self
                .header
                .attribute_classes
                .get(&EntryType::Spectrum)
                .into_iter()
                .flatten()
            {
                if attr_set.id == name || attr_set.id == "all" {
                    for mut attr in attr_set.attributes.iter().cloned() {
                        if let Some(gi) = attr.group_id {
                            if gi != last_group_id {
                                next_group_id += 1;
                                last_group_id = gi;
                            }
                            attr.group_id = Some(next_group_id);
                        }
                        spec.add_attribute(attr);
                    }
                }
            }
        }

        if matches!(self.state, ParserState::Peaks) {
            loop {
                let z = self
                    .read_next_line(&mut buf)
                    .map_err(MzSpecLibTextParseError::IOError)?;
                if z == 0 {
                    break;
                }
                let buf_trimmed = buf.trim();
                if buf_trimmed.starts_with("<") {
                    self.push_back_line(buf);
                    self.state = ParserState::Between;
                    break;
                } else if buf_trimmed.is_empty() {
                    self.state = ParserState::Between;
                    break;
                }

                let peak = self.parse_peak_line(&buf_trimmed)?;
                spec.peaks.push(peak);
            }
        }
        Ok(spec)
    }

    pub fn read_next(&mut self) -> Result<LibrarySpectrum, MzSpecLibTextParseError> {
        self.read_spectrum()
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

impl<R: io::Read> Iterator for MzSpecLibParser<R> {
    type Item = Result<LibrarySpectrum, MzSpecLibTextParseError>;

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
    pub key: IdType,
    pub index: usize,
    pub line_number: u64,
    pub name: Box<str>,
}

#[derive(Debug, Default, Clone)]
pub struct LibraryIndex {
    pub init: bool,
    primary_index: IndexMap<usize, LibraryIndexRecord>,
    key_index: HashMap<IdType, usize>,
    name_index: HashMap<Box<str>, usize>,
}

impl LibraryIndex {
    pub fn insert(&mut self, index: usize, record: LibraryIndexRecord) {
        self.key_index.insert(record.key, index);
        self.name_index.insert(record.name.clone(), index);
        self.primary_index.insert(index, record);
    }

    pub fn get_by_key(&self, key: IdType) -> Option<&LibraryIndexRecord> {
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

impl<R: io::Read + io::Seek> MzSpecLibParser<R> {
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
                let (_, rest) = buf.trim().split_once("=").unwrap();
                if rest.ends_with(">") {
                    if let Some(current_record) = current_record {
                        offset_index.insert(current_record.index, current_record);
                    }
                    current_record = Some(Default::default());
                    if let Ok(key) = rest[..rest.len() - 1].parse::<IdType>() {
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
    ) -> impl Iterator<Item = Result<LibrarySpectrum, MzSpecLibTextParseError>> {
        (0..self.len()).map(|i| self.get_spectrum_by_index(i).unwrap())
    }

    pub fn get_spectrum_by_key(
        &mut self,
        key: IdType,
    ) -> Option<Result<LibrarySpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init {
            if let Err(e) = self.build_index() {
                return Some(Err(MzSpecLibTextParseError::IOError(e)));
            }
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
    ) -> Option<Result<LibrarySpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init {
            if let Err(e) = self.build_index() {
                return Some(Err(MzSpecLibTextParseError::IOError(e)));
            }
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
    ) -> Option<Result<LibrarySpectrum, MzSpecLibTextParseError>> {
        if !self.offsets.init {
            if let Err(e) = self.build_index() {
                return Some(Err(MzSpecLibTextParseError::IOError(e)));
            }
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
