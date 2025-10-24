use std::{
    collections::{HashMap, VecDeque},
    fmt::Display,
    fs::File,
    io::{self, BufReader, prelude::*},
    path::{Path, PathBuf},
};

use context_error::{BoxedError, Context, CreateError, ErrorKind, FullErrorContent};
use indexmap::IndexMap;
use itertools::Itertools;
use mzcore::{
    ontology::CustomDatabase,
    prelude::*,
    system::{MassOverCharge, isize::Charge},
};
use mzdata::{
    curie,
    mzpeaks::prelude::PeakCollectionMut,
    params::{ParamValue, Unit, Value},
    spectrum::{Precursor, ScanEvent, ScanWindow, SelectedIon},
};

use crate::{
    fragment::Fragment,
    helper_functions::explain_number_error,
    mzspeclib::{
        Analyte, AnalyteTarget, Attribute, AttributeParseError, AttributeSet, AttributeValue,
        EntryType, Id, Interpretation, LibraryHeader, ProteinDescription, merge_attributes,
        populate_spectrum_description_from_attributes,
    },
    spectrum::{AnnotatedPeak, AnnotatedSpectrum},
    term,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum ParserState {
    Initial,
    Header,
    AttributeSet,
    Spectrum,
    Analyte,
    Interpretation,
    Peaks,
    Between,
    Eof,
}

impl Display for ParserState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

/// The kind of error that can occur when reading an mzSpecLib file
#[derive(Debug, Default, PartialEq, Eq, Clone, Copy)]
pub enum MzSpecLibErrorKind {
    /// An error concerning the underlying stream
    #[default]
    IO,
    /// An error concerning the parsing of a declaration
    Declaration,
    /// An error concerning the parsing of an attribute
    Attribute,
    /// An error concerning the parsing of a peak
    Peak,
    /// The end of the file was found but more data was expected
    Eof,
    /// The state is inconsistent (for example the charge was defined multiple times with different values)
    InconsistentState,
    /// A unit was missing on an attribute that needs one
    MissingUnit,
    /// A ProForma string was invalid
    ProForma,
    /// An mzPAF string was invalid
    MzPAF,
}

impl ErrorKind for MzSpecLibErrorKind {
    type Settings = ();
    fn descriptor(&self) -> &'static str {
        "error"
    }
    fn ignored(&self, _settings: Self::Settings) -> bool {
        false
    }
    fn is_error(&self, _settings: Self::Settings) -> bool {
        true
    }
}

/// The state for a parser for mzSpecLib.txt files
#[derive(Debug)]
pub struct MzSpecLibTextParser<'custom_database, Reader: Read> {
    inner: Reader,
    header: LibraryHeader,
    header_attribute_sets_with_context: HashMap<
        EntryType,
        HashMap<String, HashMap<Option<u32>, Vec<(Attribute, Context<'static>)>>>,
    >,
    state: ParserState,
    line_cache: VecDeque<String>,
    line_index: u32,
    path: Option<PathBuf>,
    offsets: LibraryIndex,
    /// A flag to indicate is the last spectrum failed, if so it will scan to the start of the next spectrum next time a spectrum is requested.
    last_error: bool,
    /// The last compound peptidoform ion, used when parsing the mzPAF peaks in a spectrum.
    last_compound_peptidoform: Vec<(u32, AnalyteTarget)>,
    custom_database: Option<&'custom_database CustomDatabase>,
}

impl<'a> MzSpecLibTextParser<'a, BufReader<File>> {
    /// Parse a mzSpecLib file from the given file.
    /// # Errors
    /// If the file does not contain valid mzSpecLib data.
    pub fn open_file(
        path: &Path,
        custom_database: Option<&'a CustomDatabase>,
    ) -> Result<Self, BoxedError<'static, MzSpecLibErrorKind>> {
        Self::open(
            BufReader::new(File::open(path).map_err(|e| {
                BoxedError::new(
                    MzSpecLibErrorKind::IO,
                    "Could not open mzSpecLib file",
                    e.to_string(),
                    Context::none().source(path.to_string_lossy()).to_owned(),
                )
            })?),
            Some(path.to_path_buf()),
            custom_database,
        )
    }
}

impl<'a, R: BufRead> MzSpecLibTextParser<'a, R> {
    /// Parse a mzSpecLib file from the given stream, the original filepath can be given for nicer error messages.
    /// # Errors
    /// If the file does not contain valid mzSpecLib data.
    pub fn open(
        reader: R,
        path: Option<PathBuf>,
        custom_database: Option<&'a CustomDatabase>,
    ) -> Result<Self, BoxedError<'static, MzSpecLibErrorKind>> {
        let mut this = Self {
            inner: reader,
            header: LibraryHeader::default(),
            header_attribute_sets_with_context: HashMap::default(),
            state: ParserState::Initial,
            line_cache: VecDeque::new(),
            line_index: 0,
            path,
            offsets: LibraryIndex {
                done: false,
                primary: IndexMap::new(),
                key: HashMap::new(),
                name: HashMap::new(),
                scan_number: HashMap::new(),
            },
            last_error: false,
            last_compound_peptidoform: Vec::new(),
            custom_database,
        };
        this.read_header()?;
        Ok(this)
    }

    /// Get the header of the file
    pub const fn header(&self) -> &LibraryHeader {
        &self.header
    }

    /// Get the current context, meaning the line index and file name (if known), the line itself will have to be added by the caller.
    fn current_context(&self) -> Context<'_> {
        self.path.as_ref().map_or_else(
            || Context::none().line_index(self.line_index),
            |p| {
                Context::none()
                    .line_index(self.line_index)
                    .source(p.to_string_lossy())
            },
        )
    }

    fn push_back_line(&mut self, line: String) {
        self.line_cache.push_front(line);
        self.line_index -= 1;
    }

    /// Read the next line in the reader and update the internal state.
    /// # Errors
    /// If the reader errors.
    fn read_next_line(&mut self, buf: &mut String) -> io::Result<usize> {
        buf.clear();
        if let Some(cached) = self.line_cache.pop_front() {
            *buf = cached;
            self.line_index += 1;
            Ok(buf.len())
        } else {
            let mut z = self.inner.read_line(buf)?;
            if z != 0 {
                self.line_index += 1;
            } else {
                self.state = ParserState::Eof;
                return Ok(0);
            }

            while buf.trim().is_empty() || buf.starts_with('#') {
                buf.clear();
                let z1 = if let Some(cached) = self.line_cache.pop_front() {
                    *buf = cached;
                    self.line_index += 1;
                    buf.len()
                } else {
                    let z2 = self.inner.read_line(buf)?;
                    if z2 != 0 {
                        self.line_index += 1;
                    }
                    z2
                };
                if z1 == 0 {
                    return Ok(z1);
                }
                z += z1;
            }
            Ok(z)
        }
    }

    /// Read a single attribute from the current position.
    /// # Errors
    /// If the current line is not a valid attribute.
    fn read_attribute(
        &mut self,
        buf: &mut String,
    ) -> Result<(Option<u32>, Attribute), BoxedError<'static, MzSpecLibErrorKind>> {
        let z = self.read_next_line(buf).map_err(|e| {
            BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )
        })?;
        if z == 0 {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Eof,
                "Early end of file",
                "Expected to read an attribute but the file already ended.",
                self.current_context().lines(0, &*buf).to_owned(),
            ));
        }
        Attribute::parse(buf.trim_ascii_end()).map_err(|e: AttributeParseError| {
            BoxedError::new(
                MzSpecLibErrorKind::Attribute,
                "Invalid attribute",
                e.to_string(),
                self.current_context().lines(0, &*buf).to_owned(),
            )
        })
    }

    /// Read attributes sets at the current point in the reader and store them in the header.
    /// # Errors
    /// When an attribute set is invalid or if the structure is invalid.
    fn read_attribute_sets(&mut self) -> Result<(), BoxedError<'static, MzSpecLibErrorKind>> {
        type IntermediateAttributeSet = (
            EntryType,
            String,
            HashMap<Option<u32>, Vec<(Attribute, Context<'static>)>>,
        );

        fn store<R: Read>(
            current_attribute_set: Option<IntermediateAttributeSet>,
            parser: &mut MzSpecLibTextParser<'_, R>,
        ) {
            if let Some((namespace, id, attributes)) = current_attribute_set {
                let mut attribute_vec = vec![Vec::new(); 1];
                attribute_vec[0].extend(
                    attributes
                        .get(&None)
                        .iter()
                        .flat_map(|a| a.iter())
                        .map(|(a, _)| a.clone()),
                );
                attribute_vec.extend(
                    attributes
                        .iter()
                        .filter_map(|(g, a)| g.map(|g| (g, a)))
                        .sorted_by_key(|(g, _)| *g)
                        .map(|(_, a)| a.iter().map(|(a, _)| a).cloned().collect()),
                );

                parser
                    .header
                    .attribute_classes
                    .entry(namespace)
                    .or_default()
                    .push(AttributeSet {
                        id: id.clone(),
                        namespace,
                        attributes: attribute_vec,
                    });
                let base = parser
                    .header_attribute_sets_with_context
                    .entry(namespace)
                    .or_default()
                    .entry(id)
                    .or_default();
                for (key, value) in &attributes {
                    base.entry(*key).or_default().extend_from_slice(value);
                }
            }
        }

        let mut buf = String::new();
        let mut current_attribute_set: Option<IntermediateAttributeSet> = None;

        loop {
            match self.read_attribute(&mut buf) {
                Ok((group_id, attr)) => {
                    current_attribute_set
                        .as_mut()
                        .ok_or_else(|| {
                            BoxedError::new(
                                // This should be impossible to reach because any lingering attributes should have been picked up in the header
                                MzSpecLibErrorKind::Declaration,
                                "Invalid attribute set",
                                "An attribute was defined before an attribute set was defined",
                                self.current_context().lines(0, &buf).to_owned(),
                            )
                        })?
                        .2
                        .entry(group_id)
                        .or_default()
                        .push((attr, self.current_context().lines(0, &buf).to_owned()));
                }
                Err(e) => {
                    if buf.starts_with('<') {
                        if buf.starts_with("<Spectrum=") {
                            store(current_attribute_set.take(), self);
                            self.state = ParserState::Spectrum;
                            self.push_back_line(buf);
                            break;
                        } else if buf.starts_with("<AttributeSet") {
                            self.state = ParserState::AttributeSet;
                            if let Some((_, rest)) = buf.trim().split_once(' ') {
                                if !rest.ends_with('>') {
                                    return Err(BoxedError::new(
                                        MzSpecLibErrorKind::Declaration,
                                        "Invalid attribute set",
                                        "The closing bracket '>' is missing",
                                        self.current_context().lines(0, &buf).to_owned(),
                                    ));
                                }
                                if let Some((entry_tp, id)) = rest[..rest.len() - 1].split_once('=')
                                {
                                    store(current_attribute_set.take(), self);
                                    let set_entry_tp = match entry_tp {
                                        "Spectrum" => EntryType::Spectrum,
                                        "Analyte" => EntryType::Analyte,
                                        "Interpretation" => EntryType::Interpretation,
                                        "Cluster" => EntryType::Cluster,
                                        _ => {
                                            return Err(BoxedError::new(
                                                MzSpecLibErrorKind::Declaration,
                                                "Invalid attribute set type",
                                                "Use 'Spectrum', 'Analyte', 'Interpretation', or 'Cluster' and note that this is case sensitive",
                                                self.current_context().lines(0, &buf).to_owned(),
                                            ));
                                        }
                                    };
                                    let set_id = id.to_string();
                                    current_attribute_set =
                                        Some((set_entry_tp, set_id, HashMap::new()));
                                } else {
                                    return Err(BoxedError::new(
                                        MzSpecLibErrorKind::Declaration,
                                        "Invalid attribute set",
                                        "The id for the attribute set is missing, is should look like '<AttributeSet Spectrum=all>' and the equals '=' was missing",
                                        self.current_context().lines(0, &buf).to_owned(),
                                    ));
                                }
                            } else {
                                return Err(BoxedError::new(
                                    MzSpecLibErrorKind::Declaration,
                                    "Invalid attribute set",
                                    "The id for the attribute set is missing, is should look like '<AttributeSet Spectrum=all>' and the space was missing",
                                    self.current_context().lines(0, &buf).to_owned(),
                                ));
                            }
                        } else {
                            return Err(BoxedError::new(
                                MzSpecLibErrorKind::Declaration,
                                "Invalid declaration",
                                "In the header only attribute sets can be defined '<AttributeSet Spectrum=all>' or a header can be followed by a '<Spectrum=XX>'",
                                self.current_context().lines(0, &buf).to_owned(),
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

    /// Read the mzSpecLib header and store it in this structure.
    /// # Errors
    /// If the header is invalid
    fn read_header(&mut self) -> Result<(), BoxedError<'static, MzSpecLibErrorKind>> {
        let mut buf = String::new();
        let z = self.read_next_line(&mut buf).map_err(|e| {
            BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )
        })?;
        if z == 0 {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Eof,
                "Early end of file",
                "Expected to read the header but the file already ended.",
                self.current_context().lines(0, &buf).to_owned(),
            ));
        }

        if !buf.starts_with("<mzSpecLib>") {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Declaration,
                "Invalid file start",
                "The first line in an mzSpecLib.txt file should be '<mzSpecLib>'",
                self.current_context().lines(0, &buf).to_owned(),
            ));
        }
        self.state = ParserState::Header;

        loop {
            match self.read_attribute(&mut buf) {
                Ok((group_id, attr)) => match attr.name.accession {
                    curie!(MS:1003186) => {
                        self.header.format_version = attr.value.to_string();
                    }
                    _ => {
                        let index = group_id.map_or(0, |i| i as usize + 1);
                        if self.header.attributes.len() <= index {
                            self.header.attributes.extend(std::iter::repeat_n(
                                Vec::new(),
                                index - self.header.attributes.len() + 1,
                            ));
                        }
                        self.header.attributes[index].push(attr);
                    }
                },
                Err(e) => {
                    if buf.starts_with('<') {
                        if buf.starts_with("<Spectrum=") {
                            self.state = ParserState::Spectrum;
                        } else if buf.starts_with("<AttributeSet") {
                            self.state = ParserState::AttributeSet;
                        } else {
                            return Err(BoxedError::new(
                                MzSpecLibErrorKind::Declaration,
                                "Invalid declaration",
                                "In the header only attribute sets can be defined '<AttributeSet Spectrum=all>' or a header can be followed by a '<Spectrum=XX>'",
                                self.current_context().lines(0, &buf).to_owned(),
                            ));
                        }

                        self.push_back_line(buf);
                        break;
                    } else if buf.is_empty() && e.get_kind() != MzSpecLibErrorKind::Eof {
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

        let mut first = true;
        self.header.attributes.retain(|v| {
            let retain = !v.is_empty() || first;
            first = false;
            retain
        });

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
    ) -> Result<u32, BoxedError<'static, MzSpecLibErrorKind>> {
        let id = if buf.starts_with(declaration) {
            self.state = to_state;
            if let Some((before, val)) = buf[..buf.trim_ascii_end().len() - 1].split_once('=') {
                val.trim().parse::<Id>().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Declaration,
                        "Invalid declaration",
                        format!("The ID was {}", explain_number_error(&e)),
                        self.current_context()
                            .lines(0, buf)
                            .to_owned()
                            .add_highlight((0, before.len() + 1, val.len())),
                    )
                })?
            } else {
                return Err(BoxedError::new(
                    MzSpecLibErrorKind::Declaration,
                    "Invalid declaration",
                    "The ID needs to contain an equals symbol `=` and this was missing",
                    self.current_context().lines(0, buf).to_owned(),
                ));
            }
        } else {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Declaration,
                "Invalid declaration",
                format!("The declaration `{declaration}` was expected but not found"),
                self.current_context().lines(0, buf).to_owned(),
            ));
        };
        Ok(id)
    }

    /// Parse an analyte from the stream.
    /// # Errors
    /// If the analyte
    fn read_analyte(&mut self) -> Result<Analyte, BoxedError<'static, MzSpecLibErrorKind>> {
        let mut buf = String::new();
        let z = self.read_next_line(&mut buf).map_err(|e| {
            BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )
        })?;
        if z == 0 {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Eof,
                "Early end of file",
                "Expected to read an analyte but the file already ended.",
                self.current_context().to_owned(),
            ));
        }
        let id = self.parse_open_declaration(&buf, "<Analyte=", ParserState::Analyte)?;
        let mut groups: HashMap<u32, Vec<(Attribute, Context)>> = HashMap::new();

        let mut analyte = Analyte::new(id, AnalyteTarget::default());
        let mut protein = ProteinDescription::default();
        loop {
            match self.read_attribute(&mut buf) {
                Ok((group_id, attr)) => {
                    if attr.name.accession == curie!(MS:1003270)
                        && let Value::String(value) = attr.value.scalar().as_ref()
                    {
                        let mut peptidoform_ion =
                            PeptidoformIon::pro_forma(value, self.custom_database).map_err(
                                |e| e.to_owned().convert::<MzSpecLibErrorKind, BoxedError<'static, MzSpecLibErrorKind>>(|_| MzSpecLibErrorKind::ProForma),
                            )?;
                        if peptidoform_ion.get_charge_carriers().is_none() {
                            match analyte.target {
                                AnalyteTarget::Unknown(Some(c)) => {
                                    peptidoform_ion
                                        .set_charge_carriers(Some(MolecularCharge::proton(c)));
                                }
                                AnalyteTarget::MolecularFormula(f)
                                    if peptidoform_ion.get_charge_carriers().is_none()
                                        && f.charge().value != 0 =>
                                {
                                    peptidoform_ion.set_charge_carriers(Some(
                                        MolecularCharge::proton(f.charge()),
                                    ));
                                }
                                _ => (),
                            }
                        }

                        analyte.target = AnalyteTarget::PeptidoformIon(peptidoform_ion);
                    } else if attr.name.accession == curie!(MS:1000866) {
                        let value = attr.value.scalar().to_string();
                        let mut formula = MolecularFormula::from_pro_forma(
                            &value,
                            0..value.len(),
                            false,
                            false,
                            true,
                            false,
                        )
                        .map_err(|e| e.to_owned().convert::<MzSpecLibErrorKind, BoxedError<'static, MzSpecLibErrorKind>>(|_| MzSpecLibErrorKind::ProForma))?;
                        analyte.target = match analyte.target {
                            AnalyteTarget::Unknown(c) => {
                                if let Some(c) = c {
                                    formula.set_charge(c);
                                }
                                AnalyteTarget::MolecularFormula(formula)
                            }
                            AnalyteTarget::MolecularFormula(_) => {
                                AnalyteTarget::MolecularFormula(formula) // TODO: detect double definitions?
                            }
                            AnalyteTarget::PeptidoformIon(pep) => {
                                AnalyteTarget::PeptidoformIon(pep) // TODO: detect incongruent definitions?
                            }
                        };
                    } else if attr.name.accession == curie!(MS:1000041)
                        && let Value::Int(value) = attr.value.scalar().as_ref()
                    {
                        let charge = Charge::new::<mzcore::system::e>(*value as isize);
                        analyte.target.set_charge(charge);
                    } else if [
                        curie!(MS:1000888), // Stripped peptide sequence
                        curie!(MS:1003043), // number of residues
                        curie!(MS:1003053), // theoretical monoisotopic m/z
                        curie!(MS:1003054), // theoretical average m/z
                        curie!(MS:1000224), // molecular mass (average weight MH+)
                        curie!(MS:1003243), // adduct ion mass (monoisotopic MH+)
                        curie!(MS:1001117), // theoretical neutral mass (monoisotopic M)
                    ]
                    .contains(&attr.name.accession)
                    {
                        // Ignore, can be calculated easily on the fly
                    } else if let Some(group_id) = group_id {
                        groups
                            .entry(group_id)
                            .or_default()
                            .push((attr, self.current_context().lines(0, &buf).to_owned()));
                    } else if !protein.populate_from_attribute(
                        &attr,
                        &self.current_context().lines(0, &buf).to_owned(),
                    )? {
                        analyte.params.push(attr.into());
                    }
                }
                Err(e) => {
                    if buf.starts_with('<') {
                        if !protein.is_empty() {
                            analyte.proteins.push(protein);
                        }
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

        for group in groups.values() {
            if let Some((value, _)) = group
                .iter()
                .find(|a| a.0.name.accession == curie!(MS:1003275))
                && let Some((name, _)) = group
                    .iter()
                    .find(|a| a.0.name.accession == curie!(MS:1003276))
                && group.len() == 2
            {
                analyte.params.push(match &name.value {
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
                let mut protein = ProteinDescription::default();
                let attr_sets: Vec<_> = group
                    .iter()
                    .filter(|(a, _)| a.name == term!(MS:1003212|library attribute set name))
                    .map(|(v, _)| v.value.to_string())
                    .collect();
                for (attribute, context) in group.iter().chain(
                    self.header_attribute_sets_with_context
                        .get(&EntryType::Analyte)
                        .into_iter()
                        .flatten()
                        .filter(|set| attr_sets.contains(set.0))
                        .flat_map(|a| a.1.values())
                        .flatten(),
                ) {
                    if !protein.populate_from_attribute(attribute, context)? {
                        protein.attributes.push(attribute.clone());
                    }
                }
                if group.len() == protein.attributes.len() {
                    // If zero attributes were understood as protein descriptions just store
                    // them as analyte parameters
                    analyte
                        .params
                        .extend(protein.attributes.into_iter().map(Into::into));
                } else {
                    analyte.proteins.push(protein);
                }
            }
        }

        Ok(analyte)
    }

    /// Read an interpretation from the current point in the reader. Assumes the current lines starts with `<Interpretation=`.
    /// # Errors
    /// If there is no valid interpretation at the current point or if the interpretation tag is missing.
    fn read_interpretation(
        &mut self,
    ) -> Result<Interpretation, BoxedError<'static, MzSpecLibErrorKind>> {
        let mut buf = String::new();
        let z = self.read_next_line(&mut buf).map_err(|e| {
            BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )
        })?;
        if z == 0 {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Eof,
                "Early end of file",
                "Expected to read an interpretation but the file already ended.",
                self.current_context().to_owned(),
            ));
        }
        let id =
            self.parse_open_declaration(&buf, "<Interpretation=", ParserState::Interpretation)?;
        let mut interp = Interpretation {
            id,
            attributes: vec![Vec::new(); 1],
            ..Default::default()
        };
        loop {
            match self.read_attribute(&mut buf) {
                Ok((group_id, attribute)) => {
                    if attribute.name == term!(MS:1002357|PSM-level probability) {
                        interp.probability =
                            Some(attribute.value.scalar().to_f64().map_err(|e| {
                                BoxedError::new(
                                    MzSpecLibErrorKind::Attribute,
                                    "Invalid PSM-level probability",
                                    e.to_string(),
                                    self.current_context().lines(0, &buf).to_owned(),
                                )
                            })?);
                    } else if [
                        curie!(MS:1003288), // number of unassigned peaks
                        curie!(MS:1003079), // total unassigned intensity fraction
                        curie!(MS:1003080), // top 20 peak unassigned intensity fraction
                        curie!(MS:1003290), // number of unassigned peaks among top 20 peaks
                        curie!(MS:1003289), // intensity of highest unassigned peak
                        curie!(MS:1001975), // delta m/z
                        curie!(MS:1003209), // monoisotopic m/z deviation
                        curie!(UO:0000000), // unit (for m/z deviation)
                    ]
                    .contains(&attribute.name.accession)
                    {
                        // Ignore, can be calculated easily on the fly
                    } else {
                        let index = group_id.map_or(0, |i| i as usize + 1);
                        if interp.attributes.len() <= index {
                            interp.attributes.extend(std::iter::repeat_n(
                                Vec::new(),
                                index - interp.attributes.len() + 1,
                            ));
                        }
                        interp.attributes[index].push(attribute);
                    }
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

        let attr_sets: Vec<_> = interp.attributes[0]
            .iter()
            .filter(|a| a.name == term!(MS:1003212|library attribute set name))
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
                    merge_attributes(&mut interp.attributes, &attr_set.attributes); // TODO: try to interpret the merged attributes as well.
                }
            }
        }
        let mut first = true;
        interp.attributes.retain(|v| {
            let retain = !v.is_empty() || first;
            first = false;
            retain
        });
        Ok(interp)
    }

    /// Parse an mzSpecLib peak line.
    /// # Errors
    /// If the line is not a valid mzSpecLib peak line.
    fn parse_peak_line(
        &self,
        buf: &str,
    ) -> Result<AnnotatedPeak<Fragment>, BoxedError<'static, MzSpecLibErrorKind>> {
        let mut field_offset = buf.chars().take_while(char::is_ascii_whitespace).count(); // Because only ASCII this is a bytes offset as well
        let mut it = buf[field_offset..].split('\t');

        let mz = it
            .next()
            .ok_or_else(|| {
                BoxedError::new(
                    MzSpecLibErrorKind::Peak,
                    "Peak m/z is missing",
                    "At least two columns are neccessary in the peaks data lines",
                    self.current_context().lines(0, buf).to_owned(),
                )
            })
            .and_then(|v| {
                let r = v.parse::<f64>().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Peak,
                        "Peak m/z is not a number",
                        e.to_string(),
                        self.current_context()
                            .lines(0, buf)
                            .to_owned()
                            .add_highlight((0, field_offset, v.len())),
                    )
                });
                field_offset += v.len() + 1;
                r
            })?;

        let intensity = it
            .next()
            .ok_or_else(|| {
                BoxedError::new(
                    MzSpecLibErrorKind::Peak,
                    "Peak intensity is missing",
                    "At least two columns are neccessary in the peaks data lines",
                    self.current_context().lines(0, buf).to_owned(),
                )
            })
            .and_then(|v| {
                let r = v.parse::<f32>().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Peak,
                        "Peak intensity is not a number",
                        e.to_string(),
                        self.current_context()
                            .lines(0, buf)
                            .to_owned()
                            .add_highlight((0, field_offset, v.len())),
                    )
                });
                field_offset += v.len() + 1;
                r
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
                    let annots = Fragment::mz_paf_substring(
                        &self.current_context().lines(0, buf),
                        buf,
                        field_offset..field_offset + v.len(),
                        self.custom_database,
                        &self.last_compound_peptidoform,
                    )
                    .map_err(|e| {
                        e.to_owned()
                            .convert::<MzSpecLibErrorKind, BoxedError<'static, MzSpecLibErrorKind>>(
                                |_| MzSpecLibErrorKind::MzPAF,
                            )
                    })?;
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

        // TODO: what to do with any further columns?

        Ok(peak)
    }

    /// Read the next full spectrum. This assumes it is already at the "<Spectrum=XX>" line.
    /// # Errors
    /// IF the next spectrum contains data that is not a valid spectrum.
    fn read_spectrum(
        &mut self,
    ) -> Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>> {
        let mut buf = String::new();

        let z = self.read_next_line(&mut buf).map_err(|e| {
            BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )
        })?;
        if z == 0 {
            return Err(BoxedError::new(
                MzSpecLibErrorKind::Eof,
                "Early end of file",
                "Expected to read a spectrum but the file already ended.",
                self.current_context().to_owned(),
            ));
        }
        let key = self.parse_open_declaration(&buf, "<Spectrum=", ParserState::Spectrum)?;
        let mut spec = AnnotatedSpectrum {
            key,
            attributes: vec![Vec::new(); 1],
            ..Default::default()
        };

        // Set default assumptions on mzSpecLib files
        spec.description.ms_level = 2;
        spec.description.signal_continuity = mzdata::spectrum::SignalContinuity::Centroid;
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

        let mut term_collection: HashMap<Option<u32>, Vec<(Attribute, Context<'static>)>> =
            HashMap::new();

        loop {
            match self.read_attribute(&mut buf) {
                Ok((group_id, attr)) => term_collection.entry(group_id).or_default().push((
                    attr,
                    self.current_context().to_owned().lines(0, buf.clone()),
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

        let set_names: Vec<_> = spec.attributes[0]
            .iter()
            .filter(|a| a.name == term!(MS:1003212|library attribute set name))
            .map(|v| v.value.to_string())
            .collect();

        populate_spectrum_description_from_attributes(
            term_collection.iter(),
            &mut spec.description,
            &mut spec.attributes,
        )?;
        if let Some(sets) = self
            .header_attribute_sets_with_context
            .get(&EntryType::Spectrum)
        {
            populate_spectrum_description_from_attributes(
                sets.iter()
                    .filter_map(|attr_set| {
                        (set_names.contains(attr_set.0) || attr_set.0 == "all")
                            .then_some(attr_set.1.iter())
                    })
                    .flatten(),
                &mut spec.description,
                &mut spec.attributes,
            )?;
        }

        spec.description.precursor[0]
            .activation
            ._extract_methods_from_params();

        self.last_compound_peptidoform = spec
            .analytes
            .iter()
            .map(|a| (a.id, a.target.clone()))
            .collect();

        if matches!(self.state, ParserState::Peaks) {
            loop {
                let z = self.read_next_line(&mut buf).map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::IO,
                        "IO error",
                        e.to_string(),
                        self.current_context().to_owned(),
                    )
                })?;
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
                // TODO: validate that the peaks only reference existing analytes, also make sure that
                // spectra with both formula and peptidoforms actually get assigned the correct peaks
                spec.peaks.push(peak);
            }
        }
        Ok(spec)
    }

    /// Read the next spectrum. If the previous spectrum resulted in an error this will scan
    /// forward until the next "<Spectrum=XX>" is found.
    ///
    /// # Errors
    /// If the next spectrum contains invalid data.
    pub fn read_next(
        &mut self,
    ) -> Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>> {
        self.last_compound_peptidoform = Vec::new();
        if self.last_error {
            // If the last element went awry scan the buffer until the start of the next spectrum
            // is found, otherwise it generates an error for each skipped line.

            let mut buf = String::new();

            loop {
                let z = self.read_next_line(&mut buf).map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::IO,
                        "IO error",
                        e.to_string(),
                        self.current_context().to_owned(),
                    )
                })?;
                if z == 0 {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Eof,
                        "Early end of file",
                        "Expected to read the next spectrum but the file already ended.",
                        self.current_context().to_owned(),
                    ));
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
}

impl<R: BufRead> Iterator for MzSpecLibTextParser<'_, R> {
    type Item = Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>>;

    fn next(&mut self) -> Option<Self::Item> {
        if matches!(self.state, ParserState::Eof) {
            return None;
        }
        Some(self.read_next())
    }
}

/// The start of a library item for indexed access
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd)]
pub struct LibraryIndexRecord {
    /// The byte offset into the file
    pub offset: u64,
    /// The line index
    pub line_index: u32,
}

/// An index to find the byte offsets of all records in a file
#[derive(Debug, Clone)]
pub struct LibraryIndex {
    /// Indicate if the full file is scanned yet
    done: bool,
    primary: IndexMap<usize, LibraryIndexRecord>,
    key: HashMap<Id, usize>,
    name: HashMap<Box<str>, usize>,
    scan_number: HashMap<usize, usize>,
}

impl LibraryIndex {
    /// Insert a record into the spectrum index
    fn insert(
        &mut self,
        index: usize,
        key: Id,
        name: Option<Box<str>>,
        scan_number: Option<usize>,
        record: LibraryIndexRecord,
    ) {
        self.key.insert(key, index);
        if let Some(name) = name {
            self.name.insert(name, index);
        }
        if let Some(scan_number) = scan_number {
            self.scan_number.insert(scan_number, index);
        }
        self.primary.insert(index, record);
    }

    /// Get a record by index (the 0 based index number of the order of the spectra in the file)
    pub fn get_by_index(&self, index: usize) -> Option<&LibraryIndexRecord> {
        self.primary.get(&index)
    }

    /// Get a record by ID (the record ID as defined "<Spectrum=XX>")
    pub fn get_by_key(&self, key: Id) -> Option<&LibraryIndexRecord> {
        self.primary.get(self.key.get(&key)?)
    }

    /// Get a record by name/native ID (MS:1003061)
    pub fn get_by_name(&self, name: &str) -> Option<&LibraryIndexRecord> {
        self.primary.get(self.name.get(name)?)
    }

    /// Get a record by scan number (MS:1003057)
    pub fn get_by_scan_number(&self, scan_number: usize) -> Option<&LibraryIndexRecord> {
        self.primary.get(self.scan_number.get(&scan_number)?)
    }

    /// Get the number of spectra in the index
    pub fn len(&self) -> usize {
        self.primary.len()
    }

    /// Check if the spectrum index is empty
    pub fn is_empty(&self) -> bool {
        self.primary.is_empty()
    }
}

impl<R: BufRead + Seek> MzSpecLibTextParser<'_, R> {
    /// Build the spectrum index.
    ///
    /// # Errors
    /// When the underlying stream errors.
    pub fn build_index(&mut self) -> io::Result<()> {
        type IntermediateRecord = (
            usize,              // index
            Id,                 // key
            Option<Box<str>>,   // name (if found)
            Option<usize>,      // scan number (if found)
            LibraryIndexRecord, // record
        );

        self.inner.rewind()?;
        self.line_cache.clear();
        let mut buf = String::new();
        let mut line_count = 0;
        let mut offset_index = LibraryIndex {
            done: false,
            primary: IndexMap::new(),
            key: HashMap::new(),
            name: HashMap::new(),
            scan_number: HashMap::new(),
        };
        let mut offset = 0;
        let mut current_record: Option<IntermediateRecord> = None;
        let mut index = 0;
        while let Ok(z) = self.read_next_line(&mut buf) {
            if z == 0 {
                break;
            }

            if let Some(rest) = buf.strip_prefix("<Spectrum=")
                && let Some(key) = rest.trim().strip_suffix('>')
            {
                if let Some((index, key, name, scan_number, record)) = current_record {
                    offset_index.insert(index, key, name, scan_number, record);
                }
                if let Ok(key) = key.parse::<Id>() {
                    current_record = Some((
                        index,
                        key,
                        None,
                        None,
                        LibraryIndexRecord {
                            offset,
                            line_index: line_count,
                        },
                    ));
                    index += 1;
                } else {
                    // This is an 'all' (or something like that) spectrum, so reset all info but do not increase the index number
                    current_record = None;
                }
            } else if let Some(rest) = buf.strip_prefix("MS:1003061|library spectrum name=")
                && let Some((_, _, name, _, _)) = &mut current_record
            {
                *name = Some(rest.trim().to_string().into_boxed_str());
            } else if let Some(rest) = buf.strip_prefix("MS:1003057|scan number=")
                && let Ok(scan_number) = rest.trim().parse::<usize>()
                && let Some((_, _, _, num, _)) = &mut current_record
            {
                *num = Some(scan_number);
            }
            line_count += 1;
            offset += z as u64;
        }
        if let Some((index, key, name, scan_number, record)) = current_record {
            offset_index.insert(index, key, name, scan_number, record);
        }
        offset_index.done = true;
        self.offsets = offset_index;
        self.inner.rewind()?;
        Ok(())
    }

    /// Get the total number of spectra in the spectrum index.
    pub fn len(&self) -> usize {
        self.offsets.len()
    }

    /// See if the spectrum index is empty.
    pub fn is_empty(&self) -> bool {
        self.offsets.is_empty()
    }

    /// Get a spectrum by the ID, if the index is not yet built this first builds the index.
    /// It returns None when the key does not exist and an error if the spectrum index could not
    /// be built, the correct reader could not be repositioned to the correct location, or the
    /// spectrum could not be correctly parsed.
    pub fn get_spectrum_by_key(
        &mut self,
        key: Id,
    ) -> Option<Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>>> {
        if !self.offsets.done
            && let Err(e) = self.build_index()
        {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        let rec = self.offsets.get_by_key(key)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        self.line_index = rec.line_index;
        self.line_cache.clear();
        Some(self.read_next())
    }

    /// Get a spectrum by the scan number (MS:1003057), if the index is not yet built this first
    /// builds the index. It returns None when the key does not exist and an error if the spectrum
    /// index could not be built, the correct reader could not be repositioned to the correct
    /// location, or the spectrum could not be correctly parsed.
    pub fn get_spectrum_by_scan_number(
        &mut self,
        scan_number: usize,
    ) -> Option<Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>>> {
        if !self.offsets.done
            && let Err(e) = self.build_index()
        {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        let rec = self.offsets.get_by_scan_number(scan_number)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        self.line_index = rec.line_index;
        self.line_cache.clear();
        Some(self.read_next())
    }

    /// Get a spectrum by the index, if the index is not yet built this first builds the index.
    /// It returns None when the index does not exist and an error if the spectrum index could not
    /// be built, the correct reader could not be repositioned to the correct location, or the
    /// spectrum could not be correctly parsed.
    pub fn get_spectrum_by_index(
        &mut self,
        index: usize,
    ) -> Option<Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>>> {
        if !self.offsets.done
            && let Err(e) = self.build_index()
        {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        let rec = self.offsets.get_by_index(index)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        self.line_index = rec.line_index;
        self.line_cache.clear();
        Some(self.read_next())
    }

    /// Get a spectrum by the name (MS:1003061), if the index is not yet built this first builds
    /// the index. It returns None when the name does not exist and an error if the spectrum index
    /// could not be built, the correct reader could not be repositioned to the correct location,
    /// or the spectrum could not be correctly parsed.
    pub fn get_spectrum_by_name(
        &mut self,
        name: &str,
    ) -> Option<Result<AnnotatedSpectrum, BoxedError<'static, MzSpecLibErrorKind>>> {
        if !self.offsets.done
            && let Err(e) = self.build_index()
        {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        let rec = self.offsets.get_by_name(name)?;
        if let Err(e) = self.inner.seek(io::SeekFrom::Start(rec.offset)) {
            return Some(Err(BoxedError::new(
                MzSpecLibErrorKind::IO,
                "IO error",
                e.to_string(),
                self.current_context().to_owned(),
            )));
        }
        self.line_index = rec.line_index;
        self.line_cache.clear();
        Some(self.read_next())
    }
}

#[allow(clippy::missing_panics_doc, clippy::missing_errors_doc)]
#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_indexed_reader() {
        let buf = BufReader::new(
            File::open("../data/chinese_hamster_hcd_selected_head.mzspeclib.txt").unwrap(),
        );

        let mut this = MzSpecLibTextParser::open(buf, None, None).unwrap();

        this.build_index().unwrap();

        eprintln!("{:#?}", this.offsets);

        assert_eq!(this.len(), 7);

        let spec = this.get_spectrum_by_index(3).unwrap().unwrap();
        assert_eq!(spec.key, 4);

        let spec = this.get_spectrum_by_scan_number(5538).unwrap().unwrap();
        assert_eq!(spec.key, 1);

        let spec = this
            .get_spectrum_by_name("AAAALGSHGSCSSEVEK/2_1(10,C,CAM)_50eV")
            .unwrap()
            .unwrap();
        assert_eq!(spec.key, 6);

        let spec = this.get_spectrum_by_index(6).unwrap().unwrap();
        assert_eq!(spec.key, 7);
    }

    #[test]
    fn test_header() -> Result<(), BoxedError<'static, MzSpecLibErrorKind>> {
        let buf = BufReader::new(
            File::open("../data/chinese_hamster_hcd_selected_head.mzspeclib.txt").unwrap(),
        );

        let mut this = MzSpecLibTextParser::open(buf, None, None)?;
        let header = this.header();
        eprintln!("{header:?}");
        // TODO: switch to assigning these basic attributes to header fields
        assert_eq!(header.attributes.len(), 1);

        let _spec = this.read_next()?;

        Ok(())
    }

    #[test]
    fn unknown_cv_value() {
        let text = r"<mzSpecLib>
MS:1003186|library format version=UW:0000000|text";

        let res = MzSpecLibTextParser::open(text.as_bytes(), None, None);
        assert!(res.is_err());
    }

    #[test]
    fn unknown_cv_term_hang() {
        let text = r"<mzSpecLib>
MS:0000000|unknown=a";

        let res = MzSpecLibTextParser::open(text.as_bytes(), None, None);
        assert!(res.is_err());
    }

    #[test]
    fn stack_overflow() {
        let text = b"<mzSpecLib\n";
        let bytes: Vec<u8> = text
            .iter()
            .copied()
            .chain(std::iter::repeat_n(b'\n', 10000))
            .collect();

        let res = MzSpecLibTextParser::open(bytes.as_slice(), None, None);
        assert!(res.is_err());
    }
}
