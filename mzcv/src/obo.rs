use std::{collections::HashMap, fs::File, io::BufReader, path::Path};

use chrono::FixedOffset;
use context_error::{BoxedError, Context, CreateError, ErrorKind};
use flate2::bufread::GzDecoder;

use crate::{CVVersion, hash_buf_reader::HashBufReader};

/// An Obo ontology. This can be read from a file with [`Self::from_file`] or from a raw reader
/// with [`Self::from_raw`].
#[derive(Clone, Debug, Default)]
pub struct OboOntology {
    /// The data version
    pub data_version: Option<String>,
    /// The last updated date. (format: year, month, day, hour, min)
    pub date: Option<(u16, u8, u8, u8, u8)>,
    /// The hash of the file. This defaults to [`sha2::Sha256`] with [`Self::from_file`] but can be set
    /// differently with [`Self::from_raw`].
    pub hash: Vec<u8>,
    /// The other headers of the Obo file. (tag, value, trailing modifiers, comment)
    pub headers: Vec<(String, String, Vec<(String, String)>, Option<String>)>,
    /// All enclosed objects
    pub objects: Vec<OboStanza>,
}

/// An Obo stanza.
#[derive(Clone, Debug, Default)]
pub struct OboStanza {
    /// The stanza type
    pub stanza_type: OboStanzaType,
    /// The id, split into the CV and local identifiers, the local id is stored as string as some ontologies use non numeric values
    pub id: (Option<String>, String),
    pub definition: Option<(
        String,
        Vec<(Option<String>, String)>,
        Vec<(String, String)>,
        Option<String>,
    )>,
    /// The synonyms for this stanza
    pub synonyms: Vec<OboSynonym>,
    /// The tags that are defined for this stanza (value, trialing modifiers, comment)
    pub lines: HashMap<String, Vec<(String, Vec<(String, String)>, Option<String>)>>,
    /// All property value tags parsed as the defined value type (value, trialing modifiers, comment)
    pub property_values: HashMap<String, Vec<(OboValue, Vec<(String, String)>, Option<String>)>>,
}

/// A synonym in an Obo stanza
#[derive(Clone, Debug, Default)]
pub struct OboSynonym {
    /// The synonym itself
    pub synonym: String,
    /// The type or scope of a synonym
    pub scope: SynonymScope,
    /// Optional synonym type name
    pub type_name: Option<String>,
    /// The dbxref list
    pub cross_references: Vec<(Option<String>, String)>,
    /// The trailing modifiers
    pub trailing_modifiers: Vec<(String, String)>,
    /// The comment
    pub comment: Option<String>,
}

/// The type or scope for a synonym
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum SynonymScope {
    /// An exact relation
    Exact,
    /// A broad relation
    Broad,
    /// A narrow relation
    Narrow,
    /// A related term
    #[default]
    Related,
}

impl std::str::FromStr for SynonymScope {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "EXACT" => Ok(Self::Exact),
            "BROAD" => Ok(Self::Broad),
            "NARROW" => Ok(Self::Narrow),
            "RELATED" => Ok(Self::Related),
            _ => Err(()),
        }
    }
}

/// The type for an Obo stanza
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum OboStanzaType {
    /// A Term stanza
    #[default]
    Term,
    /// A Typedef stanza
    Typedef,
    /// An instance stanza, generally not used
    Instance,
}

/// The value that a property value tag can have in Obo
#[derive(Clone, Debug)]
pub enum OboValue {
    /// A string, when the unit is not set or `xsd:string`
    String(String),
    /// A URI, when the unit is not set or `xsd:anyURI`
    URI(String),
    /// A 64 bit floating point, when the unit is `xsd:double`, `xsd:decimal` or `xsd:float`
    Float(f64),
    /// An integer, when the unit is `xsd:integer`, `xsd:int`, `xsd:nonNegativeInteger`, or `xsd:positiveInteger`
    Integer(isize),
    /// A boolean, when the unit is `xsd:boolean`
    Boolean(bool),
    /// A datetime in RFC 3339 format, when the unit is `xsd:dateTime`
    DateTime(chrono::DateTime<FixedOffset>),
}

impl std::fmt::Display for OboValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::String(s) | Self::URI(s) => write!(f, "{s}"),
            Self::Float(s) => write!(f, "{s}"),
            Self::Integer(s) => write!(f, "{s}"),
            Self::Boolean(s) => write!(f, "{s}"),
            Self::DateTime(s) => write!(f, "{}", s.to_rfc3339()),
        }
    }
}

impl OboValue {
    fn parse(
        unit: &str,
        value: &str,
        base_context: Context<'static>,
    ) -> Result<Self, BoxedError<'static, OboError>> {
        match unit {
            "xsd:string" => Ok(Self::String(value.to_string())),
            "xsd:anyURI" => Ok(Self::URI(value.to_string())),
            "xsd:double" | "xsd:float" | "xsd:decimal" => {
                if !value.starts_with('-') && value.contains('-') {
                    // Some ontologies use a range
                    Ok(Self::String(value.to_string()))
                } else {
                    Ok(Self::Float(value.parse::<f64>().map_err(|e| {
                        BoxedError::new(
                            OboError::InvalidFloat,
                            "Could not parse float",
                            e.to_string(),
                            base_context,
                        )
                    })?))
                }
            }
            "xsd:boolean" => Ok(Self::Boolean(value == "true" || value == "1")),
            "xsd:dateTime" => Ok(Self::DateTime(
                chrono::DateTime::parse_from_rfc3339(value).map_err(|e| {
                    BoxedError::new(
                        OboError::InvalidDate,
                        "Could not parse date time",
                        e.to_string(),
                        base_context,
                    )
                })?,
            )),
            "xsd:integer" | "xsd:int" | "xsd:nonNegativeInteger" | "xsd:positiveInteger" => {
                Ok(Self::Integer(value.parse::<isize>().map_err(|e| {
                    BoxedError::new(
                        OboError::InvalidInteger,
                        "Could not parse integer",
                        e.to_string(),
                        base_context,
                    )
                })?))
            }
            dt => Err(BoxedError::new(
                OboError::InvalidDataType,
                "Unrecognised data type",
                format!("'{dt}' is not a recognised datatype"),
                base_context,
            )),
        }
    }
}

impl OboOntology {
    /// Parse an [`OboOntology`] from a path.
    /// # Errors
    /// If the text contained is not valid according to the Obo format.
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self, BoxedError<'static, OboError>> {
        let base_context = Context::none()
            .source(path.as_ref().to_string_lossy())
            .to_owned();
        let file = File::open(path.as_ref()).map_err(|e| {
            BoxedError::new(
                OboError::CouldNotOpenFile,
                "Could not open file",
                e.to_string(),
                base_context.clone(),
            )
        })?;
        if path
            .as_ref()
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            Self::from_raw_internal(
                HashBufReader::<_, sha2::Sha256>::new(GzDecoder::new(BufReader::new(file))),
                &base_context,
            )
        } else {
            Self::from_raw_internal(HashBufReader::<_, sha2::Sha256>::new(file), &base_context)
        }
    }

    /// Parse an [`OboOntology`] from a raw reader.
    /// # Errors
    /// If the text contained is not valid according to the Obo format.
    pub fn from_raw(
        reader: HashBufReader<impl std::io::Read, impl sha2::Digest>,
    ) -> Result<Self, BoxedError<'static, OboError>> {
        Self::from_raw_internal(reader, &Context::none())
    }

    fn from_raw_internal(
        mut reader: HashBufReader<impl std::io::Read, impl sha2::Digest>,
        base_context: &Context<'static>,
    ) -> Result<Self, BoxedError<'static, OboError>> {
        let mut obo = Self::default();
        let mut recent_obj = None;

        for (line_index, line) in reader.lines().enumerate() {
            // let line = match reader.next_line() {
            //     Ok(Some(line)) => line,
            //     Ok(None) => break,
            //     Err(e) => {
            //         return Err(BoxedError::new(
            //             OboError::CouldNotReadLine,
            //             "Could not read line",
            //             e.to_string(),
            //             base_context.clone().line_index(line_index as u32),
            //         ));
            //     }
            // };
            let line = line
                .map_err(|e| {
                    BoxedError::new(
                        OboError::CouldNotReadLine,
                        "Could not read line",
                        e.to_string(),
                        base_context.clone().line_index(line_index as u32),
                    )
                })?
                .trim_end()
                .to_string();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                if let Some(obj) = recent_obj {
                    obo.objects.push(obj);
                }
                let stanza_type = match line[1..=line.len() - 2]
                    .trim()
                    .to_ascii_lowercase()
                    .as_str()
                {
                    "term" => OboStanzaType::Term,
                    "typedef" => OboStanzaType::Typedef,
                    "instance" => OboStanzaType::Instance,
                    _ => {
                        return Err(BoxedError::new(
                            OboError::InvalidStanzaType,
                            "Invalid Obo stanza",
                            "The stanza types has to be any of Term, Typedef, or Instance",
                            base_context
                                .clone()
                                .lines(0, line.to_string())
                                .line_index(line_index as u32)
                                .add_highlight((0, 1..=line.len() - 2)),
                        ));
                    }
                };
                recent_obj = Some(OboStanza::new(stanza_type));
            } else if let Some((id, value_line)) = line.split_once(':') {
                let (value_line, trailing_modifiers, comment) = tokenise_obo_value_line(
                    value_line,
                    base_context
                        .clone()
                        .line_index(line_index as u32)
                        .lines(0, line.clone()),
                )?;

                let value_line = value_line.trim();

                if let Some(obj) = &mut recent_obj {
                    if id == "property_value" {
                        let first_space = value_line
                            .char_indices()
                            .find_map(|(i, c)| (c == ' ').then_some(i))
                            .unwrap();
                        let last_space = value_line
                            .char_indices()
                            .rfind(|(_, c)| *c == ' ')
                            .map(|(i, _)| i)
                            .unwrap();
                        if first_space == last_space {
                            let name = value_line[..first_space].trim();
                            let value = value_line[first_space..].trim().trim_matches('"');
                            obj.property_values
                                .entry(name.to_string())
                                .or_insert(Vec::new())
                                .push((
                                    OboValue::String(value.to_string()),
                                    trailing_modifiers,
                                    comment,
                                ));
                        } else {
                            let name = value_line[..first_space].trim().trim_end_matches(':');
                            let value =
                                value_line[first_space..last_space].trim().trim_matches('"');
                            let unit = value_line[last_space..].trim();
                            obj.property_values
                                .entry(name.to_string())
                                .or_insert(Vec::new())
                                .push((
                                    OboValue::parse(
                                        unit,
                                        value,
                                        base_context
                                            .clone()
                                            .line_index(line_index as u32)
                                            .lines(0, line.clone()),
                                    )?,
                                    trailing_modifiers,
                                    comment,
                                ));
                        }
                    } else if id == "id" {
                        if let Some((cv, id)) = value_line
                            .split_once(':')
                            .or_else(|| value_line.split_once('_'))
                        {
                            obj.id = (Some(unescape(cv)), unescape(id));
                        } else {
                            obj.id = (None, unescape(id));
                        }
                    } else if id == "synonym" {
                        let parts = tokenise(value_line).map_err(|close| {
                            BoxedError::new(
                                OboError::InvalidLine,
                                "Invalid synonym line",
                                format!(
                                    "The line did not contain the closing delimiter: `{close}`",
                                ),
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            )
                        })?;
                        if parts.is_empty() || parts.len() > 4 {
                            return Err(BoxedError::new(
                                OboError::InvalidSynonym,
                                "Invalid synonym line",
                                format!(
                                    "The number of elements ({}) is not correct, one to four elements expected `\"def\" TYPE? NAME? DBXREF?`",
                                    parts.len()
                                ),
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            ));
                        }
                        if parts.len() == 4 && parts[3].0 != Some(']') {
                            return Err(BoxedError::new(
                                OboError::InvalidSynonym,
                                "Invalid synonym line",
                                "The DBXREF should be enclosed in square brackets `[]`",
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            ));
                        }
                        let synonym = unescape(parts[0].1);
                        let scope = parts.get(1).map(|(_, ty)| ty.parse::<SynonymScope>().map_err(|()| BoxedError::new(
                                OboError::InvalidSynonym,
                                "Invalid synonym line",
                                "The type is not correct, expected one of: EXACT, BROAD, NARROW, RELATED",
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            ))).transpose()?.unwrap_or_default();
                        let type_name = parts.get(2).map(|(_, n)| unescape(n));
                        let cross_references = parts
                            .get(3)
                            .map(|(_, n)| parse_dbxref(n))
                            .unwrap_or_default();
                        obj.synonyms.push(OboSynonym {
                            synonym,
                            scope,
                            type_name,
                            cross_references,
                            trailing_modifiers,
                            comment,
                        });
                    } else if id == "def" {
                        let parts = tokenise(value_line).map_err(|close| {
                            BoxedError::new(
                                OboError::InvalidLine,
                                "Invalid def line",
                                format!(
                                    "The line did not contain the closing delimiter: `{close}`",
                                ),
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            )
                        })?;
                        if parts.len() != 2 {
                            return Err(BoxedError::new(
                                OboError::InvalidLine,
                                "Invalid sef line",
                                format!(
                                    "The number of elements ({}) is not correct, two elements expected `\"def\" [DBXREF]`",
                                    parts.len()
                                ),
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            ));
                        }
                        if parts[1].0 != Some(']') {
                            return Err(BoxedError::new(
                                OboError::InvalidLine,
                                "Invalid def line",
                                "The DBXREF should be enclosed in square brackets `[]`",
                                base_context
                                    .clone()
                                    .lines(0, line.clone())
                                    .line_index(line_index as u32),
                            ));
                        }
                        let def = unescape(parts[0].1);
                        let cross_references = parse_dbxref(parts[1].1);
                        obj.definition = Some((def, cross_references, trailing_modifiers, comment));
                    } else {
                        obj.lines
                            .entry(id.trim().to_string())
                            .or_insert(Vec::new())
                            .push((unescape(value_line), trailing_modifiers, comment));
                    }
                } else if id.eq_ignore_ascii_case("data-version") {
                    obo.data_version = Some(unescape(value_line));
                } else if id.eq_ignore_ascii_case("date") {
                    let base = BoxedError::new(
                        OboError::InvalidDate,
                        "Invalid date time",
                        "",
                        base_context
                            .clone()
                            .lines(0, line.clone())
                            .line_index(line_index as u32),
                    );
                    if let Some((date, time)) = value_line.trim().split_once(' ')
                        && let Some((hour, minute)) = time.split_once(':')
                    {
                        let date_parts = date.splitn(3, ':').collect::<Vec<_>>();
                        if date_parts.len() != 3 {
                            return Err(base.clone().long_description(
                                "The date part does not follow the format of 'dd:mm:yyyy'",
                            ));
                        }

                        let d = date_parts[0].parse().map_err(|e| {
                            base.clone().long_description(format!(
                                "The day part is not a valid number '{e}'"
                            ))
                        })?;
                        let m = date_parts[1].parse().map_err(|e| {
                            base.clone().long_description(format!(
                                "The month part is not a valid number '{e}'"
                            ))
                        })?;
                        let y = date_parts[2].parse().map_err(|e| {
                            base.clone().long_description(format!(
                                "The year part is not a valid number '{e}'"
                            ))
                        })?;
                        let hour = hour.parse().map_err(|e| {
                            base.clone().long_description(format!(
                                "The hour part is not a valid number '{e}'"
                            ))
                        })?;
                        let minute = minute.parse().map_err(|e| {
                            base.clone().long_description(format!(
                                "The minute part is not a valid number '{e}'"
                            ))
                        })?;

                        obo.date = Some((y, m, d, hour, minute));
                    } else {
                        return Err(base.clone().long_description(
                            "The date time does not follow the format of 'dd:mm::yyyy hh:mm'",
                        ));
                    }
                } else {
                    obo.headers.push((
                        id.to_string(),
                        unescape(value_line),
                        trailing_modifiers,
                        comment,
                    ));
                }
            } else {
                return Err(BoxedError::new(
                    OboError::InvalidLine,
                    "Invalid Obo line",
                    "This line could not be recognised as a valid line in the Obo format",
                    base_context
                        .clone()
                        .line_index(line_index as u32)
                        .lines(0, line),
                ));
            }
        }
        if let Some(obj) = recent_obj {
            obo.objects.push(obj);
        }
        obo.hash = reader.hash().to_vec();
        Ok(obo)
    }

    /// Get the version of this Obo file as parsed from the file
    pub fn version(&self) -> CVVersion {
        CVVersion {
            last_updated: self.date,
            version: self.data_version.clone(),
            hash: self.hash.clone(),
        }
    }
}

fn unescape(value: &str) -> String {
    let mut result = String::new();
    for c in value.trim().chars() {
        if c != '\\' {
            result.push(c);
        }
    }
    result
}

impl OboStanza {
    pub(crate) fn new(stanza_type: OboStanzaType) -> Self {
        Self {
            stanza_type,
            ..Self::default()
        }
    }
}

/// All possible errors when reading an Obo file
#[derive(Clone, Copy, Debug, Default, Eq, Ord, PartialEq, PartialOrd)]
pub enum OboError {
    /// If the file could not be opened
    #[default]
    CouldNotOpenFile,
    /// If a line give an error while reading
    CouldNotReadLine,
    /// If a line is wholly invalid
    InvalidLine,
    /// If a synonym line is invalid
    InvalidSynonym,
    /// If a date was formatted incorrectly
    InvalidDate,
    /// If a float was formatted incorrectly
    InvalidFloat,
    /// If an integer was formatted incorrectly
    InvalidInteger,
    /// If an unrecognised datatype was specified
    InvalidDataType,
    /// If an invalid stanza type was specified
    InvalidStanzaType,
}

impl ErrorKind for OboError {
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

// Tokenise a value line, the `<value>` part is left untouched the trailing modifiers and comment are unescaped.
// `<value> {<trailing modifiers>} ! <comment>`
fn tokenise_obo_value_line<'a>(
    text: &'a str,
    context: Context<'static>,
) -> Result<(&'a str, Vec<(String, String)>, Option<String>), BoxedError<'static, OboError>> {
    let mut trailing_start = None;
    let mut comment_start = None;
    let mut enclosed: Option<(char, usize)> = None;
    let mut escaped = false;

    for (index, char) in text.char_indices() {
        if !escaped {
            if let Some((c, i)) = enclosed.as_ref().copied() {
                if c == char {
                    enclosed = None;
                    if char == '}' {
                        trailing_start = Some(i);
                    }
                }
            } else if char == '\"' {
                enclosed = Some(('\"', index));
            } else if char == '{' {
                enclosed = Some(('}', index));
            } else if char == '!' {
                comment_start = Some(index);
                break;
            }
        }

        escaped = char == '\\';
    }

    if let Some((c, _)) = enclosed {
        return Err(BoxedError::new(
            OboError::InvalidLine,
            "Invalid Obo line",
            format!("This line is enclosed with '{c}' but this is not closed"),
            context,
        ));
    }

    let value = &text[..trailing_start.or(comment_start).unwrap_or(text.len())];
    let mut trailing = Vec::new();
    if let Some(start) = trailing_start {
        let slice = &text[start..comment_start.unwrap_or(text.len())];
        let mut name_start = 0;
        let mut value_start = None;

        // Parse trailing modifiers while allowing both escaping and enclosing
        for (index, char) in slice.char_indices() {
            if let Some((c, _)) = enclosed.as_ref().copied() {
                if c == char && !escaped {
                    enclosed = None;
                }
            } else if char == '\"' {
                enclosed = Some(('\"', index));
            } else if char == '=' && !escaped {
                value_start = Some(index);
            } else if char == ',' && !escaped {
                let Some(value) = value_start else {
                    return Err(BoxedError::new(
                        OboError::InvalidLine,
                        "Invalid Obo trailing modifier",
                        "A trailing modifier should be a name value pair separated by '=' but this symbol was not found",
                        context,
                    ));
                };

                trailing.push((
                    unescape(&slice[name_start..value]),
                    unescape(&slice[value + 1..index]),
                ));
                name_start = index + 1;
                value_start = None;
            }

            escaped = char == '\\';
        }
    }
    let comment = comment_start.map(|s| unescape(&text[s..]));

    Ok((value.trim(), trailing, comment))
}

/// Split into parts with the enclosing characters
fn tokenise(text: &str) -> Result<Vec<(Option<char>, &str)>, char> {
    let mut parts = Vec::new();
    let mut enclosed = None;
    let mut start = None;
    let mut escaped = false;
    for (index, char) in text.char_indices() {
        if !escaped {
            if let Some(s) = start {
                if let Some(close) = enclosed {
                    if char == close {
                        parts.push((Some(close), text[s + 1..index].trim()));
                        start = None;
                        enclosed = None;
                    }
                } else if char == ' ' {
                    parts.push((None, text[s..index].trim()));
                    start = None;
                }
            } else {
                match char {
                    '\"' => {
                        enclosed = Some('\"');
                    }
                    '[' => {
                        enclosed = Some(']');
                    }
                    '{' => {
                        enclosed = Some('}');
                    }
                    ' ' => continue, // Skip double spaces
                    _ => (),
                }
                start = Some(index);
            }
        }
        escaped = char == '\\';
    }
    if let Some(s) = start {
        if let Some(close) = enclosed {
            return Err(close);
        }
        parts.push((None, text[s..].trim()));
    }
    Ok(parts)
}

/// Assumes one opening symbol is already read
fn enclosed_by(text: &str, closing: char) -> Option<usize> {
    let mut escaped = false;
    for (index, char) in text.char_indices() {
        if char == closing && !escaped {
            return Some(index);
        }
        escaped = char == '\\';
    }
    None
}

fn parse_dbxref(text: &str) -> Vec<(Option<String>, String)> {
    text.split(',')
        .filter(|s| !s.is_empty())
        .map(|id| {
            id.split_once(':').map_or_else(
                || (None, unescape(id)),
                |(r, i)| (Some(unescape(r)), unescape(i)),
            )
        })
        .collect()
}
