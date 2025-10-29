//! Methods for reading and parsing CSV files. (Internal use mostly).

use std::{
    borrow::Cow,
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader, Write},
    ops::Range,
    str::FromStr,
    sync::Arc,
};

use context_error::*;
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::helper_functions::check_extension;

/// A single line in a CSV file
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct CsvLine {
    pub(super) line_index: usize,
    pub(super) line: String,
    pub(super) fields: Vec<(Arc<String>, Range<usize>)>,
    pub(super) file: Option<Arc<Box<str>>>,
}

#[allow(dead_code)]
impl CsvLine {
    /// Get the line index (0 based)
    pub const fn line_index(&self) -> usize {
        self.line_index
    }

    /// Get the full line
    pub fn line(&self) -> &str {
        &self.line
    }

    /// Get the column headers
    pub fn headers(&self) -> impl Iterator<Item = &str> {
        self.fields.iter().map(|f| f.0.as_str())
    }

    /// Get the column values
    pub fn values(&self) -> impl Iterator<Item = (Arc<String>, &str)> {
        self.fields
            .iter()
            .map(|f| (f.0.clone(), &self.line[f.1.clone()]))
    }

    /// Get the number of columns
    pub const fn number_of_columns(&self) -> usize {
        self.fields.len()
    }

    /// Get the context applicable to the specified column
    pub fn column_context(&self, column: usize) -> Context<'_> {
        let base = Context::none()
            .line_index(self.line_index as u32)
            .lines(0, &self.line)
            .add_highlight((
                0,
                self.fields[column].1.clone(),
                self.fields[column].0.as_str(),
            ));
        if let Some(source) = &self.file {
            base.source(source.as_ref().as_ref())
        } else {
            base
        }
    }

    /// Get the context for the specified range in the original line
    pub fn range_context<'a>(
        &'a self,
        range: Range<usize>,
        comment: Option<Cow<'a, str>>,
    ) -> Context<'a> {
        let base = Context::none()
            .line_index(self.line_index as u32)
            .lines(0, &self.line);
        let base = if let Some(comment) = comment {
            base.add_highlight((0, range, comment))
        } else {
            base.add_highlight((0, range))
        };
        if let Some(source) = &self.file {
            base.source(source.as_ref().as_ref())
        } else {
            base
        }
    }

    /// Get the context for the whole line
    pub fn full_context(&self) -> Context<'_> {
        let base = Context::none()
            .line_index(self.line_index as u32)
            .lines(0, &self.line);
        if let Some(source) = &self.file {
            base.source(source.as_ref().as_ref())
        } else {
            base
        }
    }

    /// Get the range of a specified column
    pub fn range(&self, index: usize) -> &Range<usize> {
        &self.fields[index].1
    }

    /// Get the specified column, by column name
    /// # Errors
    /// If the given name is not a column header return an error
    pub fn index_column<'a>(
        &'a self,
        name: &str,
    ) -> Result<(&'a str, &'a Range<usize>), BoxedError<'a, BasicKind>> {
        self.fields
            .iter()
            .find(|f| f.0.eq_ignore_ascii_case(name))
            .map(|f| (&self.line[f.1.clone()], &f.1))
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not find given column",
                    format!("This CSV file does not contain the needed column '{name}'"),
                    self.full_context(),
                )
            })
    }

    /// Parse a column into the given format
    /// # Errors
    /// If erroneous extend the base error with the correct context and return that
    pub fn parse_column<'a, F: FromStr>(
        &'a self,
        column: usize,
        base_error: BoxedError<'a, BasicKind>,
    ) -> Result<F, BoxedError<'a, BasicKind>> {
        self[column]
            .parse()
            .map_err(|_| base_error.replace_context(self.column_context(column)))
    }

    /// Parse a column into the given format
    /// # Errors
    /// If erroneous extend the base error with the correct context and return that
    pub fn parse_column_or_empty<'a, F: FromStr>(
        &'a self,
        column: usize,
        base_error: BoxedError<'a, BasicKind>,
    ) -> Result<Option<F>, BoxedError<'a, BasicKind>> {
        let text = &self[column];
        if text.is_empty() || text == "-" {
            Ok(None)
        } else {
            Ok(Some(text.parse().map_err(|_| {
                base_error.replace_context(self.column_context(column))
            })?))
        }
    }
}

impl<Hasher: ::std::hash::BuildHasher + Default> From<&CsvLine>
    for HashMap<String, String, Hasher>
{
    fn from(value: &CsvLine) -> Self {
        value
            .fields
            .iter()
            .map(|(name, range)| (name.to_string(), value.line[range.clone()].to_string()))
            .collect()
    }
}

impl From<&CsvLine> for BTreeMap<String, String> {
    fn from(value: &CsvLine) -> Self {
        value
            .fields
            .iter()
            .map(|(name, range)| (name.to_string(), value.line[range.clone()].to_string()))
            .collect()
    }
}

impl std::ops::Index<usize> for CsvLine {
    type Output = str;
    fn index(&self, index: usize) -> &str {
        &self.line[self.fields[index].1.clone()]
    }
}

/// Parse a CSV file into an iterator with the parsed lines.
/// # Errors
/// If the file cannot be opened it returns `Err` with the error.
/// If any single line cannot be read it returns an error for that line.
pub fn parse_csv(
    path: impl AsRef<std::path::Path>,
    separator: u8,
    provided_header: Option<Vec<String>>,
) -> Result<
    Box<dyn Iterator<Item = Result<CsvLine, BoxedError<'static, BasicKind>>>>,
    BoxedError<'static, BasicKind>,
> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|e| {
        BoxedError::new(
            BasicKind::Error,
            "Could not open file",
            e.to_string(),
            Context::default().source(path.to_string_lossy()).to_owned(),
        )
    })?;
    if check_extension(path, "gz") {
        Ok(Box::new(parse_csv_raw(
            GzDecoder::new(BufReader::new(file)),
            separator,
            provided_header,
            Some(path.to_string_lossy().to_string().into_boxed_str()),
        )?))
    } else {
        Ok(Box::new(parse_csv_raw(
            file,
            separator,
            provided_header,
            Some(path.to_string_lossy().to_string().into_boxed_str()),
        )?))
    }
}

/// Parse a CSV file from a raw `BufReader`
/// # Errors
/// If no header is provided and the first line could not be read as a header line.
/// Or if the 'sep=C' uses a character that is more than 1 byte wide in utf8.
pub fn parse_csv_raw<T: std::io::Read>(
    reader: T,
    mut separator: u8,
    provided_header: Option<Vec<String>>,
    path: Option<Box<str>>,
) -> Result<CsvLineIter<T>, BoxedError<'static, BasicKind>> {
    let reader = BufReader::new(reader);
    let mut lines = reader.lines().enumerate().peekable();
    let mut skip = false;
    if let Some(sep) = lines
        .peek()
        .and_then(|(_, l)| l.as_ref().ok())
        .map(|l| l.trim_start_matches("\u{feff}"))
        .and_then(|l| l.strip_prefix("sep="))
    {
        skip = true;
        if let Some(c) = sep.chars().next() {
            if c.len_utf8() == 1 {
                separator = c as u8;
            } else {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Unicode value separators not supported",
                    "This is a character that takes more than 1 byte to represent in Unicode, this is not supported in parsing CSV files.",
                    Context::line_with_comment(Some(0), format!("sep={sep}"), 4, sep.len(), None),
                ));
            }
        }
    }
    if skip {
        // Actually consume this line
        let _unused = lines.next();
    }
    let column_headers = if let Some(header) = provided_header {
        let (_, column_headers) = lines.peek().ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Could parse csv file",
                "The file is empty",
                Context::default(),
            )
        })?;
        let header_line = column_headers
            .as_ref()
            .map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Could not read header line",
                    err.to_string(),
                    Context::default(),
                )
            })?
            .trim_start_matches("\u{feff}");
        let first_line = csv_separate(header_line, separator)
            .map_err(BoxedError::to_owned)?
            .into_iter()
            .map(|r| Arc::new(header_line[r].to_lowercase()))
            .collect_vec();
        let provided_header = header.into_iter().map(Arc::new).collect();
        if first_line == provided_header {
            drop(lines.next()); // Ignore the first line if the first line is identical to the provided header
        }
        provided_header
    } else {
        let (_, column_headers) = lines.next().ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Could parse csv file",
                "The file is empty",
                Context::none(),
            )
        })?;
        let header_line = column_headers.map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not read header line",
                err.to_string(),
                Context::none(),
            )
        })?;
        let header_line = header_line.trim_start_matches("\u{feff}");

        csv_separate(header_line, separator)
            .map_err(BoxedError::to_owned)?
            .into_iter()
            .map(|r| Arc::new(header_line[r].to_lowercase()))
            .collect()
    };

    Ok(CsvLineIter {
        lines,
        header: column_headers,
        separator,
        file: path.map(Arc::new),
    })
}

/// An iterator returning CSV lines
#[derive(Debug)]
pub struct CsvLineIter<T: std::io::Read> {
    lines: std::iter::Peekable<std::iter::Enumerate<std::io::Lines<BufReader<T>>>>,
    header: Vec<Arc<String>>,
    separator: u8,
    file: Option<Arc<Box<str>>>,
}

impl<T: std::io::Read> Iterator for CsvLineIter<T> {
    type Item = Result<CsvLine, BoxedError<'static, BasicKind>>;
    fn next(&mut self) -> Option<Self::Item> {
        self.lines.next().map(|(line_index, line)| {
            let line = line.map_err(|err|BoxedError::new(BasicKind::Error,
                    "Could not read line",
                    err.to_string(),
                    Context::default().line_index(line_index as u32),
                ))?;
            csv_separate(&line, self.separator).map_err(BoxedError::to_owned).and_then(|row| {
                if self.header.len() == row.len() {
                    Ok(CsvLine {
                        line_index,
                        line,
                        fields: self.header.iter().cloned().zip(row).collect(),
                        file: self.file.clone(),
                    })
                } else {
                    Err(BoxedError::new(BasicKind::Error,
                        "Incorrect number of columns",
                        format!("It does not have the correct number of columns. {} columns were expected but {} were found.", self.header.len(), row.len()),
                        Context::full_line(line_index as u32, line),
                    ))
                }
            })
        })
    }
}

/// # Errors
/// If the line is empty.
pub(crate) fn csv_separate(
    line: &str,
    separator: u8,
) -> Result<Vec<Range<usize>>, BoxedError<'_, BasicKind>> {
    if line.is_empty() {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Empty line",
            "The line is empty",
            Context::none(),
        ));
    }
    let mut enclosed = None;
    let mut was_enclosed = false;
    let mut was_double = false;
    let mut row = Vec::new();
    let mut start = None;
    let mut last_non_whitespace = None;
    for (index, ch) in line.bytes().enumerate() {
        match (ch, enclosed, start) {
            (b'\"' | b'\'', None, None) => {
                enclosed = Some(ch);
                start = Some(index + 1);
                was_double = false;
            }
            (c, Some(e), Some(s)) if c == e => {
                // Ignore an 'escaped' enclosing token by doubling it up
                if line.as_bytes().get(index + 1).copied() == Some(e) && !was_double {
                    // skip the next one
                    was_double = true;
                } else if was_double {
                    was_double = false;
                } else {
                    enclosed = None;
                    if c.is_ascii_whitespace() {
                        row.push(s..last_non_whitespace.unwrap_or(index));
                    } else {
                        row.push(s..index);
                    }
                    start = None;
                    last_non_whitespace = None;
                    was_enclosed = true;
                }
            }
            (sep, None, Some(s)) if sep == separator => {
                if sep.is_ascii_whitespace() {
                    row.push(s..last_non_whitespace.unwrap_or(index));
                } else {
                    row.push(s..index);
                }
                start = None;
                last_non_whitespace = None;
                was_enclosed = false;
            }
            (sep, None, None) if sep == separator => {
                if !was_enclosed {
                    // ignore any comma directly after an enclosed field
                    row.push(index..index);
                    start = None;
                    last_non_whitespace = None;
                }
                was_enclosed = false;
            }
            (c, _, _) if c.is_ascii_whitespace() => (), // ignore
            (_, _, None) => {
                start = Some(index);
                last_non_whitespace = Some(index + 1);
            }
            _ => last_non_whitespace = Some(index + 1),
        }
    }
    if let Some(s) = start {
        row.push(s..last_non_whitespace.unwrap_or(line.len()));
    } else if !was_enclosed {
        row.push(line.len()..line.len());
    }
    Ok(row)
}

impl std::fmt::Display for CsvLine {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(
            f,
            "{}",
            Context::default()
                .line_index(self.line_index as u32)
                .lines(0, &self.line)
                .add_highlights(self.fields.iter().map(|f| (0, f.1.clone())))
        )
    }
}

/// Write a CSV file. It fill empty columns with empty space, ensures the correct amount of columns
/// on each line, and auto wraps any separator containing values and headers in double quotes (").
/// It also replaces any double quotes (") in wrapped fields in single quotes (').
/// # Errors
/// If the `Write` implementation errors.
#[allow(dead_code)]
pub fn write_csv(
    mut f: impl Write,
    data: impl IntoIterator<Item = impl IntoIterator<Item = (String, String)>>,
    separator: char,
) -> Result<(), std::io::Error> {
    let mut order: Vec<String> = Vec::new();
    let sorted: Vec<Vec<String>> = data
        .into_iter()
        .map(|row| {
            let mut new_row = vec![String::new(); order.len()];
            for (mut column, mut value) in row {
                if value.contains(separator) {
                    value = format!("\"{}\"", value.replace('\"', "\'"));
                }
                if let Some(index) = order.iter().position(|i| *i == column) {
                    new_row[index] = value;
                } else {
                    if column.contains(separator) {
                        column = format!("\"{}\"", column.replace('\"', "\'"));
                    }
                    order.push(column);
                    new_row.push(value);
                }
            }
            new_row
        })
        .collect_vec();
    let separator = separator.to_string();
    writeln!(f, "{}", order.iter().join(&separator))?;
    for row in sorted {
        let len = order.len() - row.len();
        writeln!(
            f,
            "{}",
            row.into_iter()
                .chain(std::iter::repeat_n(String::new(), len))
                .join(&separator)
        )?;
    }
    Ok(())
}
