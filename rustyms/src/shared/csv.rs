//! Methods for reading and parsing CSV files. (Internal use mostly).

use std::{
    collections::{BTreeMap, HashMap},
    fmt::Debug,
    fs::File,
    io::{BufRead, BufReader, Write},
    ops::Range,
    str::FromStr,
    sync::Arc,
};

use flate2::bufread::GzDecoder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::check_extension,
};

/// A single line in a CSV file
#[derive(Clone, Eq, PartialEq, Hash, Debug, Serialize, Deserialize)]
pub struct CsvLine {
    line_index: usize,
    line: String,
    fields: Vec<(Arc<String>, Range<usize>)>,
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
    pub fn number_of_columns(&self) -> usize {
        self.fields.len()
    }
    /// Get the context applicable to the specified column
    pub fn column_context(&self, column: usize) -> crate::error::Context {
        crate::error::Context::line(
            Some(self.line_index),
            self.line.clone(),
            self.fields[column].1.start,
            self.fields[column].1.len(),
        )
    }
    /// Get the context for the specified range in the original line
    pub fn range_context(&self, range: std::ops::Range<usize>) -> crate::error::Context {
        crate::error::Context::line(
            Some(self.line_index()),
            self.line.clone(),
            range.start,
            range.len(),
        )
    }
    /// Get the context for the whole line
    pub fn full_context(&self) -> crate::error::Context {
        crate::error::Context::full_line(self.line_index, self.line.clone())
    }
    /// Get the range of a specified column
    pub fn range(&self, index: usize) -> &Range<usize> {
        &self.fields[index].1
    }
    /// Get the specified column, by column name
    /// # Errors
    /// If the given name is not a column header return an error
    pub fn index_column(&self, name: &str) -> Result<(&str, &Range<usize>), CustomError> {
        self.fields
            .iter()
            .find(|f| *f.0 == *name)
            .map(|f| (&self.line[f.1.clone()], &f.1))
            .ok_or_else(|| {
                CustomError::error(
                    "Could not find given column",
                    format!("This CSV file does not contain the needed column '{name}'"),
                    self.full_context(),
                )
            })
    }

    /// Parse a column into the given format
    /// # Errors
    /// If erroneous extend the base error with the correct context and return that
    pub fn parse_column<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<F, CustomError> {
        self[column]
            .parse()
            .map_err(|_| base_error.with_context(self.column_context(column)))
    }
    /// Parse a column into the given format
    /// # Errors
    /// If erroneous extend the base error with the correct context and return that
    pub fn parse_column_or_empty<F: FromStr>(
        &self,
        column: usize,
        base_error: &CustomError,
    ) -> Result<Option<F>, CustomError> {
        let text = &self[column];
        if text.is_empty() || text == "-" {
            Ok(None)
        } else {
            Ok(Some(text.parse().map_err(|_| {
                base_error.with_context(self.column_context(column))
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
) -> Result<Box<dyn Iterator<Item = Result<CsvLine, CustomError>>>, CustomError> {
    let file = File::open(path.as_ref()).map_err(|e| {
        CustomError::error(
            "Could not open file",
            e,
            crate::error::Context::Show {
                line: path.as_ref().to_string_lossy().to_string(),
            },
        )
    })?;
    if check_extension(path, "gz") {
        Ok(Box::new(parse_csv_raw(
            GzDecoder::new(BufReader::new(file)),
            separator,
            provided_header,
        )?))
    } else {
        Ok(Box::new(parse_csv_raw(file, separator, provided_header)?))
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
) -> Result<CsvLineIter<T>, CustomError> {
    let reader = BufReader::new(reader);
    let mut lines = reader.lines().enumerate().peekable();
    let mut skip = false;
    if let Some(sep) = lines
        .peek()
        .and_then(|(_, l)| l.as_ref().ok())
        .and_then(|l| l.strip_prefix("sep="))
    {
        skip = true;
        if let Some(c) = sep.chars().next() {
            if c.len_utf8() == 1 {
                separator = c as u8;
            } else {
                return Err(CustomError::error("Unicode value separators not supported", "This is a character that takes more than 1 byte to represent in Unicode, this is not supported in parsing CSV files.", Context::line(Some(0), format!("sep={sep}"), 4, sep.len())));
            }
        }
    }
    if skip {
        // Actually consume this line
        let _ = lines.next();
    }
    let column_headers = if let Some(header) = provided_header {
        header.into_iter().map(Arc::new).collect()
    } else {
        let (_, column_headers) = lines.next().ok_or_else(|| {
            CustomError::error("Could parse csv file", "The file is empty", Context::None)
        })?;
        let header_line = column_headers
            .map_err(|err| CustomError::error("Could not read header line", err, Context::None))?;
        csv_separate(&header_line, separator)?
            .into_iter()
            .map(|r| Arc::new(header_line[r].to_lowercase()))
            .collect()
    };

    Ok(CsvLineIter {
        lines,
        header: column_headers,
        separator,
    })
}

/// An iterator returning CSV lines
pub struct CsvLineIter<T: std::io::Read> {
    lines: std::iter::Peekable<std::iter::Enumerate<std::io::Lines<BufReader<T>>>>,
    header: Vec<Arc<String>>,
    separator: u8,
}

impl<T: std::io::Read> Iterator for CsvLineIter<T> {
    type Item = Result<CsvLine, CustomError>;
    fn next(&mut self) -> Option<Self::Item> {
        self.lines.next().map(|(line_index, line)| {
            let line = line.map_err(|err|CustomError::error(
                    "Could not read line",
                    err,
                    Context::full_line(line_index, "(failed)"),
                ))?;
            csv_separate(&line, self.separator).and_then(|row| {
                if self.header.len() == row.len() {
                    Ok(CsvLine {
                        line_index,
                        line,
                        fields: self.header.iter().cloned().zip(row).collect(),
                    })
                } else {
                    Err(CustomError::error(
                        "Incorrect number of columns",
                        format!("It does not have the correct number of columns. {} columns were expected but {} were found.", self.header.len(), row.len()),
                        Context::full_line(line_index, line),
                    ))
                }
            })
        })
    }
}

/// # Errors
/// If the line is empty.
pub(crate) fn csv_separate(line: &str, separator: u8) -> Result<Vec<Range<usize>>, CustomError> {
    if line.is_empty() {
        return Err(CustomError::error(
            "Empty line",
            "The line is empty",
            Context::None,
        ));
    }
    let mut enclosed = None;
    let mut was_enclosed = false;
    let mut row = Vec::new();
    let mut start = None;
    let mut last_non_whitespace = None;
    for (index, ch) in line.bytes().enumerate() {
        match (ch, enclosed, start) {
            (b'\"' | b'\'', None, None) => {
                enclosed = Some(ch);
                start = Some(index + 1);
            }
            (c, Some(e), Some(s)) if c == e => {
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
        let mut line = self.line.clone();
        let mut offset = 0;
        let mut display_field = |range: &Range<usize>| match range.len() {
            0 => {
                line.insert(range.start + offset, '⌧'); // Insert an empty icon to be able to indicate a field of length 0
                offset += '⌧'.len_utf8();
                "ô".to_string()
            }
            1 => "^".to_string(),
            n => {
                let full = format!("{}-{}", range.start, range.end);
                if full.len() <= n - 2 {
                    format!("└{:·^w$}┘", full, w = n - 2)
                } else {
                    format!("└{}┘", "·".repeat(n - 2))
                }
            }
        };
        let mut fields = String::new();
        let mut index = 0;
        for (_name, field) in &self.fields {
            if field.start >= index {
                fields = format!(
                    "{fields}{}{}",
                    " ".repeat(field.start - index),
                    display_field(field)
                );
                index = field.end;
            } else {
                fields = format!(
                    "{fields}\n{}{}",
                    " ".repeat(field.start),
                    display_field(field)
                );
            }
        }
        write!(
            f,
            "{}: {}\n{}{}",
            self.line_index + 1,
            line,
            " ".repeat((self.line_index + 1).to_string().len() + 2),
            fields
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
                .chain(std::iter::repeat(String::new()).take(len))
                .join(&separator)
        )?;
    }
    Ok(())
}
