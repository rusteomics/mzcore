use crate::{
    error::{Context, CustomError},
    identification::csv::CsvLine,
};
use std::{ops::Range, str::FromStr};

/// The way to set up a format family.
///
/// It starts with defining the name, complexity, peptidoform availability, all versions (parsing
/// is tried in this order), which separator is used for CSV files (character separated values),
/// and if needed the header for CSV files.
///
/// After that the columns are defined, first the required columns (present in all format versions)
/// followed by the optional columns (missing/optional in at least one version). For each column,
/// the name, type, and the lambda function to parse the text are given. Note that optional column
/// types are automatically wrapped in an `Option` so no additional `Option` should be specified.
///
/// Lastly, a post processing function can be specified.
///
/// Please check the other file formats for how the precisely use the syntax.
///
/// # Notes
/// * Do not forget to implement `MetaData` for the `<$format>Data` type, this will set up the logic
///   for an automatic `From<format> for IdentifiedPeptidoformData` implementation.
/// * Do not forget to create a `<$format>Version` enum which contains all supported version of the
///   format, which additionally should implement `IdentifiedPeptidoformVersion`.
/// * For each version a public constant should be generated that contains an instantiation of the
///   `<$format>Format` type, these are also the ones that need to be listed in the versions list.
///
macro_rules! format_family {
     ($format:ident,
     $complexity:ident, $peptidoform_availability:ident, $versions:expr, $separator:expr, $header:expr;
     required { $($(#[doc = $rdoc:expr])? $rname:ident: $rtyp:ty, $rf:expr;)* }
     optional { $($(#[doc = $odoc:expr])? $oname:ident: $otyp:ty, $of:expr;)*}
     $($post_process:item)?) => {paste::paste!{
        use super::super::common_parser::HasLocation;

        #[doc = "The type to contain the format description for " $format " files."]
        #[non_exhaustive]
        #[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
        pub struct [<$format Format>] {
            $($rname: &'static str,)*
            $($oname: crate::identification::common_parser::OptionalColumn,)*
            version: [<$format Version>]
        }

        #[doc = "The data for individual entries in " $format " files."]
        #[non_exhaustive]
        #[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
        #[allow(missing_docs)]
        pub struct [<$format Data>] {
            $($(#[doc = $rdoc])? pub $rname: $rtyp,)*
            $($(#[doc = $odoc])? pub $oname: Option<$otyp>,)*
            /// The version used to read in the data
            pub version: [<$format Version>],
            /// The stored columns if kept
            columns: Option<Vec<(std::sync::Arc<String>, String)>>,
        }

        impl IdentifiedPeptidoformSource for [<$format Data>] {
            type Source = CsvLine;
            type Format = [<$format Format>];
            type Complexity = $complexity;
            type PeptidoformAvailability = $peptidoform_availability;
            type Version = [<$format Version>];
            fn parse(source: &Self::Source, custom_database: Option<&crate::ontology::CustomDatabase>, keep_all_columns: bool) -> Result<(Self, &'static Self::Format), CustomError> {
                let mut errors = Vec::new();
                for format in $versions {
                    match Self::parse_specific(source, format, custom_database, keep_all_columns) {
                        Ok(peptide) => return Ok((peptide, format)),
                        Err(err) => errors.push(err.with_version(&format.version)),
                    }
                }
                Err(CustomError::error(
                    format!("Invalid {} line", stringify!($format)),
                    "The correct format could not be determined automatically",
                    source.full_context(),
                ).with_underlying_errors(errors))
            }
            fn parse_file(
                path: impl AsRef<std::path::Path>,
                custom_database: Option<&crate::ontology::CustomDatabase>,
                keep_all_columns: bool,
                version: Option<Self::Version>,
            ) -> Result<BoxedIdentifiedPeptideIter<Self>, CustomError> {
                let format = version.map(|v| v.format());
                parse_csv(path, $separator, $header).and_then(|lines| {
                    let mut i = Self::parse_many::<Box<dyn Iterator<Item = Result<Self::Source, CustomError>>>>(
                        Box::new(lines), custom_database, keep_all_columns, format);
                    if let Some(Err(e)) = i.peek() {
                        Err(e.clone())
                    } else {
                        Ok(i)
                    }
                })
            }
            fn parse_reader<'a>(
                reader: impl std::io::Read + 'a,
                custom_database: Option<&'a crate::ontology::CustomDatabase>,
                keep_all_columns: bool,
                version: Option<Self::Version>,
            ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, CustomError> {
                let format = version.map(|v| v.format());
                crate::identification::csv::parse_csv_raw(reader, $separator, $header).and_then(move |lines| {
                    let mut i = Self::parse_many::<Box<dyn Iterator<Item = Result<Self::Source, CustomError>>>>(
                        Box::new(lines), custom_database, keep_all_columns, format);
                    if let Some(Err(e)) = i.peek() {
                        Err(e.clone())
                    } else {
                        Ok(i)
                    }
                })
            }
            #[allow(clippy::redundant_closure_call)] // Macro magic
            fn parse_specific(source: &Self::Source, format: &[<$format Format>], custom_database: Option<&crate::ontology::CustomDatabase>, keep_all_columns: bool) -> Result<Self, CustomError> {
                #[allow(unused_imports)]
                use crate::helper_functions::InvertResult;

                let parsed = Self {
                    $($rname: $rf(source.column(format.$rname)?, custom_database)?,)*
                    $($oname: format.$oname.open_column(source).and_then(|l: Option<Location>| l.map(|value: Location| $of(value, custom_database)).invert())?,)*
                    version: format.version.clone(),
                    columns: keep_all_columns.then(|| source.values().map(|(h, v)| (h, v.to_string())).collect()),
                };
                Self::post_process(source, parsed, custom_database)
            }
            $($post_process)?
        }

        impl [<$format Data>] {
            /// Get all original columns from the CSV file, only available if the file was opened with the option `keep_all_columns` turned on.
            pub fn full_csv_line(&self) -> Option<&[(std::sync::Arc<String>, String)]> {
                self.columns.as_deref()
            }
        }

        impl From<[<$format Data>]> for IdentifiedPeptidoform<$complexity, $peptidoform_availability> {
            fn from(value: [<$format Data>]) -> Self {
                Self {
                    score: value.confidence(),
                    local_confidence: value.local_confidence().map(|v| v.to_vec()),
                    data: IdentifiedPeptidoformData::$format(value),
                    complexity_marker: PhantomData,
                    peptidoform_availability_marker: PhantomData,
                }
            }
        }
    }};
}

/// The possible options for an optional column
#[derive(
    Copy,
    Clone,
    PartialEq,
    Eq,
    PartialOrd,
    Ord,
    Hash,
    Debug,
    Default,
    serde::Serialize,
    serde::Deserialize,
)]
pub(super) enum OptionalColumn {
    /// This column is not avalable in this version
    #[default]
    NotAvailable,
    /// This column is optional in this version
    Optional(&'static str),
    /// This column is required in this version (but as it is an `OptionalColumn` not in some other version)
    Required(&'static str),
}

impl OptionalColumn {
    /// Open the column
    /// # Errors
    /// while creating the correct error messages for missing columns
    pub(super) fn open_column(self, source: &CsvLine) -> Result<Option<Location>, CustomError> {
        match self {
            Self::NotAvailable => Ok(None),
            Self::Optional(s) => Ok(source.column(s).ok()),
            Self::Required(s) => source.column(s).map(Some),
        }
    }
}

pub(super) trait HasLocation {
    /// Get the specified column.
    /// # Errors
    /// If the column does not exist.
    fn column<'a>(&'a self, name: &'a str) -> Result<Location<'a>, CustomError>;
}

impl HasLocation for CsvLine {
    /// Get the specified column
    /// # Errors
    /// If the given column does not exist
    fn column<'a>(&'a self, name: &'a str) -> Result<Location<'a>, CustomError> {
        self.index_column(name).map(|(_v, c)| Location {
            line: self,
            location: c.clone(),
            column: Some(name),
        })
    }
}

/// The base location type to keep track of the location of to be parsed pieces in the monadic parser combinators below
#[derive(Clone)]
pub(super) struct Location<'a> {
    pub(super) line: &'a CsvLine,
    pub(super) location: Range<usize>,
    pub(super) column: Option<&'a str>,
}

impl Location<'_> {
    pub(super) fn len(&self) -> usize {
        self.location.len()
    }

    pub(super) fn array(self, sep: char) -> std::vec::IntoIter<Self> {
        let mut offset = 0;
        let mut output = Vec::new();
        for part in self.as_str().split(sep) {
            output.push(Location {
                line: self.line,
                location: self.location.start + offset..self.location.start + offset + part.len(),
                column: self.column,
            });
            offset += part.len() + 1;
        }
        output.into_iter()
    }

    pub(super) fn or_empty(self) -> Option<Self> {
        let text = self.as_str();
        if text.is_empty() || text == "-" {
            None
        } else {
            Some(self)
        }
    }

    pub(super) fn ignore(self, pattern: &str) -> Option<Self> {
        let text = self.as_str();
        if text == pattern { None } else { Some(self) }
    }

    pub(super) fn skip(self, bytes: usize) -> Self {
        Self {
            line: self.line,
            location: self
                .location
                .start
                .saturating_add(bytes)
                .min(self.location.end)..self.location.end,
            column: self.column,
        }
    }

    pub(super) fn trim_end_matches(mut self, pattern: &str) -> Self {
        let trimmed = self.as_str().trim_end_matches(pattern);
        let dif = self.location.len() - trimmed.len();
        self.location = self.location.start..self.location.end - dif;
        self
    }

    pub(super) fn trim_start_matches(mut self, pattern: &str) -> Self {
        let trimmed = self.as_str().trim_start_matches(pattern);
        let dif = self.location.len() - trimmed.len();
        self.location = self.location.start + dif..self.location.end;
        self
    }

    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    pub(super) fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<T, CustomError> {
        self.as_str().trim().parse().map_err(|_| {
            CustomError::error(
                base_error.0,
                base_error.1,
                self.line
                    .range_context(self.location, self.column.map(ToString::to_string)),
            )
        })
    }

    /// # Errors
    /// If the provided parse method fails.
    pub(super) fn parse_with<T>(
        self,
        f: impl Fn(Self) -> Result<T, CustomError>,
    ) -> Result<T, CustomError> {
        f(self)
    }

    /// # Errors
    /// If the text could not be read as a valid id.
    pub(super) fn get_id(
        self,
        base_error: (&str, &str),
    ) -> Result<(Option<usize>, usize), CustomError> {
        if let Some((start, end)) = self.as_str().split_once(':') {
            Ok((
                Some(
                    Self {
                        line: self.line,
                        location: self.location.start + 1..self.location.start + start.len(),
                        column: self.column,
                    }
                    .parse(base_error)?,
                ),
                Self {
                    line: self.line,
                    location: self.location.start + start.len() + 1
                        ..self.location.start + start.len() + 1 + end.len(),
                    column: self.column,
                }
                .parse(base_error)?,
            ))
        } else {
            Ok((None, self.parse(base_error)?))
        }
    }

    pub(super) fn get_string(self) -> String {
        self.as_str().to_string()
    }

    pub(super) fn as_str(&self) -> &str {
        &self.line.line()[self.location.clone()]
    }

    pub(super) fn full_line(&self) -> &str {
        self.line.line()
    }

    pub(super) fn context(&self) -> Context {
        Context::line_range_with_comment(
            Some(self.line.line_index()),
            self.full_line(),
            self.location.clone(),
            self.column.map(ToString::to_string),
        )
    }

    pub(super) fn trim(&self) -> Self {
        let str = self.as_str();
        let length = str.len();
        let trimmed_start = length - str.trim_start().len();
        let trimmed_end = length - str.trim_end().len();

        Self {
            line: self.line,
            location: if self.location.start + trimmed_end > self.location.end - trimmed_end {
                self.location.start..self.location.start
            } else {
                self.location.start + trimmed_start..self.location.end - trimmed_end
            },
            column: self.column,
        }
    }

    pub(super) fn apply(self, f: impl FnOnce(Self) -> Self) -> Self {
        f(self)
    }

    pub(super) fn split_once(self, p: char) -> Option<(Self, Self)> {
        self.as_str().split_once(p).map(|(start, end)| {
            (
                Self {
                    line: self.line,
                    location: self.location.start..self.location.start + start.len(),
                    column: self.column,
                },
                Self {
                    line: self.line,
                    location: self.location.end - end.len()..self.location.end,
                    column: self.column,
                },
            )
        })
    }
}

#[expect(dead_code)]
pub(super) trait OptionalLocation<'a> {
    fn or_empty(self) -> Option<Location<'a>>;
    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<Option<T>, CustomError>;
    /// # Errors
    /// If the provided parse method fails.
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError>;
    /// # Errors
    /// If the text could not be read as a valid id.
    fn get_id(
        self,
        base_error: (&str, &str),
    ) -> Result<Option<(Option<usize>, usize)>, CustomError>;
    fn get_string(self) -> Option<String>;
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Option<Location<'a>>;
    type ArrayIter: Iterator<Item = Location<'a>>;
    fn array(self, sep: char) -> Self::ArrayIter;
    fn optional_array(self, sep: char) -> Option<Self::ArrayIter>;
    fn ignore(self, pattern: &str) -> Option<Location<'a>>;
}

impl<'a> OptionalLocation<'a> for Option<Location<'a>> {
    fn or_empty(self) -> Self {
        self.and_then(Location::or_empty)
    }
    fn parse<T: FromStr>(self, base_error: (&str, &str)) -> Result<Option<T>, CustomError> {
        self.map(|l| l.parse::<T>(base_error)).transpose()
    }
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, CustomError>,
    ) -> Result<Option<T>, CustomError> {
        self.map(f).transpose()
    }
    fn get_id(
        self,
        base_error: (&str, &str),
    ) -> Result<Option<(Option<usize>, usize)>, CustomError> {
        self.map(|l| l.get_id(base_error)).transpose()
    }
    fn get_string(self) -> Option<String> {
        self.map(Location::get_string)
    }
    fn apply(self, f: impl FnOnce(Location<'a>) -> Location<'a>) -> Self {
        self.map(|s| s.apply(f))
    }
    type ArrayIter = std::vec::IntoIter<Location<'a>>;
    fn array(self, sep: char) -> Self::ArrayIter {
        self.map(|l| l.array(sep)).unwrap_or_default()
    }
    fn optional_array(self, sep: char) -> Option<Self::ArrayIter> {
        self.map(|l| l.array(sep))
    }
    fn ignore(self, pattern: &str) -> Self {
        self.and_then(|s| s.ignore(pattern))
    }
}
