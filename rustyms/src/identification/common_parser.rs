use crate::{
    error::{Context, CustomError},
    identification::csv::CsvLine,
};
use std::{ops::Range, str::FromStr};

macro_rules! format_family {
    (#[doc = $format_doc:expr]
     $format:ident,
     #[doc = $data_doc:expr]
     $data:ident,
     $version:ident, $versions:expr, $separator:expr, $header:expr;
     required { $($(#[doc = $rdoc:expr])? $rname:ident: $rtyp:ty, $rf:expr;)* }
     optional { $($(#[doc = $odoc:expr])? $oname:ident: $otyp:ty, $of:expr;)*}
     $($post_process:item)?) => {
        use super::super::common_parser::HasLocation;

        #[non_exhaustive]
        #[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Default, Serialize, Deserialize)]
        #[doc = $format_doc]
        pub struct $format {
            $($rname: &'static str,)*
            $($oname: crate::identification::common_parser::OptionalColumn,)*
            version: $version
        }

        #[non_exhaustive]
        #[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
        #[doc = $data_doc]
        #[expect(missing_docs)]
        pub struct $data {
            $($(#[doc = $rdoc])? pub $rname: $rtyp,)*
            $($(#[doc = $odoc])? pub $oname: Option<$otyp>,)*
            /// The version used to read in the data
            pub version: $version,
            /// The stored columns if kept
            columns: Option<Vec<(std::sync::Arc<String>, String)>>,
        }

        impl IdentifiedPeptidoformSource for $data {
            type Source = CsvLine;
            type Format = $format;
            type Version = $version;
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
                version: Option<$version>,
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
                version: Option<$version>,
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
            fn parse_specific(source: &Self::Source, format: &$format, custom_database: Option<&crate::ontology::CustomDatabase>, keep_all_columns: bool) -> Result<Self, CustomError> {
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

        impl $data {
            /// Get all original columns from the CSV file, only available if the file was opened with the option `keep_all_columns` turned on.
            pub fn full_csv_line(&self) -> Option<&[(std::sync::Arc<String>, String)]> {
                self.columns.as_deref()
            }
        }
    };
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
    fn column<'a>(&'a self, name: &str) -> Result<Location<'a>, CustomError>;
}

impl HasLocation for CsvLine {
    /// Get the specified column
    /// # Errors
    /// If the given column does not exist
    fn column<'a>(&'a self, name: &str) -> Result<Location<'a>, CustomError> {
        self.index_column(name).map(|(_v, c)| Location {
            line: self,
            location: c.clone(),
        })
    }
}

/// The base location type to keep track of the location of to be parsed pieces in the monadic parser combinators below
#[derive(Clone)]
pub(super) struct Location<'a> {
    pub(super) line: &'a CsvLine,
    pub(super) location: Range<usize>,
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
        if text == pattern {
            None
        } else {
            Some(self)
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
                self.line.range_context(self.location),
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
                    }
                    .parse(base_error)?,
                ),
                Self {
                    line: self.line,
                    location: self.location.start + start.len() + 1
                        ..self.location.start + start.len() + 1 + end.len(),
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
        Context::line_range(
            Some(self.line.line_index()),
            self.full_line(),
            self.location.clone(),
        )
    }

    pub(super) fn trim(&self) -> Self {
        let start_trim = self.as_str().trim_start().len();
        let end_trim = self.as_str().trim_end().len();
        let length = self.as_str().len();

        Self {
            line: self.line,
            location: self.location.start + (length - start_trim)
                ..self.location.end - (length - end_trim),
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
                },
                Self {
                    line: self.line,
                    location: self.location.end - end.len()..self.location.end,
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
