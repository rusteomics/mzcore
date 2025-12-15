use std::{borrow::Cow, ops::Range, str::FromStr};

use context_error::*;

use mzcore::csv::CsvLine;

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
/// * Do not forget to implement `PSMMetaData` for the `<$format>PSM` type, this will set up the logic
///   for an automatic `From<format> for PSMData` implementation.
/// * Do not forget to create a `<$format>Version` enum which contains all supported version of the
///   format, which additionally should implement `PSMFileFormatVersion`.
/// * For each version a public constant should be generated that contains an instantiation of the
///   `<$format>Format` type, these are also the ones that need to be listed in the versions list.
///
macro_rules! format_family {
     ($(#[doc = $ddoc:expr])*
     $format:ident,
     $complexity:ident, $peptidoform_availability:ident, $versions:expr, $separator:expr, $header:expr;
     required { $($(#[doc = $rdoc:expr])? $rname:ident: $rtyp:ty, $rf:expr;)* }
     optional { $($(#[doc = $odoc:expr])? $(#[cfg(feature = $ocfg:literal)])? $oname:ident: $otyp:ty, $of:expr;)*}
     $($post_process:item)?
     $(protein {
            $accession_name:ident => $accessionf:expr;
            required { $($(#[doc = $prdoc:expr])? $prname:ident: $prtyp:ty, $prf:expr;)* }
            optional { $($(#[doc = $podoc:expr])? $(#[cfg(feature = $pocfg:literal)])? $poname:ident: $potyp:ty, $pof:expr;)*}
    })?
    ) => {paste::paste!{
        #[allow(unused_imports)] // Needed sometimes, but not all invocations of the macro
        use context_error::*;

        use super::super::common_parser::HasLocation;

        #[doc = "The type to contain the format description for " $format " files."]
        #[non_exhaustive]
        #[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
        pub struct [<$format Format>] {
            $($rname: &'static str,)*
            $($(#[cfg(feature = $ocfg)])?  $oname: crate::common_parser::OptionalColumn,)*
            $(
                $accession_name: &'static str,
                $($prname: &'static str,)*
                $($(#[cfg(feature = $pocfg)])?  $poname: crate::common_parser::OptionalColumn,)*
            )?
            version: [<$format Version>]
        }

        #[doc = "The data for individual entries in " $format " files."]
        $(#[doc = $ddoc])*
        #[non_exhaustive]
        #[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
        #[allow(missing_docs)]
        pub struct [<$format PSM>] {
            $($(#[doc = $rdoc])? pub $rname: $rtyp,)*
            $($(#[doc = $odoc])? $(#[cfg(feature = $ocfg)])?  pub $oname: Option<$otyp>,)*
            /// The version used to read in the data
            pub version: [<$format Version>],
            /// The stored columns if kept
            columns: Option<Vec<(std::sync::Arc<String>, String)>>,
            $(
            pub $accession_name: std::sync::Arc<<[<$format PSM>] as PSMMetaData>::Protein>
            )?
        }

        $(
            #[doc = "The data for proteins in " $format " files."]
            #[non_exhaustive]
            #[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
            #[allow(missing_docs)]
            pub struct [<$format Protein>] {
                pub $accession_name: String,
                $($(#[doc = $prdoc])? pub $prname: $prtyp,)*
                $($(#[doc = $podoc])? $(#[cfg(feature = $pocfg)])?  pub $poname: Option<$potyp>,)*
            }

            impl mzcore::space::Space for [<$format Protein>] {
                fn space(&self) -> mzcore::space::UsedSpace {
                    ( mzcore::space::UsedSpace::default()
                    $(+ self.$prname.space())*
                    $(+ self.$poname.space())*)
                    .set_total::<Self>()
                }
            }

            impl From<[<$format Protein>]> for crate::ProteinData {
            fn from(value: [<$format Protein>]) -> Self {
                Self::$format(value)
            }
        }
        )?

        impl PSMSource for [<$format PSM>] {
            type Source = CsvLine;
            type Format = [<$format Format>];
            type Complexity = $complexity;
            type PeptidoformAvailability = $peptidoform_availability;
            type Version = [<$format Version>];

            fn parse(source: &Self::Source, ontologies: &mzcore::ontology::Ontologies, keep_all_columns: bool, proteins: &mut std::collections::HashMap<String, std::sync::Arc<<[<$format PSM>] as PSMMetaData>::Protein>>) -> Result<(Self, &'static Self::Format), BoxedError<'static, BasicKind>> {
                let mut errors = Vec::new();
                for format in $versions {
                    match Self::parse_specific(source, format, ontologies, keep_all_columns, proteins) {
                        Ok(peptide) => return Ok((peptide, format)),
                        Err(err) => errors.push(err.version(format.version.to_string())),
                    }
                }
                Err(BoxedError::new(BasicKind::Error,
                    format!("Invalid {} line", stringify!($format)),
                    "The correct format could not be determined automatically",
                    source.full_context().to_owned(),
                ).add_underlying_errors(errors))
            }

            fn parse_file(
                path: impl AsRef<std::path::Path>,
                ontologies: &mzcore::ontology::Ontologies,
                keep_all_columns: bool,
                version: Option<Self::Version>,
            ) -> Result<BoxedIdentifiedPeptideIter<'_, Self>, BoxedError<'static, BasicKind>> {
                let format = version.map(|v| v.format());
                parse_csv(path, $separator, $header).and_then(|lines| {
                    let mut i = Self::parse_many::<Box<dyn Iterator<Item = Result<Self::Source, BoxedError<'_, BasicKind>>>>>(
                        Box::new(lines), ontologies, keep_all_columns, format);
                    if let Some(Err(e)) = i.peek() {
                        Err(e.clone())
                    } else {
                        Ok(i)
                    }
                })
            }

            fn parse_reader<'a>(
                reader: impl std::io::Read + 'a,
                ontologies: &'a mzcore::ontology::Ontologies,
                keep_all_columns: bool,
                version: Option<Self::Version>,
            ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, BoxedError<'static, BasicKind>> {
                let format = version.map(|v| v.format());
                mzcore::csv::parse_csv_raw(reader, $separator, $header, None).and_then(move |lines| {
                    let mut i = Self::parse_many::<Box<dyn Iterator<Item = Result<Self::Source, BoxedError<'_, BasicKind>>>>>(
                        Box::new(lines), ontologies, keep_all_columns, format);
                    if let Some(Err(e)) = i.peek() {
                        Err(e.clone())
                    } else {
                        Ok(i)
                    }
                })
            }

            #[allow(clippy::redundant_closure_call, unused_variables)] // Macro magic
            fn parse_specific(source: &Self::Source, format: &[<$format Format>], ontologies: &mzcore::ontology::Ontologies, keep_all_columns: bool, proteins: &mut std::collections::HashMap<String, std::sync::Arc<<[<$format PSM>] as PSMMetaData>::Protein>>) -> Result<Self, BoxedError<'static, BasicKind>> {
                #[allow(unused_imports)]
                use crate::helper_functions::InvertResult;

                /// Shadowing `Result::Ok` to inject the correct error type which otherwise leads
                /// to the compiler complaining about its absence when a field cannot fail.
                #[allow(non_snake_case, dead_code, clippy::missing_errors_doc)]
                const fn Ok<T>(value: T) -> Result<T, BoxedError<'static, BasicKind>> {
                    Result::Ok(value)
                }

                $(
                    let key = $accessionf(source.column(format.$accession_name).map_err(BoxedError::to_owned)?, ontologies)?;

                    let protein = if let Some(p) = proteins.get(&key) {
                        p.clone()
                    } else {
                        // Get around parser limitations by creating a type synonym.
                        type TempProtein = <[<$format PSM>] as PSMMetaData>::Protein;
                        let parsed = std::sync::Arc::new(TempProtein {
                            $accession_name: key.clone(),
                            $($prname: $prf(source.column(format.$prname).map_err(BoxedError::to_owned)?, ontologies)?,)*
                            $($(#[cfg(feature = $pocfg)])?  $poname: format.$poname.open_column(source).and_then(|l: Option<Location>| l.map(|value: Location| $pof(value, ontologies)).invert()).map_err(BoxedError::to_owned)?,)*
                        });
                        proteins.insert(key, parsed.clone());
                        parsed
                    };
                )?

                let parsed = Self {
                    $($rname: $rf(source.column(format.$rname).map_err(BoxedError::to_owned)?, ontologies)?,)*
                    $($(#[cfg(feature = $ocfg)])?  $oname: format.$oname.open_column(source).and_then(|l: Option<Location>| l.map(|value: Location| $of(value, ontologies)).invert()).map_err(BoxedError::to_owned)?,)*
                    version: format.version.clone(),
                    columns: keep_all_columns.then(|| source.values().map(|(h, v)| (h, v.to_string())).collect()),
                    $($accession_name: protein,)?
                };
                Self::post_process(source, parsed, ontologies)
            }
            $($post_process)?
        }

        impl [<$format PSM>] {
            /// Get all original columns from the CSV file, only available if the file was opened with the option `keep_all_columns` turned on.
            pub fn full_csv_line(&self) -> Option<&[(std::sync::Arc<String>, String)]> {
                self.columns.as_deref()
            }
        }

        impl From<[<$format PSM>]> for PSM<$complexity, $peptidoform_availability> {
            fn from(value: [<$format PSM>]) -> Self {
                Self {
                    score: value.confidence(),
                    local_confidence: value.local_confidence().map(|v| v.to_vec()),
                    data: PSMData::$format(value),
                    complexity_marker: PhantomData,
                    peptidoform_availability_marker: PhantomData,
                }
            }
        }

        impl mzcore::space::Space for [<$format PSM>] {
            fn space(&self) -> mzcore::space::UsedSpace {
                ( mzcore::space::UsedSpace::default()
                $(+ self.$rname.space())*
                $(+ self.$oname.space())*)
                .set_total::<Self>()
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
    pub(super) fn open_column(
        self,
        source: &CsvLine,
    ) -> Result<Option<Location<'_>>, BoxedError<'static, BasicKind>> {
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
    fn column<'a>(&'a self, name: &'a str) -> Result<Location<'a>, BoxedError<'static, BasicKind>>;
}

impl HasLocation for CsvLine {
    /// Get the specified column
    /// # Errors
    /// If the given column does not exist
    fn column<'a>(&'a self, name: &'a str) -> Result<Location<'a>, BoxedError<'static, BasicKind>> {
        self.index_column(name)
            .map(|(_v, c)| Location {
                line: self,
                location: c.clone(),
                column: Some(name),
            })
            .map_err(BoxedError::to_owned)
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

    pub(super) fn is_empty(&self) -> bool {
        self.location.is_empty()
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

    pub(super) fn get_string(self) -> String {
        self.as_str().to_string()
    }

    pub(super) fn get_boxed_str(self) -> Box<str> {
        self.as_str().into()
    }

    pub(super) fn as_str(&self) -> &str {
        &self.line.line()[self.location.clone()]
    }

    pub(super) fn full_line(&self) -> &str {
        self.line.line()
    }
}

impl<'a> Location<'a> {
    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    pub(super) fn parse<T: FromStr>(
        self,
        base_error: (&'static str, &'static str),
    ) -> Result<T, BoxedError<'static, BasicKind>> {
        self.as_str().trim().parse().map_err(|_| {
            BoxedError::new(
                BasicKind::Error,
                base_error.0,
                base_error.1,
                self.line
                    .range_context(self.location, self.column.map(Cow::Borrowed))
                    .to_owned(),
            )
        })
    }

    /// # Errors
    /// If the provided parse method fails.
    pub(super) fn parse_with<T>(
        self,
        f: impl Fn(Self) -> Result<T, BoxedError<'static, BasicKind>>,
    ) -> Result<T, BoxedError<'static, BasicKind>> {
        f(self)
    }

    pub(super) fn context(&'a self) -> Context<'a> {
        let base = Context::none()
            .line_index(self.line.line_index as u32)
            .lines(0, self.full_line());
        let base = if let Some(comment) = self.column {
            base.add_highlight((0, self.location.clone(), comment))
        } else {
            base.add_highlight((0, self.location.clone()))
        };
        if let Some(source) = &self.line.file {
            base.source(source.as_ref().as_ref())
        } else {
            base
        }
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

    // fn apply(self, f: impl FnOnce(Self) -> Self) -> Self {
    //     f(self)
    // }

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

    /// Split twice on the character, split on the first and last occurrence of the given character.
    /// So any additional occurrences of the characters are in the middle segment.
    pub(super) fn split_twice(self, p: char) -> Option<(Self, Self, Self)> {
        let (start, after) = self.as_str().split_once(p)?;
        let (middle, end) = after.rsplit_once(p)?;
        let start_middle = self.location.start + start.len() + p.len_utf8();
        Some((
            Self {
                line: self.line,
                location: self.location.start..self.location.start + start.len(),
                column: self.column,
            },
            Self {
                line: self.line,
                location: start_middle..start_middle + middle.len(),
                column: self.column,
            },
            Self {
                line: self.line,
                location: self.location.end - end.len()..self.location.end,
                column: self.column,
            },
        ))
    }
}

pub(super) trait OptionalLocation<'a> {
    /// # Errors
    /// If the parse method fails. See [`FromStr::parse`].
    fn parse<T: FromStr>(
        self,
        base_error: (&'static str, &'static str),
    ) -> Result<Option<T>, BoxedError<'static, BasicKind>>;
    /// # Errors
    /// If the provided parse method fails.
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, BoxedError<'static, BasicKind>>,
    ) -> Result<Option<T>, BoxedError<'static, BasicKind>>;
    fn get_string(self) -> Option<String>;
    type ArrayIter: Iterator<Item = Location<'a>>;
    fn array(self, sep: char) -> Self::ArrayIter;
    fn optional_array(self, sep: char) -> Option<Self::ArrayIter>;
    fn ignore(self, pattern: &str) -> Option<Location<'a>>;
}

impl<'a> OptionalLocation<'a> for Option<Location<'a>> {
    fn parse<T: FromStr>(
        self,
        base_error: (&'static str, &'static str),
    ) -> Result<Option<T>, BoxedError<'static, BasicKind>> {
        self.map(|l| l.parse::<T>(base_error)).transpose()
    }
    fn parse_with<T>(
        self,
        f: impl Fn(Location<'a>) -> Result<T, BoxedError<'static, BasicKind>>,
    ) -> Result<Option<T>, BoxedError<'static, BasicKind>> {
        self.map(f).transpose()
    }
    fn get_string(self) -> Option<String> {
        self.map(Location::get_string)
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
