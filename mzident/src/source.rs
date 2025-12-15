use std::fmt::Display;

use context_error::{BasicKind, BoxedError};

use crate::{MaybePeptidoform, PSM, PSMMetaData};
use mzcore::sequence::{AtLeast, Linked};

/// A version for an PSM file format
pub trait PSMFileFormatVersion<Format>: Copy + Display {
    /// The format for this version
    fn format(self) -> Format;
    /// The name for this version
    fn name(self) -> &'static str;
}

/// The required methods for any source of PSMs
pub trait PSMSource
where
    Self: Sized + PSMMetaData,
{
    /// The source data where the peptides are parsed form
    type Source;

    /// The format type
    type Format: Clone;

    /// The complexity marker type
    type Complexity;

    /// The peptidoform availability marker type
    type PeptidoformAvailability;

    /// The version type
    type Version: Display + PSMFileFormatVersion<Self::Format>;

    /// Parse a single PSM from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(
        source: &Self::Source,
        ontologies: &mzcore::ontology::Ontologies,
        keep_all_columns: bool,
    ) -> Result<(Self, &'static Self::Format), BoxedError<'static, BasicKind>>;

    /// Parse a single PSM with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &Self::Format,
        ontologies: &mzcore::ontology::Ontologies,
        keep_all_columns: bool,
    ) -> Result<Self, BoxedError<'static, BasicKind>>;

    /// Parse a source of multiple peptides using the given format or automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Result<Self::Source, BoxedError<'static, BasicKind>>>>(
        iter: I,
        ontologies: &mzcore::ontology::Ontologies,
        keep_all_columns: bool,
        format: Option<Self::Format>,
    ) -> PSMIter<'_, Self, I> {
        PSMIter {
            iter: Box::new(iter),
            format,
            ontologies,
            keep_all_columns,
            peek: None,
        }
    }

    /// Parse a file with PSMs.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
        ontologies: &mzcore::ontology::Ontologies,
        keep_all_columns: bool,
        version: Option<Self::Version>,
    ) -> Result<BoxedIdentifiedPeptideIter<'_, Self>, BoxedError<'static, BasicKind>>;

    /// Parse a reader with PSMs.
    /// # Errors
    /// When the file is empty or no headers are present.
    fn parse_reader<'a>(
        reader: impl std::io::Read + 'a,
        ontologies: &'a mzcore::ontology::Ontologies,
        keep_all_columns: bool,
        version: Option<Self::Version>,
    ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, BoxedError<'static, BasicKind>>;

    /// Allow post processing of the peptide
    /// # Errors
    /// On errors in the post processing, format specific
    #[expect(unused_variables)]
    fn post_process(
        source: &Self::Source,
        parsed: Self,
        ontologies: &mzcore::ontology::Ontologies,
    ) -> Result<Self, BoxedError<'static, BasicKind>> {
        Ok(parsed)
    }
}

/// A general generic PSM iterator from any source format
pub type GeneralPSMs<'lifetime> = Box<
    dyn Iterator<Item = Result<PSM<Linked, MaybePeptidoform>, BoxedError<'static, BasicKind>>>
        + 'lifetime,
>;

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = PSMIter<
    'lifetime,
    T,
    Box<
        dyn Iterator<Item = Result<<T as PSMSource>::Source, BoxedError<'static, BasicKind>>>
            + 'lifetime,
    >,
>;

/// An iterator returning parsed PSMs
#[derive(Debug)]
pub struct PSMIter<
    'lifetime,
    Source: PSMSource,
    Iter: Iterator<Item = Result<Source::Source, BoxedError<'static, BasicKind>>>,
> {
    iter: Box<Iter>,
    format: Option<Source::Format>,
    ontologies: &'lifetime mzcore::ontology::Ontologies,
    keep_all_columns: bool,
    peek: Option<Result<Source, BoxedError<'static, BasicKind>>>,
}

impl<R: PSMSource + Clone, I: Iterator<Item = Result<R::Source, BoxedError<'static, BasicKind>>>>
    PSMIter<'_, R, I>
where
    R::Format: 'static,
{
    /// Peek at the next item in the iterator
    pub fn peek(&mut self) -> Option<Result<R, BoxedError<'static, BasicKind>>> {
        if self.peek.is_some() {
            return self.peek.clone();
        }

        let peek = if let Some(format) = &self.format {
            self.iter.next().map(|source| {
                R::parse_specific(&source?, format, self.ontologies, self.keep_all_columns)
                    .map_err(BoxedError::to_owned)
            })
        } else {
            match self.iter.next().map(|source| {
                R::parse(&source?, self.ontologies, self.keep_all_columns)
                    .map_err(BoxedError::to_owned)
            }) {
                None => None,
                Some(Ok((pep, format))) => {
                    self.format = Some(format.clone());
                    Some(Ok(pep))
                }
                Some(Err(e)) => Some(Err(e)),
            }
        };
        self.peek.clone_from(&peek);
        peek
    }
}

impl<R: PSMSource, I: Iterator<Item = Result<R::Source, BoxedError<'static, BasicKind>>>> Iterator
    for PSMIter<'_, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, BoxedError<'static, BasicKind>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek.is_some() {
            return self.peek.take();
        }

        if let Some(format) = &self.format {
            self.iter.next().map(|source| {
                R::parse_specific(&source?, format, self.ontologies, self.keep_all_columns)
                    .map_err(BoxedError::to_owned)
            })
        } else {
            match self.iter.next().map(|source| {
                R::parse(&source?, self.ontologies, self.keep_all_columns)
                    .map_err(BoxedError::to_owned)
            }) {
                None => None,
                Some(Ok((pep, format))) => {
                    self.format = Some(format.clone());
                    Some(Ok(pep))
                }
                Some(Err(e)) => Some(Err(e)),
            }
        }
    }
}

impl<'lifetime, Source, Iter> PSMIter<'lifetime, Source, Iter>
where
    Source: PSMSource + Into<PSM<Source::Complexity, Source::PeptidoformAvailability>> + 'lifetime,
    Iter: Iterator<Item = Result<Source::Source, BoxedError<'static, BasicKind>>> + 'lifetime,
    Source::Format: 'static,
    MaybePeptidoform: From<Source::PeptidoformAvailability>,
{
    /// Make this into a generic boxed iterator for merging with other PSM file formats
    pub fn into_box<Complexity: AtLeast<Source::Complexity>>(
        self,
    ) -> Box<
        dyn Iterator<
                Item = Result<PSM<Complexity, MaybePeptidoform>, BoxedError<'static, BasicKind>>,
            > + 'lifetime,
    > {
        Box::new(self.map(
            |p: Result<Source, BoxedError<'static, BasicKind>>| match p {
                Ok(p) => Ok(p.into().cast()),
                Err(e) => Err(e),
            },
        ))
    }
}
