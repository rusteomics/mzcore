use std::fmt::Display;

use custom_error::BoxedError;

use crate::{
    identification::{IdentifiedPeptidoform, MaybePeptidoform, MetaData},
    ontology::CustomDatabase,
    sequence::{AtLeast, Linked},
};

/// A version for an identified peptide version
pub trait IdentifiedPeptidoformVersion<Format>: Copy + Display {
    /// The format for this version
    fn format(self) -> Format;
    /// The name for this version
    fn name(self) -> &'static str;
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptidoformSource
where
    Self: Sized + MetaData,
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
    type Version: Display + IdentifiedPeptidoformVersion<Self::Format>;

    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(
        source: &Self::Source,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
    ) -> Result<(Self, &'static Self::Format), BoxedError<'static>>;

    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &Self::Format,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
    ) -> Result<Self, BoxedError<'static>>;

    /// Parse a source of multiple peptides using the given format or automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Result<Self::Source, BoxedError<'static>>>>(
        iter: I,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
        format: Option<Self::Format>,
    ) -> IdentifiedPeptidoformIter<'_, Self, I> {
        IdentifiedPeptidoformIter {
            iter: Box::new(iter),
            format,
            custom_database,
            keep_all_columns,
            peek: None,
        }
    }

    /// Parse a file with identified peptides.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
        version: Option<Self::Version>,
    ) -> Result<BoxedIdentifiedPeptideIter<'_, Self>, BoxedError<'static>>;

    /// Parse a reader with identified peptides.
    /// # Errors
    /// When the file is empty or no headers are present.
    fn parse_reader<'a>(
        reader: impl std::io::Read + 'a,
        custom_database: Option<&'a CustomDatabase>,
        keep_all_columns: bool,
        version: Option<Self::Version>,
    ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, BoxedError<'static>>;

    /// Allow post processing of the peptide
    /// # Errors
    /// On errors in the post processing, format specific
    #[expect(unused_variables)]
    fn post_process(
        source: &Self::Source,
        parsed: Self,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, BoxedError<'static>> {
        Ok(parsed)
    }
}

/// A general generic identified peptidoform iterator from any source format
pub type GeneralIdentifiedPeptidoforms<'lifetime> = Box<
    dyn Iterator<
            Item = Result<IdentifiedPeptidoform<Linked, MaybePeptidoform>, BoxedError<'static>>,
        > + 'lifetime,
>;

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = IdentifiedPeptidoformIter<
    'lifetime,
    T,
    Box<
        dyn Iterator<Item = Result<<T as IdentifiedPeptidoformSource>::Source, BoxedError<'static>>>
            + 'lifetime,
    >,
>;

/// An iterator returning parsed identified peptides
#[derive(Debug)]
pub struct IdentifiedPeptidoformIter<
    'lifetime,
    Source: IdentifiedPeptidoformSource,
    Iter: Iterator<Item = Result<Source::Source, BoxedError<'static>>>,
> {
    iter: Box<Iter>,
    format: Option<Source::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
    keep_all_columns: bool,
    peek: Option<Result<Source, BoxedError<'static>>>,
}

impl<
    R: IdentifiedPeptidoformSource + Clone,
    I: Iterator<Item = Result<R::Source, BoxedError<'static>>>,
> IdentifiedPeptidoformIter<'_, R, I>
where
    R::Format: 'static,
{
    /// Peek at the next item in the iterator
    pub fn peek(&mut self) -> Option<Result<R, BoxedError<'static>>> {
        if self.peek.is_some() {
            return self.peek.clone();
        }

        let peek = if let Some(format) = &self.format {
            self.iter.next().map(|source| {
                R::parse_specific(
                    &source?,
                    format,
                    self.custom_database,
                    self.keep_all_columns,
                )
                .map_err(BoxedError::to_owned)
            })
        } else {
            match self.iter.next().map(|source| {
                R::parse(&source?, self.custom_database, self.keep_all_columns)
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

impl<R: IdentifiedPeptidoformSource, I: Iterator<Item = Result<R::Source, BoxedError<'static>>>>
    Iterator for IdentifiedPeptidoformIter<'_, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, BoxedError<'static>>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek.is_some() {
            return self.peek.take();
        }

        if let Some(format) = &self.format {
            self.iter.next().map(|source| {
                R::parse_specific(
                    &source?,
                    format,
                    self.custom_database,
                    self.keep_all_columns,
                )
                .map_err(BoxedError::to_owned)
            })
        } else {
            match self.iter.next().map(|source| {
                R::parse(&source?, self.custom_database, self.keep_all_columns)
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

impl<'lifetime, Source, Iter> IdentifiedPeptidoformIter<'lifetime, Source, Iter>
where
    Source: IdentifiedPeptidoformSource
        + Into<IdentifiedPeptidoform<Source::Complexity, Source::PeptidoformAvailability>>
        + 'lifetime,
    Iter: Iterator<Item = Result<Source::Source, BoxedError<'static>>> + 'lifetime,
    Source::Format: 'static,
    MaybePeptidoform: From<Source::PeptidoformAvailability>,
{
    /// Make this into a generic boxed iterator for merging with other identified peptidoform formats
    pub fn into_box<Complexity: AtLeast<Source::Complexity>>(
        self,
    ) -> Box<
        dyn Iterator<
                Item = Result<
                    IdentifiedPeptidoform<Complexity, MaybePeptidoform>,
                    BoxedError<'static>,
                >,
            > + 'lifetime,
    > {
        Box::new(self.map(|p: Result<Source, BoxedError>| match p {
            Ok(p) => Ok(p.into().cast()),
            Err(e) => Err(e),
        }))
    }
}
