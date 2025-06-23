use std::fmt::Display;

use crate::{
    error::CustomError, identification::IdentifiedPeptidoform, ontology::CustomDatabase,
    sequence::AtLeast,
};

/// A version for an identified peptide version
pub trait IdentifiedPeptidoformVersion<Format>: Copy {
    /// The format for this version
    fn format(self) -> Format;
    /// The name for this version
    fn name(self) -> &'static str;
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptidoformSource
where
    Self: Sized,
{
    /// The source data where the peptides are parsed form
    type Source;

    /// The format type
    type Format: Clone;

    /// The complexity marker type
    type Complexity;

    /// The version type
    type Version: Display + IdentifiedPeptidoformVersion<Self::Format>;

    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(
        source: &Self::Source,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
    ) -> Result<(Self, &'static Self::Format), CustomError>;

    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &Self::Format,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
    ) -> Result<Self, CustomError>;

    /// Parse a source of multiple peptides using the given format or automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Result<Self::Source, CustomError>>>(
        iter: I,
        custom_database: Option<&CustomDatabase>,
        keep_all_columns: bool,
        format: Option<Self::Format>,
    ) -> IdentifiedPeptidoformIter<Self, I> {
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
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, CustomError>;

    /// Parse a reader with identified peptides.
    /// # Errors
    /// When the file is empty or no headers are present.
    fn parse_reader<'a>(
        reader: impl std::io::Read + 'a,
        custom_database: Option<&'a CustomDatabase>,
        keep_all_columns: bool,
        version: Option<Self::Version>,
    ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, CustomError>;

    /// Allow post processing of the peptide
    /// # Errors
    /// On errors in the post processing, format specific
    #[expect(unused_variables)]
    fn post_process(
        source: &Self::Source,
        parsed: Self,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        Ok(parsed)
    }
}

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = IdentifiedPeptidoformIter<
    'lifetime,
    T,
    Box<
        dyn Iterator<Item = Result<<T as IdentifiedPeptidoformSource>::Source, CustomError>>
            + 'lifetime,
    >,
>;

/// An iterator returning parsed identified peptides
#[derive(Debug)]
pub struct IdentifiedPeptidoformIter<
    'lifetime,
    Source: IdentifiedPeptidoformSource,
    Iter: Iterator<Item = Result<Source::Source, CustomError>>,
> {
    iter: Box<Iter>,
    format: Option<Source::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
    keep_all_columns: bool,
    peek: Option<Result<Source, CustomError>>,
}

impl<R: IdentifiedPeptidoformSource + Clone, I: Iterator<Item = Result<R::Source, CustomError>>>
    IdentifiedPeptidoformIter<'_, R, I>
where
    R::Format: 'static,
{
    /// Peek at the next item in the iterator
    pub fn peek(&mut self) -> Option<Result<R, CustomError>> {
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
            })
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source?, self.custom_database, self.keep_all_columns))
            {
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

impl<R: IdentifiedPeptidoformSource, I: Iterator<Item = Result<R::Source, CustomError>>> Iterator
    for IdentifiedPeptidoformIter<'_, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, CustomError>;
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
            })
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source?, self.custom_database, self.keep_all_columns))
            {
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
    Source:
        IdentifiedPeptidoformSource + Into<IdentifiedPeptidoform<Source::Complexity>> + 'lifetime,
    Iter: Iterator<Item = Result<Source::Source, CustomError>> + 'lifetime,
    Source::Format: 'static,
{
    /// Make this into a generic boxed iterator for merging with other identified peptidoform formats
    pub fn into_box<Complexity: AtLeast<Source::Complexity>>(
        self,
    ) -> Box<dyn Iterator<Item = Result<IdentifiedPeptidoform<Complexity>, CustomError>> + 'lifetime>
    {
        Box::new(self.map(|p: Result<Source, CustomError>| match p {
            Ok(p) => Ok(p.into().cast()),
            Err(e) => Err(e),
        }))
    }
}
