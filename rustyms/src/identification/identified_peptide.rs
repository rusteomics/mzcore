use serde::{Deserialize, Serialize};

use super::fasta::FastaData;
use super::novor::NovorData;
use super::opair::OpairData;
use super::peaks::PeaksData;
use super::{MSFraggerData, MaxQuantData, SageData};
use crate::error::CustomError;
use crate::ontologies::CustomDatabase;
use crate::peptide::VerySimple;
use crate::system::usize::Charge;
use crate::LinearPeptide;

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub struct IdentifiedPeptide {
    /// The peptide sequence
    pub peptide: LinearPeptide<VerySimple>,
    /// The local confidence of this peptide (same length as the peptide)
    pub local_confidence: Option<Vec<f64>>,
    /// The score -1.0..=1.0 if available in the original format
    pub score: Option<f64>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, PartialEq, Debug, Default, Serialize, Deserialize)]
pub enum MetaData {
    /// No metadata
    #[default]
    None,
    /// Peaks metadata
    Peaks(PeaksData),
    /// Novor metadata
    Novor(NovorData),
    /// OPair metadata
    Opair(OpairData),
    /// Fasta metadata
    Fasta(FastaData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
    /// Sage metadata
    Sage(SageData),
    /// Sage metadata
    MSFragger(MSFraggerData),
}

impl MetaData {
    /// The charge of the precursor, if known
    pub const fn charge(&self) -> Option<Charge> {
        match self {
            Self::Peaks(PeaksData { z, .. })
            | Self::Novor(NovorData { z, .. })
            | Self::Opair(OpairData { z, .. })
            | Self::Sage(SageData { z, .. })
            | Self::MSFragger(MSFraggerData { z, .. })
            | Self::MaxQuant(MaxQuantData { z, .. }) => Some(*z),
            Self::Fasta(_) | Self::None => None,
        }
    }

    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match self {
            Self::Peaks(PeaksData { mode, .. }) => Some(mode),
            _ => None,
        }
    }

    /// The scan number of the spectrum for this identified peptide, if known.
    // TODO: Allow multiple scan numbers to be returned, but think about merging spectra then as well
    pub fn scan_number(&self) -> Option<usize> {
        match self {
            Self::Peaks(PeaksData { scan, .. }) => {
                scan.first().and_then(|i| i.scans.first().copied())
            }
            Self::Novor(NovorData { scan, .. }) | Self::Opair(OpairData { scan, .. }) => {
                Some(*scan)
            }
            Self::MaxQuant(MaxQuantData { scan_number, .. }) => scan_number.first().copied(),
            Self::Sage(SageData { scan_nr, .. }) => Some(scan_nr.2),
            Self::MSFragger(MSFraggerData { spectrum, .. }) => Some(spectrum.scan.0),
            Self::Fasta(_) | Self::None => None,
        }
    }

    /// Get the file name for the rawfile
    pub fn raw_file(&self) -> Option<&str> {
        match self {
            Self::Peaks(PeaksData { source_file, .. }) => source_file.as_ref().map(String::as_str),
            Self::Opair(OpairData { file_name, .. }) => Some(file_name),
            Self::MaxQuant(MaxQuantData { raw_file, .. }) => Some(raw_file),
            Self::Sage(SageData { filename, .. }) => Some(filename),
            Self::MSFragger(MSFraggerData { spectrum, .. }) => Some(&spectrum.file),
            Self::Novor(_) | Self::Fasta(_) | Self::None => None,
        }
    }
}

/// The required methods for any source of identified peptides
pub trait IdentifiedPeptideSource
where
    Self: std::marker::Sized,
{
    /// The source data where the peptides are parsed form
    type Source;
    /// The format type
    type Format: Clone;
    /// Parse a single identified peptide from its source and return the detected format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse(
        source: &Self::Source,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<(Self, &'static Self::Format), CustomError>;
    /// Parse a single identified peptide with the given format
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_specific(
        source: &Self::Source,
        format: &Self::Format,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError>;
    /// Parse a source of multiple peptides automatically determining the format to use by the first item
    /// # Errors
    /// When the source is not a valid peptide
    fn parse_many<I: Iterator<Item = Self::Source>>(
        iter: I,
        custom_database: Option<&CustomDatabase>,
    ) -> IdentifiedPeptideIter<Self, I> {
        IdentifiedPeptideIter {
            iter: Box::new(iter),
            format: None,
            custom_database,
        }
    }
    /// Parse a file with identified peptides.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, CustomError>;
}

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = IdentifiedPeptideIter<
    'lifetime,
    T,
    Box<dyn Iterator<Item = <T as IdentifiedPeptideSource>::Source>>,
>;

/// An iterator returning parsed identified peptides
pub struct IdentifiedPeptideIter<
    'lifetime,
    R: IdentifiedPeptideSource,
    I: Iterator<Item = R::Source>,
> {
    iter: Box<I>,
    format: Option<R::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
}

impl<'lifetime, R: IdentifiedPeptideSource, I: Iterator<Item = R::Source>> Iterator
    for IdentifiedPeptideIter<'lifetime, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, CustomError>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Some(format) = &self.format {
            self.iter
                .next()
                .map(|source| R::parse_specific(&source, format, self.custom_database))
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source, self.custom_database))
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
