use std::{borrow::Cow, fmt::Display, path::PathBuf};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    formula::MultiChemical,
    identification::{
        deepnovofamily::DeepNovoFamilyData, fasta::FastaData, instanovo::InstaNovoData,
        novor::NovorData, opair::OpairData, peaks::PeaksData, plink::PLinkData,
        system::MassOverCharge, MSFraggerData, MZTabData, MaxQuantData, SageData,
    },
    ontologies::CustomDatabase,
    peptide::SemiAmbiguous,
    system::usize::Charge,
    system::Time,
    LinearPeptide, Peptidoform,
};

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct IdentifiedPeptide {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[allow(clippy::large_enum_variant)]
pub enum MetaData {
    /// DeepNovo/PointNovo/PGPointNovo metadata
    DeepNovoFamily(DeepNovoFamilyData),
    /// Fasta metadata
    Fasta(FastaData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
    /// InstaNovo metadata
    InstaNovo(InstaNovoData),
    /// MSFragger metadata
    MSFragger(MSFraggerData),
    /// mzTab metadata
    MZTab(MZTabData),
    /// Novor metadata
    Novor(NovorData),
    /// OPair metadata
    Opair(OpairData),
    /// Peaks metadata
    Peaks(PeaksData),
    /// pLink metadata
    PLink(PLinkData),
    /// Sage metadata
    Sage(SageData),
}

/// A peptide as stored in a identified peptide file, either a simple linear one or a cross-linked peptidoform
#[derive(Debug, Clone, Copy)]
pub enum ReturnedPeptide<'a> {
    /// A simple linear peptide
    Linear(&'a LinearPeptide<SemiAmbiguous>),
    /// A (potentially) cross-linked peptidoform
    Peptidoform(&'a Peptidoform),
}

impl<'a> MultiChemical for ReturnedPeptide<'a> {
    fn formulas_inner(
        &self,
        _sequence_index: super::SequencePosition,
        _peptide_index: usize,
    ) -> super::Multi<super::MolecularFormula> {
        match self {
            Self::Linear(p) => p.formulas(),
            Self::Peptidoform(p) => p.formulas(),
        }
    }
}

impl<'a> std::fmt::Display for ReturnedPeptide<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Linear(p) => write!(f, "{p}"),
            Self::Peptidoform(p) => write!(f, "{p}"),
        }
    }
}

#[allow(dead_code)]
impl<'a> ReturnedPeptide<'a> {
    /// Get the underlying peptide, or None if the underlying result was a peptidoform
    pub const fn peptide(self) -> Option<&'a LinearPeptide<SemiAmbiguous>> {
        match self {
            Self::Linear(p) => Some(p),
            Self::Peptidoform(_) => None,
        }
    }
    /// Get the underlying result as a peptidoform, if it was a peptide make a new peptidoform from it
    pub fn peptidoform(self) -> Cow<'a, Peptidoform> {
        match self {
            Self::Linear(p) => Cow::Owned(p.clone().into()),
            Self::Peptidoform(p) => Cow::Borrowed(p),
        }
    }
    /// Display this peptidoform.
    /// `specification_compliant` Displays this peptidoform either normalised to the internal representation or as fully spec compliant ProForma
    /// (no glycan structure or custom modifications).
    /// # Panics
    /// When some peptides do not have the same global isotope modifications.
    /// # Errors
    /// If the underlying formatter errors.
    pub fn display(
        &self,
        f: &mut impl std::fmt::Write,
        show_global_mods: bool,
        specification_compliant: bool,
    ) -> std::fmt::Result {
        match self {
            Self::Linear(p) => p.display(f, show_global_mods, specification_compliant),
            Self::Peptidoform(p) => p.display(f, show_global_mods, specification_compliant),
        }
    }
}

impl IdentifiedPeptide {
    /// Get the peptide, as pLink can have cross-linked peptides the return type is either a simple peptide or a cross-linked peptidoform
    pub fn peptide(&self) -> Option<ReturnedPeptide<'_>> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { peptide, .. })
            | MetaData::Novor(NovorData { peptide, .. })
            | MetaData::Opair(OpairData { peptide, .. })
            | MetaData::InstaNovo(InstaNovoData { peptide, .. })
            | MetaData::Sage(SageData { peptide, .. })
            | MetaData::MZTab(MZTabData { peptide, .. }) => Some(ReturnedPeptide::Linear(peptide)),
            MetaData::MSFragger(MSFraggerData { peptide, .. })
            | MetaData::MaxQuant(MaxQuantData { peptide, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(ReturnedPeptide::Linear)
            }
            MetaData::Fasta(f) => Some(ReturnedPeptide::Linear(f.peptide())),
            MetaData::PLink(PLinkData { peptidoform, .. }) => {
                Some(ReturnedPeptide::Peptidoform(peptidoform))
            }
        }
    }

    /// Get the name of the format
    pub const fn format_name(&self) -> &'static str {
        match &self.metadata {
            MetaData::DeepNovoFamily(_) => "DeepNovo Family",
            MetaData::Fasta(_) => "Fasta",
            MetaData::InstaNovo(_) => "InstaNovo",
            MetaData::MaxQuant(_) => "MaxQuant",
            MetaData::MSFragger(_) => "MSFragger",
            MetaData::MZTab(_) => "mzTab",
            MetaData::Novor(_) => "Novor",
            MetaData::Opair(_) => "OPair",
            MetaData::Peaks(_) => "PEAKS",
            MetaData::PLink(_) => "pLink",
            MetaData::Sage(_) => "Sage",
        }
    }

    /// Get the format version detected
    pub fn format_version(&self) -> String {
        match &self.metadata {
            MetaData::Fasta(_) => "Fasta".to_string(),
            MetaData::MZTab(_) => "mzTab 1.0".to_string(),
            MetaData::MaxQuant(MaxQuantData { version, .. }) => version.to_string(),
            MetaData::DeepNovoFamily(DeepNovoFamilyData { version, .. }) => version.to_string(),
            MetaData::MSFragger(MSFraggerData { version, .. }) => version.to_string(),
            MetaData::Novor(NovorData { version, .. }) => version.to_string(),
            MetaData::Opair(OpairData { version, .. }) => version.to_string(),
            MetaData::Peaks(PeaksData { version, .. }) => version.to_string(),
            MetaData::PLink(PLinkData { version, .. }) => version.to_string(),
            MetaData::InstaNovo(InstaNovoData { version, .. }) => version.to_string(),
            MetaData::Sage(SageData { version, .. }) => version.to_string(),
        }
    }

    /// Get the identifier
    pub fn id(&self) -> String {
        match &self.metadata {
            MetaData::Peaks(PeaksData { scan, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => scan.iter().join(";"),
            MetaData::Novor(NovorData { id, scan, .. }) => id.unwrap_or(*scan).to_string(),
            MetaData::Opair(OpairData { scan, .. })
            | MetaData::InstaNovo(InstaNovoData { scan, .. }) => scan.to_string(),
            MetaData::Sage(SageData { id, .. }) | MetaData::MZTab(MZTabData { id, .. }) => {
                id.to_string()
            }
            MetaData::Fasta(f) => f.identifier().accession().to_string(),
            MetaData::MSFragger(MSFraggerData { scan, .. }) => scan.to_string(),
            MetaData::PLink(PLinkData { order, .. }) => order.to_string(),
            MetaData::MaxQuant(MaxQuantData { id, scan, .. }) => {
                id.map_or_else(|| scan.iter().join(";"), |id| id.to_string())
            }
        }
    }

    /// Get the local confidence, it is the same length as the peptide with a local score in 0..=1
    pub fn local_confidence(&self) -> Option<&[f64]> {
        match &self.metadata {
            MetaData::Peaks(PeaksData {
                local_confidence, ..
            })
            | MetaData::InstaNovo(InstaNovoData {
                local_confidence, ..
            }) => Some(local_confidence),

            MetaData::DeepNovoFamily(DeepNovoFamilyData {
                local_confidence, ..
            })
            | MetaData::Novor(NovorData {
                local_confidence, ..
            })
            | MetaData::MZTab(MZTabData {
                local_confidence, ..
            }) => local_confidence.as_deref(),
            _ => None,
        }
    }

    /// The charge of the precursor, if known
    pub fn charge(&self) -> Option<Charge> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { z, .. })
            | MetaData::Novor(NovorData { z, .. })
            | MetaData::Opair(OpairData { z, .. })
            | MetaData::Sage(SageData { z, .. })
            | MetaData::MSFragger(MSFraggerData { z, .. })
            | MetaData::MaxQuant(MaxQuantData { z, .. })
            | MetaData::PLink(PLinkData { z, .. })
            | MetaData::InstaNovo(InstaNovoData { z, .. })
            | MetaData::MZTab(MZTabData { z, .. }) => Some(*z),
            MetaData::DeepNovoFamily(DeepNovoFamilyData { z, .. }) => *z,
            MetaData::Fasta(_) => None,
        }
    }

    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mode, .. }) => Some(mode),
            MetaData::MaxQuant(MaxQuantData { fragmentation, .. }) => fragmentation.as_deref(),
            _ => None,
        }
    }

    /// The retention time, if known
    pub fn retention_time(&self) -> Option<Time> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { rt, .. })
            | MetaData::Opair(OpairData { rt, .. })
            | MetaData::Sage(SageData { rt, .. })
            | MetaData::MSFragger(MSFraggerData { rt, .. }) => Some(*rt),
            MetaData::MaxQuant(MaxQuantData { rt, .. })
            | MetaData::Novor(NovorData { rt, .. })
            | MetaData::MZTab(MZTabData { rt, .. }) => *rt,
            MetaData::DeepNovoFamily(_)
            | MetaData::InstaNovo(_)
            | MetaData::Fasta(_)
            | MetaData::PLink(_) => None,
        }
    }

    /// The scans per rawfile that are at the basis for this identified peptide, if the rawfile is unknown there will be one
    pub fn scans(&self) -> SpectrumIds {
        match &self.metadata {
            MetaData::Peaks(PeaksData { raw_file, scan, .. }) => raw_file.clone().map_or_else(
                || {
                    SpectrumIds::FileNotKnown(
                        scan.iter()
                            .flat_map(|s| s.scans.clone())
                            .map(SpectrumId::Index)
                            .collect(),
                    )
                },
                |raw_file| {
                    SpectrumIds::FileKnown(vec![(
                        raw_file,
                        scan.iter()
                            .flat_map(|s| s.scans.clone())
                            .map(SpectrumId::Index)
                            .collect(),
                    )])
                },
            ),
            MetaData::Novor(NovorData { scan, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)])
            }
            MetaData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => SpectrumIds::FileNotKnown(
                scan.iter()
                    .flat_map(|s| s.scans.clone())
                    .map(SpectrumId::Index)
                    .collect(),
            ),

            MetaData::Opair(OpairData { raw_file, scan, .. })
            | MetaData::InstaNovo(InstaNovoData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(raw_file.clone(), vec![SpectrumId::Index(*scan)])])
            }
            MetaData::MaxQuant(MaxQuantData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(
                    raw_file.clone(),
                    scan.iter().copied().map(SpectrumId::Index).collect(),
                )])
            }
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => spectra_ref.clone(),
            MetaData::MSFragger(MSFraggerData { raw_file, scan, .. }) => {
                raw_file.clone().map_or_else(
                    || SpectrumIds::FileNotKnown(vec![scan.clone()]),
                    |raw_file| SpectrumIds::FileKnown(vec![(raw_file, vec![scan.clone()])]),
                )
            }
            MetaData::PLink(PLinkData {
                raw_file,
                scan,
                title,
                ..
            }) => scan.map_or_else(
                || SpectrumIds::FileNotKnown(vec![SpectrumId::Native(title.clone())]),
                |scan| {
                    raw_file.clone().map_or_else(
                        || SpectrumIds::FileNotKnown(vec![SpectrumId::Index(scan)]),
                        |raw_file| {
                            SpectrumIds::FileKnown(vec![(raw_file, vec![SpectrumId::Index(scan)])])
                        },
                    )
                },
            ),
            MetaData::Sage(SageData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(raw_file.clone(), vec![scan.clone()])])
            }
            MetaData::Fasta(_) => SpectrumIds::None,
        }
    }

    /// Get the mz as experimentally determined
    pub fn experimental_mz(&self) -> Option<MassOverCharge> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mz, .. })
            | MetaData::Novor(NovorData { mz, .. })
            | MetaData::Opair(OpairData { mz, .. })
            | MetaData::InstaNovo(InstaNovoData { mz, .. })
            | MetaData::MSFragger(MSFraggerData { mz, .. }) => Some(*mz),
            MetaData::MZTab(MZTabData { mz, .. }) | MetaData::MaxQuant(MaxQuantData { mz, .. }) => {
                *mz
            }
            MetaData::Sage(SageData { mass, z, .. })
            | MetaData::PLink(PLinkData { mass, z, .. }) => {
                Some(MassOverCharge::new::<crate::system::mz>(
                    mass.value / (z.value as f64),
                ))
            }
            MetaData::DeepNovoFamily(_) | MetaData::Fasta(_) => None,
        }
    }

    /// Get the mass as experimentally determined
    pub fn experimental_mass(&self) -> Option<crate::system::Mass> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mass, .. })
            | MetaData::Novor(NovorData { mass, .. })
            | MetaData::Opair(OpairData { mass, .. })
            | MetaData::MSFragger(MSFraggerData { mass, .. })
            | MetaData::PLink(PLinkData { mass, .. })
            | MetaData::Sage(SageData { mass, .. }) => Some(*mass),
            MetaData::MaxQuant(MaxQuantData { mass, .. }) => *mass,
            MetaData::MZTab(MZTabData { mz, z, .. }) => mz.map(|mz| mz * z.to_float()),
            MetaData::InstaNovo(InstaNovoData { mz, z, .. }) => Some(*mz * z.to_float()),
            MetaData::DeepNovoFamily(DeepNovoFamilyData { mz, z, .. }) => {
                mz.and_then(|mz| z.map(|z| (mz, z)).map(|(mz, z)| mz * z.to_float()))
            }
            MetaData::Fasta(_) => None,
        }
    }

    /// Get the absolute ppm error between the experimental and theoretical precursor mass
    pub fn ppm_error(&self) -> Option<crate::system::Ratio> {
        let exp_mass = self.experimental_mass()?;
        let theo_mass = self
            .peptide()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;

        Some(theo_mass.ppm(exp_mass))
    }

    /// Get the absolute mass error between the experimental and theoretical precursor mass
    pub fn mass_error(&self) -> Option<crate::system::Mass> {
        let exp_mass = self.experimental_mass()?;
        let theo_mass = self
            .peptide()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;

        Some((exp_mass - theo_mass).abs())
    }
}

/// Multiple spectrum identifiers
#[derive(Clone, Default, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum SpectrumIds {
    /// When no spectra references are known at all
    #[default]
    None,
    /// When the source file is now known
    FileNotKnown(Vec<SpectrumId>),
    /// When the source file is known, grouped per file
    FileKnown(Vec<(PathBuf, Vec<SpectrumId>)>),
}

/// A spectrum identifier
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
pub enum SpectrumId {
    /// A native id, the format differs between vendors
    Native(String),
    /// A spectrum index
    Index(usize),
}

impl Default for SpectrumId {
    fn default() -> Self {
        Self::Index(0)
    }
}

impl std::fmt::Display for SpectrumId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Index(i) => write!(f, "{i}"),
            Self::Native(n) => write!(f, "{n}"),
        }
    }
}

impl SpectrumId {
    /// Get the index if this is an index
    pub const fn index(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::Native(_) => None,
        }
    }

    /// Get the native ID if this is a native ID
    pub fn native(&self) -> Option<&str> {
        match self {
            Self::Index(_) => None,
            Self::Native(n) => Some(n),
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

    /// The version type
    type Version;

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
    fn parse_many<I: Iterator<Item = Result<Self::Source, CustomError>>>(
        iter: I,
        custom_database: Option<&CustomDatabase>,
    ) -> IdentifiedPeptideIter<Self, I> {
        IdentifiedPeptideIter {
            iter: Box::new(iter),
            format: None,
            custom_database,
            peek: None,
        }
    }

    /// Parse a file with identified peptides.
    /// # Errors
    /// Returns Err when the file could not be opened
    fn parse_file(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<BoxedIdentifiedPeptideIter<Self>, CustomError>;

    /// Parse a reader with identified peptides.
    /// # Errors
    /// When the file is empty or no headers are present.
    fn parse_reader<'a>(
        reader: impl std::io::Read + 'a,
        custom_database: Option<&'a CustomDatabase>,
    ) -> Result<BoxedIdentifiedPeptideIter<'a, Self>, CustomError>;

    /// Allow post processing of the peptide
    /// # Errors
    /// On errors in the post processing, format specific
    #[allow(unused_variables)]
    fn post_process(
        source: &Self::Source,
        parsed: Self,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        Ok(parsed)
    }
}

/// Convenience type to not have to type out long iterator types
pub type BoxedIdentifiedPeptideIter<'lifetime, T> = IdentifiedPeptideIter<
    'lifetime,
    T,
    Box<
        dyn Iterator<Item = Result<<T as IdentifiedPeptideSource>::Source, CustomError>>
            + 'lifetime,
    >,
>;

/// An iterator returning parsed identified peptides
pub struct IdentifiedPeptideIter<
    'lifetime,
    R: IdentifiedPeptideSource,
    I: Iterator<Item = Result<R::Source, CustomError>>,
> {
    iter: Box<I>,
    format: Option<R::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
    peek: Option<Result<R, CustomError>>,
}

impl<
        'lifetime,
        R: IdentifiedPeptideSource + Clone,
        I: Iterator<Item = Result<R::Source, CustomError>>,
    > IdentifiedPeptideIter<'lifetime, R, I>
where
    R::Format: 'static,
{
    /// Peek at the next item in the iterator
    pub fn peek(&mut self) -> Option<Result<R, CustomError>> {
        if self.peek.is_some() {
            return self.peek.clone();
        }

        let peek = if let Some(format) = &self.format {
            self.iter
                .next()
                .map(|source| R::parse_specific(&source?, format, self.custom_database))
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source?, self.custom_database))
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

impl<'lifetime, R: IdentifiedPeptideSource, I: Iterator<Item = Result<R::Source, CustomError>>>
    Iterator for IdentifiedPeptideIter<'lifetime, R, I>
where
    R::Format: 'static,
{
    type Item = Result<R, CustomError>;
    fn next(&mut self) -> Option<Self::Item> {
        if self.peek.is_some() {
            return self.peek.take();
        }

        if let Some(format) = &self.format {
            self.iter
                .next()
                .map(|source| R::parse_specific(&source?, format, self.custom_database))
        } else {
            match self
                .iter
                .next()
                .map(|source| R::parse(&source?, self.custom_database))
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

impl<'lifetime, R, I> IdentifiedPeptideIter<'lifetime, R, I>
where
    R: IdentifiedPeptideSource + Into<IdentifiedPeptide> + 'lifetime,
    I: Iterator<Item = Result<R::Source, CustomError>> + 'lifetime,
    R::Format: 'static,
{
    pub(super) fn into_box(
        self,
    ) -> Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'lifetime> {
        Box::new(self.map(|p: Result<R, CustomError>| match p {
            Ok(p) => Ok(p.into()),
            Err(e) => Err(e),
        }))
    }
}

/// Test a dataset for common errors in identified peptide parsing
#[cfg(test)]
pub(crate) fn test_format<T: IdentifiedPeptideSource + Into<IdentifiedPeptide>>(
    reader: impl std::io::Read,
    custom_database: Option<&CustomDatabase>,
    allow_mass_mods: bool,
    expect_lc: bool,
    format: Option<T::Version>,
) -> Result<usize, String>
where
    T::Format: 'static,
    T::Version: std::fmt::Display,
{
    let mut number = 0;
    for peptide in T::parse_reader(reader, custom_database).map_err(|e| e.to_string())? {
        let peptide: IdentifiedPeptide = peptide.map_err(|e| e.to_string())?.into();
        number += 1;

        if peptide
            .peptide()
            .and_then(|p| p.peptide())
            .map(LinearPeptide::len)
            != peptide.local_confidence().map(<[f64]>::len)
        {
            if expect_lc && peptide.local_confidence().is_none() {
                return Err(format!(
                    "No local confidence was provided for peptide {}",
                    peptide.id()
                ));
            } else if peptide.local_confidence().is_some() {
                return Err(format!("The local confidence ({}) does not have the same number of elements as the peptide ({}) for peptide {}", peptide.local_confidence().map_or(0, <[f64]>::len), peptide.peptide().and_then(|p| p.peptide()).map_or(0, LinearPeptide::len), peptide.id()));
            }
        }
        if peptide.score.is_some_and(|s| !(-1.0..=1.0).contains(&s)) {
            return Err(format!(
                "The score {} for peptide {} is outside of range",
                peptide.score.unwrap(),
                peptide.id()
            ));
        }
        if peptide
            .local_confidence()
            .is_some_and(|s| s.iter().any(|s| !(-1.0..=1.0).contains(s)))
        {
            return Err(format!(
                "The local score {} for peptide {} is outside of range",
                peptide.local_confidence().unwrap().iter().join(","),
                peptide.id()
            ));
        }
        if !allow_mass_mods
            && peptide.peptide().is_some_and(|p| {
                p.peptidoform().peptides().iter().any(|p| {
                    p.sequence().iter().any(|s| {
                        s.modifications.iter().any(|m| {
                            matches!(
                                m,
                                crate::Modification::Simple(
                                    crate::modification::SimpleModification::Mass(_)
                                )
                            )
                        })
                    })
                })
            })
        {
            return Err(format!(
                "Peptide {} contains mass modifications, sequence {}",
                peptide.id(),
                peptide.peptide().unwrap(),
            ));
        }
        if format
            .as_ref()
            .is_some_and(|f| f.to_string() != peptide.format_version())
        {
            return Err(format!(
                "Peptide {} was detected as the wrong version ({} instead of {})",
                peptide.id(),
                peptide.format_version(),
                format.unwrap(),
            ));
        }
    }
    Ok(number)
}

/// The scans identifier for a peaks identification
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub struct PeaksFamilyId {
    /// The file, if defined
    pub file: Option<usize>,
    /// The scan(s)
    pub scans: Vec<usize>,
}

impl Display for PeaksFamilyId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self.file.map_or(String::new(), |f| format!("F{f}:")),
            self.scans.iter().join(",")
        )
    }
}

impl std::str::FromStr for PeaksFamilyId {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((start, end)) = s.split_once(':') {
            if start.is_empty() || end.is_empty() {
                Err(())
            } else {
                Ok(Self {
                    file: Some(start[1..].parse().map_err(|_| ())?),
                    scans: end
                        .split(' ')
                        .map(str::parse)
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|_| ())?,
                })
            }
        } else {
            Ok(Self {
                file: None,
                scans: vec![s.parse().map_err(|_| ())?],
            })
        }
    }
}
