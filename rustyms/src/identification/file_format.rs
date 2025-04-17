use std::{fmt::Display, path::Path};

use serde::{Deserialize, Serialize};

use crate::{error::CustomError, identification::*, ontologies::CustomDatabase};

/// A file format that is fully known
#[derive(Clone, PartialEq, Eq, Debug, Copy, Serialize, Deserialize)]
#[allow(clippy::upper_case_acronyms, missing_docs)]
pub enum KnownFileFormat {
    BasicCSV(BasicCSVVersion),
    DeepNovoFamily(DeepNovoFamilyVersion),
    Fasta,
    InstaNovo(InstaNovoVersion),
    MSFragger(MSFraggerVersion),
    MaxQuant(MaxQuantVersion),
    MZTab,
    NovoB(NovoBVersion),
    Novor(NovorVersion),
    Opair(OpairVersion),
    PepNet(PepNetVersion),
    PLGS(PLGSVersion),
    PLink(PLinkVersion),
    Peaks(PeaksVersion),
    PowerNovo(PowerNovoVersion),
    Sage(SageVersion),
    SpectrumSequenceList(SpectrumSequenceListVersion),
}

impl KnownFileFormat {
    /// Get the name of the format
    pub const fn name(self) -> &'static str {
        match self {
            Self::BasicCSV(_) => "CSV",
            Self::SpectrumSequenceList(_) => "SpectrumSequenceList",
            Self::DeepNovoFamily(_) => "DeepNovo Family",
            Self::Fasta => "Fasta",
            Self::InstaNovo(_) => "InstaNovo",
            Self::MaxQuant(_) => "MaxQuant",
            Self::MSFragger(_) => "MSFragger",
            Self::MZTab => "mzTab",
            Self::NovoB(_) => "NovoB",
            Self::Novor(_) => "Novor",
            Self::Opair(_) => "OPair",
            Self::Peaks(_) => "PEAKS",
            Self::PepNet(_) => "PepNet",
            Self::PLGS(_) => "ProteinLynx Global Server",
            Self::PLink(_) => "pLink",
            Self::PowerNovo(_) => "PowerNovo",
            Self::Sage(_) => "Sage",
        }
    }

    /// Get the format version
    pub fn version(self) -> Option<String> {
        match self {
            Self::SpectrumSequenceList(version) => Some(version.to_string()),
            Self::DeepNovoFamily(version) => Some(version.to_string()),
            Self::BasicCSV(version) => Some(version.to_string()),
            Self::Fasta => None,
            Self::InstaNovo(version) => Some(version.to_string()),
            Self::MaxQuant(version) => Some(version.to_string()),
            Self::MSFragger(version) => Some(version.to_string()),
            Self::MZTab => Some("1.0".to_string()),
            Self::NovoB(version) => Some(version.to_string()),
            Self::Novor(version) => Some(version.to_string()),
            Self::Opair(version) => Some(version.to_string()),
            Self::Peaks(version) => Some(version.to_string()),
            Self::PepNet(version) => Some(version.to_string()),
            Self::PLGS(version) => Some(version.to_string()),
            Self::PLink(version) => Some(version.to_string()),
            Self::PowerNovo(version) => Some(version.to_string()),
            Self::Sage(version) => Some(version.to_string()),
        }
    }
}

impl Display for KnownFileFormat {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}{}",
            self.name(),
            self.version().map_or(String::new(), |v| format!(" - {v}"))
        )
    }
}

/// A file format that might not be (fully) known
#[derive(Clone, PartialEq, Eq, Debug, Copy, Serialize, Deserialize)]
#[allow(clippy::upper_case_acronyms, missing_docs)]
pub enum FileFormat {
    Undefined,
    BasicCSV(Option<BasicCSVVersion>),
    DeepNovoFamily(Option<DeepNovoFamilyVersion>),
    Fasta,
    InstaNovo(Option<InstaNovoVersion>),
    MSFragger(Option<MSFraggerVersion>),
    MaxQuant(Option<MaxQuantVersion>),
    MZTab,
    NovoB(Option<NovoBVersion>),
    Novor(Option<NovorVersion>),
    Opair(Option<OpairVersion>),
    PepNet(Option<PepNetVersion>),
    PLGS(Option<PLGSVersion>),
    PLink(Option<PLinkVersion>),
    Peaks(Option<PeaksVersion>),
    PowerNovo(Option<PowerNovoVersion>),
    Sage(Option<SageVersion>),
    SSL(Option<SpectrumSequenceListVersion>),
}

impl FileFormat {
    /// Open a file with this file format
    /// # Errors
    /// If the file is not valid according to the file format chosen
    pub fn open<'a>(
        self,
        path: &Path,
        custom_database: Option<&'a CustomDatabase>,
    ) -> Result<
        Box<dyn Iterator<Item = Result<IdentifiedPeptidoform, CustomError>> + 'a>,
        CustomError,
    > {
        match self {
            Self::Undefined => open_identified_peptides_file(path, custom_database, false),
            Self::BasicCSV(version) => {
                BasicCSVData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::DeepNovoFamily(version) => {
                DeepNovoFamilyData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Fasta => FastaData::parse_file(path).map(|sequences| {
                let b: Box<dyn Iterator<Item = Result<IdentifiedPeptidoform, CustomError>>> =
                    Box::new(sequences.into_iter().map(|p| Ok(p.into())));
                b
            }),
            Self::InstaNovo(version) => {
                InstaNovoData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::MSFragger(version) => {
                MSFraggerData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::MaxQuant(version) => {
                MaxQuantData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::MZTab => MZTabData::parse_file(path, custom_database).map(|sequences| {
                let b: Box<dyn Iterator<Item = Result<IdentifiedPeptidoform, CustomError>>> =
                    Box::new(sequences.into_iter().map(|p| p.map(Into::into)));
                b
            }),
            Self::NovoB(version) => NovoBData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::Novor(version) => NovorData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::Opair(version) => OpairData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::PLGS(version) => PLGSData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::PLink(version) => PLinkData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::Peaks(version) => PeaksData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::PepNet(version) => PepNetData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::PowerNovo(version) => {
                PowerNovoData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Sage(version) => SageData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::SSL(version) => {
                SpectrumSequenceListData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
        }
    }
}
