use std::{fmt::Display, path::Path};

use custom_error::BoxedError;
use serde::{Deserialize, Serialize};

use crate::{
    identification::*,
    ontology::CustomDatabase,
    sequence::{Linked, SemiAmbiguous},
};

/// A file format that is fully known
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[allow(clippy::upper_case_acronyms, missing_docs)]
pub enum KnownFileFormat {
    BasicCSV(BasicCSVVersion),
    DeepNovoFamily(DeepNovoFamilyVersion),
    Fasta,
    InstaNovo(InstaNovoVersion),
    MaxQuant(MaxQuantVersion),
    MZTab,
    NovoB(NovoBVersion),
    Novor(NovorVersion),
    Opair(OpairVersion),
    Peaks(PeaksVersion),
    PepNet(PepNetVersion),
    PLGS(PLGSVersion),
    PLink(PLinkVersion),
    PowerNovo(PowerNovoVersion),
    Proteoscape(ProteoscapeVersion),
    Sage(SageVersion),
    MSFragger(MSFraggerVersion),
    SpectrumSequenceList(SpectrumSequenceListVersion),
}

impl KnownFileFormat {
    /// Get the name of the format
    pub const fn name(self) -> &'static str {
        match self {
            Self::BasicCSV(_) => "CSV",
            Self::DeepNovoFamily(_) => "DeepNovo Family",
            Self::Fasta => "Fasta",
            Self::InstaNovo(_) => "InstaNovo",
            Self::MaxQuant(_) => "MaxQuant",
            Self::MZTab => "mzTab",
            Self::NovoB(_) => "NovoB",
            Self::Novor(_) => "Novor",
            Self::Opair(_) => "OPair",
            Self::Peaks(_) => "PEAKS",
            Self::PepNet(_) => "PepNet",
            Self::PLGS(_) => "ProteinLynx Global Server",
            Self::PLink(_) => "pLink",
            Self::PowerNovo(_) => "PowerNovo",
            Self::Proteoscape(_) => "Proteoscape",
            Self::Sage(_) => "Sage",
            Self::MSFragger(_) => "MSFragger",
            Self::SpectrumSequenceList(_) => "SpectrumSequenceList",
        }
    }

    /// Get the format version
    pub fn version(self) -> Option<String> {
        match self {
            Self::BasicCSV(version) => Some(version.to_string()),
            Self::DeepNovoFamily(version) => Some(version.to_string()),
            Self::Fasta => None,
            Self::InstaNovo(version) => Some(version.to_string()),
            Self::MaxQuant(version) => Some(version.to_string()),
            Self::MZTab => Some("1.0".to_string()),
            Self::NovoB(version) => Some(version.to_string()),
            Self::Novor(version) => Some(version.to_string()),
            Self::Opair(version) => Some(version.to_string()),
            Self::Peaks(version) => Some(version.to_string()),
            Self::PepNet(version) => Some(version.to_string()),
            Self::PLGS(version) => Some(version.to_string()),
            Self::PLink(version) => Some(version.to_string()),
            Self::PowerNovo(version) => Some(version.to_string()),
            Self::Proteoscape(version) => Some(version.to_string()),
            Self::Sage(version) => Some(version.to_string()),
            Self::MSFragger(version) => Some(version.to_string()),
            Self::SpectrumSequenceList(version) => Some(version.to_string()),
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

impl From<KnownFileFormat> for FileFormat {
    fn from(value: KnownFileFormat) -> Self {
        match value {
            KnownFileFormat::BasicCSV(version) => Self::BasicCSV(Some(version)),
            KnownFileFormat::DeepNovoFamily(version) => Self::DeepNovoFamily(Some(version)),
            KnownFileFormat::Fasta => Self::Fasta,
            KnownFileFormat::InstaNovo(version) => Self::InstaNovo(Some(version)),
            KnownFileFormat::MaxQuant(version) => Self::MaxQuant(Some(version)),
            KnownFileFormat::MZTab => Self::MZTab,
            KnownFileFormat::NovoB(version) => Self::NovoB(Some(version)),
            KnownFileFormat::Novor(version) => Self::Novor(Some(version)),
            KnownFileFormat::Opair(version) => Self::Opair(Some(version)),
            KnownFileFormat::Peaks(version) => Self::Peaks(Some(version)),
            KnownFileFormat::PepNet(version) => Self::PepNet(Some(version)),
            KnownFileFormat::PLGS(version) => Self::PLGS(Some(version)),
            KnownFileFormat::PLink(version) => Self::PLink(Some(version)),
            KnownFileFormat::PowerNovo(version) => Self::PowerNovo(Some(version)),
            KnownFileFormat::Proteoscape(version) => Self::Proteoscape(Some(version)),
            KnownFileFormat::Sage(version) => Self::Sage(Some(version)),
            KnownFileFormat::MSFragger(version) => Self::MSFragger(Some(version)),
            KnownFileFormat::SpectrumSequenceList(version) => {
                Self::SpectrumSequenceList(Some(version))
            }
        }
    }
}

/// A file format that might not be (fully) known
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
#[allow(clippy::upper_case_acronyms, missing_docs)]
pub enum FileFormat {
    BasicCSV(Option<BasicCSVVersion>),
    DeepNovoFamily(Option<DeepNovoFamilyVersion>),
    Fasta,
    InstaNovo(Option<InstaNovoVersion>),
    MaxQuant(Option<MaxQuantVersion>),
    MZTab,
    NovoB(Option<NovoBVersion>),
    Novor(Option<NovorVersion>),
    Opair(Option<OpairVersion>),
    Peaks(Option<PeaksVersion>),
    PepNet(Option<PepNetVersion>),
    PLGS(Option<PLGSVersion>),
    PLink(Option<PLinkVersion>),
    PowerNovo(Option<PowerNovoVersion>),
    Proteoscape(Option<ProteoscapeVersion>),
    Sage(Option<SageVersion>),
    MSFragger(Option<MSFraggerVersion>),
    SpectrumSequenceList(Option<SpectrumSequenceListVersion>),
    #[default]
    Undefined,
}

impl FileFormat {
    /// Open a file with this file format. If the file format is [`Self::Undefined`] it uses
    /// [`open_identified_peptides_file`] to automatically try all possible options based on the
    /// extension. If the file format is specific format but without a defined version, for example
    /// `Self::Peaks(None)`, all known versions are tried and the first successful one is assumed.
    /// # Errors
    /// If the file is not valid according to the file format chosen
    pub fn open<'a>(
        self,
        path: &Path,
        custom_database: Option<&'a CustomDatabase>,
    ) -> Result<GeneralIdentifiedPeptidoforms<'a>, BoxedError<'static>> {
        match self {
            Self::BasicCSV(version) => {
                BasicCSVData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::DeepNovoFamily(version) => {
                DeepNovoFamilyData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Fasta => FastaData::parse_file(path).map(|sequences| {
                let b: Box<
                    dyn Iterator<
                        Item = Result<IdentifiedPeptidoform<Linked, MaybePeptidoform>, BoxedError>,
                    >,
                > = Box::new(sequences.into_iter().map(|p| {
                    Ok(IdentifiedPeptidoform::<SemiAmbiguous, PeptidoformPresent>::from(p).cast())
                }));
                b
            }),
            Self::InstaNovo(version) => {
                InstaNovoData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::MaxQuant(version) => {
                MaxQuantData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::MZTab => MZTabData::parse_file(path, custom_database).map(|sequences| {
                let b: Box<
                    dyn Iterator<
                        Item = Result<IdentifiedPeptidoform<Linked, MaybePeptidoform>, BoxedError>,
                    >,
                > = Box::new(sequences.into_iter().map(|p| {
                    p.map(|p| {
                        IdentifiedPeptidoform::<SemiAmbiguous, MaybePeptidoform>::from(p).cast()
                    })
                }));
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
            Self::MSFragger(version) => {
                MSFraggerData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::PowerNovo(version) => {
                PowerNovoData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Proteoscape(version) => {
                ProteoscapeData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Sage(version) => SageData::parse_file(path, custom_database, false, version)
                .map(IdentifiedPeptidoformIter::into_box),
            Self::SpectrumSequenceList(version) => {
                SpectrumSequenceListData::parse_file(path, custom_database, false, version)
                    .map(IdentifiedPeptidoformIter::into_box)
            }
            Self::Undefined => open_identified_peptidoforms_file(path, custom_database, false),
        }
    }
}
