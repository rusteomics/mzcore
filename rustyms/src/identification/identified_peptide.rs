use std::{
    borrow::Cow,
    fmt::Display,
    ops::{Range, RangeInclusive},
    path::PathBuf,
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    formula::MultiChemical,
    identification::*,
    ontologies::CustomDatabase,
    peptidoform::{SemiAmbiguous, SimpleLinear},
    system::{usize::Charge, MassOverCharge, OrderedTime, Time},
    Peptidoform, PeptidoformIon,
};

use super::{BasicCSVData, CompoundPeptidoformIon};

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
pub struct IdentifiedPeptidoform {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The local confidence, if available, in range -1.0..=1.0
    pub local_confidence: Option<Vec<f64>>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, PartialEq, Debug, Serialize, Deserialize)]
#[expect(clippy::upper_case_acronyms)]
pub enum MetaData {
    /// A basic csv format
    BasicCSV(BasicCSVData),
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
    /// NovoB metadata
    NovoB(NovoBData),
    /// Novor metadata
    Novor(NovorData),
    /// OPair metadata
    Opair(OpairData),
    /// Peaks metadata
    Peaks(PeaksData),
    /// PepNet metadata
    PepNet(PepNetData),
    /// PLGS metadata
    PLGS(PLGSData),
    /// pLink metadata
    PLink(PLinkData),
    /// PowerNovo metadata
    PowerNovo(PowerNovoData),
    /// Sage metadata
    Sage(SageData),
    /// SpectrumSequenceList metadata
    SpectrumSequenceList(SpectrumSequenceListData),
}

/// A peptide as stored in a identified peptide file, either a simple linear one or a cross-linked peptidoform
#[derive(Debug, Clone)]
pub enum ReturnedPeptidoform<'a> {
    /// A semi ambiguous linear peptide
    PeptidoformSemiAmbiguous(&'a Peptidoform<SemiAmbiguous>),
    /// A simple linear peptide
    PeptidoformSimpleLinear(&'a Peptidoform<SimpleLinear>),
    /// A (potentially) cross-linked peptidoform
    PeptidoformIon(&'a PeptidoformIon),
    /// A (potentially) cross-linked chimeric set of peptidoforms
    CompoundPeptidoformIon(Cow<'a, CompoundPeptidoformIon>),
}

impl MultiChemical for ReturnedPeptidoform<'_> {
    fn formulas_inner(
        &self,
        _sequence_index: super::SequencePosition,
        _peptidoform_index: usize,
    ) -> super::Multi<super::MolecularFormula> {
        match self {
            Self::PeptidoformSemiAmbiguous(p) => p.formulas(),
            Self::PeptidoformSimpleLinear(p) => p.formulas(),
            Self::PeptidoformIon(p) => p.formulas(),
            Self::CompoundPeptidoformIon(p) => p.formulas(),
        }
    }
}

impl std::fmt::Display for ReturnedPeptidoform<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::PeptidoformSemiAmbiguous(p) => write!(f, "{p}"),
            Self::PeptidoformSimpleLinear(p) => write!(f, "{p}"),
            Self::PeptidoformIon(p) => write!(f, "{p}"),
            Self::CompoundPeptidoformIon(p) => write!(f, "{p}"),
        }
    }
}

impl<'a> ReturnedPeptidoform<'a> {
    /// Get the underlying peptide, or None if the underlying result was a peptidoform
    pub fn peptidoform(self) -> Option<Cow<'a, Peptidoform<SimpleLinear>>> {
        match self {
            Self::PeptidoformSemiAmbiguous(p) => Some(Cow::Owned(p.clone().into())),
            Self::PeptidoformSimpleLinear(p) => Some(Cow::Borrowed(p)),
            Self::PeptidoformIon(_) | Self::CompoundPeptidoformIon(_) => None,
        }
    }
    /// Get the underlying result as a peptidoform, if it was a peptide make a new peptidoform from it
    pub fn peptidoform_ion(self) -> Option<Cow<'a, PeptidoformIon>> {
        match self {
            Self::PeptidoformSemiAmbiguous(p) => Some(Cow::Owned(p.clone().into())),
            Self::PeptidoformSimpleLinear(p) => Some(Cow::Owned(p.clone().into())),
            Self::PeptidoformIon(p) => Some(Cow::Borrowed(p)),
            Self::CompoundPeptidoformIon(_) => None,
        }
    }
    /// Get the underlying result as a compound peptidoform, if it was a peptide make a new peptidoform from it
    pub fn compound_peptidoform_ion(self) -> Cow<'a, CompoundPeptidoformIon> {
        match self {
            Self::PeptidoformSemiAmbiguous(p) => Cow::Owned(p.clone().into()),
            Self::PeptidoformSimpleLinear(p) => Cow::Owned(p.clone().into()),
            Self::PeptidoformIon(p) => Cow::Owned(p.clone().into()),
            Self::CompoundPeptidoformIon(p) => p,
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
            Self::PeptidoformSemiAmbiguous(p) => {
                p.display(f, show_global_mods, specification_compliant)
            }
            Self::PeptidoformSimpleLinear(p) => {
                p.display(f, show_global_mods, specification_compliant)
            }
            Self::PeptidoformIon(p) => p.display(f, show_global_mods, specification_compliant),
            Self::CompoundPeptidoformIon(p) => p.display(f, specification_compliant),
        }
    }
}

impl IdentifiedPeptidoform {
    /// Get the peptide, as pLink can have cross-linked peptides the return type is either a simple peptide or a cross-linked peptidoform
    pub fn peptidoform(&self) -> Option<ReturnedPeptidoform<'_>> {
        match &self.metadata {
            MetaData::Novor(NovorData { peptide, .. })
            | MetaData::InstaNovo(InstaNovoData { peptide, .. })
            | MetaData::Opair(OpairData { peptide, .. })
            | MetaData::PepNet(PepNetData { peptide, .. })
            | MetaData::PowerNovo(PowerNovoData { peptide, .. })
            | MetaData::Sage(SageData { peptide, .. }) => {
                Some(ReturnedPeptidoform::PeptidoformSemiAmbiguous(peptide))
            }
            MetaData::PLGS(PLGSData { peptide, .. }) => {
                Some(ReturnedPeptidoform::PeptidoformSimpleLinear(peptide))
            }
            MetaData::Peaks(PeaksData { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(ReturnedPeptidoform::PeptidoformSemiAmbiguous(&peptide.1[0]))
                } else {
                    Some(ReturnedPeptidoform::CompoundPeptidoformIon(Cow::Owned(
                        peptide.1.clone().into(),
                    )))
                }
            }
            MetaData::BasicCSV(BasicCSVData { sequence, .. }) => Some(
                ReturnedPeptidoform::CompoundPeptidoformIon(Cow::Borrowed(sequence)),
            ),
            MetaData::MSFragger(MSFraggerData { peptide, .. })
            | MetaData::SpectrumSequenceList(SpectrumSequenceListData { peptide, .. })
            | MetaData::MaxQuant(MaxQuantData { peptide, .. })
            | MetaData::MZTab(MZTabData { peptide, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => peptide
                .as_ref()
                .map(ReturnedPeptidoform::PeptidoformSemiAmbiguous),
            MetaData::Fasta(f) => Some(ReturnedPeptidoform::PeptidoformSemiAmbiguous(f.peptide())),
            MetaData::PLink(PLinkData { peptidoform, .. }) => {
                Some(ReturnedPeptidoform::PeptidoformIon(peptidoform))
            }
            MetaData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => {
                if score_forward >= score_reverse {
                    peptide_forward
                        .as_ref()
                        .map(ReturnedPeptidoform::PeptidoformSemiAmbiguous)
                } else {
                    peptide_reverse
                        .as_ref()
                        .map(ReturnedPeptidoform::PeptidoformSemiAmbiguous)
                }
            }
        }
    }

    /// Get the format and version for this peptidoform
    pub const fn format(&self) -> KnownFileFormat {
        match &self.metadata {
            MetaData::BasicCSV(BasicCSVData { version, .. }) => KnownFileFormat::BasicCSV(*version),
            MetaData::SpectrumSequenceList(SpectrumSequenceListData { version, .. }) => {
                KnownFileFormat::SpectrumSequenceList(*version)
            }
            MetaData::DeepNovoFamily(DeepNovoFamilyData { version, .. }) => {
                KnownFileFormat::DeepNovoFamily(*version)
            }
            MetaData::Fasta(_) => KnownFileFormat::Fasta,
            MetaData::InstaNovo(InstaNovoData { version, .. }) => {
                KnownFileFormat::InstaNovo(*version)
            }
            MetaData::MaxQuant(MaxQuantData { version, .. }) => KnownFileFormat::MaxQuant(*version),
            MetaData::MSFragger(MSFraggerData { version, .. }) => {
                KnownFileFormat::MSFragger(*version)
            }
            MetaData::MZTab(_) => KnownFileFormat::MZTab,
            MetaData::NovoB(NovoBData { version, .. }) => KnownFileFormat::NovoB(*version),
            MetaData::Novor(NovorData { version, .. }) => KnownFileFormat::Novor(*version),
            MetaData::Opair(OpairData { version, .. }) => KnownFileFormat::Opair(*version),
            MetaData::Peaks(PeaksData { version, .. }) => KnownFileFormat::Peaks(*version),
            MetaData::PepNet(PepNetData { version, .. }) => KnownFileFormat::PepNet(*version),
            MetaData::PLGS(PLGSData { version, .. }) => KnownFileFormat::PLGS(*version),
            MetaData::PLink(PLinkData { version, .. }) => KnownFileFormat::PLink(*version),
            MetaData::PowerNovo(PowerNovoData { version, .. }) => {
                KnownFileFormat::PowerNovo(*version)
            }
            MetaData::Sage(SageData { version, .. }) => KnownFileFormat::Sage(*version),
        }
    }

    /// Get the identifier
    pub fn id(&self) -> String {
        match &self.metadata {
            MetaData::Peaks(PeaksData {
                id,
                scan_number,
                feature,
                ..
            }) => id.map_or(
                scan_number.as_ref().map_or(
                    feature
                        .as_ref()
                        .map_or("-".to_string(), ToString::to_string),
                    |s| s.iter().join(";"),
                ),
                |i| i.to_string(),
            ),
            MetaData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => scan.iter().join(";"),
            MetaData::Novor(NovorData {
                id, scan_number, ..
            }) => id.unwrap_or(*scan_number).to_string(),
            MetaData::Opair(OpairData {
                scan_number: scan, ..
            })
            | MetaData::NovoB(NovoBData { scan, .. })
            | MetaData::SpectrumSequenceList(SpectrumSequenceListData { scan, .. })
            | MetaData::InstaNovo(InstaNovoData {
                scan_number: scan, ..
            })
            | MetaData::BasicCSV(BasicCSVData {
                scan_index: scan, ..
            }) => scan.to_string(),
            MetaData::Sage(SageData { id, .. }) | MetaData::MZTab(MZTabData { id, .. }) => {
                id.to_string()
            }
            MetaData::Fasta(f) => f.identifier().accession().to_string(),
            MetaData::MSFragger(MSFraggerData { scan, .. }) => scan.to_string(),
            MetaData::PLink(PLinkData { order, .. }) => order.to_string(),
            MetaData::MaxQuant(MaxQuantData {
                id, scan_number, ..
            }) => id.map_or_else(|| scan_number.iter().join(";"), |id| id.to_string()),
            MetaData::PowerNovo(PowerNovoData { scan, .. }) => {
                scan.as_ref().map_or("-".to_string(), ToString::to_string)
            }
            MetaData::PepNet(_) => "-".to_string(),
            MetaData::PLGS(PLGSData {
                peptide_component_id,
                ..
            }) => peptide_component_id.to_string(),
        }
    }

    /// Get the original local confidence, it is the same length as the peptide with a local score
    pub fn local_confidence(&self) -> Option<&[f64]> {
        match &self.metadata {
            MetaData::InstaNovo(InstaNovoData {
                local_confidence, ..
            })
            | MetaData::PowerNovo(PowerNovoData {
                local_confidence, ..
            })
            | MetaData::PepNet(PepNetData {
                local_confidence, ..
            }) => Some(local_confidence),

            MetaData::Peaks(PeaksData {
                local_confidence, ..
            })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData {
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
            MetaData::Novor(NovorData { z, .. })
            | MetaData::Opair(OpairData { z, .. })
            | MetaData::Sage(SageData { z, .. })
            | MetaData::MSFragger(MSFraggerData { z, .. })
            | MetaData::MaxQuant(MaxQuantData { z, .. })
            | MetaData::NovoB(NovoBData { z, .. })
            | MetaData::PLGS(PLGSData { precursor_z: z, .. })
            | MetaData::PLink(PLinkData { z, .. })
            | MetaData::InstaNovo(InstaNovoData { z, .. })
            | MetaData::MZTab(MZTabData { z, .. })
            | MetaData::BasicCSV(BasicCSVData { z, .. }) => Some(*z),
            MetaData::Peaks(PeaksData { z, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { z, .. }) => *z,
            MetaData::SpectrumSequenceList(SpectrumSequenceListData { z, .. }) => {
                (z.value >= 0).then_some(Charge::new::<crate::system::charge::e>(z.value as usize))
            }
            MetaData::Fasta(_) | MetaData::PowerNovo(_) | MetaData::PepNet(_) => None,
        }
    }

    /// Which fragmentation mode was used, if known
    pub fn mode(&self) -> Option<&str> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mode, .. })
            | MetaData::BasicCSV(BasicCSVData { mode, .. })
            | MetaData::MaxQuant(MaxQuantData { mode, .. }) => mode.as_deref(),
            _ => None,
        }
    }

    /// The retention time, if known
    pub fn retention_time(&self) -> Option<Time> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { rt, .. })
            | MetaData::Opair(OpairData { rt, .. })
            | MetaData::Sage(SageData { rt, .. })
            | MetaData::PLGS(PLGSData {
                precursor_rt: rt, ..
            })
            | MetaData::MSFragger(MSFraggerData { rt, .. }) => Some(*rt),
            MetaData::MaxQuant(MaxQuantData { rt, .. })
            | MetaData::Novor(NovorData { rt, .. })
            | MetaData::SpectrumSequenceList(SpectrumSequenceListData { rt, .. })
            | MetaData::MZTab(MZTabData { rt, .. }) => *rt,
            MetaData::DeepNovoFamily(_)
            | MetaData::InstaNovo(_)
            | MetaData::Fasta(_)
            | MetaData::NovoB(_)
            | MetaData::PowerNovo(_)
            | MetaData::PepNet(_)
            | MetaData::PLink(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    /// The scans per rawfile that are at the basis for this identified peptide, if the rawfile is unknown there will be one
    pub fn scans(&self) -> SpectrumIds {
        match &self.metadata {
            MetaData::Peaks(PeaksData {
                raw_file,
                scan_number,
                ..
            }) => scan_number
                .as_ref()
                .map_or(SpectrumIds::None, |scan_number| {
                    raw_file.clone().map_or_else(
                        || {
                            SpectrumIds::FileNotKnown(
                                scan_number
                                    .iter()
                                    .flat_map(|s| s.scans.clone())
                                    .map(SpectrumId::Number)
                                    .collect(),
                            )
                        },
                        |raw_file| {
                            SpectrumIds::FileKnown(vec![(
                                raw_file,
                                scan_number
                                    .iter()
                                    .flat_map(|s| s.scans.clone())
                                    .map(SpectrumId::Number)
                                    .collect(),
                            )])
                        },
                    )
                }),
            MetaData::Novor(NovorData { scan_number, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Number(*scan_number)])
            }
            MetaData::NovoB(NovoBData { scan, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)])
            }
            MetaData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => SpectrumIds::FileNotKnown(
                scan.iter()
                    .flat_map(|s| s.scans.clone())
                    .map(SpectrumId::Index)
                    .collect(),
            ),

            MetaData::Opair(OpairData {
                raw_file,
                scan_number,
                ..
            })
            | MetaData::InstaNovo(InstaNovoData {
                raw_file,
                scan_number,
                ..
            }) => SpectrumIds::FileKnown(vec![(
                raw_file.clone(),
                vec![SpectrumId::Number(*scan_number)],
            )]),
            MetaData::SpectrumSequenceList(SpectrumSequenceListData {
                raw_file,
                scan: scan_index,
                ..
            })
            | MetaData::BasicCSV(BasicCSVData {
                raw_file,
                scan_index,
                ..
            }) => SpectrumIds::FileKnown(vec![(
                raw_file.clone(),
                vec![SpectrumId::Index(*scan_index)],
            )]),

            MetaData::PowerNovo(PowerNovoData { raw_file, scan, .. }) => {
                scan.as_ref().map_or(SpectrumIds::None, |scan| {
                    raw_file.clone().map_or_else(
                        || SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)]),
                        |raw_file| {
                            SpectrumIds::FileKnown(vec![(raw_file, vec![SpectrumId::Index(*scan)])])
                        },
                    )
                })
            }

            MetaData::MaxQuant(MaxQuantData {
                raw_file,
                scan_number,
                ..
            }) => raw_file.as_ref().map_or_else(
                || {
                    SpectrumIds::FileNotKnown(
                        scan_number
                            .iter()
                            .copied()
                            .map(SpectrumId::Number)
                            .collect(),
                    )
                },
                |raw_file| {
                    SpectrumIds::FileKnown(vec![(
                        raw_file.clone(),
                        scan_number
                            .iter()
                            .copied()
                            .map(SpectrumId::Number)
                            .collect(),
                    )])
                },
            ),
            MetaData::MZTab(MZTabData { spectra_ref, .. }) => spectra_ref.clone(),
            MetaData::MSFragger(MSFraggerData { raw_file, scan, .. }) => {
                raw_file.clone().map_or_else(
                    || SpectrumIds::FileNotKnown(vec![scan.clone()]),
                    |raw_file| SpectrumIds::FileKnown(vec![(raw_file, vec![scan.clone()])]),
                )
            }
            MetaData::PLink(PLinkData {
                raw_file,
                scan_number: scan,
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
            MetaData::PLGS(PLGSData {
                precursor_lift_off_rt,
                precursor_touch_down_rt,
                ..
            }) => SpectrumIds::FileNotKnown(vec![SpectrumId::RetentionTime(
                OrderedTime::from(*precursor_lift_off_rt)
                    ..=OrderedTime::from(*precursor_touch_down_rt),
            )]),
            MetaData::Fasta(_) | MetaData::PepNet(_) => SpectrumIds::None,
        }
    }

    /// Get the mz as experimentally determined
    pub fn experimental_mz(&self) -> Option<MassOverCharge> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mz, .. })
            | MetaData::Novor(NovorData { mz, .. })
            | MetaData::Opair(OpairData { mz, .. })
            | MetaData::InstaNovo(InstaNovoData { mz, .. })
            | MetaData::PLGS(PLGSData {
                precursor_mz: mz, ..
            })
            | MetaData::MSFragger(MSFraggerData { mz, .. }) => Some(*mz),
            MetaData::MZTab(MZTabData { mz, .. })
            | MetaData::MaxQuant(MaxQuantData { mz, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { mz, .. }) => *mz,
            MetaData::Sage(SageData { mass, z, .. })
            | MetaData::NovoB(NovoBData { mass, z, .. })
            | MetaData::PLink(PLinkData { mass, z, .. }) => {
                Some(MassOverCharge::new::<crate::system::mz>(
                    mass.value / (z.value as f64),
                ))
            }
            MetaData::Fasta(_)
            | MetaData::SpectrumSequenceList(_)
            | MetaData::PowerNovo(_)
            | MetaData::PepNet(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    /// Get the mass as experimentally determined
    pub fn experimental_mass(&self) -> Option<crate::system::Mass> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { mass, mz, z, .. }) => {
                mass.map_or(z.map_or(None, |z| Some(*mz * z.to_float())), Some)
            }
            MetaData::Novor(NovorData { mass, .. })
            | MetaData::Opair(OpairData { mass, .. })
            | MetaData::PLGS(PLGSData {
                precursor_mass: mass,
                ..
            })
            | MetaData::NovoB(NovoBData { mass, .. })
            | MetaData::MSFragger(MSFraggerData { mass, .. })
            | MetaData::PLink(PLinkData { mass, .. })
            | MetaData::Sage(SageData { mass, .. }) => Some(*mass),
            MetaData::MaxQuant(MaxQuantData { mass, .. }) => *mass,
            MetaData::MZTab(MZTabData { mz, z, .. }) => mz.map(|mz| mz * z.to_float()),
            MetaData::InstaNovo(InstaNovoData { mz, z, .. }) => Some(*mz * z.to_float()),
            MetaData::DeepNovoFamily(DeepNovoFamilyData { mz, z, .. }) => {
                mz.and_then(|mz| z.map(|z| (mz, z)).map(|(mz, z)| mz * z.to_float()))
            }
            MetaData::Fasta(_)
            | MetaData::PowerNovo(_)
            | MetaData::SpectrumSequenceList(_)
            | MetaData::PepNet(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    /// Get the absolute ppm error between the experimental and theoretical precursor mass
    pub fn ppm_error(&self) -> Option<crate::system::Ratio> {
        match &self.metadata {
            MetaData::PepNet(p) => Some(p.ppm_diff),
            MetaData::NovoB(p) => Some(if p.score_forward >= p.score_reverse {
                p.ppm_diff_forward
            } else {
                p.ppm_diff_reverse
            }),
            _ => {
                let exp_mass = self.experimental_mass()?;
                let theo_mass = self
                    .peptidoform()
                    .and_then(|p| p.formulas().to_vec().pop())
                    .map(|f| f.monoisotopic_mass())?;

                Some(theo_mass.ppm(exp_mass))
            }
        }
    }

    /// Get the absolute mass error between the experimental and theoretical precursor mass
    pub fn mass_error(&self) -> Option<crate::system::Mass> {
        let exp_mass = self.experimental_mass()?;
        let theo_mass = self
            .peptidoform()
            .and_then(|p| p.formulas().to_vec().pop())
            .map(|f| f.monoisotopic_mass())?;

        Some((exp_mass - theo_mass).abs())
    }

    /// Get the protein name if this was database matched data
    pub fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        match &self.metadata {
            MetaData::Peaks(PeaksData {
                protein_accession, ..
            }) => protein_accession.clone(),
            MetaData::Opair(OpairData { protein_name, .. }) => Some(protein_name.clone()),
            MetaData::PLGS(PLGSData {
                protein_description,
                ..
            }) => Some(protein_description.clone()),
            MetaData::MSFragger(MSFraggerData { protein, .. }) => Some(protein.clone()),
            MetaData::MZTab(MZTabData { accession, .. }) => accession
                .as_ref()
                .map(|a| FastaIdentifier::Undefined(a.clone())),
            MetaData::NovoB(_)
            | MetaData::MaxQuant(_)
            | MetaData::Sage(_)
            | MetaData::PLink(_)
            | MetaData::Novor(_)
            | MetaData::Fasta(_)
            | MetaData::DeepNovoFamily(_)
            | MetaData::InstaNovo(_)
            | MetaData::PowerNovo(_)
            | MetaData::SpectrumSequenceList(_)
            | MetaData::PepNet(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    /// Get the protein id if this was database matched data
    pub const fn protein_id(&self) -> Option<usize> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { protein_id, .. }) => *protein_id,
            MetaData::Novor(NovorData { protein, .. }) => *protein,
            MetaData::PLGS(PLGSData { protein_id, .. }) => Some(*protein_id),
            MetaData::MSFragger(_)
            | MetaData::MZTab(_)
            | MetaData::MaxQuant(_)
            | MetaData::Sage(_)
            | MetaData::PLink(_)
            | MetaData::NovoB(_)
            | MetaData::Opair(_)
            | MetaData::Fasta(_)
            | MetaData::PowerNovo(_)
            | MetaData::DeepNovoFamily(_)
            | MetaData::SpectrumSequenceList(_)
            | MetaData::InstaNovo(_)
            | MetaData::PepNet(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    /// Get the protein location if this was database matched data
    pub fn protein_location(&self) -> Option<Range<usize>> {
        match &self.metadata {
            MetaData::Peaks(PeaksData { start, end, .. }) => start.and_then(|s| end.map(|e| s..e)),
            MetaData::Novor(NovorData {
                protein_start,
                peptide,
                ..
            }) => protein_start.map(|s| s..s + peptide.len()),
            MetaData::Opair(OpairData {
                protein_location, ..
            }) => Some(protein_location.clone()),
            MetaData::PLGS(PLGSData {
                peptide_start,
                peptide,
                ..
            }) => Some(*peptide_start..*peptide_start + peptide.len()),
            MetaData::MSFragger(MSFraggerData {
                protein_start,
                protein_end,
                ..
            }) => Some(*protein_start..*protein_end),
            MetaData::MZTab(MZTabData { start, end, .. }) => start.and_then(|s| end.map(|e| s..e)),
            MetaData::InstaNovo(_)
            | MetaData::DeepNovoFamily(_)
            | MetaData::MaxQuant(_)
            | MetaData::Sage(_)
            | MetaData::PLink(_)
            | MetaData::NovoB(_)
            | MetaData::Fasta(_)
            | MetaData::PowerNovo(_)
            | MetaData::SpectrumSequenceList(_)
            | MetaData::PepNet(_)
            | MetaData::BasicCSV(_) => None,
        }
    }

    // Get the matched fragments, potentially with m/z and intensity
    // #[doc(hidden)]
    // pub fn matched_fragments(
    //     &self,
    // ) -> Option<Vec<(Option<MassOverCharge>, Option<f64>, Fragment)>> {
    //     // OPair, MaxQuant, PLGS
    //     None
    // }
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
    /// A scan number, unless there is a better alternative should be interpreted as index+1
    Number(usize),
    /// Time range, assumes all MS2 spectra within this range are selected
    RetentionTime(RangeInclusive<OrderedTime>),
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
            Self::Number(i) => write!(f, "\x23{i}"),
            Self::Native(n) => write!(f, "{n}"),
            Self::RetentionTime(n) => {
                write!(f, "{:.3} — {:.3} min", n.start().value, n.end().value)
            }
        }
    }
}

impl SpectrumId {
    /// Get the index if this is an index or scan number
    pub fn index(&self) -> Option<usize> {
        match self {
            Self::Index(i) => Some(*i),
            Self::Number(i) => Some(i - 1),
            Self::Native(_) | Self::RetentionTime(_) => None,
        }
    }

    /// Get the native ID if this is a native ID
    pub fn native(&self) -> Option<&str> {
        match self {
            Self::Native(n) => Some(n),
            Self::Index(_) | Self::RetentionTime(_) | Self::Number(_) => None,
        }
    }
}

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
    Self: std::marker::Sized,
{
    /// The source data where the peptides are parsed form
    type Source;

    /// The format type
    type Format: Clone;

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
pub struct IdentifiedPeptidoformIter<
    'lifetime,
    R: IdentifiedPeptidoformSource,
    I: Iterator<Item = Result<R::Source, CustomError>>,
> {
    iter: Box<I>,
    format: Option<R::Format>,
    custom_database: Option<&'lifetime CustomDatabase>,
    keep_all_columns: bool,
    peek: Option<Result<R, CustomError>>,
}

impl<
        R: IdentifiedPeptidoformSource + Clone,
        I: Iterator<Item = Result<R::Source, CustomError>>,
    > IdentifiedPeptidoformIter<'_, R, I>
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

impl<'lifetime, R, I> IdentifiedPeptidoformIter<'lifetime, R, I>
where
    R: IdentifiedPeptidoformSource + Into<IdentifiedPeptidoform> + 'lifetime,
    I: Iterator<Item = Result<R::Source, CustomError>> + 'lifetime,
    R::Format: 'static,
{
    /// Make this into a generic boxed iterator for merging with other identified peptidoform formats
    pub fn into_box(
        self,
    ) -> Box<dyn Iterator<Item = Result<IdentifiedPeptidoform, CustomError>> + 'lifetime> {
        Box::new(self.map(|p: Result<R, CustomError>| match p {
            Ok(p) => Ok(p.into()),
            Err(e) => Err(e),
        }))
    }
}

/// Test a dataset for common errors in identified peptide parsing
/// # Errors
/// * If the peptide was not identified as the correct version of the format (see parameters).
/// * See errors at``test_identified_peptide()``
#[expect(clippy::missing_panics_doc)]
#[cfg(test)]
pub fn test_format<T: IdentifiedPeptidoformSource + Into<IdentifiedPeptidoform>>(
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
    for peptide in
        T::parse_reader(reader, custom_database, false, None).map_err(|e| e.to_string())?
    {
        let peptide: IdentifiedPeptidoform = peptide.map_err(|e| e.to_string())?.into();
        number += 1;

        test_identified_peptide(&peptide, allow_mass_mods, expect_lc)?;

        if format.as_ref().is_some_and(|f| {
            peptide
                .format()
                .version()
                .is_some_and(|version| f.to_string() != version)
        }) {
            return Err(format!(
                "Peptide {} was detected as the wrong version ({} instead of {})",
                peptide.id(),
                peptide
                    .format()
                    .version()
                    .unwrap_or_else(|| "-".to_string()),
                format.unwrap(),
            ));
        }
    }
    Ok(number)
}

/// Test a peptide for common errors in identified peptide parsing
/// # Errors
/// * If the local confidence has to be there and is not there (see parameter).
/// * If the local confidence is not the same length as the peptide.
/// * If the score of the peptide is outside of the range -1.0..=1.0.
/// * If any of the local scores is outside of range -1.0..=1.0.
/// * If the peptide contains mass modifications (see parameters).
#[expect(clippy::missing_panics_doc)]
#[cfg(test)]
pub fn test_identified_peptide(
    peptide: &IdentifiedPeptidoform,
    allow_mass_mods: bool,
    expect_lc: bool,
) -> Result<(), String> {
    if peptide
        .peptidoform()
        .and_then(ReturnedPeptidoform::peptidoform)
        .map(|p| p.len())
        != peptide.local_confidence().map(<[f64]>::len)
    {
        if expect_lc && peptide.local_confidence.is_none() {
            return Err(format!(
                "No local confidence was provided for peptide {}",
                peptide.id()
            ));
        } else if peptide.local_confidence.is_some() {
            return Err(format!("The local confidence ({}) does not have the same number of elements as the peptide ({}) for peptide {}", peptide.local_confidence().map_or(0, <[f64]>::len), peptide.peptidoform().and_then(ReturnedPeptidoform::peptidoform).map_or(0,|p| p.len()), peptide.id()));
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
        .local_confidence
        .as_ref()
        .is_some_and(|s| s.iter().any(|s| !(-1.0..=1.0).contains(s)))
    {
        return Err(format!(
            "The local score {} for peptide {} is outside of range",
            peptide.local_confidence().unwrap().iter().join(","),
            peptide.id()
        ));
    }
    if !allow_mass_mods
        && peptide.peptidoform().is_some_and(|p| {
            p.compound_peptidoform_ion().peptidoforms().any(|p| {
                p.sequence().iter().any(|s| {
                    s.modifications.iter().any(|m| {
                        m.simple().is_some_and(|m| {
                            matches!(**m, crate::modification::SimpleModificationInner::Mass(_))
                        })
                    })
                })
            })
        })
    {
        return Err(format!(
            "Peptide {} contains mass modifications, sequence {}",
            peptide.id(),
            peptide.peptidoform().unwrap(),
        ));
    }
    if let Err(err) = peptide.peptidoform().map_or(Ok(()), |p| {
        p.compound_peptidoform_ion()
            .peptidoforms()
            .try_for_each(Peptidoform::enforce_modification_rules)
    }) {
        return Err(format!(
            "Peptide {} contains misplaced modifications, sequence {}\n{}",
            peptide.id(),
            peptide.peptidoform().unwrap(),
            err
        ));
    }
    Ok(())
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
                scans: s
                    .split(' ')
                    .map(str::parse)
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|_| ())?,
            })
        }
    }
}
