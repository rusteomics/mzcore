use std::{borrow::Cow, ops::Range};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{MolecularFormula, MultiChemical},
    identification::*,
    sequence::{
        CompoundPeptidoformIon, Peptidoform, PeptidoformIon, SemiAmbiguous, SequencePosition,
        SimpleLinear,
    },
    system::{MassOverCharge, OrderedTime, Time, isize::Charge},
};

/// A peptide that is identified by a de novo or database matching program
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct IdentifiedPeptidoform {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The local confidence, if available, in range -1.0..=1.0
    pub local_confidence: Option<Vec<f64>>,
    /// The full metadata of this peptide
    pub metadata: MetaData,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[expect(clippy::upper_case_acronyms)]
pub enum MetaData {
    /// A basic CSV format
    BasicCSV(BasicCSVData),
    /// DeepNovo/PointNovo/PGPointNovo metadata
    DeepNovoFamily(DeepNovoFamilyData),
    /// Fasta metadata
    Fasta(FastaData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
    /// InstaNovo metadata
    InstaNovo(InstaNovoData),
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
    /// MSFragger metadata
    MSFragger(MSFraggerData),
    /// SpectrumSequenceList metadata
    SpectrumSequenceList(SpectrumSequenceListData),
}

/// A peptide as stored in a identified peptide file, either a simple linear one or a cross-linked peptidoform
#[derive(Clone, Debug)]
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
        _sequence_index: SequencePosition,
        _peptidoform_index: usize,
        _peptidoform_ion_index: usize,
    ) -> crate::quantities::Multi<MolecularFormula> {
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
            MetaData::MSFragger(MSFraggerData { peptide, .. })
            | MetaData::PLGS(PLGSData { peptide, .. }) => {
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
            MetaData::SpectrumSequenceList(SpectrumSequenceListData { peptide, .. })
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
            | MetaData::MaxQuant(MaxQuantData { z, .. })
            | MetaData::NovoB(NovoBData { z, .. })
            | MetaData::PLGS(PLGSData { precursor_z: z, .. })
            | MetaData::PLink(PLinkData { z, .. })
            | MetaData::MSFragger(MSFraggerData { z, .. })
            | MetaData::InstaNovo(InstaNovoData { z, .. })
            | MetaData::MZTab(MZTabData { z, .. })
            | MetaData::BasicCSV(BasicCSVData { z, .. })
            | MetaData::SpectrumSequenceList(SpectrumSequenceListData { z, .. }) => Some(*z),
            MetaData::Peaks(PeaksData { z, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { z, .. }) => *z,
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
            | MetaData::MSFragger(MSFraggerData { rt, .. })
            | MetaData::PLGS(PLGSData {
                precursor_rt: rt, ..
            }) => Some(*rt),
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
                    .map(SpectrumId::Number)
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
            }) => Some(*mz),
            MetaData::MZTab(MZTabData { mz, .. })
            | MetaData::MaxQuant(MaxQuantData { mz, .. })
            | MetaData::DeepNovoFamily(DeepNovoFamilyData { mz, .. })
            | MetaData::MSFragger(MSFraggerData { mz, .. }) => *mz,
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
            | MetaData::PLink(PLinkData { mass, .. })
            | MetaData::MSFragger(MSFraggerData { mass, .. })
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
            MetaData::MZTab(_)
            | MetaData::MaxQuant(_)
            | MetaData::Sage(_)
            | MetaData::PLink(_)
            | MetaData::MSFragger(_)
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
            }) => protein_start.and_then(|start| protein_end.map(|end| start..end)),
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
