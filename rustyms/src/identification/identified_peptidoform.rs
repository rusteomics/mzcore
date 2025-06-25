use std::{borrow::Cow, marker::PhantomData, ops::Range};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    identification::*,
    sequence::{
        AtLeast, CompoundPeptidoformIon, HasPeptidoformImpl, Linear, Peptidoform, SemiAmbiguous,
        SimpleLinear, UnAmbiguous,
    },
    system::{MassOverCharge, OrderedTime, Time, isize::Charge},
};

/// A peptide that is identified by a _de novo_ or database matching program
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct IdentifiedPeptidoform<Complexity, PeptidoformAvailability> {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The local confidence, if available, in range -1.0..=1.0
    pub local_confidence: Option<Vec<f64>>,
    /// The full metadata of this peptide
    pub metadata: IdentifiedPeptidoformData,
    /// The marker for the complexity, Linked means full [`CompoundPeptidoformIon`] anything below means [`Peptidoform`], see [Complexity](crate::sequence::Complexity)
    pub(super) complexity_marker: PhantomData<Complexity>,
    /// The marker for availability of the peptidoform, see [`PeptidoformAvailability`]
    pub(super) peptidoform_availability_marker: PhantomData<PeptidoformAvailability>,
}

/// The definition of all special metadata for all types of identified peptides that can be read
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
#[expect(clippy::upper_case_acronyms)]
pub enum IdentifiedPeptidoformData {
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

impl<PeptidoformAvailability> IdentifiedPeptidoform<Linear, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a linear peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<Linear>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Sage(SageData { peptide, .. }) => Some(peptide.as_ref()),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { peptide, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { peptide, .. }) => Some(peptide.as_ref()),
            IdentifiedPeptidoformData::Peaks(PeaksData { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(peptide.1[0].as_ref())
                } else {
                    None
                }
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                peptide,
                ..
            })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { peptide, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            IdentifiedPeptidoformData::Fasta(f) => Some(f.peptide().as_ref()),
            IdentifiedPeptidoformData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => if score_forward >= score_reverse {
                peptide_forward.as_ref()
            } else {
                peptide_reverse.as_ref()
            }
            .map(AsRef::as_ref),
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { sequence, .. }) => sequence
                .singular_peptidoform_ref()
                .and_then(|p| p.as_linear()),
            IdentifiedPeptidoformData::PLink(PLinkData { peptidoform, .. }) => {
                peptidoform.singular_ref().and_then(|p| p.as_linear())
            }
        }
    }
}

impl<PeptidoformAvailability> IdentifiedPeptidoform<SimpleLinear, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a simple linear peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<SimpleLinear>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Sage(SageData { peptide, .. }) => Some(peptide.as_ref()),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { peptide, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { peptide, .. }) => Some(peptide),
            IdentifiedPeptidoformData::Peaks(PeaksData { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(peptide.1[0].as_ref())
                } else {
                    None
                }
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                peptide,
                ..
            })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { peptide, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            IdentifiedPeptidoformData::Fasta(f) => Some(f.peptide().as_ref()),
            IdentifiedPeptidoformData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => if score_forward >= score_reverse {
                peptide_forward.as_ref()
            } else {
                peptide_reverse.as_ref()
            }
            .map(AsRef::as_ref),
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { sequence, .. }) => sequence
                .singular_peptidoform_ref()
                .and_then(|p| p.as_simple_linear()),
            IdentifiedPeptidoformData::PLink(PLinkData { peptidoform, .. }) => peptidoform
                .singular_ref()
                .and_then(|p| p.as_simple_linear()),
        }
    }
}

impl<PeptidoformAvailability> IdentifiedPeptidoform<SemiAmbiguous, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a semi ambiguous peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<SemiAmbiguous>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Sage(SageData { peptide, .. }) => Some(peptide),
            IdentifiedPeptidoformData::Peaks(PeaksData { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(&peptide.1[0])
                } else {
                    None
                }
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                peptide,
                ..
            })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { peptide, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref()
            }
            IdentifiedPeptidoformData::Fasta(f) => Some(f.peptide()),
            IdentifiedPeptidoformData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => {
                if score_forward >= score_reverse {
                    peptide_forward.as_ref()
                } else {
                    peptide_reverse.as_ref()
                }
            }
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { peptide, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { peptide, .. }) => {
                peptide.as_semi_ambiguous()
            }
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { sequence, .. }) => sequence
                .singular_peptidoform_ref()
                .and_then(|p| p.as_semi_ambiguous()),
            IdentifiedPeptidoformData::PLink(PLinkData { peptidoform, .. }) => peptidoform
                .singular_ref()
                .and_then(|p| p.as_semi_ambiguous()),
        }
    }
}

impl<PeptidoformAvailability> IdentifiedPeptidoform<UnAmbiguous, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as an unambiguous peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<UnAmbiguous>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Sage(SageData { peptide, .. }) => peptide.as_unambiguous(),
            IdentifiedPeptidoformData::Peaks(PeaksData { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    peptide.1[0].as_unambiguous()
                } else {
                    None
                }
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                peptide,
                ..
            })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { peptide, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().and_then(|p| p.as_unambiguous())
            }
            IdentifiedPeptidoformData::Fasta(f) => f.peptide().as_unambiguous(),
            IdentifiedPeptidoformData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => {
                if score_forward >= score_reverse {
                    peptide_forward.as_ref().and_then(|p| p.as_unambiguous())
                } else {
                    peptide_reverse.as_ref().and_then(|p| p.as_unambiguous())
                }
            }
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { peptide, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { peptide, .. }) => peptide.as_unambiguous(),
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { sequence, .. }) => sequence
                .singular_peptidoform_ref()
                .and_then(|p| p.as_unambiguous()),
            IdentifiedPeptidoformData::PLink(PLinkData { peptidoform, .. }) => {
                peptidoform.singular_ref().and_then(|p| p.as_unambiguous())
            }
        }
    }
}

impl<Complexity, PeptidoformAvailability>
    IdentifiedPeptidoform<Complexity, PeptidoformAvailability>
{
    /// Check if this identified peptidoform is linear and contains a peptide
    pub fn into_linear(self) -> Option<IdentifiedPeptidoform<Linear, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| {
                p.singular_peptidoform_ref()
                    .is_some_and(Peptidoform::is_linear)
            })
            .then(|| self.mark())
    }

    /// Check if this identified peptidoform is simple linear and contains a peptide
    pub fn into_simple_linear(
        self,
    ) -> Option<IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| {
                p.singular_peptidoform_ref()
                    .is_some_and(Peptidoform::is_simple_linear)
            })
            .then(|| self.mark())
    }

    /// Check if this identified peptidoform is semi ambiguous and contains a peptide
    pub fn into_semi_ambiguous(
        self,
    ) -> Option<IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| {
                p.singular_peptidoform_ref()
                    .is_some_and(Peptidoform::is_semi_ambiguous)
            })
            .then(|| self.mark())
    }

    /// Check if this identified peptidoform is unambiguous and contains a peptide
    pub fn into_unambiguous(
        self,
    ) -> Option<IdentifiedPeptidoform<UnAmbiguous, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| {
                p.singular_peptidoform_ref()
                    .is_some_and(Peptidoform::is_unambiguous)
            })
            .then(|| self.mark())
    }
}

impl HasPeptidoformImpl for IdentifiedPeptidoform<Linear, PeptidoformPresent> {
    type Complexity = Linear;
    fn peptidoform(&self) -> &Peptidoform<Linear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &IdentifiedPeptidoform<Linear, PeptidoformPresent> {
    type Complexity = Linear;
    fn peptidoform(&self) -> &Peptidoform<Linear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent> {
    type Complexity = SimpleLinear;
    fn peptidoform(&self) -> &Peptidoform<SimpleLinear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent> {
    type Complexity = SimpleLinear;
    fn peptidoform(&self) -> &Peptidoform<SimpleLinear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent> {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<SemiAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent> {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<SemiAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for IdentifiedPeptidoform<UnAmbiguous, PeptidoformPresent> {
    type Complexity = UnAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<UnAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &IdentifiedPeptidoform<UnAmbiguous, PeptidoformPresent> {
    type Complexity = UnAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<UnAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl IdentifiedPeptidoform<Linear, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a linear peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<Linear>> {
        self.inner_peptidoform()
    }
}

impl IdentifiedPeptidoform<SimpleLinear, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a simple linear peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<SimpleLinear>> {
        self.inner_peptidoform()
    }
}

impl IdentifiedPeptidoform<SemiAmbiguous, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a semi ambiguous peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<SemiAmbiguous>> {
        self.inner_peptidoform()
    }
}

impl IdentifiedPeptidoform<UnAmbiguous, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as an unambiguous peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<UnAmbiguous>> {
        self.inner_peptidoform()
    }
}

impl<Complexity, PeptidoformAvailability>
    IdentifiedPeptidoform<Complexity, PeptidoformAvailability>
{
    /// Mark this with the following complexity, be warned that the complexity level is not checked.
    fn mark<C, A>(self) -> IdentifiedPeptidoform<C, A> {
        IdentifiedPeptidoform {
            score: self.score,
            local_confidence: self.local_confidence,
            metadata: self.metadata,
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }

    /// Cast this identified peptidoform to a higher complexity level. This does not change the
    /// content of the peptidoform. It only allows to pass this as higher complexity if needed.
    pub fn cast<
        NewComplexity: AtLeast<Complexity>,
        NewAvailability: From<PeptidoformAvailability>,
    >(
        self,
    ) -> IdentifiedPeptidoform<NewComplexity, NewAvailability> {
        self.mark()
    }
}

impl<Complexity, PeptidoformAvailability> MetaData
    for IdentifiedPeptidoform<Complexity, PeptidoformAvailability>
{
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Sage(SageData { peptide, .. }) => {
                Some(Cow::Owned(peptide.clone().into()))
            }
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { peptide, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { peptide, .. }) => {
                Some(Cow::Owned(peptide.clone().into()))
            }
            IdentifiedPeptidoformData::Peaks(PeaksData { peptide, .. }) => {
                Some(Cow::Owned(peptide.1.clone().into()))
            }
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { sequence, .. }) => {
                Some(Cow::Borrowed(sequence))
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                peptide,
                ..
            })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { peptide, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
            }
            IdentifiedPeptidoformData::Fasta(f) => Some(Cow::Owned(f.peptide().clone().into())),
            IdentifiedPeptidoformData::PLink(PLinkData { peptidoform, .. }) => {
                Some(Cow::Owned(peptidoform.clone().into()))
            }
            IdentifiedPeptidoformData::NovoB(NovoBData {
                score_forward,
                score_reverse,
                peptide_forward,
                peptide_reverse,
                ..
            }) => if score_forward >= score_reverse {
                peptide_forward.as_ref()
            } else {
                peptide_reverse.as_ref()
            }
            .map(|p| Cow::Owned(p.clone().into())),
        }
    }

    /// Get the format and version for this peptidoform
    fn format(&self) -> KnownFileFormat {
        match &self.metadata {
            IdentifiedPeptidoformData::BasicCSV(BasicCSVData { version, .. }) => {
                KnownFileFormat::BasicCSV(*version)
            }
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                version,
                ..
            }) => KnownFileFormat::SpectrumSequenceList(*version),
            IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { version, .. }) => {
                KnownFileFormat::DeepNovoFamily(*version)
            }
            IdentifiedPeptidoformData::Fasta(_) => KnownFileFormat::Fasta,
            IdentifiedPeptidoformData::InstaNovo(InstaNovoData { version, .. }) => {
                KnownFileFormat::InstaNovo(*version)
            }
            IdentifiedPeptidoformData::MaxQuant(MaxQuantData { version, .. }) => {
                KnownFileFormat::MaxQuant(*version)
            }
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { version, .. }) => {
                KnownFileFormat::MSFragger(*version)
            }
            IdentifiedPeptidoformData::MZTab(_) => KnownFileFormat::MZTab,
            IdentifiedPeptidoformData::NovoB(NovoBData { version, .. }) => {
                KnownFileFormat::NovoB(*version)
            }
            IdentifiedPeptidoformData::Novor(NovorData { version, .. }) => {
                KnownFileFormat::Novor(*version)
            }
            IdentifiedPeptidoformData::Opair(OpairData { version, .. }) => {
                KnownFileFormat::Opair(*version)
            }
            IdentifiedPeptidoformData::Peaks(PeaksData { version, .. }) => {
                KnownFileFormat::Peaks(*version)
            }
            IdentifiedPeptidoformData::PepNet(PepNetData { version, .. }) => {
                KnownFileFormat::PepNet(*version)
            }
            IdentifiedPeptidoformData::PLGS(PLGSData { version, .. }) => {
                KnownFileFormat::PLGS(*version)
            }
            IdentifiedPeptidoformData::PLink(PLinkData { version, .. }) => {
                KnownFileFormat::PLink(*version)
            }
            IdentifiedPeptidoformData::PowerNovo(PowerNovoData { version, .. }) => {
                KnownFileFormat::PowerNovo(*version)
            }
            IdentifiedPeptidoformData::Sage(SageData { version, .. }) => {
                KnownFileFormat::Sage(*version)
            }
        }
    }

    /// Get the identifier
    fn id(&self) -> String {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData {
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
            IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => {
                scan.iter().join(";")
            }
            IdentifiedPeptidoformData::Novor(NovorData {
                id, scan_number, ..
            }) => id.unwrap_or(*scan_number).to_string(),
            IdentifiedPeptidoformData::Opair(OpairData {
                scan_number: scan, ..
            })
            | IdentifiedPeptidoformData::NovoB(NovoBData { scan, .. })
            | IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                scan,
                ..
            })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData {
                scan_number: scan, ..
            })
            | IdentifiedPeptidoformData::BasicCSV(BasicCSVData {
                scan_index: scan, ..
            }) => scan.to_string(),
            IdentifiedPeptidoformData::Sage(SageData { id, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { id, .. }) => id.to_string(),
            IdentifiedPeptidoformData::Fasta(f) => (*f.identifier().accession()).to_string(),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { scan, .. }) => scan.to_string(),
            IdentifiedPeptidoformData::PLink(PLinkData { order, .. }) => order.to_string(),
            IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
                id, scan_number, ..
            }) => id.map_or_else(|| scan_number.iter().join(";"), |id| id.to_string()),
            IdentifiedPeptidoformData::PowerNovo(PowerNovoData { scan, .. }) => {
                scan.as_ref().map_or("-".to_string(), ToString::to_string)
            }
            IdentifiedPeptidoformData::PepNet(_) => "-".to_string(),
            IdentifiedPeptidoformData::PLGS(PLGSData {
                peptide_component_id,
                ..
            }) => peptide_component_id.to_string(),
        }
    }

    fn confidence(&self) -> Option<f64> {
        self.score
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| Cow::Borrowed(lc.as_slice()))
    }

    fn original_confidence(&self) -> Option<f64> {
        todo!()
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        match &self.metadata {
            IdentifiedPeptidoformData::InstaNovo(InstaNovoData {
                local_confidence, ..
            })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData {
                local_confidence, ..
            })
            | IdentifiedPeptidoformData::PepNet(PepNetData {
                local_confidence, ..
            }) => Some(local_confidence),

            IdentifiedPeptidoformData::Peaks(PeaksData {
                local_confidence, ..
            })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData {
                local_confidence,
                ..
            })
            | IdentifiedPeptidoformData::Novor(NovorData {
                local_confidence, ..
            })
            | IdentifiedPeptidoformData::MZTab(MZTabData {
                local_confidence, ..
            }) => local_confidence.as_deref(),
            _ => None,
        }
    }

    /// The charge of the precursor, if known
    fn charge(&self) -> Option<Charge> {
        match &self.metadata {
            IdentifiedPeptidoformData::Novor(NovorData { z, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { z, .. })
            | IdentifiedPeptidoformData::Sage(SageData { z, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { z, .. })
            | IdentifiedPeptidoformData::NovoB(NovoBData { z, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData { precursor_z: z, .. })
            | IdentifiedPeptidoformData::PLink(PLinkData { z, .. })
            | IdentifiedPeptidoformData::MSFragger(MSFraggerData { z, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { z, .. })
            | IdentifiedPeptidoformData::MZTab(MZTabData { z, .. })
            | IdentifiedPeptidoformData::BasicCSV(BasicCSVData { z, .. })
            | IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                z, ..
            }) => Some(*z),
            IdentifiedPeptidoformData::Peaks(PeaksData { z, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { z, .. }) => *z,
            IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::PepNet(_) => None,
        }
    }

    /// Which fragmentation mode was used, if known
    fn mode(&self) -> Option<&str> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { mode, .. })
            | IdentifiedPeptidoformData::BasicCSV(BasicCSVData { mode, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { mode, .. }) => mode.as_deref(),
            _ => None,
        }
    }

    /// The retention time, if known
    fn retention_time(&self) -> Option<Time> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { rt, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { rt, .. })
            | IdentifiedPeptidoformData::Sage(SageData { rt, .. })
            | IdentifiedPeptidoformData::MSFragger(MSFraggerData { rt, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData {
                precursor_rt: rt, ..
            }) => Some(*rt),
            IdentifiedPeptidoformData::MaxQuant(MaxQuantData { rt, .. })
            | IdentifiedPeptidoformData::Novor(NovorData { rt, .. })
            | IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                rt,
                ..
            })
            | IdentifiedPeptidoformData::MZTab(MZTabData { rt, .. }) => *rt,
            IdentifiedPeptidoformData::DeepNovoFamily(_)
            | IdentifiedPeptidoformData::InstaNovo(_)
            | IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::NovoB(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::PLink(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
        }
    }

    /// The scans per rawfile that are at the basis for this identified peptide, if the rawfile is unknown there will be one
    fn scans(&self) -> SpectrumIds {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData {
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
            IdentifiedPeptidoformData::Novor(NovorData { scan_number, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Number(*scan_number)])
            }
            IdentifiedPeptidoformData::NovoB(NovoBData { scan, .. }) => {
                SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)])
            }
            IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { scan, .. }) => {
                SpectrumIds::FileNotKnown(
                    scan.iter()
                        .flat_map(|s| s.scans.clone())
                        .map(SpectrumId::Number)
                        .collect(),
                )
            }

            IdentifiedPeptidoformData::Opair(OpairData {
                raw_file,
                scan_number,
                ..
            })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData {
                raw_file,
                scan_number,
                ..
            }) => SpectrumIds::FileKnown(vec![(
                raw_file.clone(),
                vec![SpectrumId::Number(*scan_number)],
            )]),
            IdentifiedPeptidoformData::SpectrumSequenceList(SpectrumSequenceListData {
                raw_file,
                scan: scan_index,
                ..
            })
            | IdentifiedPeptidoformData::BasicCSV(BasicCSVData {
                raw_file,
                scan_index,
                ..
            }) => SpectrumIds::FileKnown(vec![(
                raw_file.clone(),
                vec![SpectrumId::Index(*scan_index)],
            )]),

            IdentifiedPeptidoformData::PowerNovo(PowerNovoData { raw_file, scan, .. }) => {
                scan.as_ref().map_or(SpectrumIds::None, |scan| {
                    raw_file.clone().map_or_else(
                        || SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)]),
                        |raw_file| {
                            SpectrumIds::FileKnown(vec![(raw_file, vec![SpectrumId::Index(*scan)])])
                        },
                    )
                })
            }

            IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
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
            IdentifiedPeptidoformData::MZTab(MZTabData { spectra_ref, .. }) => spectra_ref.clone(),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { raw_file, scan, .. }) => {
                raw_file.clone().map_or_else(
                    || SpectrumIds::FileNotKnown(vec![scan.clone()]),
                    |raw_file| SpectrumIds::FileKnown(vec![(raw_file, vec![scan.clone()])]),
                )
            }
            IdentifiedPeptidoformData::PLink(PLinkData {
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
            IdentifiedPeptidoformData::Sage(SageData { raw_file, scan, .. }) => {
                SpectrumIds::FileKnown(vec![(raw_file.clone(), vec![scan.clone()])])
            }
            IdentifiedPeptidoformData::PLGS(PLGSData {
                precursor_lift_off_rt,
                precursor_touch_down_rt,
                ..
            }) => SpectrumIds::FileNotKnown(vec![SpectrumId::RetentionTime(
                OrderedTime::from(*precursor_lift_off_rt)
                    ..=OrderedTime::from(*precursor_touch_down_rt),
            )]),
            IdentifiedPeptidoformData::Fasta(_) | IdentifiedPeptidoformData::PepNet(_) => {
                SpectrumIds::None
            }
        }
    }

    /// Get the mz as experimentally determined
    fn experimental_mz(&self) -> Option<MassOverCharge> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { mz, .. })
            | IdentifiedPeptidoformData::Novor(NovorData { mz, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { mz, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { mz, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData {
                precursor_mz: mz, ..
            }) => Some(*mz),
            IdentifiedPeptidoformData::MZTab(MZTabData { mz, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData { mz, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { mz, .. })
            | IdentifiedPeptidoformData::MSFragger(MSFraggerData { mz, .. }) => *mz,
            IdentifiedPeptidoformData::Sage(SageData { mass, z, .. })
            | IdentifiedPeptidoformData::NovoB(NovoBData { mass, z, .. })
            | IdentifiedPeptidoformData::PLink(PLinkData { mass, z, .. }) => Some(
                MassOverCharge::new::<crate::system::mz>(mass.value / (z.value as f64)),
            ),
            IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::SpectrumSequenceList(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
        }
    }

    /// Get the mass as experimentally determined
    fn experimental_mass(&self) -> Option<crate::system::Mass> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { mass, mz, z, .. }) => {
                mass.map_or(z.map_or(None, |z| Some(*mz * z.to_float())), Some)
            }
            IdentifiedPeptidoformData::Novor(NovorData { mass, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { mass, .. })
            | IdentifiedPeptidoformData::PLGS(PLGSData {
                precursor_mass: mass,
                ..
            })
            | IdentifiedPeptidoformData::NovoB(NovoBData { mass, .. })
            | IdentifiedPeptidoformData::PLink(PLinkData { mass, .. })
            | IdentifiedPeptidoformData::MSFragger(MSFraggerData { mass, .. })
            | IdentifiedPeptidoformData::Sage(SageData { mass, .. }) => Some(*mass),
            IdentifiedPeptidoformData::MaxQuant(MaxQuantData { mass, .. }) => *mass,
            IdentifiedPeptidoformData::MZTab(MZTabData { mz, z, .. }) => {
                mz.map(|mz| mz * z.to_float())
            }
            IdentifiedPeptidoformData::InstaNovo(InstaNovoData { mz, z, .. }) => {
                Some(*mz * z.to_float())
            }
            IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { mz, z, .. }) => {
                mz.and_then(|mz| z.map(|z| (mz, z)).map(|(mz, z)| mz * z.to_float()))
            }
            IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::SpectrumSequenceList(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
        }
    }

    /// Get the protein name if this was database matched data
    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData {
                protein_accession, ..
            }) => protein_accession.clone(),
            IdentifiedPeptidoformData::Opair(OpairData { protein_name, .. }) => {
                Some(protein_name.clone())
            }
            IdentifiedPeptidoformData::PLGS(PLGSData {
                protein_description,
                ..
            }) => Some(protein_description.clone()),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData { protein, .. }) => {
                Some(protein.clone())
            }
            IdentifiedPeptidoformData::MZTab(MZTabData { accession, .. }) => accession
                .as_ref()
                .map(|a| FastaIdentifier::Undefined(a.clone())),
            IdentifiedPeptidoformData::NovoB(_)
            | IdentifiedPeptidoformData::MaxQuant(_)
            | IdentifiedPeptidoformData::Sage(_)
            | IdentifiedPeptidoformData::PLink(_)
            | IdentifiedPeptidoformData::Novor(_)
            | IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::DeepNovoFamily(_)
            | IdentifiedPeptidoformData::InstaNovo(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::SpectrumSequenceList(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
        }
    }

    /// Get the protein id if this was database matched data
    fn protein_id(&self) -> Option<usize> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { protein_id, .. }) => *protein_id,
            IdentifiedPeptidoformData::Novor(NovorData { protein, .. }) => *protein,
            IdentifiedPeptidoformData::PLGS(PLGSData { protein_id, .. }) => Some(*protein_id),
            IdentifiedPeptidoformData::MZTab(_)
            | IdentifiedPeptidoformData::MaxQuant(_)
            | IdentifiedPeptidoformData::Sage(_)
            | IdentifiedPeptidoformData::PLink(_)
            | IdentifiedPeptidoformData::MSFragger(_)
            | IdentifiedPeptidoformData::NovoB(_)
            | IdentifiedPeptidoformData::Opair(_)
            | IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::DeepNovoFamily(_)
            | IdentifiedPeptidoformData::SpectrumSequenceList(_)
            | IdentifiedPeptidoformData::InstaNovo(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
        }
    }

    /// Get the protein location if this was database matched data
    fn protein_location(&self) -> Option<Range<usize>> {
        match &self.metadata {
            IdentifiedPeptidoformData::Peaks(PeaksData { start, end, .. }) => {
                start.and_then(|s| end.map(|e| s..e))
            }
            IdentifiedPeptidoformData::Novor(NovorData {
                protein_start,
                peptide,
                ..
            }) => protein_start.map(|s| s..s + peptide.len()),
            IdentifiedPeptidoformData::Opair(OpairData {
                protein_location, ..
            }) => Some(protein_location.clone()),
            IdentifiedPeptidoformData::PLGS(PLGSData {
                peptide_start,
                peptide,
                ..
            }) => Some(*peptide_start..*peptide_start + peptide.len()),
            IdentifiedPeptidoformData::MSFragger(MSFraggerData {
                protein_start,
                protein_end,
                ..
            }) => protein_start.and_then(|start| protein_end.map(|end| start..end)),
            IdentifiedPeptidoformData::MZTab(MZTabData { start, end, .. }) => {
                start.and_then(|s| end.map(|e| s..e))
            }
            IdentifiedPeptidoformData::InstaNovo(_)
            | IdentifiedPeptidoformData::DeepNovoFamily(_)
            | IdentifiedPeptidoformData::MaxQuant(_)
            | IdentifiedPeptidoformData::Sage(_)
            | IdentifiedPeptidoformData::PLink(_)
            | IdentifiedPeptidoformData::NovoB(_)
            | IdentifiedPeptidoformData::Fasta(_)
            | IdentifiedPeptidoformData::PowerNovo(_)
            | IdentifiedPeptidoformData::SpectrumSequenceList(_)
            | IdentifiedPeptidoformData::PepNet(_)
            | IdentifiedPeptidoformData::BasicCSV(_) => None,
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
