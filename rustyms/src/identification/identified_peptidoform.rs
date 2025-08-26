use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    identification::*,
    sequence::{
        AtLeast, CompoundPeptidoformIon, HasPeptidoformImpl, Linear, Linked, Peptidoform,
        SemiAmbiguous, SimpleLinear, UnAmbiguous,
    },
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

/// A peptide that is identified by a _de novo_ or database matching program
#[derive(Clone, Debug, Deserialize, PartialEq, Serialize)]
pub struct IdentifiedPeptidoform<Complexity, PeptidoformAvailability> {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The local confidence, if available, in range -1.0..=1.0
    pub local_confidence: Option<Vec<f64>>,
    /// The full metadata of this peptide
    pub data: IdentifiedPeptidoformData,
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
    /// InstaNovo metadata
    InstaNovo(InstaNovoData),
    /// MaxQuant metadata
    MaxQuant(MaxQuantData),
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
    /// pi-HelixNovo metadata
    PiHelixNovo(PiHelixNovoData),
    /// pi-PrimeNovo metadata
    PiPrimeNovo(PiPrimeNovoData),
    /// PLGS metadata
    PLGS(PLGSData),
    /// pLink metadata
    PLink(PLinkData),
    /// PowerNovo metadata
    PowerNovo(PowerNovoData),
    /// Proteoscape metadata
    Proteoscape(ProteoscapeData),
    /// Sage metadata
    Sage(SageData),
    /// SpectrumSequenceList metadata
    SpectrumSequenceList(SpectrumSequenceListData),
}

impl<PeptidoformAvailability> IdentifiedPeptidoform<Linear, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a linear peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<Linear>> {
        match &self.data {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PiHelixNovo(PiHelixNovoData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Proteoscape(ProteoscapeData {
                peptide: (_, peptide, _),
                ..
            })
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
            | IdentifiedPeptidoformData::PiPrimeNovo(PiPrimeNovoData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            IdentifiedPeptidoformData::MZTab(MZTabData { peptidoform, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().map(AsRef::as_ref),
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
        match &self.data {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::PiHelixNovo(PiHelixNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Proteoscape(ProteoscapeData {
                peptide: (_, peptide, _),
                ..
            })
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
            | IdentifiedPeptidoformData::PiPrimeNovo(PiPrimeNovoData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            IdentifiedPeptidoformData::MZTab(MZTabData { peptidoform, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref(),
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
        match &self.data {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PiHelixNovo(PiHelixNovoData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Proteoscape(ProteoscapeData {
                peptide: (_, peptide, _),
                ..
            })
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
            | IdentifiedPeptidoformData::PiPrimeNovo(PiPrimeNovoData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref()
            }
            IdentifiedPeptidoformData::MZTab(MZTabData { peptidoform, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().and_then(|p| p.as_semi_ambiguous()),
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
        match &self.data {
            IdentifiedPeptidoformData::Novor(NovorData { peptide, .. })
            | IdentifiedPeptidoformData::InstaNovo(InstaNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Opair(OpairData { peptide, .. })
            | IdentifiedPeptidoformData::PiHelixNovo(PiHelixNovoData { peptide, .. })
            | IdentifiedPeptidoformData::PepNet(PepNetData { peptide, .. })
            | IdentifiedPeptidoformData::PowerNovo(PowerNovoData { peptide, .. })
            | IdentifiedPeptidoformData::Proteoscape(ProteoscapeData {
                peptide: (_, peptide, _),
                ..
            })
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
            | IdentifiedPeptidoformData::PiPrimeNovo(PiPrimeNovoData { peptide, .. })
            | IdentifiedPeptidoformData::DeepNovoFamily(DeepNovoFamilyData { peptide, .. }) => {
                peptide.as_ref().and_then(|p| p.as_unambiguous())
            }
            IdentifiedPeptidoformData::MZTab(MZTabData { peptidoform, .. })
            | IdentifiedPeptidoformData::MaxQuant(MaxQuantData {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().and_then(|p| p.as_unambiguous()),
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
    fn check<T>(
        self,
        f: impl Fn(&Peptidoform<Linked>) -> bool,
    ) -> Option<IdentifiedPeptidoform<T, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| p.singular_peptidoform_ref().is_some_and(f))
            .then(|| self.mark())
    }

    /// Check if this identified peptidoform is linear and contains a peptide
    pub fn into_linear(self) -> Option<IdentifiedPeptidoform<Linear, PeptidoformPresent>> {
        self.check(Peptidoform::is_linear)
    }

    /// Check if this identified peptidoform is simple linear and contains a peptide
    pub fn into_simple_linear(
        self,
    ) -> Option<IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent>> {
        self.check(Peptidoform::is_simple_linear)
    }

    /// Check if this identified peptidoform is semi ambiguous and contains a peptide
    pub fn into_semi_ambiguous(
        self,
    ) -> Option<IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>> {
        self.check(Peptidoform::is_semi_ambiguous)
    }

    /// Check if this identified peptidoform is unambiguous and contains a peptide
    pub fn into_unambiguous(
        self,
    ) -> Option<IdentifiedPeptidoform<UnAmbiguous, PeptidoformPresent>> {
        self.check(Peptidoform::is_unambiguous)
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
            data: self.data,
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

/// Implement the [`MetaData`] trait for [`IdentifiedPeptidoform`] without having to type everything
/// out. Needs some macro fudging to allow for the proper syntax of just specifying a list of formats
/// and a list of functions to implement.
macro_rules! impl_metadata {
    (formats: $format:tt; functions: {$(fn $function:ident(&self) -> $t:ty);+;}) => {
        impl<Complexity, PeptidoformAvailability> MetaData for IdentifiedPeptidoform<Complexity, PeptidoformAvailability> {
            /// Reuse the cached normalised confidence
            fn confidence(&self) -> Option<f64> {
                self.score
            }

            /// Reuse the cached normalised local confidence
            fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
                self.local_confidence
                    .as_ref()
                    .map(|lc| Cow::Borrowed(lc.as_slice()))
            }

            $(impl_metadata!(inner: formats: $format; function: $function -> $t);)+
        }
    };
    (inner: formats: {$($format:ident),*}; function: $function:ident -> $t:ty) => {
        fn $function(&self) -> $t {
            match &self.data {
                $(IdentifiedPeptidoformData::$format(d) => d.$function()),*
            }
        }
    };
}

impl_metadata!(
    formats: {BasicCSV,DeepNovoFamily,Fasta,MaxQuant,InstaNovo,MZTab,NovoB,Novor,Opair,Peaks,PepNet,PiHelixNovo,PiPrimeNovo,PLGS,PLink,PowerNovo,Proteoscape,Sage,MSFragger,SpectrumSequenceList};
    functions: {
        fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>>;
        fn format(&self) -> KnownFileFormat;
        fn id(&self) -> String;
        fn original_confidence(&self) -> Option<f64>;
        fn original_local_confidence(&self) -> Option<&[f64]>;
        fn charge(&self) -> Option<Charge>;
        fn mode(&self) -> Option<&str>;
        fn retention_time(&self) -> Option<Time>;
        fn scans(&self) -> SpectrumIds;
        fn experimental_mz(&self) -> Option<MassOverCharge>;
        fn experimental_mass(&self) -> Option<Mass>;
        fn ppm_error(&self) -> Option<Ratio>;
        fn mass_error(&self) -> Option<Mass>;
        fn protein_name(&self) -> Option<FastaIdentifier<String>>;
        fn protein_id(&self) -> Option<usize>;
        fn protein_location(&self) -> Option<Range<u16>>;
        fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence);
        fn database(&self) -> Option<(&str, Option<&str>)>;
    }
);
