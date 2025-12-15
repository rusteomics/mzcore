use std::{borrow::Cow, marker::PhantomData, ops::Range};

#[cfg(feature = "mzannotate")]
use mzannotate::prelude::AnnotatedSpectrum;
#[cfg(not(feature = "mzannotate"))]
use serde::{Deserialize, Serialize};

use crate::*;
use mzcore::{
    sequence::{
        AtLeast, CompoundPeptidoformIon, FlankingSequence, HasPeptidoformImpl, Linear, Linked,
        Peptidoform, SemiAmbiguous, SimpleLinear, UnAmbiguous,
    },
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

/// A peptidoform that is identified by a _de novo_ or database matching program as matching to a spectrum
#[cfg_attr(not(feature = "mzannotate"), derive(Deserialize, Serialize))]
#[derive(Clone, Debug)]
pub struct PSM<Complexity, PeptidoformAvailability> {
    /// The score -1.0..=1.0 if a score was available in the original format
    pub score: Option<f64>,
    /// The local confidence, if available, in range -1.0..=1.0
    pub local_confidence: Option<Vec<f64>>,
    /// The full metadata of this peptide
    pub data: PSMData,
    /// The marker for the complexity, Linked means full [`CompoundPeptidoformIon`] anything below means [`Peptidoform`], see [Complexity](crate::sequence::Complexity)
    pub(super) complexity_marker: PhantomData<Complexity>,
    /// The marker for availability of the peptidoform, see [`PeptidoformAvailability`]
    pub(super) peptidoform_availability_marker: PhantomData<PeptidoformAvailability>,
}

/// The definition of all special metadata for all types of PSM that can be read
#[cfg_attr(not(feature = "mzannotate"), derive(Deserialize, Serialize))]
#[derive(Clone, Debug)]
#[expect(clippy::upper_case_acronyms)]
pub enum PSMData {
    /// A basic CSV format
    BasicCSV(BasicCSVPSM),
    /// DeepNovo/PointNovo/PGPointNovo metadata
    DeepNovoFamily(DeepNovoFamilyPSM),
    /// Fasta metadata
    Fasta(FastaData),
    /// InstaNovo metadata
    InstaNovo(InstaNovoPSM),
    /// MaxQuant metadata
    MaxQuant(MaxQuantPSM),
    /// MaxQuant metadata
    MetaMorpheus(MetaMorpheusPSM),
    /// MSFragger metadata
    MSFragger(MSFraggerPSM),
    /// mzTab metadata
    MzTab(MzTabPSM),
    /// NovoB metadata
    NovoB(NovoBPSM),
    /// Novor metadata
    Novor(NovorPSM),
    /// OPair metadata
    Opair(OpairPSM),
    /// Peaks metadata
    Peaks(PeaksPSM),
    /// PepNet metadata
    PepNet(PepNetPSM),
    /// pi-HelixNovo metadata
    PiHelixNovo(PiHelixNovoPSM),
    /// pi-PrimeNovo metadata
    PiPrimeNovo(PiPrimeNovoPSM),
    /// PLGS metadata
    PLGS(PLGSPSM),
    /// pLink metadata
    PLink(PLinkPSM),
    /// PowerNovo metadata
    PowerNovo(PowerNovoPSM),
    /// Proteoscape metadata
    Proteoscape(ProteoscapePSM),
    /// pUniFind metadata
    PUniFind(PUniFindPSM),
    /// Sage metadata
    Sage(SagePSM),
    /// SpectrumSequenceList metadata
    SpectrumSequenceList(SpectrumSequenceListPSM),
    /// An mzSpecLib spectrum (only available if the feature `mzannotate` is turned on)
    #[cfg(feature = "mzannotate")]
    AnnotatedSpectrum(AnnotatedSpectrum),
}

impl<PeptidoformAvailability> PSM<Linear, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a linear peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<Linear>> {
        match &self.data {
            PSMData::Novor(NovorPSM { peptide, .. })
            | PSMData::InstaNovo(InstaNovoPSM { peptide, .. })
            | PSMData::Opair(OpairPSM { peptide, .. })
            | PSMData::PiHelixNovo(PiHelixNovoPSM { peptide, .. })
            | PSMData::PepNet(PepNetPSM { peptide, .. })
            | PSMData::PowerNovo(PowerNovoPSM { peptide, .. })
            | PSMData::PUniFind(PUniFindPSM {
                peptidoform: peptide,
                ..
            })
            | PSMData::Proteoscape(ProteoscapePSM {
                peptide: (_, peptide, _),
                ..
            })
            | PSMData::Sage(SagePSM { peptide, .. }) => Some(peptide.as_ref()),
            PSMData::MSFragger(MSFraggerPSM { peptide, .. })
            | PSMData::PLGS(PLGSPSM { peptide, .. }) => Some(peptide.as_ref()),
            PSMData::Peaks(PeaksPSM { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(peptide.1[0].as_ref())
                } else {
                    None
                }
            }
            PSMData::SpectrumSequenceList(SpectrumSequenceListPSM { peptide, .. })
            | PSMData::PiPrimeNovo(PiPrimeNovoPSM { peptide, .. })
            | PSMData::DeepNovoFamily(DeepNovoFamilyPSM { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            PSMData::MzTab(MzTabPSM { peptidoform, .. })
            | PSMData::MaxQuant(MaxQuantPSM {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().map(AsRef::as_ref),
            PSMData::Fasta(f) => Some(f.peptide().as_ref()),
            PSMData::NovoB(NovoBPSM {
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
            PSMData::MetaMorpheus(MetaMorpheusPSM { peptide, .. })
            | PSMData::BasicCSV(BasicCSVPSM {
                sequence: peptide, ..
            }) => peptide
                .singular_peptidoform_ref()
                .and_then(|p| p.as_linear()),
            PSMData::PLink(PLinkPSM { peptidoform, .. }) => {
                peptidoform.singular_ref().and_then(|p| p.as_linear())
            }
            #[cfg(feature = "mzannotate")]
            PSMData::AnnotatedSpectrum(spectrum) => {
                use itertools::Itertools;
                use mzannotate::mzspeclib::AnalyteTarget;

                let ion = spectrum
                    .analytes
                    .iter()
                    .filter_map(|a| match &a.target {
                        AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                        _ => None,
                    })
                    .exactly_one()
                    .ok()?;
                ion.singular_ref().and_then(|p| p.as_linear())
            }
        }
    }
}

impl<PeptidoformAvailability> PSM<SimpleLinear, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a simple linear peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<SimpleLinear>> {
        match &self.data {
            PSMData::Novor(NovorPSM { peptide, .. })
            | PSMData::InstaNovo(InstaNovoPSM { peptide, .. })
            | PSMData::PiHelixNovo(PiHelixNovoPSM { peptide, .. })
            | PSMData::Opair(OpairPSM { peptide, .. })
            | PSMData::PepNet(PepNetPSM { peptide, .. })
            | PSMData::PowerNovo(PowerNovoPSM { peptide, .. })
            | PSMData::PUniFind(PUniFindPSM {
                peptidoform: peptide,
                ..
            })
            | PSMData::Proteoscape(ProteoscapePSM {
                peptide: (_, peptide, _),
                ..
            })
            | PSMData::Sage(SagePSM { peptide, .. }) => Some(peptide.as_ref()),
            PSMData::MSFragger(MSFraggerPSM { peptide, .. })
            | PSMData::PLGS(PLGSPSM { peptide, .. }) => Some(peptide),
            PSMData::Peaks(PeaksPSM { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(peptide.1[0].as_ref())
                } else {
                    None
                }
            }
            PSMData::SpectrumSequenceList(SpectrumSequenceListPSM { peptide, .. })
            | PSMData::PiPrimeNovo(PiPrimeNovoPSM { peptide, .. })
            | PSMData::DeepNovoFamily(DeepNovoFamilyPSM { peptide, .. }) => {
                peptide.as_ref().map(AsRef::as_ref)
            }
            PSMData::MzTab(MzTabPSM { peptidoform, .. })
            | PSMData::MaxQuant(MaxQuantPSM {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref(),
            PSMData::Fasta(f) => Some(f.peptide().as_ref()),
            PSMData::NovoB(NovoBPSM {
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
            PSMData::MetaMorpheus(MetaMorpheusPSM { peptide, .. })
            | PSMData::BasicCSV(BasicCSVPSM {
                sequence: peptide, ..
            }) => peptide
                .singular_peptidoform_ref()
                .and_then(|p| p.as_simple_linear()),
            PSMData::PLink(PLinkPSM { peptidoform, .. }) => peptidoform
                .singular_ref()
                .and_then(|p| p.as_simple_linear()),
            #[cfg(feature = "mzannotate")]
            PSMData::AnnotatedSpectrum(spectrum) => {
                use itertools::Itertools;
                use mzannotate::mzspeclib::AnalyteTarget;

                let ion = spectrum
                    .analytes
                    .iter()
                    .filter_map(|a| match &a.target {
                        AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                        _ => None,
                    })
                    .exactly_one()
                    .ok()?;
                ion.singular_ref().and_then(|p| p.as_simple_linear())
            }
        }
    }
}

impl<PeptidoformAvailability> PSM<SemiAmbiguous, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as a semi ambiguous peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<SemiAmbiguous>> {
        match &self.data {
            PSMData::Novor(NovorPSM { peptide, .. })
            | PSMData::InstaNovo(InstaNovoPSM { peptide, .. })
            | PSMData::Opair(OpairPSM { peptide, .. })
            | PSMData::PiHelixNovo(PiHelixNovoPSM { peptide, .. })
            | PSMData::PepNet(PepNetPSM { peptide, .. })
            | PSMData::PowerNovo(PowerNovoPSM { peptide, .. })
            | PSMData::PUniFind(PUniFindPSM {
                peptidoform: peptide,
                ..
            })
            | PSMData::Proteoscape(ProteoscapePSM {
                peptide: (_, peptide, _),
                ..
            })
            | PSMData::Sage(SagePSM { peptide, .. }) => Some(peptide),
            PSMData::Peaks(PeaksPSM { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    Some(&peptide.1[0])
                } else {
                    None
                }
            }
            PSMData::SpectrumSequenceList(SpectrumSequenceListPSM { peptide, .. })
            | PSMData::PiPrimeNovo(PiPrimeNovoPSM { peptide, .. })
            | PSMData::DeepNovoFamily(DeepNovoFamilyPSM { peptide, .. }) => peptide.as_ref(),
            PSMData::MzTab(MzTabPSM { peptidoform, .. })
            | PSMData::MaxQuant(MaxQuantPSM {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().and_then(|p| p.as_semi_ambiguous()),
            PSMData::Fasta(f) => Some(f.peptide()),
            PSMData::NovoB(NovoBPSM {
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
            PSMData::MSFragger(MSFraggerPSM { peptide, .. })
            | PSMData::PLGS(PLGSPSM { peptide, .. }) => peptide.as_semi_ambiguous(),
            PSMData::MetaMorpheus(MetaMorpheusPSM { peptide, .. })
            | PSMData::BasicCSV(BasicCSVPSM {
                sequence: peptide, ..
            }) => peptide
                .singular_peptidoform_ref()
                .and_then(|p| p.as_semi_ambiguous()),
            PSMData::PLink(PLinkPSM { peptidoform, .. }) => peptidoform
                .singular_ref()
                .and_then(|p| p.as_semi_ambiguous()),
            #[cfg(feature = "mzannotate")]
            PSMData::AnnotatedSpectrum(spectrum) => {
                use itertools::Itertools;

                let ion = spectrum
                    .analytes
                    .iter()
                    .filter_map(|a| match &a.target {
                        mzannotate::mzspeclib::AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                        _ => None,
                    })
                    .exactly_one()
                    .ok()?;
                ion.singular_ref().and_then(|p| p.as_semi_ambiguous())
            }
        }
    }
}

impl<PeptidoformAvailability> PSM<UnAmbiguous, PeptidoformAvailability> {
    /// If this peptidoform contains a peptidoform that is valid as an unambiguous peptidoform get a reference to the peptidoform.
    fn inner_peptidoform(&self) -> Option<&Peptidoform<UnAmbiguous>> {
        match &self.data {
            PSMData::Novor(NovorPSM { peptide, .. })
            | PSMData::InstaNovo(InstaNovoPSM { peptide, .. })
            | PSMData::Opair(OpairPSM { peptide, .. })
            | PSMData::PiHelixNovo(PiHelixNovoPSM { peptide, .. })
            | PSMData::PepNet(PepNetPSM { peptide, .. })
            | PSMData::PowerNovo(PowerNovoPSM { peptide, .. })
            | PSMData::PUniFind(PUniFindPSM {
                peptidoform: peptide,
                ..
            })
            | PSMData::Proteoscape(ProteoscapePSM {
                peptide: (_, peptide, _),
                ..
            })
            | PSMData::Sage(SagePSM { peptide, .. }) => peptide.as_unambiguous(),
            PSMData::Peaks(PeaksPSM { peptide, .. }) => {
                if peptide.1.len() == 1 {
                    peptide.1[0].as_unambiguous()
                } else {
                    None
                }
            }
            PSMData::SpectrumSequenceList(SpectrumSequenceListPSM { peptide, .. })
            | PSMData::PiPrimeNovo(PiPrimeNovoPSM { peptide, .. })
            | PSMData::DeepNovoFamily(DeepNovoFamilyPSM { peptide, .. }) => {
                peptide.as_ref().and_then(|p| p.as_unambiguous())
            }
            PSMData::MzTab(MzTabPSM { peptidoform, .. })
            | PSMData::MaxQuant(MaxQuantPSM {
                peptide: peptidoform,
                ..
            }) => peptidoform.as_ref().and_then(|p| p.as_unambiguous()),
            PSMData::Fasta(f) => f.peptide().as_unambiguous(),
            PSMData::NovoB(NovoBPSM {
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
            PSMData::MSFragger(MSFraggerPSM { peptide, .. })
            | PSMData::PLGS(PLGSPSM { peptide, .. }) => peptide.as_unambiguous(),
            PSMData::MetaMorpheus(MetaMorpheusPSM { peptide, .. })
            | PSMData::BasicCSV(BasicCSVPSM {
                sequence: peptide, ..
            }) => peptide
                .singular_peptidoform_ref()
                .and_then(|p| p.as_unambiguous()),
            PSMData::PLink(PLinkPSM { peptidoform, .. }) => {
                peptidoform.singular_ref().and_then(|p| p.as_unambiguous())
            }
            #[cfg(feature = "mzannotate")]
            PSMData::AnnotatedSpectrum(spectrum) => {
                use itertools::Itertools;

                let ion = spectrum
                    .analytes
                    .iter()
                    .filter_map(|a| match &a.target {
                        mzannotate::mzspeclib::AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                        _ => None,
                    })
                    .exactly_one()
                    .ok()?;
                ion.singular_ref().and_then(|p| p.as_unambiguous())
            }
        }
    }
}

impl<Complexity, PeptidoformAvailability> PSM<Complexity, PeptidoformAvailability> {
    fn check<T>(
        self,
        f: impl Fn(&Peptidoform<Linked>) -> bool,
    ) -> Option<PSM<T, PeptidoformPresent>> {
        self.compound_peptidoform_ion()
            .is_some_and(|p| p.singular_peptidoform_ref().is_some_and(f))
            .then(|| self.mark())
    }

    /// Check if this PSM is linear and contains a peptide
    pub fn into_linear(self) -> Option<PSM<Linear, PeptidoformPresent>> {
        self.check(Peptidoform::is_linear)
    }

    /// Check if this PSM is simple linear and contains a peptide
    pub fn into_simple_linear(self) -> Option<PSM<SimpleLinear, PeptidoformPresent>> {
        self.check(Peptidoform::is_simple_linear)
    }

    /// Check if this PSM is semi ambiguous and contains a peptide
    pub fn into_semi_ambiguous(self) -> Option<PSM<SemiAmbiguous, PeptidoformPresent>> {
        self.check(Peptidoform::is_semi_ambiguous)
    }

    /// Check if this PSM is unambiguous and contains a peptide
    pub fn into_unambiguous(self) -> Option<PSM<UnAmbiguous, PeptidoformPresent>> {
        self.check(Peptidoform::is_unambiguous)
    }
}

impl HasPeptidoformImpl for PSM<Linear, PeptidoformPresent> {
    type Complexity = Linear;
    fn peptidoform(&self) -> &Peptidoform<Linear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &PSM<Linear, PeptidoformPresent> {
    type Complexity = Linear;
    fn peptidoform(&self) -> &Peptidoform<Linear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for PSM<SimpleLinear, PeptidoformPresent> {
    type Complexity = SimpleLinear;
    fn peptidoform(&self) -> &Peptidoform<SimpleLinear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &PSM<SimpleLinear, PeptidoformPresent> {
    type Complexity = SimpleLinear;
    fn peptidoform(&self) -> &Peptidoform<SimpleLinear> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for PSM<SemiAmbiguous, PeptidoformPresent> {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<SemiAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &PSM<SemiAmbiguous, PeptidoformPresent> {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<SemiAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for PSM<UnAmbiguous, PeptidoformPresent> {
    type Complexity = UnAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<UnAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl HasPeptidoformImpl for &PSM<UnAmbiguous, PeptidoformPresent> {
    type Complexity = UnAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<UnAmbiguous> {
        self.inner_peptidoform()
            .expect("Identified peptidoform incorrectly marked as containing a peptidoform")
    }
}

impl PSM<Linear, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a linear peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<Linear>> {
        self.inner_peptidoform()
    }
}

impl PSM<SimpleLinear, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a simple linear peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<SimpleLinear>> {
        self.inner_peptidoform()
    }
}

impl PSM<SemiAmbiguous, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as a semi ambiguous peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<SemiAmbiguous>> {
        self.inner_peptidoform()
    }
}

impl PSM<UnAmbiguous, MaybePeptidoform> {
    /// If this peptidoform contains a peptidoform that is valid as an unambiguous peptidoform get a reference to the peptidoform.
    pub fn peptidoform(&self) -> Option<&Peptidoform<UnAmbiguous>> {
        self.inner_peptidoform()
    }
}

impl<Complexity, PeptidoformAvailability> PSM<Complexity, PeptidoformAvailability> {
    /// Mark this with the following complexity, be warned that the complexity level is not checked.
    fn mark<C, A>(self) -> PSM<C, A> {
        PSM {
            score: self.score,
            local_confidence: self.local_confidence,
            data: self.data,
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }

    /// Cast this PSM to a higher complexity level. This does not change the
    /// content of the peptidoform. It only allows to pass this as higher complexity if needed.
    pub fn cast<
        NewComplexity: AtLeast<Complexity>,
        NewAvailability: From<PeptidoformAvailability>,
    >(
        self,
    ) -> PSM<NewComplexity, NewAvailability> {
        self.mark()
    }
}

/// Implement the [`PSMMetaData`] trait for [`PSM`] without having to type everything
/// out. Needs some macro fudging to allow for the proper syntax of just specifying a list of formats
/// and a list of functions to implement.
macro_rules! impl_metadata {
    (formats: $format:tt; functions: {$($(#[cfg($cfg:expr)])?fn $function:ident(&self) -> $t:ty);+;}) => {
        impl<Complexity, PeptidoformAvailability> PSMMetaData for PSM<Complexity, PeptidoformAvailability> {
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

            type Protein = ProteinData;

            fn proteins(&self) -> Cow<'_, [Self::Protein]> {
                Cow::Owned(impl_metadata!(match: self, formats: $format; function: proteins;))
            }

            $(impl_metadata!(inner: formats: $format; function: $function -> $t);)+
        }
    };
    (inner: formats: {$($format:ident),*}; function: $(#[cfg($cfg:expr)])?$function:ident -> $t:ty) => {
        $(#[cfg($cfg)])?
        fn $function(&self) -> $t {
            match &self.data {
                $(PSMData::$format(d) => PSMMetaData::$function(d)),*
            }
        }
    };
    (match: $self:ident, formats: {$($format:ident),*}; function: $function:ident;) => {
         match &$self.data {
            $(PSMData::$format(d) => PSMMetaData::$function(d).iter().map(|p| p.clone().into()).collect::<Vec<_>>()),*
        }
    }
}
#[cfg(not(feature = "mzannotate"))]
impl_metadata!(
    formats: {BasicCSV,DeepNovoFamily,Fasta,MaxQuant,MetaMorpheus,InstaNovo,MzTab,NovoB,Novor,Opair,Peaks,PepNet,PiHelixNovo,PiPrimeNovo,PLGS,PLink,PowerNovo,Proteoscape,PUniFind,Sage,MSFragger,SpectrumSequenceList};
    functions: {
        fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>>;
        fn format(&self) -> KnownFileFormat;
        fn numerical_id(&self) -> Option<usize>;
        fn id(&self) -> String;
        fn search_engine(&self) -> Option<mzcv::Term>;
        fn original_confidence(&self) -> Option<(f64, mzcv::Term)>;
        fn original_local_confidence(&self) -> Option<&[f64]>;
        fn charge(&self) -> Option<Charge>;
        fn mode(&self) -> Option<Cow<'_, str>>;
        fn retention_time(&self) -> Option<Time>;
        fn scans(&self) -> SpectrumIds;
        fn experimental_mz(&self) -> Option<MassOverCharge>;
        fn experimental_mass(&self) -> Option<Mass>;
        fn ppm_error(&self) -> Option<Ratio>;
        fn mass_error(&self) -> Option<Mass>;
        fn protein_location(&self) -> Option<Range<u16>>;
        fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence);
        fn database(&self) -> Option<(&str, Option<&str>)>;
        fn unique(&self) -> Option<bool>;
        fn reliability(&self) -> Option<Reliability>;
        fn uri(&self) -> Option<String>;
    }
);

#[cfg(feature = "mzannotate")]
impl_metadata!(
    formats: {BasicCSV,DeepNovoFamily,Fasta,MaxQuant,MetaMorpheus,InstaNovo,MzTab,NovoB,Novor,Opair,Peaks,PepNet,PiHelixNovo,PiPrimeNovo,PLGS,PLink,PowerNovo,Proteoscape,PUniFind,Sage,MSFragger,SpectrumSequenceList,AnnotatedSpectrum};
    functions: {
        fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>>;
        fn format(&self) -> KnownFileFormat;
        fn numerical_id(&self) -> Option<usize>;
        fn id(&self) -> String;
        fn search_engine(&self) -> Option<mzcv::Term>;
        fn original_confidence(&self) -> Option<(f64, mzcv::Term)>;
        fn original_local_confidence(&self) -> Option<&[f64]>;
        fn charge(&self) -> Option<Charge>;
        fn mode(&self) -> Option<Cow<'_, str>>;
        fn fragmentation_model(&self) -> Option<mzannotate::annotation::model::BuiltInFragmentationModel>;
        fn retention_time(&self) -> Option<Time>;
        fn scans(&self) -> SpectrumIds;
        fn experimental_mz(&self) -> Option<MassOverCharge>;
        fn experimental_mass(&self) -> Option<Mass>;
        fn ppm_error(&self) -> Option<Ratio>;
        fn mass_error(&self) -> Option<Mass>;
        fn protein_location(&self) -> Option<Range<u16>>;
        fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence);
        fn database(&self) -> Option<(&str, Option<&str>)>;
        fn unique(&self) -> Option<bool>;
        fn reliability(&self) -> Option<Reliability>;
        fn uri(&self) -> Option<String>;
        fn annotated_spectrum(&self) -> Option<Cow<'_, AnnotatedSpectrum>>;
        fn has_annotated_spectrum(&self) -> bool;
    }
);

impl<C, P> mzcore::space::Space for PSM<C, P> {
    fn space(&self) -> mzcore::space::UsedSpace {
        self.score.space() + self.local_confidence.space() + self.data.space()
    }
}

impl mzcore::space::Space for PSMData {
    fn space(&self) -> mzcore::space::UsedSpace {
        match self {
            Self::BasicCSV(data) => data.space(),
            Self::DeepNovoFamily(data) => data.space(),
            Self::Fasta(data) => data.space(),
            Self::InstaNovo(data) => data.space(),
            Self::MaxQuant(data) => data.space(),
            Self::MetaMorpheus(data) => data.space(),
            Self::MSFragger(data) => data.space(),
            Self::MzTab(data) => data.space(),
            Self::NovoB(data) => data.space(),
            Self::Novor(data) => data.space(),
            Self::Opair(data) => data.space(),
            Self::Peaks(data) => data.space(),
            Self::PepNet(data) => data.space(),
            Self::PiHelixNovo(data) => data.space(),
            Self::PiPrimeNovo(data) => data.space(),
            Self::PLGS(data) => data.space(),
            Self::PLink(data) => data.space(),
            Self::PowerNovo(data) => data.space(),
            Self::Proteoscape(data) => data.space(),
            Self::PUniFind(data) => data.space(),
            Self::Sage(data) => data.space(),
            Self::SpectrumSequenceList(data) => data.space(),
            #[cfg(feature = "mzannotate")]
            Self::AnnotatedSpectrum(data) => data.space(),
        }
        .set_total::<Self>()
    }
}

/// The definition of all special metadata for all types of Proteins that can be read
#[derive(Clone, Debug, serde::Deserialize, serde::Serialize)]
pub enum ProteinData {
    /// No protein
    NoProtein(NoProtein),
    /// A single fasta identifier
    FastaId(FastaIdentifier<String>),
    /// A Fasta entry
    Fasta(FastaData),
    /// An mzTab protein
    MzTab(MzTabProtein),
    /// A protein from PLGS files
    PLGS(PLGSProtein),
    /// A protein from MSFragger files
    MSFragger(MSFraggerProtein),
    /// A protein from Opair files
    Opair(OpairProtein),
    /// A protein from MetaMorpheus files
    MetaMorpheus(MetaMorpheusProtein),
}

impl Default for ProteinData {
    fn default() -> Self {
        Self::NoProtein(NoProtein::default())
    }
}

impl From<NoProtein> for ProteinData {
    fn from(value: NoProtein) -> Self {
        Self::NoProtein(value)
    }
}

macro_rules! impl_protein_metadata {
    (formats: $format:tt; functions: {$($(#[cfg($cfg:expr)])?fn $function:ident(&self) -> $t:ty);+;}) => {
        impl ProteinMetaData for ProteinData {
            $(impl_protein_metadata!(inner: formats: $format; function: $function -> $t);)+
        }
    };
    (inner: formats: {$($format:ident),*}; function: $(#[cfg($cfg:expr)])?$function:ident -> $t:ty) => {
        $(#[cfg($cfg)])?
        fn $function(&self) -> $t {
            match &self {
                $(ProteinData::$format(d) => ProteinMetaData::$function(d)),*
            }
        }
    };
}

impl_protein_metadata!(
    formats: {NoProtein, FastaId, Fasta, MzTab, PLGS, MSFragger, Opair, MetaMorpheus};
    functions: {
        fn sequence(&self) -> Option<Cow<'_, Peptidoform<Linear>>>;
        fn numerical_id(&self) -> Option<usize>;
        fn id(&self) -> FastaIdentifier<&str>;
        fn description(&self) -> Option<&str>;
        fn species(&self) -> Option<mzcv::Curie>;
        fn species_name(&self) -> Option<&str>;
        fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)];
        fn ambiguity_members(&self) -> &[String];
        fn database(&self) -> Option<(&str, Option<&str>)>;
        fn modifications(&self) -> &[(Vec<(mzcore::sequence::SequencePosition, Option<f64>)>, mzcore::sequence::SimpleModification)];
        fn coverage(&self) -> Option<f64>;
        fn gene_ontology(&self) -> &[mzcv::Curie];
        fn reliability(&self) -> Option<Reliability>;
        fn uri(&self) -> Option<&str>;
    }
);

impl mzcore::space::Space for ProteinData {
    fn space(&self) -> mzcore::space::UsedSpace {
        match self {
            Self::NoProtein(data) => data.space(),
            Self::FastaId(data) => data.space(),
            Self::Fasta(data) => data.space(),
            Self::MzTab(data) => data.space(),
            Self::PLGS(data) => data.space(),
            Self::MSFragger(data) => data.space(),
            Self::Opair(data) => data.space(),
            Self::MetaMorpheus(data) => data.space(),
        }
        .set_total::<Self>()
    }
}
