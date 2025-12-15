use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, KnownFileFormat, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, PeptidoformPresent, SpectrumId, SpectrumIds,
    common_parser::Location,
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Proteoscape line",
    "This column is not a number but it is required to be a number in this Proteoscape format",
);

format_family!(
    Proteoscape,
    SemiAmbiguous, PeptidoformPresent, [&V2025B], b'\t', None;
    required {
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        /// Up to 3 leading amino acids (if present), the peptidoform itself, and up to 3 tailing amino acids (if present)
        peptide: (FlankingSequence, Peptidoform<SemiAmbiguous>, FlankingSequence), |location: Location, ontologies: &Ontologies| {
            location.clone().split_twice('.').ok_or_else(|| BoxedError::new(BasicKind::Error,"Invalid Proteoscape line", "The peptide columns should contain the previous amino acids, the peptide, and the following amino acids separated by dots.", location.context().to_owned())).and_then(|(before, peptide, after)| {
                let before = before.trim_start_matches("-");
                let after = after.trim_end_matches("-");
                Ok((
                    (!before.is_empty()).then(|| Peptidoform::sloppy_pro_forma(
                        before.full_line(),
                        before.location.clone(),
                        ontologies,
                        &SloppyParsingParameters::default(),
                    ).map_err(BoxedError::to_owned)).transpose()?
                    .map_or(FlankingSequence::Terminal, |s| FlankingSequence::Sequence(Box::new(s))),
                    Peptidoform::sloppy_pro_forma(
                        peptide.full_line(),
                        peptide.location.clone(),
                        ontologies,
                        &SloppyParsingParameters::default(),
                    ).map_err(BoxedError::to_owned)?,
                    (!after.is_empty()).then(|| Peptidoform::sloppy_pro_forma(
                        after.full_line(),
                        after.location.clone(),
                        ontologies,
                        &SloppyParsingParameters::default(),
                    ).map_err(BoxedError::to_owned)).transpose()?
                    .map_or(FlankingSequence::Terminal, |s| FlankingSequence::Sequence(Box::new(s))),
                ))
            })
        };
        xcorr_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::s>);
        corrected_o_over_k0: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        protein: FastaIdentifier<String>, |location: Location, _|  location.parse(NUMBER_ERROR);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        delta_cn_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        confidence_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        matched_ions: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        predicted_o_over_k0: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        tim_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        is_unique: bool, |location: Location, _| location.parse_with(|l| match l.as_str().to_ascii_lowercase().as_str() {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(BoxedError::new(BasicKind::Error,
                "Invalid Proteoscape line",
                "This column (Is Unique) is not a boolean but it is required to be a boolean ('true' or 'false') in this Proteoscape format",
                l.context().to_owned(),
            ))
        });
    }
    optional { }
);

/// All available Novor versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum ProteoscapeVersion {
    /// The current version
    #[default]
    V2025b,
}
impl std::fmt::Display for ProteoscapeVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<ProteoscapeFormat> for ProteoscapeVersion {
    fn format(self) -> ProteoscapeFormat {
        match self {
            Self::V2025b => V2025B,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V2025b => "2025b",
        }
    }
}

/// Version 2025b
pub const V2025B: ProteoscapeFormat = ProteoscapeFormat {
    version: ProteoscapeVersion::V2025b,
    scan_number: "ms2 id",
    peptide: "peptide sequence",
    xcorr_score: "xcorr score",
    mz: "precursor mz",
    rt: "rt",
    corrected_o_over_k0: "corrected ook0",
    protein: "protein group name",
    z: "charge",
    delta_cn_score: "delta cn score",
    confidence_score: "confidence score",
    matched_ions: "matched ions",
    predicted_o_over_k0: "predicted ook0",
    tim_score: "timscore",
    is_unique: "is unique",
};

impl PSMMetaData for ProteoscapePSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.1.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Proteoscape(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.scan_number)
    }

    fn id(&self) -> String {
        self.scan_number.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        Some((self.confidence_score / 100.0).clamp(-1.0, 1.0))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.confidence_score,
            mzcv::term!(MS:1001153|search engine specific score),
        ))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        Some(self.rt)
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Number(self.scan_number)]) // TODO: timsTOF so works differently maybe
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mz * self.z.to_float())
    }

    type Protein = FastaIdentifier<String>;
    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(std::slice::from_ref(&self.protein))
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&self.peptide.0, &self.peptide.2)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn unique(&self) -> Option<bool> {
        Some(self.is_unique)
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}
