use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, KnownFileFormat, PSM, PSMData, PSMFileFormatVersion, PSMMetaData,
    PSMSource, PeptidoformPresent, SpectrumId, SpectrumIds, common_parser::Location,
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid π-HelixNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PiHelixNovo,
    SemiAmbiguous, PeptidoformPresent, [&PIHELIXNOVO_V1_1], b'\t', Some(vec!["title".to_string(),"peptide".to_string(),"score".to_string()]);
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies| Peptidoform::sloppy_pro_forma_inner(
            &location.base_context(),
            location.full_line(),
            location.range.clone(),
            ontologies,
            &SloppyParsingParameters {
                allow_unwrapped_modifications: true,
                ..Default::default()
            },
        ).map_err(BoxedError::to_owned);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
    optional { }
);

/// The only known version of π-HelixNovo
pub const PIHELIXNOVO_V1_1: PiHelixNovoFormat = PiHelixNovoFormat {
    version: PiHelixNovoVersion::V1_1,
    title: "title",
    peptide: "peptide",
    score: "score",
};

/// All possible π-HelixNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PiHelixNovoVersion {
    #[default]
    /// π-HelixNovo version 1.1
    V1_1,
}

impl std::fmt::Display for PiHelixNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<PiHelixNovoFormat> for PiHelixNovoVersion {
    fn format(self) -> PiHelixNovoFormat {
        match self {
            Self::V1_1 => PIHELIXNOVO_V1_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_1 => "v1.1",
        }
    }
}

impl PSMMetaData for PiHelixNovoPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PiHelixNovo(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> String {
        "-".to_string() // TODO: best would be to use the scan index in some way shape or form
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.score,
            mzcv::term!(MS:1001153|search engine specific score),
        ))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        None
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Native(self.title.clone())])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
    }

    fn ppm_error(&self) -> Option<Ratio> {
        None
    }

    type Protein = crate::NoProtein;

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn unique(&self) -> Option<bool> {
        None
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}
