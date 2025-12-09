use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, KnownFileFormat, MetaData,
    PeptidoformPresent, SpectrumIds, common_parser::Location,
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
    "Invalid PepNet line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PepNet,
    SemiAmbiguous, PeptidoformPresent, [&PEPNET_V1_0], b'\t', None;
    required {
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            ontologies,
            &SloppyParsingParameters::default(),
        ).map_err(BoxedError::to_owned);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location
            .trim_start_matches("[").trim_end_matches("]")
            .array(',')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
        ppm_diff: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<mzcore::system::ratio::ppm>);
    }
    optional { }
);

/// The only known version of PepNet
pub const PEPNET_V1_0: PepNetFormat = PepNetFormat {
    version: PepNetVersion::V1_0,
    peptide: "denovo",
    score: "score",
    local_confidence: "positional score",
    ppm_diff: "ppm difference",
};

/// All possible PepNet versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PepNetVersion {
    #[default]
    /// PepNet version 1.0
    V1_0,
}

impl std::fmt::Display for PepNetVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PepNetFormat> for PepNetVersion {
    fn format(self) -> PepNetFormat {
        match self {
            Self::V1_0 => PEPNET_V1_0,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_0 => "v1.0",
        }
    }
}

impl MetaData for PepNetData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PepNet(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> String {
        "-".to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        Some(Cow::Borrowed(self.local_confidence.as_slice()))
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.score,
            mzcv::term!(MS:1001153|search engine specific score),
        ))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        Some(self.local_confidence.as_slice())
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
        SpectrumIds::None
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
    }

    fn ppm_error(&self) -> Option<Ratio> {
        Some(self.ppm_diff)
    }

    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

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
