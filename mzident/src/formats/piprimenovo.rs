use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, PSMData,
    PSMSource, PSMFileFormatVersion, KnownFileFormat, MaybePeptidoform,
    PSM, PSMMetaData, SpectrumId, SpectrumIds, common_parser::Location,
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
    "Invalid π-PrimeNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PiPrimeNovo,
    SemiAmbiguous, MaybePeptidoform, [&PIPRIMENOVO_V0_1], b'\t', None;
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, ontologies: &Ontologies| location.or_empty().map(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            ontologies,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned)).transpose();
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        z: Charge, |location: Location, _| location
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<mzcore::system::e>);
    }
    optional { }
);

/// The only known version of π-PrimeNovo
pub const PIPRIMENOVO_V0_1: PiPrimeNovoFormat = PiPrimeNovoFormat {
    version: PiPrimeNovoVersion::V0_1,
    title: "label",
    peptide: "prediction",
    score: "score",
    z: "charge",
};

/// All possible π-PrimeNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PiPrimeNovoVersion {
    #[default]
    /// π-PrimeNovo version 0.1
    V0_1,
}

impl std::fmt::Display for PiPrimeNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<PiPrimeNovoFormat> for PiPrimeNovoVersion {
    fn format(self) -> PiPrimeNovoFormat {
        match self {
            Self::V0_1 => PIPRIMENOVO_V0_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V0_1 => "v0.1",
        }
    }
}

impl PSMMetaData for PiPrimeNovoPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PiPrimeNovo(self.version)
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
        Some(self.z)
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

    type Protein<'a> = crate::NoProtein;
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
