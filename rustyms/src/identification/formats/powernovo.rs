use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
    sync::LazyLock,
};

use serde::{Deserialize, Serialize};

use crate::{
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::CompoundPeptidoformIon,
    sequence::{Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid PowerNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PowerNovo,
    SemiAmbiguous, PeptidoformPresent, [&POWERNOVO_V1_0_1], b',', None;
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default(),
        ).map_err(BoxedError::to_owned);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location.array(' ')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
    }
    optional {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, BoxedError<'static>> {
        if let Some(m) = IDENTIFER_REGEX
            .captures(&parsed.title)
        {
            parsed.raw_file = Some(PathBuf::from(m.get(1).unwrap().as_str()));
            parsed.scan = Some(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
        }
        Ok(parsed)
    }
);

/// The Regex to match against PowerNovo scan fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"^(.*):index=(\d+)$").unwrap());

/// The only known version of PowerNovo
pub const POWERNOVO_V1_0_1: PowerNovoFormat = PowerNovoFormat {
    version: PowerNovoVersion::V1_0_1,
    scan: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::NotAvailable,
    title: "spectrum name",
    peptide: "powernovo peptides",
    score: "powernovo score",
    local_confidence: "powernovo aascore",
};

/// All possible PowerNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PowerNovoVersion {
    #[default]
    /// PowerNovo version 1.0.1
    V1_0_1,
}

impl std::fmt::Display for PowerNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PowerNovoFormat> for PowerNovoVersion {
    fn format(self) -> PowerNovoFormat {
        match self {
            Self::V1_0_1 => POWERNOVO_V1_0_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_0_1 => "v1.0.1",
        }
    }
}

impl MetaData for PowerNovoData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PowerNovo(self.version)
    }

    fn id(&self) -> String {
        self.scan
            .as_ref()
            .map_or("-".to_string(), ToString::to_string)
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        Some(Cow::Borrowed(self.local_confidence.as_slice()))
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        Some(self.local_confidence.as_slice())
    }

    fn charge(&self) -> Option<Charge> {
        None
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        self.scan.as_ref().map_or(SpectrumIds::None, |scan| {
            self.raw_file.clone().map_or_else(
                || SpectrumIds::FileNotKnown(vec![SpectrumId::Index(*scan)]),
                |raw_file| SpectrumIds::FileKnown(vec![(raw_file, vec![SpectrumId::Index(*scan)])]),
            )
        })
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
    }

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
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
}
