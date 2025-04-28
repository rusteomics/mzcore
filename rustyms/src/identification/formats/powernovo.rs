use std::{
    path::{Path, PathBuf},
    sync::LazyLock,
};

use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    identification::{
        BoxedIdentifiedPeptideIter, IdentifiedPeptidoform, IdentifiedPeptidoformSource,
        IdentifiedPeptidoformVersion, MetaData,
        common_parser::Location,
        common_parser::OptionalColumn,
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    sequence::{Peptidoform, SemiAmbiguous, SloppyParsingParameters},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid PowerNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    /// The format for any PowerNovo file
    PowerNovoFormat,
    /// The data from any PowerNovo file
    PowerNovoData,
    PowerNovoVersion, [&POWERNOVO_V1_0_1], b',', None;
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default(),
        );
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location.array(' ')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
    }
    optional {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
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

impl From<PowerNovoData> for IdentifiedPeptidoform {
    fn from(value: PowerNovoData) -> Self {
        Self {
            score: Some(value.score),
            local_confidence: Some(value.local_confidence.clone()),
            metadata: MetaData::PowerNovo(value),
        }
    }
}

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
