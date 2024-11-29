use crate::{
    error::CustomError,
    identification::{
        common_parser::OptionalColumn, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
    },
    ontologies::CustomDatabase,
    LinearPeptide, SemiAmbiguous, SloppyParsingParameters,
};

use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter,
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
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
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
            .get_or_init(|| regex::Regex::new(r"^(.*):index=(\d+)$").unwrap())
            .captures(&parsed.title)
        {
            parsed.raw_file = Some(PathBuf::from(m.get(1).unwrap().as_str()));
            parsed.scan = Some(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
        }
        Ok(parsed)
    }
);

/// The Regex to match against PowerNovo scan fields
static IDENTIFER_REGEX: std::sync::OnceLock<regex::Regex> = std::sync::OnceLock::new();

impl From<PowerNovoData> for IdentifiedPeptide {
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
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum PowerNovoVersion {
    #[default]
    /// PowerNovo version 1.0.1
    V1_0_1,
}

impl std::fmt::Display for PowerNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V1_0_1 => "v1.0.1",
            }
        )
    }
}
