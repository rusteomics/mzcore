use crate::{
    error::CustomError,
    identification::{IdentifiedPeptide, IdentifiedPeptideSource, MetaData},
    modification::{Ontology, SimpleModification},
    ontologies::CustomDatabase,
    system::{usize::Charge, MassOverCharge},
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
    "Invalid InstaNovo line",
    "This column is not a number but it is required to be a number in this format",
);

static BUILT_IN_MODIFICATIONS: std::sync::OnceLock<SloppyParsingParameters> =
    std::sync::OnceLock::new();

format_family!(
    /// The format for any InstaNovo file
    InstaNovoFormat,
    /// The data from any InstaNovo file
    InstaNovoData,
    InstaNovoVersion, [&INSTANOVO_V1_0_0], b',', None;
    required {
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        peptide: LinearPeptide<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            BUILT_IN_MODIFICATIONS.get_or_init(|| SloppyParsingParameters {
                replace_mass_modifications: Some(
                    vec![Ontology::Unimod.find_id(35, None).unwrap(), Ontology::Unimod.find_id(21, None).unwrap(), Ontology::Unimod.find_id(4, None).unwrap()]
                ),
                ..Default::default()
            }
        ));

        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location
            .trim_start_matches("[").trim_end_matches("]")
            .array(',')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
    }
    optional { }
);

impl From<InstaNovoData> for IdentifiedPeptide {
    fn from(value: InstaNovoData) -> Self {
        Self {
            score: Some(2.0 / (1.0 + 1.01_f64.powf(-value.score))),
            local_confidence: Some(
                value
                    .local_confidence
                    .iter()
                    .map(|v| 2.0 / (1.0 + 1.25_f64.powf(-v)))
                    .collect(),
            ),
            metadata: MetaData::InstaNovo(value),
        }
    }
}

/// The only known version of InstaNovo
pub const INSTANOVO_V1_0_0: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::V1_0_0,
    scan: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "preds",
    score: "log_probs",
    local_confidence: "token_log_probs",
};

/// All possible InstaNovo versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum InstaNovoVersion {
    #[default]
    /// InstaNovo version 1.0.0
    V1_0_0,
}

impl std::fmt::Display for InstaNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V1_0_0 => "v1.0.0",
            }
        )
    }
}
