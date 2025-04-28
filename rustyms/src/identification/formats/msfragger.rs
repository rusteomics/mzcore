use std::{
    path::{Path, PathBuf},
    sync::LazyLock,
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::explain_number_error,
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform,
        IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, MetaData, SpectrumId,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    sequence::{Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Time, usize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid MSFragger line",
    "This column is not a number but it is required to be a number in this MSFragger format",
);
static IDENTIFIER_ERROR: (&str, &str) = (
    "Invalid MSFragger line",
    "This column is not a fasta identifier but is required to be one in this MSFragger format",
);
static BOOL_ERROR: (&str, &str) = (
    "Invalid MSFragger line",
    "This column is not a boolean but it is required to be a boolean ('true' or 'false') in this MSFragger format",
);

format_family!(
    /// The format for any MSFragger file
    MSFraggerFormat,
    /// The data from any MSFragger file
    MSFraggerData,
    MSFraggerVersion, [&V21, &V22], b'\t', None;
    required {
        scan: SpectrumId, |location: Location, _| Ok(SpectrumId::Native(location.get_string()));
        spectrum_file: PathBuf, |location: Location, _| Ok(location.get_string().into());
        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().parse_with(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
        ));
        extended_peptide: Box<[Option<Peptidoform<SemiAmbiguous>>; 3]>, |location: Location, custom_database: Option<&CustomDatabase>| {
            let peptides = location.clone().array('.').map(|l| l.or_empty().parse_with(|location| Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                custom_database,
                &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
            ))).collect::<Result<Vec<_>,_>>()?;
            if peptides.len() == 3 {
                Ok(Box::new([peptides[0].clone(), peptides[1].clone(), peptides[2].clone()]))
            } else {
                Err(CustomError::error("Invalid extened peptide", "The extended peptide should contain the prefix.peptide.suffix for all peptides.", location.context()))
            }
        };
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::s>);
        /// Experimental mass
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        calibrated_experimental_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        /// Experimental mz
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        calibrated_experimental_mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        expectation: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        hyperscore: f64, |location: Location, _| location.parse(NUMBER_ERROR).map(|s: f64| s / 100.0);
        next_score: f64, |location: Location, _| location.parse(NUMBER_ERROR).map(|s: f64| s / 100.0);
        peptide_prophet_probability: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        enzymatic_termini: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        missed_cleavages: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_start: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_end: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        assigned_modifications: String, |location: Location, _| Ok(location.get_string());
        purity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        is_unique: bool, |location: Location, _| location.parse(BOOL_ERROR);
        protein: FastaIdentifier<String>, |location: Location, _| location.parse(IDENTIFIER_ERROR);
        protein_id: String, |location: Location, _| Ok(location.get_string());
        entry_name: String, |location: Location, _| Ok(location.get_string());
        gene: String, |location: Location, _| Ok(location.get_string());
        protein_description: String, |location: Location, _| Ok(location.get_string());
        mapped_genes: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect_vec());
        mapped_proteins: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect_vec());
    }
    optional {
        raw_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
        condition: String, |location: Location, _| Ok(Some(location.get_string()));
        group: String, |location: Location, _| Ok(Some(location.get_string()));
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        if let SpectrumId::Native(native) = &parsed.scan {
            if let Some(m) = IDENTIFER_REGEX
                .captures(native)
            {
                parsed.raw_file = Some(m.get(1).unwrap().as_str().into());
                parsed.scan =
                    SpectrumId::Number(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
            }
        }
        if parsed.peptide.is_none() {
            parsed.peptide.clone_from(&parsed.extended_peptide[1]);
        }
        Ok(parsed)
    }
);

/// The Regex to match against MSFragger scan fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"([^/]+)\.(\d+)\.\d+.\d+").unwrap());

impl From<MSFraggerData> for IdentifiedPeptidoform {
    fn from(value: MSFraggerData) -> Self {
        Self {
            score: Some(value.hyperscore),
            local_confidence: None,
            metadata: MetaData::MSFragger(value),
        }
    }
}

/// All possible MSFragger versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum MSFraggerVersion {
    /// Version 21
    #[default]
    V21,
    /// Version 22
    V22,
}

impl std::fmt::Display for MSFraggerVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<MSFraggerFormat> for MSFraggerVersion {
    fn format(self) -> MSFraggerFormat {
        match self {
            Self::V21 => V21,
            Self::V22 => V22,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V21 => "v21",
            Self::V22 => "v22",
        }
    }
}

/// v21
pub const V21: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::V21,
    scan: "spectrum",
    raw_file: OptionalColumn::NotAvailable,
    spectrum_file: "spectrum file",
    peptide: "modified peptide",
    extended_peptide: "extended peptide",
    z: "charge",
    rt: "retention",
    mass: "observed mass",
    calibrated_experimental_mass: "calibrated observed mass",
    mz: "observed m/z",
    calibrated_experimental_mz: "calibrated observed m/z",
    expectation: "expectation",
    hyperscore: "hyperscore",
    next_score: "nextscore",
    peptide_prophet_probability: "peptideprophet probability",
    enzymatic_termini: "number of enzymatic termini",
    missed_cleavages: "number of missed cleavages",
    protein_start: "protein start",
    protein_end: "protein end",
    intensity: "intensity",
    assigned_modifications: "assigned modifications",
    purity: "purity",
    is_unique: "is unique",
    protein: "protein",
    protein_id: "protein id",
    entry_name: "entry name",
    gene: "gene",
    protein_description: "protein description",
    mapped_genes: "mapped genes",
    mapped_proteins: "mapped proteins",
    condition: OptionalColumn::Optional("condition"),
    group: OptionalColumn::Optional("group"),
};

/// v22
pub const V22: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::V22,
    scan: "spectrum",
    raw_file: OptionalColumn::NotAvailable,
    spectrum_file: "spectrum file",
    peptide: "modified peptide",
    extended_peptide: "extended peptide",
    z: "charge",
    rt: "retention",
    mass: "observed mass",
    calibrated_experimental_mass: "calibrated observed mass",
    mz: "observed m/z",
    calibrated_experimental_mz: "calibrated observed m/z",
    expectation: "expectation",
    hyperscore: "hyperscore",
    next_score: "nextscore",
    peptide_prophet_probability: "probability",
    enzymatic_termini: "number of enzymatic termini",
    missed_cleavages: "number of missed cleavages",
    protein_start: "protein start",
    protein_end: "protein end",
    intensity: "intensity",
    assigned_modifications: "assigned modifications",
    purity: "purity",
    is_unique: "is unique",
    protein: "protein",
    protein_id: "protein id",
    entry_name: "entry name",
    gene: "gene",
    protein_description: "protein description",
    mapped_genes: "mapped genes",
    mapped_proteins: "mapped proteins",
    condition: OptionalColumn::Optional("condition"),
    group: OptionalColumn::Optional("group"),
};

/// The scans identifier for a MSFragger identification
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
pub struct MSFraggerID {
    /// The file, if defined
    pub file: PathBuf,
    /// The scan number triplet
    pub scan: (usize, usize, usize),
}

impl std::fmt::Display for MSFraggerID {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}.{}.{}.{}",
            self.file.to_string_lossy(),
            self.scan.0,
            self.scan.1,
            self.scan.2,
        )
    }
}

impl std::str::FromStr for MSFraggerID {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let split = s.rsplitn(4, '.').collect_vec();
        if split.len() == 4 {
            Ok(Self {
                file: Path::new(&split[3]).to_owned(),
                scan: (
                    split[2].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid MSFragger ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                    split[1].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid MSFragger ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                    split[0].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid MSFragger ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                ),
            })
        } else {
            Err(CustomError::error(
                "Invalid MSFragger ID",
                "An MSFragger ID should consist of 4 parts separated by dots (file.scan.scan.scan)",
                Context::None,
            ))
        }
    }
}
