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
    sequence::{
        Modification, Peptidoform, SemiAmbiguous, SequencePosition, SimpleModification,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, usize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid FragPipe line",
    "This column is not a number but it is required to be a number in this FragPipe format",
);
static IDENTIFIER_ERROR: (&str, &str) = (
    "Invalid FragPipe line",
    "This column is not a fasta identifier but is required to be one in this FragPipe format",
);

// TODO: it might be nice to be able to parse the glycan composition and put it on the peptide in the right location
format_family!(
    /// The format for any FragPipe file
    FragPipeFormat,
    /// The data from any FragPipe file
    FragpipeData,
    FragPipeVersion, [&V20_OR_21, &V22, &PHILOSOPHER], b'\t', None;
    required {
        scan: SpectrumId, |location: Location, _| Ok(SpectrumId::Native(location.get_string()));
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| location.parse_with(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
        ));
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
        intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        assigned_modifications: Vec<(SequencePosition, SimpleModification)>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().array(',').map(|m| if let Some((head, tail)) = m.clone().split_once('(') {
            let head_trim = head.as_str().trim();
            Ok((
            if head_trim.eq_ignore_ascii_case("N-term") {
                SequencePosition::NTerm
            } else if head_trim.eq_ignore_ascii_case("C-term") {
                SequencePosition::CTerm
            } else {
                // Format: `14M` so take only the numeric part
                head.as_str()[..head.len()-1].trim().parse::<usize>().map(|i| SequencePosition::Index(i-1)).map_err(|err| CustomError::error(
                    "Invalid FragPipe modification location",
                    format!("The location number {}", explain_number_error(&err)),
                    head.context(),
                ))?
            },
            Modification::sloppy_modification(tail.full_line(), tail.location.clone(), None, custom_database)?
        ))
        } else {
            Err(CustomError::error(
                "Invalid FragPipe modification",
                "The format `location(modification)` could not be recognised",
                m.context(),
            ))
        }).collect::<Result<Vec<_>, _>>();
        observed_modifications: String, |location: Location, _| Ok(location.get_string());
        is_unique: bool, |location: Location, _| location.parse_with(|l| match l.as_str().to_ascii_lowercase().as_str() {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(CustomError::error(
                "Invalid FragPipe line",
                "This column (Is Unique) is not a boolean but it is required to be a boolean ('true' or 'false') in this FragPipe format",
                l.context(),
            ))
        });
        protein: FastaIdentifier<String>, |location: Location, _| location.parse(IDENTIFIER_ERROR);
        protein_id: String, |location: Location, _| Ok(location.get_string());
        entry_name: String, |location: Location, _| Ok(location.get_string());
        gene: String, |location: Location, _| Ok(location.get_string());
        protein_description: String, |location: Location, _| Ok(location.get_string());
        mapped_genes: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect_vec());
        mapped_proteins: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect_vec());
    }
    optional {
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
        condition: String, |location: Location, _| Ok(Some(location.get_string()));
        ions_best_position: usize, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
        ion_mobility: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        glycan_q_value: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        glycan_score: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        group: String, |location: Location, _| Ok(Some(location.get_string()));
        msfragger_localisation: String, |location: Location, _| Ok(location.or_empty().get_string());
        position_scores: Vec<f64>, |location: Location, _| {
            let data = location.array(')').filter_map(|l| (l.len() > 2).then(|| l.skip(2).parse::<f64>((
                "Invalid FragPipe line",
                "This position score is not a number but it is required to be a number in this FragPipe format",
            )))).collect::<Result<Vec<_>, _>>()?;
            if data.is_empty() {
                Ok(None)
            } else {
                Ok(Some(data))
            }
        };
        protein_end: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        protein_start: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        purity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
        spectrum_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
        total_glycan_composition: String, |location: Location, _| Ok(location.or_empty().get_string());
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
        for (location, modification) in &parsed.assigned_modifications {
            parsed.peptide.add_simple_modification(*location, modification.clone());
        }
        Ok(parsed)
    }
);

/// The Regex to match against FragPipe scan fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"([^/]+)\.(\d+)\.\d+.\d+").unwrap());

impl From<FragpipeData> for IdentifiedPeptidoform {
    fn from(value: FragpipeData) -> Self {
        Self {
            score: Some(value.hyperscore),
            local_confidence: None,
            metadata: MetaData::FragPipe(value),
        }
    }
}

/// All possible FragPipe versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum FragPipeVersion {
    /// Version 20 or 21
    #[default]
    V20Or21,
    /// Version 22
    V22,
    /// Philpsopher
    Philosopher,
}

impl std::fmt::Display for FragPipeVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<FragPipeFormat> for FragPipeVersion {
    fn format(self) -> FragPipeFormat {
        match self {
            Self::V20Or21 => V20_OR_21,
            Self::V22 => V22,
            Self::Philosopher => PHILOSOPHER,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V20Or21 => "v20/v21",
            Self::V22 => "v22",
            Self::Philosopher => "Philosopher",
        }
    }
}

/// Philosopher
pub const PHILOSOPHER: FragPipeFormat = FragPipeFormat {
    version: FragPipeVersion::Philosopher,
    assigned_modifications: "assigned modifications",
    calibrated_experimental_mass: "calibrated observed mass",
    calibrated_experimental_mz: "calibrated observed m/z",
    condition: OptionalColumn::Optional("condition"),
    entry_name: "entry name",
    enzymatic_termini: "number of enzymatic termini",
    expectation: "expectation",
    extended_peptide: OptionalColumn::Optional("extended peptide"),
    gene: "gene",
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyperscore: "hyperscore",
    intensity: "intensity",
    ions_best_position: OptionalColumn::NotAvailable,
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    is_unique: "is unique",
    mapped_genes: "mapped genes",
    mapped_proteins: "mapped proteins",
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    msfragger_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: "observed m/z",
    next_score: "nextscore",
    observed_modifications: "observed modifications",
    peptide_prophet_probability: "peptideprophet probability",
    peptide: "peptide",
    position_scores: OptionalColumn::NotAvailable,
    protein_description: "protein description",
    protein_end: OptionalColumn::NotAvailable,
    protein_id: "protein id",
    protein_start: OptionalColumn::NotAvailable,
    protein: "protein",
    purity: OptionalColumn::Optional("purity"),
    raw_file: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    spectrum_file: OptionalColumn::Optional("spectrum file"),
    total_glycan_composition: OptionalColumn::NotAvailable,
    z: "charge",
};

/// v20 or v21
pub const V20_OR_21: FragPipeFormat = FragPipeFormat {
    version: FragPipeVersion::V20Or21,
    assigned_modifications: "assigned modifications",
    calibrated_experimental_mass: "calibrated observed mass",
    calibrated_experimental_mz: "calibrated observed m/z",
    condition: OptionalColumn::Optional("condition"),
    entry_name: "entry name",
    enzymatic_termini: "number of enzymatic termini",
    expectation: "expectation",
    extended_peptide: OptionalColumn::Optional("extended peptide"),
    gene: "gene",
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyperscore: "hyperscore",
    intensity: "intensity",
    ions_best_position: OptionalColumn::NotAvailable,
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    is_unique: "is unique",
    mapped_genes: "mapped genes",
    mapped_proteins: "mapped proteins",
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    msfragger_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: "observed m/z",
    next_score: "nextscore",
    observed_modifications: "observed modifications",
    peptide_prophet_probability: "peptideprophet probability",
    peptide: "peptide",
    position_scores: OptionalColumn::NotAvailable,
    protein_description: "protein description",
    protein_end: OptionalColumn::Required("protein end"),
    protein_id: "protein id",
    protein_start: OptionalColumn::Required("protein start"),
    protein: "protein",
    purity: OptionalColumn::Required("purity"),
    raw_file: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    spectrum_file: OptionalColumn::Required("spectrum file"),
    total_glycan_composition: OptionalColumn::NotAvailable,
    z: "charge",
};

/// v22
pub const V22: FragPipeFormat = FragPipeFormat {
    version: FragPipeVersion::V22,
    assigned_modifications: "assigned modifications",
    calibrated_experimental_mass: "calibrated observed mass",
    calibrated_experimental_mz: "calibrated observed m/z",
    condition: OptionalColumn::Optional("condition"),
    entry_name: "entry name",
    enzymatic_termini: "number of enzymatic termini",
    expectation: "expectation",
    extended_peptide: OptionalColumn::Required("extended peptide"),
    gene: "gene",
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyperscore: "hyperscore",
    intensity: "intensity",
    ions_best_position: OptionalColumn::Optional("ions best position"),
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    is_unique: "is unique",
    mapped_genes: "mapped genes",
    mapped_proteins: "mapped proteins",
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    msfragger_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: "observed m/z",
    next_score: "nextscore",
    observed_modifications: "observed modifications",
    peptide_prophet_probability: "probability",
    peptide: "peptide",
    position_scores: OptionalColumn::Optional("position scores"),
    protein_description: "protein description",
    protein_end: OptionalColumn::Required("protein end"),
    protein_id: "protein id",
    protein_start: OptionalColumn::Required("protein start"),
    protein: "protein",
    purity: OptionalColumn::Required("purity"),
    raw_file: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    spectrum_file: OptionalColumn::Required("spectrum file"),
    total_glycan_composition: OptionalColumn::Optional("total glycan composition"),
    z: "charge",
};

/// The scans identifier for a MSFragger identification
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct FragPipeID {
    /// The file, if defined
    pub file: PathBuf,
    /// The scan number triplet
    pub scan: (usize, usize, usize),
}

impl std::fmt::Display for FragPipeID {
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

impl std::str::FromStr for FragPipeID {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let split = s.rsplitn(4, '.').collect_vec();
        if split.len() == 4 {
            Ok(Self {
                file: Path::new(&split[3]).to_owned(),
                scan: (
                    split[2].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid FragPipe ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                    split[1].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid FragPipe ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                    split[0].parse().map_err(|err| {
                        CustomError::error(
                            "Invalid FragPipe ID",
                            format!("The scan number {}", explain_number_error(&err)),
                            Context::None,
                        )
                    })?,
                ),
            })
        } else {
            Err(CustomError::error(
                "Invalid FragPipe ID",
                "An FragPipe ID should consist of 4 parts separated by dots (file.scan.scan.scan)",
                Context::None,
            ))
        }
    }
}
