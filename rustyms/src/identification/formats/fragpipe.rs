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
    sequence::{Peptidoform, SemiAmbiguous, SequencePosition, SloppyParsingParameters},
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
static BOOL_ERROR: (&str, &str) = (
    "Invalid FragPipe line",
    "This column is not a boolean but it is required to be a boolean ('true' or 'false') in this FragPipe format",
);

// TODO: it might be nice to be able to parse the glycan composition and put it on the peptide in the right location
format_family!(
    /// The format for any FragPipe file
    FragPipeFormat,
    /// The data from any FragPipe file
    FragpipeData,
    FragPipeVersion, [&V20_OR_21, &V22], b'\t', None;
    required {
        scan: SpectrumId, |location: Location, _| Ok(SpectrumId::Native(location.get_string()));
        spectrum_file: PathBuf, |location: Location, _| Ok(location.get_string().into());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| location.parse_with(|location| Peptidoform::sloppy_pro_forma(
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
        assigned_modifications: Vec<(SequencePosition, f64)>, |location: Location, _| location.or_empty().array(',').map(|m| if let Some((head, tail)) = m.clone().split_once('(') {
            let tail = tail.trim_end_matches(")");
            Ok((
            if head.as_str() == "N-term" {
                SequencePosition::NTerm
            } else if head.as_str() == "C-term" {
                SequencePosition::CTerm
            } else {
                head.as_str()[..head.len()-1].parse::<usize>().map(|i| SequencePosition::Index(i-1)).map_err(|err| CustomError::error(
                    "Invalid FragPipe modification location",
                    format!("The location number {}", explain_number_error(&err)),
                    head.context(),
                ))?
            },
            tail.parse::<f64>(NUMBER_ERROR)?
        ))
        } else {
            Err(CustomError::error(
                "Invalid FragPipe modification",
                "The format `location(modification)` could not be recognised",
                m.context(),
            ))
        }).collect::<Result<Vec<_>, _>>();
        observed_modifications: String, |location: Location, _| Ok(location.get_string());
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
        glycan_score: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        glycan_q_value: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        total_glycan_composition: String, |location: Location, _| Ok(location.or_empty().get_string());
        msfragger_localisation: String, |location: Location, _| Ok(location.or_empty().get_string());
        position_scores: Vec<f64>, |location: Location, _| {
            let data = location.array(')').filter_map(|l| (l.len() > 2).then(|| l.skip(2).parse::<f64>(NUMBER_ERROR))).collect::<Result<Vec<_>, _>>()?;
            if data.is_empty() {
                Ok(None)
            } else {
                Ok(Some(data))
            }
        };
        ions_best_position: usize, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
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
        for (location, mass) in &parsed.assigned_modifications {
            let mass = crate::sequence::SimpleModificationInner::Mass(Mass::new::<crate::system::dalton>(*mass).into());
            match location {
                SequencePosition::NTerm => parsed.peptide.add_simple_n_term(mass.into()),
                SequencePosition::Index(_) => parsed.peptide[*location].modifications.push(mass.into()),
                SequencePosition::CTerm => parsed.peptide.add_simple_c_term(mass.into()),
            }

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
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V20Or21 => "v20/v21",
            Self::V22 => "v22",
        }
    }
}

/// v20 or v21
pub const V20_OR_21: FragPipeFormat = FragPipeFormat {
    version: FragPipeVersion::V20Or21,
    scan: "spectrum",
    raw_file: OptionalColumn::NotAvailable,
    spectrum_file: "spectrum file",
    peptide: "peptide",
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
    observed_modifications: "observed modifications",
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
    glycan_score: OptionalColumn::Optional("glycan score"),
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    total_glycan_composition: OptionalColumn::NotAvailable,
    msfragger_localisation: OptionalColumn::Optional("msfragger localization"),
    position_scores: OptionalColumn::NotAvailable,
    ions_best_position: OptionalColumn::NotAvailable,
};

/// v22
pub const V22: FragPipeFormat = FragPipeFormat {
    version: FragPipeVersion::V22,
    scan: "spectrum",
    raw_file: OptionalColumn::NotAvailable,
    spectrum_file: "spectrum file",
    peptide: "peptide",
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
    observed_modifications: "observed modifications",
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
    glycan_score: OptionalColumn::Optional("glycan score"),
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    total_glycan_composition: OptionalColumn::Optional("total glycan composition"),
    msfragger_localisation: OptionalColumn::Required("msfragger localization"),
    position_scores: OptionalColumn::Required("position scores"),
    ions_best_position: OptionalColumn::Required("ions best position"),
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
