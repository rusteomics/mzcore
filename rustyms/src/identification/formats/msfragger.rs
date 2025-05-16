use crate::{
    error::CustomError,
    helper_functions::explain_number_error,
    identification::{
        BoxedIdentifiedPeptideIter, IdentifiedPeptidoform, IdentifiedPeptidoformSource,
        IdentifiedPeptidoformVersion, MetaData,
        common_parser::{Location, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::SequencePosition,
    sequence::{
        Modification, Peptidoform, SemiAmbiguous, SimpleModification, SloppyParsingParameters,
    },
    system::{Mass, Time, usize::Charge},
};
use serde::{Deserialize, Serialize};

use super::FastaIdentifier;

static NUMBER_ERROR: (&str, &str) = (
    "Invalid MSFragger line",
    "This column is not a number but it is required to be a number in this MSFragger format",
);
static IDENTIFIER_ERROR: (&str, &str) = (
    "Invalid MSFragger line",
    "This column is not a fasta identifier but is required to be one in this MSFragger format",
);

format_family!(
    /// The format for MSFragger data
    MSFraggerFormat,
    /// The data for MSFragger data
    MSFraggerData,
    MSFraggerVersion, [&VERSION_V4_2], b'\t', None;
    required {
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        ion_mobility: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        rank: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| location.parse_with(|location| {
            Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                custom_database,
                &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
            )});
        protein: FastaIdentifier<String>, |location: Location, _| location.parse(IDENTIFIER_ERROR);
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        delta_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        missed_cleavages: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        hyperscore: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        nextscore: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        expectscore: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        score_without_delta_mass: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        best_score_with_delta_mass: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        second_best_score_with_delta_mass: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        delta_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        num_matched_ions: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        tot_num_ions: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        num_tol_term: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        modification_info: Vec<(SequencePosition, SimpleModification)>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().array(',').map(|m| if let Some((head, tail)) = m.clone().split_once('(') {
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

        // modification_info,  best_locs, localization_scores
    }
    optional { }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        for (location, modification) in &parsed.modification_info {
            parsed.peptide.add_simple_modification(*location, modification.clone());
        }
        Ok(parsed)
    }
);

impl From<MSFraggerData> for IdentifiedPeptidoform {
    fn from(value: MSFraggerData) -> Self {
        Self {
            score: Some(value.hyperscore / 100.0),
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
    /// Current MSFragger version
    #[default]
    V4_2,
}

impl std::fmt::Display for MSFraggerVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<MSFraggerFormat> for MSFraggerVersion {
    fn format(self) -> MSFraggerFormat {
        match self {
            Self::V4_2 => VERSION_V4_2,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V4_2 => "v4.2",
        }
    }
}

/// The only supported format for MSFragger data
pub const VERSION_V4_2: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::V4_2,
    scan_number: "scannum",
    mass: "precursor_neutral_mass",
    rt: "retention_time",
    z: "charge",
    ion_mobility: "ion_mobility",
    rank: "hit_rank",
    peptide: "peptide",
    protein: "proteins",
    theoretical_mass: "calc_neutral_pep_mass",
    delta_mass: "massdiff",
    missed_cleavages: "num_missed_cleavages",
    hyperscore: "hyperscore",
    nextscore: "nextscore",
    expectscore: "expectscore",
    score_without_delta_mass: "score_without_delta_mass",
    best_score_with_delta_mass: "best_score_with_delta_mass",
    second_best_score_with_delta_mass: "second_best_score_with_delta_mass",
    delta_score: "delta_score",
    num_matched_ions: "num_matched_ions",
    tot_num_ions: "tot_num_ions",
    num_tol_term: "num_tol_term",
    modification_info: "modification_info",
};
