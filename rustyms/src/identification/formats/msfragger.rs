use std::{marker::PhantomData, path::PathBuf, sync::LazyLock};

use crate::{
    error::CustomError,
    glycan::MonoSaccharide,
    helper_functions::explain_number_error,
    identification::{
        BoxedIdentifiedPeptideIter, IdentifiedPeptidoform, IdentifiedPeptidoformSource,
        IdentifiedPeptidoformVersion, MetaData, PeptidoformPresent, SpectrumId,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::{Chemical, SequencePosition},
    quantities::{Tolerance, WithinTolerance},
    sequence::{
        Modification, Peptidoform, SemiAmbiguous, SimpleLinear, SimpleModification,
        SimpleModificationInner, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};
use ordered_float::OrderedFloat;
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

// TODO: it might be nice to be able to put the open modification on the peptide in the right location
// TODO: in the localisation the lowercase character(s) indicate the position of the opensearch (observed modifications).
// It would be best to use this to place the mod properly (also with the position scoring if present and the mod is ambiguous).
format_family!(
    /// The format for MSFragger data
    MSFraggerFormat,
    /// The data for MSFragger data
    MSFraggerData,
    SimpleLinear, PeptidoformPresent, MSFraggerVersion, [&VERSION_V4_2, &FRAGPIPE_V20_OR_21, &FRAGPIPE_V22, &PHILOSOPHER], b'\t', None;
    required {
        expectation_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        hyper_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        missed_cleavages: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        next_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide: Peptidoform<SimpleLinear>, |location: Location, custom_database: Option<&CustomDatabase>| location.parse_with(|location| {
            Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                custom_database,
                &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
        ).map(Into::into)});
        protein: FastaIdentifier<String>, |location: Location, _| location.parse(IDENTIFIER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        scan: SpectrumId, |location: Location, _| Ok(SpectrumId::Native(location.get_string()));
        modifications: Vec<(SequencePosition, SimpleModification)>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().array(',').map(|m| if let Some((head, tail)) = m.clone().split_once('(') {
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
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
    }
    optional {
        best_score_with_delta_mass: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        calibrated_experimental_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        calibrated_experimental_mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        condition: String, |location: Location, _| Ok(Some(location.get_string()));
        delta_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        entry_name: String, |location: Location, _| Ok(location.get_string());
        enzymatic_termini: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
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
        glycan_q_value: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        glycan_score: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        group: String, |location: Location, _| Ok(Some(location.get_string()));
        intensity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ion_mobility: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        ions_best_position: usize, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
        is_unique: bool, |location: Location, _| location.parse_with(|l| match l.as_str().to_ascii_lowercase().as_str() {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(CustomError::error(
                "Invalid FragPipe line",
                "This column (Is Unique) is not a boolean but it is required to be a boolean ('true' or 'false') in this FragPipe format",
                l.context(),
            ))
        });
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        gene: String, |location: Location, _| Ok(location.get_string());
        mapped_genes: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect::<Vec<_>>());
        mapped_proteins: Vec<String>, |location: Location, _| Ok(location.get_string().split(',').map(|s| s.trim().to_string()).collect::<Vec<_>>());
        num_matched_ions: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        num_tol_term: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        open_search_localisation: String, |location: Location, _| Ok(location.or_empty().get_string());
        /// Either glycan 'HexNAc(4)Hex(5)Fuc(1)NeuAc(2) % 2350.8304' or mod 'Mod1: First isotopic peak (Theoretical: 1.0024)' 'Mod1: Deamidation (PeakApex: 0.9836, Theoretical: 0.9840)'
        open_search_modification: MSFraggerOpenModification, |location: Location, _|
            location.or_empty().parse_with(|location| location.as_str().find('%').map_or_else(
                || Ok(MSFraggerOpenModification::Other(location.as_str().to_owned())),
                |end| MonoSaccharide::from_byonic_composition(&location.as_str()[..end]).map(MSFraggerOpenModification::Glycan)));
        open_search_position_scores: Vec<f64>, |location: Location, _| {
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
        peptide_prophet_probability: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        protein_description: String, |location: Location, _| Ok(location.get_string());
        protein_end: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        protein_id: String, |location: Location, _| Ok(location.get_string());
        protein_start: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        purity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        rank: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        score_without_delta_mass: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        second_best_score_with_delta_mass: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
        total_glycan_composition: Vec<(MonoSaccharide, isize)>, |location: Location, _| location.or_empty().parse_with(|location| location.as_str().find('%').map_or_else(
                || MonoSaccharide::from_byonic_composition(location.as_str()),
                |end| MonoSaccharide::from_byonic_composition(&location.as_str()[..end])));
        tot_num_ions: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        // Parse the scan identifier to retrieve the file and number separately
        if let SpectrumId::Native(native) = &parsed.scan {
            if let Some(m) = IDENTIFER_REGEX
                .captures(native)
            {
                parsed.raw_file = Some(m.get(1).unwrap().as_str().into());
                parsed.scan =
                    SpectrumId::Number(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
            }
        }
        // Normal modifications
        for (location, modification) in &parsed.modifications {
            parsed.peptide.add_simple_modification(*location, modification.clone());
        }
        // Open modification - check for a glycan composition mod
        if let Some(glycan) = match parsed.open_search_modification.as_ref() {Some(MSFraggerOpenModification::Glycan(glycan)) => Some(glycan), _ => None}.or(parsed.total_glycan_composition.as_ref()) {
            let modification: SimpleModification = SimpleModificationInner::Glycan(glycan.clone()).into();
            let target_mass = Mass::new::<crate::system::dalton>(glycan.iter().fold(0.0, |acc, (sugar, amount)| sugar.formula().monoisotopic_mass().value.mul_add(*amount as f64, acc)));
            let mut placed = false;

            // Check if the location field contains any lowercase chars indication the location
            if let Some(location) = parsed.open_search_localisation.as_ref().map(|l| l.chars().enumerate().filter_map(|(i, c)| c.is_ascii_lowercase().then_some(i)).collect::<Vec<_>>()) {
                match location.len() {
                    0 => (),
                    1 => {
                        // Check if the modification is already placed as mass modification (this is not necessarily always present )
                        let mut index = None;
                        for (i, m) in parsed.peptide[location[0]].modifications.iter().enumerate() {
                            if let Some(SimpleModificationInner::Mass(mass)) = m.simple().map(AsRef::as_ref) {
                                if Tolerance::Absolute(Mass::new::<crate::system::dalton>(1.0)).within(&mass.into_inner(), &target_mass) {
                                    index = Some(i);
                                    break;
                                }
                            }
                        }
                        if let Some(index) = index {
                            parsed.peptide[location[0]].modifications[index] = modification.clone().into();
                        } else {
                            parsed.peptide[location[0]].modifications.push(modification.clone().into());
                        }
                        placed = true;
                    }
                    _ => {
                        // If there are multiple locations add this as an ambiguous modification
                        let positions = parsed.open_search_position_scores.as_ref().map_or_else(|| location.iter().map(|i| (SequencePosition::Index(*i), None)).collect::<Vec<(_, _)>>(), |scores| location.iter().map(|i| SequencePosition::Index(*i)).zip(scores.iter().map(|s| Some(OrderedFloat(*s)))).collect::<Vec<(_, _)>>());
                        let _ = parsed.peptide.add_ambiguous_modification(modification.clone(), None, &positions, None, None, false);
                        placed = true;
                    }
                }
            }
            // If the location field does not exist or does not contain lowercase characters try to find the modification as a mass modification
            if !placed {
                let mut index = None;
                for seq in parsed.peptide.sequence_mut() {
                    for (i, m) in seq.modifications.iter().enumerate() {
                        if let Some(SimpleModificationInner::Mass(mass)) = m.simple().map(AsRef::as_ref) {
                            if Tolerance::Absolute(Mass::new::<crate::system::dalton>(1.0)).within(&mass.into_inner(), &target_mass) {
                                index = Some(i);
                                break;
                            }
                        }
                    }
                    if let Some(index) = index {
                        seq.modifications[index] = modification.into();
                        break;
                    }
                }
            }
        }
        Ok(parsed)
    }
);

/// The Regex to match against MSFragger scan fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"([^/]+)\.(\d+)\.\d+.\d+").unwrap());

impl From<MSFraggerData> for IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent> {
    fn from(value: MSFraggerData) -> Self {
        Self {
            score: Some(value.hyper_score / 100.0),
            local_confidence: None,
            metadata: MetaData::MSFragger(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

/// A MSFragger open search modification
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum MSFraggerOpenModification {
    /// A glycan composition
    Glycan(Vec<(MonoSaccharide, isize)>),
    /// Any other modification
    Other(String),
}

/// All possible MSFragger versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum MSFraggerVersion {
    /// Current MSFragger version
    #[default]
    V4_2,
    /// Version 20 or 21
    FragPipeV20Or21,
    /// Version 22
    FragPipeV22,
    /// Philpsopher
    Philosopher,
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
            Self::FragPipeV20Or21 => FRAGPIPE_V20_OR_21,
            Self::FragPipeV22 => FRAGPIPE_V22,
            Self::Philosopher => PHILOSOPHER,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V4_2 => "v4.2",
            Self::FragPipeV20Or21 => "FragPipe v20/v21",
            Self::FragPipeV22 => "FragPipe v22",
            Self::Philosopher => "Philosopher",
        }
    }
}

/// The only supported format for MSFragger data
pub const VERSION_V4_2: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::V4_2,
    best_score_with_delta_mass: OptionalColumn::Required("best_score_with_delta_mass"),
    calibrated_experimental_mass: OptionalColumn::NotAvailable,
    calibrated_experimental_mz: OptionalColumn::NotAvailable,
    condition: OptionalColumn::NotAvailable,
    delta_score: OptionalColumn::Required("delta_score"),
    entry_name: OptionalColumn::NotAvailable,
    enzymatic_termini: OptionalColumn::NotAvailable,
    expectation_score: "expectscore",
    extended_peptide: OptionalColumn::NotAvailable,
    gene: OptionalColumn::NotAvailable,
    glycan_q_value: OptionalColumn::NotAvailable,
    glycan_score: OptionalColumn::NotAvailable,
    group: OptionalColumn::NotAvailable,
    hyper_score: "hyperscore",
    intensity: OptionalColumn::NotAvailable,
    ion_mobility: OptionalColumn::Required("ion_mobility"),
    ions_best_position: OptionalColumn::NotAvailable,
    is_unique: OptionalColumn::NotAvailable,
    mapped_genes: OptionalColumn::NotAvailable,
    mapped_proteins: OptionalColumn::NotAvailable,
    mass: "precursor_neutral_mass", // Is this experimental M, if so this must be changed to MH+
    missed_cleavages: "num_missed_cleavages",
    modifications: "modification_info",
    open_search_localisation: OptionalColumn::Required("best_locs"),
    mz: OptionalColumn::NotAvailable,
    next_score: "nextscore",
    num_matched_ions: OptionalColumn::Required("num_matched_ions"),
    num_tol_term: OptionalColumn::Required("num_tol_term"),
    open_search_modification: OptionalColumn::NotAvailable,
    peptide_prophet_probability: OptionalColumn::NotAvailable,
    peptide: "peptide",
    open_search_position_scores: OptionalColumn::Required("localization_scores"),
    protein_description: OptionalColumn::NotAvailable,
    protein_end: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_start: OptionalColumn::NotAvailable,
    protein: "proteins",
    purity: OptionalColumn::NotAvailable,
    rank: OptionalColumn::Required("hit_rank"),
    rt: "retention_time",
    scan: "scannum",
    score_without_delta_mass: OptionalColumn::Required("score_without_delta_mass"),
    second_best_score_with_delta_mass: OptionalColumn::Required(
        "second_best_score_with_delta_mass",
    ),
    raw_file: OptionalColumn::NotAvailable,
    tot_num_ions: OptionalColumn::Required("tot_num_ions"),
    total_glycan_composition: OptionalColumn::NotAvailable,
    z: "charge",
};

/// Philosopher
pub const PHILOSOPHER: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::Philosopher,
    best_score_with_delta_mass: OptionalColumn::NotAvailable,
    calibrated_experimental_mass: OptionalColumn::Required("calibrated observed mass"),
    calibrated_experimental_mz: OptionalColumn::Required("calibrated observed m/z"),
    condition: OptionalColumn::Optional("condition"),
    delta_score: OptionalColumn::NotAvailable,
    entry_name: OptionalColumn::Required("entry name"),
    enzymatic_termini: OptionalColumn::Required("number of enzymatic termini"),
    expectation_score: "expectation",
    extended_peptide: OptionalColumn::Optional("extended peptide"),
    gene: OptionalColumn::Required("gene"),
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyper_score: "hyperscore",
    intensity: OptionalColumn::Required("intensity"),
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    ions_best_position: OptionalColumn::NotAvailable,
    is_unique: OptionalColumn::Required("is unique"),
    mapped_genes: OptionalColumn::Required("mapped genes"),
    mapped_proteins: OptionalColumn::Required("mapped proteins"),
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    modifications: "assigned modifications",
    open_search_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: OptionalColumn::Required("observed m/z"),
    next_score: "nextscore",
    num_matched_ions: OptionalColumn::NotAvailable,
    num_tol_term: OptionalColumn::NotAvailable,
    open_search_modification: OptionalColumn::Required("observed modifications"),
    peptide_prophet_probability: OptionalColumn::Required("peptideprophet probability"),
    peptide: "peptide",
    open_search_position_scores: OptionalColumn::NotAvailable,
    protein_description: OptionalColumn::Required("protein description"),
    protein_end: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::Required("protein id"),
    protein_start: OptionalColumn::NotAvailable,
    protein: "protein",
    purity: OptionalColumn::Optional("purity"),
    rank: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    score_without_delta_mass: OptionalColumn::NotAvailable,
    second_best_score_with_delta_mass: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Optional("spectrum file"),
    tot_num_ions: OptionalColumn::NotAvailable,
    total_glycan_composition: OptionalColumn::NotAvailable,
    z: "charge",
};

/// v20 or v21
pub const FRAGPIPE_V20_OR_21: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::FragPipeV20Or21,
    best_score_with_delta_mass: OptionalColumn::NotAvailable,
    calibrated_experimental_mass: OptionalColumn::Required("calibrated observed mass"),
    calibrated_experimental_mz: OptionalColumn::Required("calibrated observed m/z"),
    condition: OptionalColumn::Optional("condition"),
    delta_score: OptionalColumn::NotAvailable,
    entry_name: OptionalColumn::Required("entry name"),
    enzymatic_termini: OptionalColumn::Required("number of enzymatic termini"),
    expectation_score: "expectation",
    extended_peptide: OptionalColumn::Optional("extended peptide"),
    gene: OptionalColumn::Required("gene"),
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyper_score: "hyperscore",
    intensity: OptionalColumn::Required("intensity"),
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    ions_best_position: OptionalColumn::NotAvailable,
    is_unique: OptionalColumn::Required("is unique"),
    mapped_genes: OptionalColumn::Required("mapped genes"),
    mapped_proteins: OptionalColumn::Required("mapped proteins"),
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    modifications: "assigned modifications",
    open_search_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: OptionalColumn::Required("observed m/z"),
    next_score: "nextscore",
    num_matched_ions: OptionalColumn::NotAvailable,
    num_tol_term: OptionalColumn::NotAvailable,
    open_search_modification: OptionalColumn::Required("observed modifications"),
    peptide_prophet_probability: OptionalColumn::Required("peptideprophet probability"),
    peptide: "peptide",
    open_search_position_scores: OptionalColumn::NotAvailable,
    protein_description: OptionalColumn::Required("protein description"),
    protein_end: OptionalColumn::Required("protein end"),
    protein_id: OptionalColumn::Required("protein id"),
    protein_start: OptionalColumn::Required("protein start"),
    protein: "protein",
    purity: OptionalColumn::Required("purity"),
    rank: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    score_without_delta_mass: OptionalColumn::NotAvailable,
    second_best_score_with_delta_mass: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Required("spectrum file"),
    tot_num_ions: OptionalColumn::NotAvailable,
    total_glycan_composition: OptionalColumn::NotAvailable,
    z: "charge",
};

/// v22
pub const FRAGPIPE_V22: MSFraggerFormat = MSFraggerFormat {
    version: MSFraggerVersion::FragPipeV22,
    best_score_with_delta_mass: OptionalColumn::NotAvailable,
    calibrated_experimental_mass: OptionalColumn::Required("calibrated observed mass"),
    calibrated_experimental_mz: OptionalColumn::Required("calibrated observed m/z"),
    condition: OptionalColumn::Optional("condition"),
    delta_score: OptionalColumn::NotAvailable,
    entry_name: OptionalColumn::Required("entry name"),
    enzymatic_termini: OptionalColumn::Required("number of enzymatic termini"),
    expectation_score: "expectation",
    extended_peptide: OptionalColumn::Required("extended peptide"),
    gene: OptionalColumn::Required("gene"),
    glycan_q_value: OptionalColumn::Optional("glycan q-value"),
    glycan_score: OptionalColumn::Optional("glycan score"),
    group: OptionalColumn::Optional("group"),
    hyper_score: "hyperscore",
    intensity: OptionalColumn::Required("intensity"),
    ion_mobility: OptionalColumn::Optional("ion mobility"),
    ions_best_position: OptionalColumn::Optional("ions best position"),
    is_unique: OptionalColumn::Required("is unique"),
    mapped_genes: OptionalColumn::Required("mapped genes"),
    mapped_proteins: OptionalColumn::Required("mapped proteins"),
    mass: "observed mass",
    missed_cleavages: "number of missed cleavages",
    modifications: "assigned modifications",
    open_search_localisation: OptionalColumn::Optional("msfragger localization"),
    mz: OptionalColumn::Required("observed m/z"),
    next_score: "nextscore",
    num_matched_ions: OptionalColumn::NotAvailable,
    num_tol_term: OptionalColumn::NotAvailable,
    open_search_modification: OptionalColumn::Required("observed modifications"),
    peptide_prophet_probability: OptionalColumn::Required("probability"),
    peptide: "peptide",
    open_search_position_scores: OptionalColumn::Optional("position scores"),
    protein_description: OptionalColumn::Required("protein description"),
    protein_end: OptionalColumn::Required("protein end"),
    protein_id: OptionalColumn::Required("protein id"),
    protein_start: OptionalColumn::Required("protein start"),
    protein: "protein",
    purity: OptionalColumn::Required("purity"),
    rank: OptionalColumn::NotAvailable,
    rt: "retention",
    scan: "spectrum",
    score_without_delta_mass: OptionalColumn::NotAvailable,
    second_best_score_with_delta_mass: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Required("spectrum file"),
    tot_num_ions: OptionalColumn::NotAvailable,
    total_glycan_composition: OptionalColumn::Optional("total glycan composition"),
    z: "charge",
};
