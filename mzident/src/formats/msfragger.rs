use std::{borrow::Cow, marker::PhantomData, ops::Range, path::PathBuf, sync::LazyLock};

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    BoxedIdentifiedPeptideIter, CVTerm, FastaIdentifier, KnownFileFormat, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, PeptidoformPresent, ProteinMetaData, Reliability,
    SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
    helper_functions::explain_number_error,
};
use mzcore::{
    chemistry::Chemical,
    csv::{CsvLine, parse_csv},
    glycan::MonoSaccharide,
    ontology::Ontologies,
    quantities::{Tolerance, WithinTolerance},
    sequence::{
        CompoundPeptidoformIon, FlankingSequence, Modification, Peptidoform, SemiAmbiguous,
        SequencePosition, SimpleLinear, SimpleModification, SimpleModificationInner,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

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
    MSFragger,
    SimpleLinear, PeptidoformPresent, [&VERSION_V4_2, &FRAGPIPE_V20_OR_21, &FRAGPIPE_V22, &PHILOSOPHER], b'\t', None;
    required {
        expectation_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        hyper_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        missed_cleavages: u8, |location: Location, _| location.parse(NUMBER_ERROR);
        next_score: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide: Peptidoform<SimpleLinear>, |location: Location, ontologies: &Ontologies| location.parse_with(|location| {
            Peptidoform::sloppy_pro_forma_inner(
                &location.base_context(),
                location.full_line(),
                location.location.clone(),
                ontologies,
                &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
        ).map(Into::into).map_err(BoxedError::to_owned)});
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        scan: SpectrumId, |location: Location, _| Ok(location.clone().parse::<usize>(NUMBER_ERROR).map_or_else(|_| SpectrumId::Native(location.get_string()), SpectrumId::Number));
        modifications: ThinVec<(SequencePosition, SimpleModification)>, |location: Location, ontologies: &Ontologies| location.or_empty().array(',').map(|m| if let Some((head, tail)) = m.clone().split_once('(') {
            let head_trim = head.as_str().trim();
            Ok((
                if head_trim.eq_ignore_ascii_case("N-term") {
                    SequencePosition::NTerm
                } else if head_trim.eq_ignore_ascii_case("C-term") {
                    SequencePosition::CTerm
                } else {
                    // Format: `14M` so take only the numeric part
                    head.as_str()[..head.len()-1].trim().parse::<usize>().map(|i| SequencePosition::Index(i-1)).map_err(|err| BoxedError::new(BasicKind::Error,
                        "Invalid FragPipe modification location",
                        format!("The location number {}", explain_number_error(&err)),
                        head.context(),
                    )).map_err(BoxedError::to_owned)?
                },
                Modification::sloppy_modification(tail.full_line(), tail.location.clone(), None, ontologies).map_err(BoxedError::to_owned)?
            ))
        } else {
            Err(BoxedError::new(BasicKind::Error,
                "Invalid FragPipe modification",
                "The format `location(modification)` could not be recognised",
                m.context().to_owned(),
            ))
        }).collect::<Result<ThinVec<_>, _>>();
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
    }
    optional {
        best_score_with_delta_mass: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        calibrated_experimental_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        calibrated_experimental_mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        condition: Box<str>, |location: Location, _| Ok(Some(location.get_boxed_str()));
        delta_score: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        entry_name: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        enzymatic_termini: u8, |location: Location, _| location.parse::<u8>(NUMBER_ERROR);
        extended_peptide: (FlankingSequence, Option<Peptidoform<SemiAmbiguous>>, FlankingSequence), |location: Location, ontologies: &Ontologies| {
            let mut peptides = location.clone().array('.').map(|l| l.or_empty().parse_with(|location| Peptidoform::sloppy_pro_forma_inner(
                &location.base_context(),
                location.full_line(),
                location.location.clone(),
                ontologies,
                &SloppyParsingParameters {ignore_prefix_lowercase_n: true, ..Default::default()},
            ).map_err(BoxedError::to_owned))).collect::<Result<Vec<_>,_>>()?;
            if peptides.len() == 3 {
                let suffix = peptides.pop().unwrap().map_or(FlankingSequence::Terminal, |s| FlankingSequence::Sequence(Box::new(s)));
                let peptide = peptides.pop().unwrap();
                let prefix = peptides.pop().unwrap().map_or(FlankingSequence::Terminal, |s| FlankingSequence::Sequence(Box::new(s)));
                Ok((prefix, peptide, suffix))
            } else {
                Err(BoxedError::new(BasicKind::Error,"Invalid extened peptide", "The extended peptide should contain the prefix.peptide.suffix for all peptides.", location.context().to_owned()))
            }
        };
        glycan_q_value: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        glycan_score: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        group: Box<str>, |location: Location, _| Ok(Some(location.get_boxed_str()));
        intensity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        ion_mobility: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        ions_best_position: usize, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
        is_unique: bool, |location: Location, _| location.parse_with(|l| match l.as_str().to_ascii_lowercase().as_str() {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(BoxedError::new(BasicKind::Error,
                "Invalid FragPipe line",
                "This column (Is Unique) is not a boolean but it is required to be a boolean ('true' or 'false') in this FragPipe format",
                l.context().to_owned(),
            ))
        });
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        num_matched_ions: u16, |location: Location, _| location.parse::<u16>(NUMBER_ERROR);
        num_tol_term: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        open_search_localisation: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        /// Either glycan 'HexNAc(4)Hex(5)Fuc(1)NeuAc(2) % 2350.8304' or mod 'Mod1: First isotopic peak (Theoretical: 1.0024)' 'Mod1: Deamidation (PeakApex: 0.9836, Theoretical: 0.9840)'
        open_search_modification: MSFraggerOpenModification, |location: Location, _|
            location.or_empty().parse_with(|location| location.as_str().find('%').map_or_else(
                || Ok(MSFraggerOpenModification::Other(location.as_str().into())),
                |end| MonoSaccharide::byonic_composition(&location.as_str()[..end]).map(|g| MSFraggerOpenModification::Glycan(g.into())).map_err(BoxedError::to_owned)));
        open_search_position_scores: ThinVec<f64>, |location: Location, _| {
            let data = location.array(')').filter_map(|l| (l.len() > 2).then(|| l.skip(2).parse::<f64>((
                "Invalid FragPipe line",
                "This position score is not a number but it is required to be a number in this FragPipe format",
            )))).collect::<Result<ThinVec<_>, _>>()?;
            if data.is_empty() {
                Ok(None)
            } else {
                Ok(Some(data))
            }
        };
        peptide_prophet_probability: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        protein_start: u16, |location: Location, _| location.parse::<u16>(NUMBER_ERROR);
        protein_end: u16, |location: Location, _| location.parse::<u16>(NUMBER_ERROR);
        purity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        rank: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        score_without_delta_mass: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        second_best_score_with_delta_mass: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| Ok(Some(location.get_string().into()));
        total_glycan_composition: ThinVec<(MonoSaccharide, isize)>, |location: Location, _| location.or_empty().parse_with(|location| location.as_str().find('%').map_or_else(
                || MonoSaccharide::byonic_composition(location.as_str()).map(ThinVec::from).map_err(BoxedError::to_owned),
                |end| MonoSaccharide::byonic_composition(&location.as_str()[..end]).map(ThinVec::from).map_err(BoxedError::to_owned)));
        tot_num_ions: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        // Parse the scan identifier to retrieve the file and number separately
        if let SpectrumId::Native(native) = &parsed.scan
            && let Some(m) = IDENTIFER_REGEX
                .captures(native)
            {
                parsed.raw_file = Some(m.get(1).unwrap().as_str().into());
                parsed.scan =
                    SpectrumId::Number(m.get(2).unwrap().as_str().parse::<usize>().unwrap());
            }
        // Normal modifications
        for (location, modification) in &parsed.modifications {
            parsed.peptide.add_simple_modification(*location, modification.clone());
        }
        // Open modification - check for a glycan composition mod
        if let Some(glycan) = match parsed.open_search_modification.as_ref() {Some(MSFraggerOpenModification::Glycan(glycan)) => Some(glycan), _ => None}.or(parsed.total_glycan_composition.as_ref()) {
            let modification: SimpleModification = SimpleModificationInner::Glycan(glycan.to_vec()).into();
            let target_mass = Mass::new::<mzcore::system::dalton>(glycan.iter().fold(0.0, |acc, (sugar, amount)| sugar.formula().monoisotopic_mass().value.mul_add(*amount as f64, acc)));
            let mut placed = false;

            // Check if the location field contains any lowercase chars indication the location
            if let Some(location) = parsed.open_search_localisation.as_ref().map(|l| l.chars().enumerate().filter_map(|(i, c)| c.is_ascii_lowercase().then_some(i)).collect::<Vec<_>>()) {
                match location.len() {
                    0 => (),
                    1 => {
                        // Check if the modification is already placed as mass modification (this is not necessarily always present )
                        let mut index = None;
                        for (i, m) in parsed.peptide[location[0]].modifications.iter().enumerate() {
                            if let Some(SimpleModificationInner::Mass(_, mass, _)) = m.simple().map(AsRef::as_ref)
                                && Tolerance::Absolute(Mass::new::<mzcore::system::dalton>(1.0)).within(&mass.into_inner(), &target_mass) {
                                    index = Some(i);
                                    break;
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
                        if let Some(SimpleModificationInner::Mass(_, mass, _)) = m.simple().map(AsRef::as_ref) && Tolerance::Absolute(Mass::new::<mzcore::system::dalton>(1.0)).within(&mass.into_inner(), &target_mass) {
                            index = Some(i);
                            break;
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

    protein {
        protein_key => (|location: Location, _| Ok(location.get_string()));

        required {
            /// The full header of the fasta entry, split into the id and the rest of the line (without the separating space)
            protein: (FastaIdentifier<Box<str>>, String), |location: Location, _| location.clone().split_once(' ').map_or_else(
            || location.parse(IDENTIFIER_ERROR).map(|id| (id, String::new())),
            |(id, tail)| id.parse(IDENTIFIER_ERROR).map(|id| (id, tail.get_string())));
        }
        optional {
            protein_id: String, |location: Location, _| Ok(location.get_string());
            protein_description: String, |location: Location, _| Ok(location.get_string());
            gene: String, |location: Location, _| Ok(location.get_string());
            mapped_genes: ThinVec<String>, |location: Location, _| Ok(location.or_empty().map(|l| l.get_string().split(',').map(|s| s.trim().to_string()).collect::<ThinVec<_>>()));
            mapped_proteins: ThinVec<String>, |location: Location, _| Ok(location.or_empty().map(|l| l.get_string().split(',').map(|s| s.trim().to_string()).collect::<ThinVec<_>>()));
        }
    }
);

/// The Regex to match against MSFragger scan fields
static IDENTIFER_REGEX: LazyLock<regex::Regex> =
    LazyLock::new(|| regex::Regex::new(r"([^/]+)\.(\d+)\.\d+.\d+").unwrap());

/// A MSFragger open search modification
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum MSFraggerOpenModification {
    /// A glycan composition
    Glycan(ThinVec<(MonoSaccharide, isize)>),
    /// Any other modification
    Other(Box<str>),
}

impl mzcore::space::Space for MSFraggerOpenModification {
    fn space(&self) -> mzcore::space::UsedSpace {
        match self {
            Self::Glycan(m) => m.space(),
            Self::Other(s) => s.space(),
        }
        .set_total::<Self>()
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

impl PSMFileFormatVersion<MSFraggerFormat> for MSFraggerVersion {
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
    protein_key: "proteins",
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
    protein_key: "protein",
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
    protein_key: "protein",
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
    protein_key: "protein",
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

impl PSMMetaData for MSFraggerPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::MSFragger(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        match self.scan {
            SpectrumId::Index(i) | SpectrumId::Number(i) => Some(i),
            _ => None,
        }
    }

    fn id(&self) -> String {
        self.scan.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        Some(mzcv::term!(MS:1003014|MSFragger))
    }

    fn confidence(&self) -> Option<f64> {
        Some(f64::from(self.hyper_score) / 100.0)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.hyper_score.into(),
            mzcv::term!(MS:1001331|X!Tandem:hyperscore),
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
        Some(self.rt)
    }

    fn scans(&self) -> SpectrumIds {
        self.raw_file.clone().map_or_else(
            || SpectrumIds::FileNotKnown(vec![self.scan.clone()]),
            |raw_file| SpectrumIds::FileKnown(vec![(raw_file, vec![self.scan.clone()])]),
        )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.mz
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mass)
    }

    type Protein = MSFraggerProtein;
    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(std::slice::from_ref(&self.protein_key))
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        self.protein_start
            .and_then(|start| self.protein_end.map(|end| start..end))
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        self.extended_peptide.as_ref().map_or(
            (&FlankingSequence::Unknown, &FlankingSequence::Unknown),
            |extended| (&extended.0, &extended.2),
        )
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn unique(&self) -> Option<bool> {
        self.is_unique
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}

impl ProteinMetaData for MSFraggerProtein {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<mzcore::sequence::Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.protein.0.as_str()
    }

    fn description(&self) -> Option<&str> {
        Some(self.protein.1.as_ref())
    }

    fn species(&self) -> Option<mzcv::Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        None
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[]
    }

    fn ambiguity_members(&self) -> &[String] {
        self.mapped_genes
            .as_ref()
            .map(ThinVec::as_slice)
            .unwrap_or_default()
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[mzcv::Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}
