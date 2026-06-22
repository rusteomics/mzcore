use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
    sync::OnceLock,
};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, KnownFileFormat, PSM, PSMData, PSMFileFormatVersion, PSMMetaData,
    PSMSource, PeptidoformPresent, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn},
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        FlankingSequence, Peptidoform, PeptidoformIonSet, SemiAmbiguous, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid InstaNovo line",
    "This column is not a number but it is required to be a number in this format",
);

static BUILT_IN_MODIFICATIONS: OnceLock<SloppyParsingParameters> = OnceLock::new();

format_family!(
    InstaNovo,
    SemiAmbiguous, PeptidoformPresent, [
        &INSTANOVO_COMBINED_V1_2_2,
        &INSTANOVO_V1_2_2,
        &INSTANOVOPLUS_V1_2_2,
        &INSTANOVO_V1_1_0,
        &INSTANOVO_V1_1_4,
        &INSTANOVOPLUS_V1_1_4,
        &INSTANOVO_V1_0_0,
    ], b',', None;
    required {
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies| Peptidoform::sloppy_pro_forma_inner(
            &location.base_context(),
            location.full_line(),
            location.range.clone(),
            ontologies,
            BUILT_IN_MODIFICATIONS.get_or_init(|| SloppyParsingParameters {
                replace_mass_modifications: Some(vec![
                ontologies.unimod().get_by_index(&mzcv::AccessionCode::Numeric(35)).unwrap(),
                ontologies.unimod().get_by_index(&mzcv::AccessionCode::Numeric(21)).unwrap(),
                ontologies.unimod().get_by_index(&mzcv::AccessionCode::Numeric(4)).unwrap(),
                ]),
                ..Default::default()
            })).map_err(BoxedError::to_owned);

        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
    optional {
        local_confidence: Vec<f64>, |location: Location, _| {
            let location = location.trim_start_matches("[").trim_end_matches("]");
            location.or_empty().map_or(Ok(Vec::new()), |location| {
                location
                    .array(',')
                    .map(|l| l.parse::<f64>(NUMBER_ERROR))
                    .collect::<Result<Vec<_>, _>>()
            })
        };
        used_model: UsedModel, |location: Location, _| location.parse::<UsedModel>(("Invalid InstaNovo line", "The selected model has to be 'diffusion' or 'transformer'."));
     }

     fn post_process(source: &CsvLine, mut parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        validate_instanovo_schema(source, &parsed)?;

        if parsed.local_confidence.as_ref().is_some_and(Vec::is_empty) {
            parsed.local_confidence = None;
        }
        // Only keep the parsed local_confidence is the `UsedModel == Transformer`
        if let Some(used_model) = parsed.used_model && used_model == UsedModel::Diffusion {
            parsed.local_confidence = None;
        }
        if let Some(local_confidence) = parsed.local_confidence.as_mut() && !parsed.peptide.get_n_term().is_empty() {
            let offset = parsed.peptide.get_n_term().len();
            if local_confidence.len() >= offset {
                *local_confidence = local_confidence[offset..].to_vec();
            }
        }
        Ok(parsed)
    }
);

/// InstaNovo version 1.0.0
pub const INSTANOVO_V1_0_0: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::V1_0_0,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "preds",
    score: "log_probs",
    local_confidence: OptionalColumn::Required("token_log_probs"),
    used_model: OptionalColumn::NotAvailable,
};

/// InstaNovo version 1.1.0
pub const INSTANOVO_V1_1_0: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::V1_1_0,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "predictions",
    score: "log_probabilities",
    local_confidence: OptionalColumn::Required("token_log_probabilities"),
    used_model: OptionalColumn::NotAvailable,
};

/// InstaNovo version 1.1.4
pub const INSTANOVO_V1_1_4: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::V1_1_4,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "preds",
    score: "log_probs",
    local_confidence: OptionalColumn::Required("token_log_probs"),
    used_model: OptionalColumn::NotAvailable,
};

/// The known InstaNovoPlus 1.1.4 output schema
pub const INSTANOVOPLUS_V1_1_4: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::PlusV1_1_4,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "final_prediction",
    score: "final_log_probabilities",
    local_confidence: OptionalColumn::Optional("transformer_token_log_probabilities"),
    used_model: OptionalColumn::Required("selected_model"),
};

/// InstaNovo version 1.2.2 transformer output
pub const INSTANOVO_V1_2_2: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::V1_2_2,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "predictions",
    score: "log_probs",
    local_confidence: OptionalColumn::Required("token_log_probs"),
    used_model: OptionalColumn::NotAvailable,
};

/// InstaNovoPlus version 1.2.2 standalone output
pub const INSTANOVOPLUS_V1_2_2: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::PlusV1_2_2,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "predictions",
    score: "log_probs",
    local_confidence: OptionalColumn::Required("token_log_probs"),
    used_model: OptionalColumn::NotAvailable,
};

/// InstaNovo version 1.2.2 combined transformer and InstaNovoPlus refined output
pub const INSTANOVO_COMBINED_V1_2_2: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::CombinedV1_2_2,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "predictions",
    score: "log_probs",
    local_confidence: OptionalColumn::Required("token_log_probs"),
    used_model: OptionalColumn::NotAvailable,
};

/// All possible InstaNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum InstaNovoVersion {
    #[default]
    /// InstaNovo version 1.0.0
    V1_0_0,
    /// InstaNovo version 1.1.0
    V1_1_0,
    /// InstaNovo version 1.1.4
    V1_1_4,
    /// InstaNovoPlus version 1.1.4 using refinement
    PlusV1_1_4,
    /// InstaNovo version 1.2.2
    V1_2_2,
    /// InstaNovoPlus version 1.2.2 standalone predictions
    PlusV1_2_2,
    /// InstaNovo version 1.2.2 combined transformer and InstaNovoPlus refined predictions
    CombinedV1_2_2,
}

impl std::fmt::Display for InstaNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<InstaNovoFormat> for InstaNovoVersion {
    fn format(self) -> InstaNovoFormat {
        match self {
            Self::V1_0_0 => INSTANOVO_V1_0_0,
            Self::V1_1_0 => INSTANOVO_V1_1_0,
            Self::V1_1_4 => INSTANOVO_V1_1_4,
            Self::PlusV1_1_4 => INSTANOVOPLUS_V1_1_4,
            Self::V1_2_2 => INSTANOVO_V1_2_2,
            Self::PlusV1_2_2 => INSTANOVOPLUS_V1_2_2,
            Self::CombinedV1_2_2 => INSTANOVO_COMBINED_V1_2_2,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_0_0 => "v1.0.0",
            Self::V1_1_0 => "v1.1.0",
            Self::V1_1_4 => "v1.1.4",
            Self::PlusV1_1_4 => "Plus v1.1.4",
            Self::V1_2_2 => "v1.2.2",
            Self::PlusV1_2_2 => "Plus v1.2.2",
            Self::CombinedV1_2_2 => "Combined v1.2.2",
        }
    }
}

fn validate_instanovo_schema(
    source: &CsvLine,
    parsed: &InstaNovoPSM,
) -> Result<(), BoxedError<'static, BasicKind>> {
    match parsed.version {
        InstaNovoVersion::V1_1_0 | InstaNovoVersion::V1_1_4 => {
            if !has_column(source, "delta_mass_ppm") {
                return Err(instanovo_schema_error(
                    source,
                    "This InstaNovo version requires the 'delta_mass_ppm' column",
                ));
            }
        }
        InstaNovoVersion::V1_2_2 => {
            if has_column(source, "instanovoplus_predictions") {
                return Err(instanovo_schema_error(
                    source,
                    "This is an InstaNovo combined output, not a transformer-only output",
                ));
            }
            if parsed.local_confidence.as_ref().is_some_and(Vec::is_empty) {
                return Err(instanovo_schema_error(
                    source,
                    "This InstaNovo transformer output requires token log probabilities",
                ));
            }
        }
        InstaNovoVersion::PlusV1_2_2 => {
            if has_column(source, "instanovoplus_predictions") {
                return Err(instanovo_schema_error(
                    source,
                    "This is an InstaNovo combined output, not a standalone InstaNovoPlus output",
                ));
            }
            if parsed
                .local_confidence
                .as_ref()
                .is_some_and(|local_confidence| !local_confidence.is_empty())
            {
                return Err(instanovo_schema_error(
                    source,
                    "This is an InstaNovo transformer output, not a standalone InstaNovoPlus output",
                ));
            }
        }
        InstaNovoVersion::CombinedV1_2_2 => {
            if !has_column(source, "instanovoplus_predictions") {
                return Err(instanovo_schema_error(
                    source,
                    "This InstaNovo version requires the 'instanovoplus_predictions' column",
                ));
            }
        }
        InstaNovoVersion::V1_0_0 | InstaNovoVersion::PlusV1_1_4 => {}
    }
    Ok(())
}

fn has_column(source: &CsvLine, column: &str) -> bool {
    source
        .fields
        .iter()
        .any(|f| f.0.eq_ignore_ascii_case(column))
}

fn instanovo_schema_error(
    source: &CsvLine,
    message: &'static str,
) -> BoxedError<'static, BasicKind> {
    BoxedError::new(
        BasicKind::Error,
        "Invalid InstaNovo line",
        message,
        source.full_context().to_owned(),
    )
}

/// The model that produced the final prediction for an InstaNovoPlus
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum UsedModel {
    /// The diffusion model
    Diffusion,
    /// The transformer model
    Transformer,
}

impl mzcore::space::Space for UsedModel {
    fn space(&self) -> mzcore::space::UsedSpace {
        mzcore::space::UsedSpace::stack(size_of::<Self>())
    }
}

impl std::fmt::Display for UsedModel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Diffusion => "diffusion",
                Self::Transformer => "transformer",
            }
        )
    }
}

impl std::str::FromStr for UsedModel {
    type Err = ();
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        if value.eq_ignore_ascii_case("diffusion") {
            Ok(Self::Diffusion)
        } else if value.eq_ignore_ascii_case("transformer") {
            Ok(Self::Transformer)
        } else {
            Err(())
        }
    }
}

impl PSMMetaData for InstaNovoPSM {
    fn peptidoform_ion_set(&self) -> Option<Cow<'_, PeptidoformIonSet>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::InstaNovo(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.scan_number)
    }

    fn id(&self) -> String {
        self.scan_number.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        Some(match self.version {
            InstaNovoVersion::V1_0_0
            | InstaNovoVersion::V1_1_0
            | InstaNovoVersion::V1_1_4
            | InstaNovoVersion::V1_2_2 => mzcv::term!(MS:1003612|InstaNovo),
            InstaNovoVersion::PlusV1_1_4
            | InstaNovoVersion::PlusV1_2_2
            | InstaNovoVersion::CombinedV1_2_2 => mzcv::term!(MS:1003613|InstaNovo+),
        })
    }

    fn confidence(&self) -> Option<f64> {
        Some(2.0 / (1.0 + 1.01_f64.powf(-self.score)))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| lc.iter().map(|v| 2.0 / (1.0 + 1.25_f64.powf(-v))).collect())
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.score,
            mzcv::term!(MS:1001153|search engine specific score),
        ))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        self.local_confidence.as_deref()
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileKnown(vec![(
            self.raw_file.clone(),
            vec![SpectrumId::Number(self.scan_number)],
        )])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mz * self.z.to_float())
    }

    type Protein = crate::NoProtein;

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn unique(&self) -> Option<bool> {
        None
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}
