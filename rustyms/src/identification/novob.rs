use crate::{
    error::CustomError,
    identification::{IdentifiedPeptide, IdentifiedPeptideSource, MetaData},
    modification::Ontology,
    ontologies::CustomDatabase,
    system::Ratio,
    system::{usize::Charge, Mass},
    AminoAcid, LinearPeptide, SemiAmbiguous, SloppyParsingParameters,
};

use serde::{Deserialize, Serialize};

use super::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, SequenceElement,
};

use std::sync::OnceLock;

static NUMBER_ERROR: (&str, &str) = (
    "Invalid NovoB line",
    "This column is not a number but it is required to be a number in this format",
);

static PARAMETERS_LOCK: OnceLock<SloppyParsingParameters> = OnceLock::new();

/// Global parsing parameters
#[allow(clippy::missing_panics_doc)] // These modifications exist
fn parameters() -> &'static SloppyParsingParameters {
    PARAMETERS_LOCK.get_or_init(|| SloppyParsingParameters {
        custom_alphabet: vec![
            (
                b's',
                SequenceElement::new(AminoAcid::Serine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(21, None).unwrap()),
            ),
            (
                b't',
                SequenceElement::new(AminoAcid::Tyrosine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(21, None).unwrap()),
            ),
            (
                b'y',
                SequenceElement::new(AminoAcid::Threonine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(21, None).unwrap()),
            ),
            (
                b'n',
                SequenceElement::new(AminoAcid::Asparagine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(7, None).unwrap()),
            ),
            (
                b'q',
                SequenceElement::new(AminoAcid::Glutamine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(7, None).unwrap()),
            ),
            (
                b'C',
                SequenceElement::new(AminoAcid::Cysteine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(6, None).unwrap()),
            ),
            (
                b'm',
                SequenceElement::new(AminoAcid::Methionine.into(), None)
                    .with_simple_modification(Ontology::Unimod.find_id(35, None).unwrap()),
            ),
        ],
        ..Default::default()
    })
}

format_family!(
    /// The format for any NovoB file
    NovoBFormat,
    /// The data from any NovoB file
    NovoBData,
    NovoBVersion, [&NOVOB_V0_0_1], b'\t', Some(vec![
        "mcount".to_string(),
        "charge".to_string(),
        "pepmass".to_string(),
        "senten".to_string(),
        "delta_mass".to_string(),
        "prob".to_string(),
        "senten_reverse".to_string(),
        "delta_mass_reverse".to_string(),
        "prob_reverse".to_string()
    ]);

    required {
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);

        score_forward: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ppm_diff_forward: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        peptide_forward: Option<LinearPeptide<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.trim_start_matches("['").trim_end_matches("']").or_empty().map(|location| LinearPeptide::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                custom_database,
                parameters()
            )).transpose();

        score_reverse: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ppm_diff_reverse: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        peptide_reverse: Option<LinearPeptide<SemiAmbiguous>>, | location: Location, custom_database: Option<&CustomDatabase>|
            location.trim_start_matches("['").trim_end_matches("']").or_empty().map(|location| LinearPeptide::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                custom_database,
                parameters(),
            )).transpose();
    }
    optional { }
);

impl From<NovoBData> for IdentifiedPeptide {
    fn from(value: NovoBData) -> Self {
        Self {
            score: Some(value.score_forward.max(value.score_reverse)),
            local_confidence: None,
            metadata: MetaData::NovoB(value),
        }
    }
}

/// The only known version of NovoB
pub const NOVOB_V0_0_1: NovoBFormat = NovoBFormat {
    version: NovoBVersion::V0_0_1,
    scan: "mcount",
    z: "charge",
    mass: "pepmass",

    score_forward: "prob",
    peptide_forward: "senten",
    ppm_diff_forward: "delta_mass",

    score_reverse: "prob_reverse",
    peptide_reverse: "senten_reverse",
    ppm_diff_reverse: "delta_mass_reverse",
};

/// All possible NovoB versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum NovoBVersion {
    #[default]
    /// NovoB version 0.0.1
    V0_0_1,
}

impl std::fmt::Display for NovoBVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V0_0_1 => "v0.0.1",
            }
        )
    }
}
