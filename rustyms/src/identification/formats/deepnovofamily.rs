use std::sync::LazyLock;

use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    identification::{
        BoxedIdentifiedPeptideIter, IdentifiedPeptidoform, IdentifiedPeptidoformSource,
        IdentifiedPeptidoformVersion, MetaData, PeaksFamilyId,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::{CustomDatabase, Ontology},
    sequence::{AminoAcid, Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{MassOverCharge, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid DeepNovoFamily line",
    "This column is not a number but it is required to be a number in this format",
);
static ID_ERROR: (&str, &str) = (
    "Invalid DeepNovoFamily line",
    "This column is not a valid ID but it is required to be in this peaks format\nExamples of valid IDs: '1234' & 'F2:1234'",
);

static PARAMETERS: LazyLock<SloppyParsingParameters> = LazyLock::new(|| SloppyParsingParameters {
    mod_indications: (
        Some("mod"),
        vec![
            (
                AminoAcid::Asparagine,
                Ontology::Unimod.find_id(7, None).unwrap(),
            ),
            (
                AminoAcid::Glutamine,
                Ontology::Unimod.find_id(7, None).unwrap(),
            ),
            (
                AminoAcid::Cysteine,
                Ontology::Unimod.find_id(6, None).unwrap(),
            ),
            (
                AminoAcid::Methionine,
                Ontology::Unimod.find_id(35, None).unwrap(),
            ),
        ],
    ),
    ..Default::default()
});

format_family!(
    /// The format for any DeepNovoFamily file
    DeepNovoFamilyFormat,
    /// The data from any DeepNovoFamily file
    DeepNovoFamilyData,
    DeepNovoFamilyVersion, [&DEEPNOVO_V0_0_1, &POINTNOVOFAMILY], b'\t', None;
    required {
        scan: Vec<PeaksFamilyId>, |location: Location, _| location.or_empty()
            .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());

        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>|
                location.or_empty().map(|location| Peptidoform::sloppy_pro_forma(
                    location.full_line(),
                    location.location.clone(),
                    custom_database,
                    &PARAMETERS
                )).transpose();
        score: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        local_confidence: Option<Vec<f64>>, |location: Location, _| location.or_empty()
                .optional_array(',').map(|array| array.map(|l| l.parse::<f64>(NUMBER_ERROR)).collect::<Result<Vec<_>, _>>())
                .transpose();
    }
    optional {
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<crate::system::e>);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        if parsed.local_confidence.as_ref().map(Vec::len)
            != parsed.peptide.as_ref().map(Peptidoform::len)
        {
            parsed.local_confidence = parsed.local_confidence.map(interpolate_lc);
        }
        Ok(parsed)
    }
);

impl From<DeepNovoFamilyData> for IdentifiedPeptidoform {
    fn from(value: DeepNovoFamilyData) -> Self {
        Self {
            score: value.score.map(|score| (2.0 / (1.0 + (-score).exp()))),
            local_confidence: value
                .local_confidence
                .as_ref()
                .map(|lc| lc.iter().map(|v| 2.0 / (1.0 + (-v).exp())).collect()),
            metadata: MetaData::DeepNovoFamily(value),
        }
    }
}

/// Interpolate the local confidence when the confidence between AAs is used instead of the confidence of a single AA
#[expect(clippy::needless_pass_by_value)] // The return value will replace the given value, so moving is fine
fn interpolate_lc(local_confidence: Vec<f64>) -> Vec<f64> {
    let mut reinterpolated = Vec::with_capacity(local_confidence.len() + 1);

    for i in 0..local_confidence.len() {
        if i == 0 {
            reinterpolated.push(local_confidence[i]);
        } else {
            let average = (local_confidence[i - 1] + local_confidence[i]) / 2.0;
            reinterpolated.push(average);
        }
    }
    reinterpolated.push(local_confidence[local_confidence.len() - 1]);

    reinterpolated
}

/// The only known version of DeepNovo
pub const DEEPNOVO_V0_0_1: DeepNovoFamilyFormat = DeepNovoFamilyFormat {
    version: DeepNovoFamilyVersion::DeepNovoV0_0_1,
    scan: "scan",
    peptide: "predicted_sequence",
    score: "predicted_score",
    local_confidence: "predicted_position_score",
    mz: OptionalColumn::NotAvailable,
    z: OptionalColumn::NotAvailable,
};

/// The only known version of the PointNovo Family
pub const POINTNOVOFAMILY: DeepNovoFamilyFormat = DeepNovoFamilyFormat {
    version: DeepNovoFamilyVersion::PointNovoFamily,
    scan: "scan_list_original",
    peptide: "predicted_sequence",
    score: "predicted_score",
    local_confidence: "predicted_position_score",
    mz: OptionalColumn::Required("precursor_mz"),
    z: OptionalColumn::Required("precursor_charge"),
};

/// All possible DeepNovoFamily versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum DeepNovoFamilyVersion {
    #[default]
    /// DeepNovo version 0.0.1
    DeepNovoV0_0_1,
    /// PointNovo version 0.0.1 & PGPointNovo version 1.0.6 & BiatNovo version 0.1
    PointNovoFamily,
}

impl std::fmt::Display for DeepNovoFamilyVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<DeepNovoFamilyFormat> for DeepNovoFamilyVersion {
    fn format(self) -> DeepNovoFamilyFormat {
        match self {
            Self::DeepNovoV0_0_1 => DEEPNOVO_V0_0_1,
            Self::PointNovoFamily => POINTNOVOFAMILY,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::DeepNovoV0_0_1 => "DeepNovo v0.0.1",
            Self::PointNovoFamily => "PointNovo v0.0.1 / PGPointNovo v1.0.6 / BiatNovo v0.1",
        }
    }
}
