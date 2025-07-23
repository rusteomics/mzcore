use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
    sync::LazyLock,
};

use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn},
        csv::{CsvLine, parse_csv},
    },
    ontology::{CustomDatabase, Ontology},
    prelude::CompoundPeptidoformIon,
    sequence::{Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid InstaNovo line",
    "This column is not a number but it is required to be a number in this format",
);

static BUILT_IN_MODIFICATIONS: LazyLock<SloppyParsingParameters> =
    LazyLock::new(|| SloppyParsingParameters {
        replace_mass_modifications: Some(vec![
            Ontology::Unimod.find_id(35, None).unwrap(),
            Ontology::Unimod.find_id(21, None).unwrap(),
            Ontology::Unimod.find_id(4, None).unwrap(),
        ]),
        ..Default::default()
    });

format_family!(
    InstaNovo,
    SemiAmbiguous, PeptidoformPresent, [&INSTANOVO_V1_0_0, &INSTANOVO_V1_1_4, &INSTANOVOPLUS_V1_1_4], b',', None;
    required {
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &BUILT_IN_MODIFICATIONS);

        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
    optional {
        local_confidence: Vec<f64>, |location: Location, _| location
            .trim_start_matches("[").trim_end_matches("]")
            .array(',')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
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
};

// TODO: Potentially add the other columns here as well?
/// The only known version of InstaNovoPlus
pub const INSTANOVOPLUS_V1_1_4: InstaNovoFormat = InstaNovoFormat {
    version: InstaNovoVersion::PlusV1_1_4,
    scan_number: "scan_number",
    mz: "precursor_mz",
    z: "precursor_charge",
    raw_file: "experiment_name",
    peptide: "final_prediction",
    score: "final_log_probabilities",
    local_confidence: OptionalColumn::NotAvailable,
};

/// All possible InstaNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum InstaNovoVersion {
    #[default]
    /// InstaNovo version 1.0.0
    V1_0_0,
    /// InstaNovo version 1.1.4 without refinement
    V1_1_4,
    /// InstaNovoPlus version 1.1.4 using refinement
    PlusV1_1_4,
}

impl std::fmt::Display for InstaNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<InstaNovoFormat> for InstaNovoVersion {
    fn format(self) -> InstaNovoFormat {
        match self {
            Self::V1_0_0 => INSTANOVO_V1_0_0,
            Self::V1_1_4 => INSTANOVO_V1_1_4,
            Self::PlusV1_1_4 => INSTANOVOPLUS_V1_1_4,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_0_0 => "v1.0.0",
            Self::V1_1_4 => "v1.1.4",
            Self::PlusV1_1_4 => "Plus v1.1.4",
        }
    }
}

impl MetaData for InstaNovoData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::InstaNovo(self.version)
    }

    fn id(&self) -> String {
        self.scan_number.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some(2.0 / (1.0 + 1.01_f64.powf(-self.score)))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| lc.iter().map(|v| 2.0 / (1.0 + 1.25_f64.powf(-v))).collect())
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        self.local_confidence.as_deref()
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<&str> {
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

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        None
    }
}
