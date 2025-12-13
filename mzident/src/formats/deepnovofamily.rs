use std::{borrow::Cow, marker::PhantomData, ops::Range, sync::OnceLock};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, KnownFileFormat, MaybePeptidoform,
    PSMMetaData, PeaksFamilyId, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid DeepNovoFamily line",
    "This column is not a number but it is required to be a number in this format",
);
static ID_ERROR: (&str, &str) = (
    "Invalid DeepNovoFamily line",
    "This column is not a valid ID but it is required to be in this peaks format\nExamples of valid IDs: '1234' & 'F2:1234'",
);

static PARAMETERS: OnceLock<SloppyParsingParameters> = OnceLock::new();

format_family!(
    DeepNovoFamily,
    SemiAmbiguous, MaybePeptidoform, [&DEEPNOVO_V0_0_1, &POINTNOVOFAMILY], b'\t', None;
    required {
        scan: Vec<PeaksFamilyId>, |location: Location, _| location.or_empty()
            .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());

        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, ontologies: &Ontologies|
                location.or_empty().map(|location| Peptidoform::sloppy_pro_forma(
                    location.full_line(),
                    location.location.clone(),
                    ontologies,
                    PARAMETERS.get_or_init(|| SloppyParsingParameters {
                        mod_indications: (
                            Some("mod"),
                            vec![
                                (
                                    AminoAcid::Asparagine,
                                    ontologies.unimod().get_by_index(&7).unwrap(),
                                ),
                                (
                                    AminoAcid::Glutamine,
                                    ontologies.unimod().get_by_index(&7).unwrap(),
                                ),
                                (
                                    AminoAcid::Cysteine,
                                    ontologies.unimod().get_by_index(&6).unwrap(),
                                ),
                                (
                                    AminoAcid::Methionine,
                                    ontologies.unimod().get_by_index(&3).unwrap(),
                                ),
                            ],
                        ),
                        ..Default::default()
                    })
                ).map_err(BoxedError::to_owned)).transpose();
        score: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        local_confidence: Option<Vec<f64>>, |location: Location, _| location.or_empty()
                .optional_array(',').map(|array| array.map(|l| l.parse::<f64>(NUMBER_ERROR)).collect::<Result<Vec<_>, _>>())
                .transpose();
    }
    optional {
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<mzcore::system::e>);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        if parsed.local_confidence.as_ref().map(Vec::len)
            != parsed.peptide.as_ref().map(Peptidoform::len)
        {
            parsed.local_confidence = parsed.local_confidence.map(interpolate_lc);
        }
        Ok(parsed)
    }
);

/// Interpolate the local confidence when the confidence between AAs is used instead of the confidence of a single AA
#[expect(clippy::needless_pass_by_value)] // The return value will replace the given value, so moving is fine
fn interpolate_lc(local_confidence: Vec<f64>) -> Vec<f64> {
    let mut reinterpolated = Vec::with_capacity(local_confidence.len() + 1);

    for i in 0..local_confidence.len() {
        if i == 0 {
            reinterpolated.push(local_confidence[i]);
        } else {
            let average = f64::midpoint(local_confidence[i - 1], local_confidence[i]);
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

impl PSMMetaData for DeepNovoFamilyData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::DeepNovoFamily(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        self.scan.first().and_then(|s| s.scans.first().cloned())
    }

    fn id(&self) -> String {
        self.scan.iter().join(";")
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        self.score.map(|score| 2.0 / (1.0 + (-score).exp()))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| lc.iter().map(|v| 2.0 / (1.0 + (-v).exp())).collect())
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        self.score
            .map(|v| (v, mzcv::term!(MS:1001153|search engine specific score)))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        self.local_confidence.as_deref()
    }

    fn charge(&self) -> Option<Charge> {
        self.z
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(
            self.scan
                .iter()
                .flat_map(|s| s.scans.clone())
                .map(SpectrumId::Number)
                .collect(),
        )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.mz
    }

    fn experimental_mass(&self) -> Option<Mass> {
        self.mz
            .and_then(|mz| self.z.map(|z| (mz, z)).map(|(mz, z)| mz * z.to_float()))
    }

    type Protein<'a> = crate::NoProtein;
    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

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
