use std::{borrow::Cow, marker::PhantomData, ops::Range, sync::OnceLock};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, KnownFileFormat, MaybePeptidoform,
    MetaData, SpectrumId, SpectrumIds, common_parser::Location,
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SequenceElement, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid NovoB line",
    "This column is not a number but it is required to be a number in this format",
);

/// Global parsing parameters
static PARAMETERS: OnceLock<SloppyParsingParameters> = OnceLock::new();

#[expect(
    clippy::missing_panics_doc,
    reason = "These mods are assumed to be in the ontology"
)]
fn parameters(ontologies: &Ontologies) -> &SloppyParsingParameters {
    PARAMETERS.get_or_init(|| SloppyParsingParameters {
        custom_alphabet: vec![
            (
                b's',
                SequenceElement::new(AminoAcid::Serine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&21).unwrap()),
            ),
            (
                b't',
                SequenceElement::new(AminoAcid::Tyrosine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&21).unwrap()),
            ),
            (
                b'y',
                SequenceElement::new(AminoAcid::Threonine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&21).unwrap()),
            ),
            (
                b'n',
                SequenceElement::new(AminoAcid::Asparagine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&7).unwrap()),
            ),
            (
                b'q',
                SequenceElement::new(AminoAcid::Glutamine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&7).unwrap()),
            ),
            (
                b'C',
                SequenceElement::new(AminoAcid::Cysteine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&6).unwrap()),
            ),
            (
                b'm',
                SequenceElement::new(AminoAcid::Methionine.into(), None)
                    .with_simple_modification(ontologies.unimod().get_by_index(&35).unwrap()),
            ),
        ],
        ..Default::default()
    })
}

format_family!(
    NovoB,
    SemiAmbiguous, MaybePeptidoform, [&NOVOB_V0_0_1], b'\t', Some(vec![
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
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);

        score_forward: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ppm_diff_forward: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<mzcore::system::ratio::ppm>);
        peptide_forward: Option<Peptidoform<SemiAmbiguous>>, |location: Location, ontologies: &Ontologies|
            location.trim_start_matches("['").trim_end_matches("']").or_empty().map(|location| Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                ontologies,
                parameters(ontologies)
            ).map_err(BoxedError::to_owned)).transpose();

        score_reverse: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ppm_diff_reverse: Ratio, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Ratio::new::<mzcore::system::ratio::ppm>);
        peptide_reverse: Option<Peptidoform<SemiAmbiguous>>, | location: Location, ontologies: &Ontologies|
            location.trim_start_matches("['").trim_end_matches("']").or_empty().map(|location| Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                ontologies,
                parameters(ontologies),
            ).map_err(BoxedError::to_owned)).transpose();
    }
    optional { }
);

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
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum NovoBVersion {
    #[default]
    /// NovoB version 0.0.1
    V0_0_1,
}

impl std::fmt::Display for NovoBVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<NovoBFormat> for NovoBVersion {
    fn format(self) -> NovoBFormat {
        match self {
            Self::V0_0_1 => NOVOB_V0_0_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V0_0_1 => "v0.0.1",
        }
    }
}

impl MetaData for NovoBData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        if self.score_forward >= self.score_reverse {
            self.peptide_forward.as_ref()
        } else {
            self.peptide_reverse.as_ref()
        }
        .map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::NovoB(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.scan)
    }

    fn id(&self) -> String {
        self.scan.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score_forward.max(self.score_reverse))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((
            self.score_forward.max(self.score_reverse),
            mzcv::term!(MS:1001153|search engine specific score),
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
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Index(self.scan)])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(MassOverCharge::new::<mzcore::system::thomson>(
            self.mass.value / self.z.to_float().value,
        ))
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mass)
    }

    fn ppm_error(&self) -> Option<Ratio> {
        Some(if self.score_forward >= self.score_reverse {
            self.ppm_diff_forward
        } else {
            self.ppm_diff_reverse
        })
    }

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
