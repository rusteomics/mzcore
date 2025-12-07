use std::{borrow::Cow, marker::PhantomData, ops::Range, str::FromStr};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, KnownFileFormat, MaybePeptidoform,
    MetaData, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalLocation},
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    prelude::{AminoAcid, SequencePosition},
    sequence::{
        CompoundPeptidoformIon, FlankingSequence, Modification, Peptidoform, PlacementRule,
        Position, SemiAmbiguous, SimpleModification, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid pUniFind line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PUniFind,
    SemiAmbiguous, MaybePeptidoform, [&PUNIFIND_V0_1], b',', Some(vec![
        "spectrum_title".to_string(),
        "score".to_string(),
        "cos_similarity".to_string(),
        "rt".to_string(),
        "missing fragment ion site".to_string(),
        "mass_error".to_string(),
        "peptide".to_string(),
        "peptidoform".to_string(),
        "modifications".to_string()
    ]);

    required {
        title: String, |location: Location, _| Ok(location.get_string());
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        cos_similarity: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::s>);
        mass_error: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::mass::dalton>);
        peptidoform: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies|
            location.trim_start_matches("['").trim_end_matches("']").parse_with(|location| Peptidoform::sloppy_pro_forma(
                location.full_line(),
                location.location.clone(),
                ontologies,
                &SloppyParsingParameters::default(),
            ).map_err(BoxedError::to_owned));

        modifications: Vec<(usize, SimpleModification, PlacementRule)>, |location: Location, ontologies: &Ontologies| location.or_empty().array(';').filter_map(
            |location| location.clone().or_empty().map(|l| l.split_once(',').and_then(|(loc, modification)| modification.split_once('[').map(|(a, b)| (loc, a, b))).map(|(loc, modification, rule)|
            Ok((
                loc.parse::<usize>(NUMBER_ERROR).map_err(BoxedError::to_owned)?,
                Modification::sloppy_modification(modification.line.line(), modification.location, None, ontologies).map_err(BoxedError::to_owned)?,
                rule.trim_end_matches("]").parse_with(|rule|
                    if let Result::Ok(position) = Position::from_str(rule.as_str()) && position != Position::Anywhere {
                        Ok(PlacementRule::Terminal(position))
                    } else if let Result::Ok(aa) = AminoAcid::from_str(rule.as_str()) {
                        Ok(PlacementRule::AminoAcid(vec![aa].into(), Position::Anywhere))
                    } else {
                        Err(BoxedError::new(BasicKind::Error, "Invalid pUniFind line", "The modification location should be either an amino acid or a terminal location", rule.context().to_owned()))
                    }
                ).map_err(BoxedError::to_owned)?,
            )))
            .ok_or_else(|| BoxedError::new(BasicKind::Error, "Invalid pUniFind line", "A pUniFind modification should be formatted like so: 'index,modification[location]'", location.context().to_owned()))?))
            .collect::<Result<Vec<(usize, SimpleModification, PlacementRule)>, BoxedError<'static, BasicKind>>>();
    }
    optional { }

    fn post_process(source: &CsvLine, mut parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        for (index, m, rule) in &parsed.modifications {
            match rule {
                PlacementRule::Terminal(Position::AnyNTerm | Position::ProteinNTerm) => parsed.peptidoform.add_simple_n_term(m.clone()),
                PlacementRule::Terminal(Position::AnyCTerm | Position::ProteinCTerm) => parsed.peptidoform.add_simple_c_term(m.clone()),
                PlacementRule::AminoAcid(list, _) => {
                    if list.len() == 1 && parsed.peptidoform.sequence().get(index - 1).is_none_or(|v| v.aminoacid.aminoacid() != list[0]) {
                        return Err(BoxedError::new(BasicKind::Error, "Invalid pUniFind line", format!("The modification location does not correspond with the listed aminoacid, for modification '{m}'"), source.full_context().to_owned()))
                    }
                    parsed.peptidoform.add_simple_modification(SequencePosition::Index(index-1), m.clone());
                },
                _ => unreachable!()
            }
        }
        Ok(parsed)
    }
);

/// The only version of pUniFind
pub const PUNIFIND_V0_1: PUniFindFormat = PUniFindFormat {
    version: PUniFindVersion::V0_1,
    title: "spectrum_title",
    score: "score",
    cos_similarity: "cos_similarity",
    rt: "rt",
    mass_error: "mass_error",
    peptidoform: "peptide",
    modifications: "modifications",
};

/// All possible pUniFind versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PUniFindVersion {
    #[default]
    /// pUniFind version 0.1 (compatible with 0.2)
    V0_1,
}

impl std::fmt::Display for PUniFindVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PUniFindFormat> for PUniFindVersion {
    fn format(self) -> PUniFindFormat {
        match self {
            Self::V0_1 => PUNIFIND_V0_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V0_1 => "v0.1",
        }
    }
}

impl MetaData for PUniFindData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptidoform.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PUniFind(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> String {
        self.title.clone()
    }

    fn confidence(&self) -> Option<f64> {
        // Some(self.score) // TODO: recalibrate
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        None
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        Some(self.rt)
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Native(self.title.clone())])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
    }

    fn ppm_error(&self) -> Option<Ratio> {
        None
    }

    fn mass_error(&self) -> Option<Mass> {
        Some(self.mass_error)
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
}
