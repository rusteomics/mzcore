use std::{borrow::Cow, marker::PhantomData, ops::Range};

use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{FlankingSequence, Linear, Peptidoform, PeptidoformIonSet, SemiAmbiguous},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};
use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, KnownFileFormat, MaybePeptidoform, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, SpectrumId, SpectrumIds, common_parser::Location,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid psm_utils line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PsmUtils,
    SemiAmbiguous, MaybePeptidoform, [&PSM_UTILS_V1_5], b'\t', None;

    required {
        peptidoform: Peptidoform<Linear>, |location: Location, ontologies: &Ontologies|
            Peptidoform::pro_forma_inner(
                &location.base_context(),
                location.full_line(),
                location.range.clone(),
                ontologies,
            ).map_err(|errs| BoxedError::new(BasicKind::Error, "Invalid psm_utils line", "This ProForma definition contains cross-links, which are not allowed here", location.base_context().add_highlight((0, location.range.clone()))).add_underlying_errors(errs).to_owned())
                .and_then(|(psm, _)| psm.into_linear().ok_or_else(|| BoxedError::new(BasicKind::Error, "Invalid psm_utils line", "This ProForma definition contains cross-links, which are not allowed here", location.base_context().add_highlight((0, location.range.clone())).to_owned())));
        spectrum_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        run: String, |location: Location, _| Ok(location.get_string());
        collection: String, |location: Location, _| Ok(location.get_string());
        is_decoy: bool, |location: Location, _| Ok(location.as_str().eq_ignore_ascii_case("true"));
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        qvalue: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        pep: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::s>);
        ion_mobility: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        protein_list: Vec<String>, |location: Location, _| Ok(location.trim_start_matches("['").trim_end_matches("']").array(',').map(|l| l.trim().trim_start_matches("\'").trim_end_matches("\'").get_string()).collect::<Vec<String>>());
        rank: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        source: String, |location: Location, _| Ok(location.get_string());
    }
    optional { }
);

/// Version 1.5 of `psm_utils`
pub const PSM_UTILS_V1_5: PsmUtilsFormat = PsmUtilsFormat {
    version: PsmUtilsVersion::V1_5,
    peptidoform: "peptidoform",
    spectrum_id: "spectrum_id",
    run: "run",
    collection: "collection",
    is_decoy: "is_decoy",
    score: "score",
    qvalue: "qvalue",
    pep: "pep",
    mz: "precursor_mz",
    rt: "retention_time",
    ion_mobility: "ion_mobility",
    protein_list: "protein_list",
    rank: "rank",
    source: "source",
};

/// All possible pUniFind versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PsmUtilsVersion {
    #[default]
    /// `psm_utils` version 1.5
    V1_5,
}

impl std::fmt::Display for PsmUtilsVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<PsmUtilsFormat> for PsmUtilsVersion {
    fn format(self) -> PsmUtilsFormat {
        match self {
            Self::V1_5 => PSM_UTILS_V1_5,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::V1_5 => "v1.5",
        }
    }
}

impl PSMMetaData for PsmUtilsPSM {
    type Protein = crate::NoProtein;
    #[cfg(feature = "mzannotate")]
    type SpectrumOutputMode = mzcore::chemistry::OutputMolecularFormula;

    fn peptidoform_ion_set(&self) -> Option<Cow<'_, PeptidoformIonSet>> {
        Some(Cow::Owned(self.peptidoform.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PsmUtils(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.spectrum_id)
    }

    fn id(&self) -> String {
        self.spectrum_id.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None // TODO: maybe detect based on the source column
    }

    fn confidence(&self) -> Option<f64> {
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((self.score, mzcv::term!(MS:1001492|percolator:score)))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn other_scores(&self) -> Option<Cow<'_, [(f64, mzcv::Term)]>> {
        Some(
            vec![
                (self.qvalue, mzcv::term!(MS:1001491|percolator:Q value)),
                (self.pep, mzcv::term!(MS:1001493|percolator:PEP)),
            ]
            .into(),
        )
    }

    fn charge(&self) -> Option<Charge> {
        self.peptidoform
            .get_charge_carriers()
            .map(mzcore::chemistry::MolecularCharge::charge)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        Some(self.rt)
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::Number(self.spectrum_id)])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
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
