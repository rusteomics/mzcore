use std::{borrow::Cow, marker::PhantomData, ops::Range, path::PathBuf};

use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{FlankingSequence, Linear, Peptidoform, PeptidoformIonSet, SemiAmbiguous},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};
use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, KnownFileFormat, MaybePeptidoform, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid psm_utils line",
    "This column is not a number but it is required to be a number in this format",
);

static FASTA_ID_ERROR: (&str, &str) = (
    "Invalid psm_utils line",
    "This column is not a valid fasta identifier, some part that needs to be a number is not a number",
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
            ).map_err(|mut errs| if errs.len() == 1 {errs.pop().unwrap().to_owned()} else {BoxedError::new(BasicKind::Error, "Invalid psm_utils line", "This ProForma definition is invalid, see ", location.base_context().add_highlight((0, location.range.clone()))).add_underlying_errors(errs).to_owned()})
                .and_then(|(psm, _)| psm.into_linear().ok_or_else(|| BoxedError::new(BasicKind::Error, "Invalid psm_utils line", "This ProForma definition contains cross-links, which are not allowed here", location.base_context().add_highlight((0, location.range.clone())).to_owned())));
        spectrum_id: SpectrumId, |location: Location, _| Ok(location.clone().parse::<usize>(NUMBER_ERROR).map_or_else(|_| SpectrumId::Native(location.get_string()), SpectrumId::Number));
        run: String, |location: Location, _| Ok(location.get_string());
        collection: String, |location: Location, _| Ok(location.get_string());
        is_decoy: bool, |location: Location, _| Ok(location.as_str().eq_ignore_ascii_case("true"));
        score: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        qvalue: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        pep: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        mz: Option<MassOverCharge>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|v| v.map(MassOverCharge::new::<mzcore::system::thomson>));
        rt: Option<Time>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|v| v.map(Time::new::<mzcore::system::time::s>));
        ion_mobility: Option<f64>, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        proteins: Vec<FastaIdentifier<String>>, |location: Location, _| location.trim_start_matches("['").trim_end_matches("']").array(',').map(|l| l.trim().trim_start_matches("\'").trim_end_matches("\'").parse(FASTA_ID_ERROR)).collect::<Result<Vec<_>, _>>();
        rank: Option<usize>, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
        /// The source column containing a tag for the filetype detected by `psm_utils`, if recognised this also gives the PSI-MS term for this software
        source: (String, Option<mzcv::Term>), |location: Location, _| Ok((location.as_str().to_string(), match location.as_str().to_ascii_lowercase().as_str() {
            "comet" => Some(mzcv::term!(MS:1002251|Comet)),
            "fragpipe" => Some(mzcv::term!(MS:1003429|FragPipe)),
            "msms" => Some(mzcv::term!(MS:1001583|MaxQuant)),
            "percolator" => Some(mzcv::term!(MS:1001490|Percolator)),
            "proteome_discoverer" => Some(mzcv::term!(MS:1000650|Proteome Discoverer)),
            "proteoscape" => Some(mzcv::term!(MS:1004009|ProteoScape)),
            "sage" => Some(mzcv::term!(MS:1004007|Sage)),
            "x!tandem" => Some(mzcv::term!(MS:1001476|X!Tandem)),
            _ => None, // MSAmanda, DiaNN, AlphaDIA, & ionbot miss a term in PSI-MS
        }));
    }
    optional {
        rawfile: PathBuf, |location: Location, _| Ok(PathBuf::from(location.get_string()));
    }
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
    proteins: "protein_list",
    rank: "rank",
    source: "source",
    rawfile: OptionalColumn::Pattern("provenance:", "_filename"),
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
    type Protein = FastaIdentifier<String>;
    #[cfg(feature = "mzannotate")]
    type SpectrumOutputMode = mzcore::chemistry::OutputMolecularFormula;

    fn peptidoform_ion_set(&self) -> Option<Cow<'_, PeptidoformIonSet>> {
        Some(Cow::Owned(self.peptidoform.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PsmUtils(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        if let SpectrumId::Number(n) = self.spectrum_id {
            Some(n)
        } else {
            None
        }
    }

    fn id(&self) -> String {
        self.spectrum_id.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        self.source.1.clone()
    }

    fn confidence(&self) -> Option<f64> {
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        self.score.map(|s| (s, mzcv::term!(MS:1001492|percolator:score)))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn other_scores(&self) -> Option<Cow<'_, [(f64, mzcv::Term)]>> {
        let mut scores = Vec::new();
        if let Some(s) = self.qvalue {
            scores.push((s, mzcv::term!(MS:1002354|PSM-level q-value)))
        }
        if let Some(s) = self.pep {
            scores.push((s, mzcv::term!(MS:1001493|percolator:PEP)))
        }
        if scores.is_empty() {
            None
        } else {
            Some(scores.into())
        }
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
        self.rt
    }

    fn scans(&self) -> SpectrumIds {
        self.rawfile.clone().map_or_else(
            || SpectrumIds::FileNotKnown(vec![self.spectrum_id.clone()]),
            |rawfile| SpectrumIds::FileKnown(vec![(rawfile, vec![self.spectrum_id.clone()])]),
        )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.mz
    }

    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(&self.proteins)
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

    fn is_decoy(&self) -> Option<bool> {
        Some(self.is_decoy)
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}
