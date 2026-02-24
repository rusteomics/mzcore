use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, KnownFileFormat, MaybePeptidoform, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid SpectrumSequenceList line",
    "This column is not a number but it is required to be a number in this format",
);

// The format for any [SSL file](https://skyline.ms/wiki/home/software/BiblioSpec/page.view?name=BiblioSpec%20input%20and%20output%20file%20formats).
format_family!(
    SpectrumSequenceList,
    SemiAmbiguous, MaybePeptidoform, [&SSL], b'\t', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(location.as_str()).to_owned());
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<mzcore::system::e>);
    }
    optional {
        start_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        end_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies| Peptidoform::pro_forma_inner(&location.context(), location.full_line(), location.range.clone(), ontologies).map_err(|errs| BoxedError::new(BasicKind::Error, "Invalid ProForma definition", "The string could not be parsed as a ProForma definition", location.context()).add_underlying_errors(errs)).map_err(BoxedError::to_owned).map(|(p, _)| p.into_semi_ambiguous().unwrap());
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score_type: String, |location: Location, _| Ok(location.get_string());
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        adduct: String, |location: Location, _| Ok(location.get_string());
        precursormz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        moleculename: String, |location: Location, _| Ok(location.get_string());
        inchikey: String, |location: Location, _| Ok(location.get_string());
        otherkeys: String, |location: Location, _| Ok(location.or_empty().get_string());
        ion_mobility: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        ion_mobility_units: String, |location: Location, _| Ok(location.get_string());
        ccs: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
);

/// General type of SSL files
pub const SSL: SpectrumSequenceListFormat = SpectrumSequenceListFormat {
    version: SpectrumSequenceListVersion::SSL,
    raw_file: "file",
    scan: "scan",
    z: "charge",
    start_time: OptionalColumn::Optional("start-time"),
    end_time: OptionalColumn::Optional("end-time"),
    peptide: OptionalColumn::Optional("sequence"),
    score: OptionalColumn::Optional("score"),
    score_type: OptionalColumn::Optional("score-type"),
    rt: OptionalColumn::Optional("retention-time"),
    adduct: OptionalColumn::Optional("adduct"),
    precursormz: OptionalColumn::Optional("precursorMZ"),
    moleculename: OptionalColumn::Optional("moleculename"),
    inchikey: OptionalColumn::Optional("inchikey"),
    otherkeys: OptionalColumn::Optional("otherkeys"),
    ion_mobility: OptionalColumn::Optional("ion-mobility"),
    ion_mobility_units: OptionalColumn::Optional("ion-mobility-units"),
    ccs: OptionalColumn::Optional("ccs"),
};

/// All possible SpectrumSequenceList versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[expect(clippy::upper_case_acronyms)]
pub enum SpectrumSequenceListVersion {
    #[default]
    /// SSL file format
    SSL,
}

impl std::fmt::Display for SpectrumSequenceListVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl PSMFileFormatVersion<SpectrumSequenceListFormat> for SpectrumSequenceListVersion {
    fn format(self) -> SpectrumSequenceListFormat {
        match self {
            Self::SSL => SSL,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::SSL => "",
        }
    }
}

impl PSMMetaData for SpectrumSequenceListPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::SpectrumSequenceList(self.version)
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
        self.score
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        self.score
            .map(|v| (v, mzcv::term!(MS:1001153|search engine specific score)))
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
        self.rt
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileKnown(vec![(
            self.raw_file.clone(),
            vec![SpectrumId::Index(self.scan)],
        )])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
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
