use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    error::CustomError,
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MaybePeptidoform, MetaData, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::CompoundPeptidoformIon,
    sequence::{Peptidoform, SemiAmbiguous},
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
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        z: Charge, |location: Location, _| location
            .trim_end_matches(".0")
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<crate::system::e>);
    }
    optional {
        start_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        end_time: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::pro_forma(location.as_str(), custom_database).map(|p|p.into_semi_ambiguous().unwrap());
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        score_type: String, |location: Location, _| Ok(location.get_string());
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        adduct: String, |location: Location, _| Ok(location.get_string());
        precursormz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
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

impl IdentifiedPeptidoformVersion<SpectrumSequenceListFormat> for SpectrumSequenceListVersion {
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

impl MetaData for SpectrumSequenceListData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::SpectrumSequenceList(self.version)
    }

    fn id(&self) -> String {
        self.scan.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        self.score
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        self.score
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<&str> {
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

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<usize>> {
        None
    }
}
