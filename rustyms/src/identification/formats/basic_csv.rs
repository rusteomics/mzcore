use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use serde::{Deserialize, Serialize};

use crate::{
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    sequence::{CompoundPeptidoformIon, Linked},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid CSV line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    BasicCSV,
    Linked, PeptidoformPresent, [&BASIC], b',', None;
    required {
        scan_index: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        sequence: CompoundPeptidoformIon, |location: Location, custom_database: Option<&CustomDatabase>|location.parse_with(|location| CompoundPeptidoformIon::pro_forma(
            location.as_str(),
            custom_database,
        ).map_err(BoxedError::to_owned));
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
    }
    optional {
        mode: String, |location: Location, _| Ok(location.get_string());
    }
);

/// All possible basic CSV versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum BasicCSVVersion {
    /// msms.txt
    #[default]
    Basic,
}

impl IdentifiedPeptidoformVersion<BasicCSVFormat> for BasicCSVVersion {
    fn format(self) -> BasicCSVFormat {
        match self {
            Self::Basic => BASIC,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::Basic => "Basic",
        }
    }
}

impl std::fmt::Display for BasicCSVVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

/// msms.txt
pub const BASIC: BasicCSVFormat = BasicCSVFormat {
    version: BasicCSVVersion::Basic,
    scan_index: "scan_index",
    sequence: "sequence",
    raw_file: "raw_file",
    z: "z",
    mode: OptionalColumn::Optional("mode"),
};

impl MetaData for BasicCSVData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Borrowed(&self.sequence))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::BasicCSV(self.version)
    }

    fn id(&self) -> String {
        self.scan_index.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        None
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<&str> {
        self.mode.as_deref()
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileKnown(vec![(
            self.raw_file.clone(),
            vec![SpectrumId::Index(self.scan_index)],
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
