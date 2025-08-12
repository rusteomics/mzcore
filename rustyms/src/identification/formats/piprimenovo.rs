use std::{borrow::Cow, marker::PhantomData, ops::Range};

use crate::{
    error::CustomError,
    identification::{
        FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
        IdentifiedPeptidoformSource, KnownFileFormat, MaybePeptidoform, MetaData, SpectrumId,
        SpectrumIds,
    },
    ontology::CustomDatabase,
    prelude::CompoundPeptidoformIon,
    sequence::{Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Ratio, Time, isize::Charge},
};

use serde::{Deserialize, Serialize};

use crate::identification::{
    BoxedIdentifiedPeptideIter, IdentifiedPeptidoformVersion,
    common_parser::Location,
    csv::{CsvLine, parse_csv},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid pi-PrimeNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PiPrimeNovo,
    SemiAmbiguous, MaybePeptidoform, [&PIPRIMENOVO_V0_1], b'\t', None;
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().map(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default()
        )).transpose();
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        z: Charge, |location: Location, _| location
            .parse::<isize>(NUMBER_ERROR)
            .map(Charge::new::<crate::system::e>);
    }
    optional { }
);

/// The only known version of pi-PrimeNovo
pub const PIPRIMENOVO_V0_1: PiPrimeNovoFormat = PiPrimeNovoFormat {
    version: PiPrimeNovoVersion::V0_1,
    title: "label",
    peptide: "prediction",
    score: "score",
    z: "charge",
};

/// All possible pi-PrimeNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PiPrimeNovoVersion {
    #[default]
    /// pi-PrimeNovo version 1.1
    V0_1,
}

impl std::fmt::Display for PiPrimeNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PiPrimeNovoFormat> for PiPrimeNovoVersion {
    fn format(self) -> PiPrimeNovoFormat {
        match self {
            Self::V0_1 => PIPRIMENOVO_V0_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V0_1 => "v0.1",
        }
    }
}

impl MetaData for PiPrimeNovoData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PiPrimeNovo(self.version)
    }

    fn id(&self) -> String {
        "-".to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score)
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
        Some(self.z)
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        match (self.title.find("File:\""), self.title.find("NativeID:\"")) {
            (None, None) => SpectrumIds::None,
            (None, Some(native_id)) => {
                let native_content_start = native_id + 10; // Skip 'NativeID:"'
                self.title[native_content_start..].find('"').map_or(
                    SpectrumIds::None,
                    |native_end| {
                        let native_id: SpectrumId = SpectrumId::Native(
                            self.title[native_content_start..native_content_start + native_end]
                                .to_string(),
                        );

                        SpectrumIds::FileNotKnown(vec![native_id])
                    },
                )
            }
            (Some(file_path), None) => {
                let file_content_start = file_path + 6; // Skip 'File:"'

                self.title[file_content_start..]
                    .find('"')
                    .map_or(SpectrumIds::None, |file_end| {
                        let raw_file = self.title
                            [file_content_start..file_content_start + file_end]
                            .to_string();

                        SpectrumIds::FileKnown(vec![(raw_file.into(), vec![])])
                    })
            }
            (Some(file_path), Some(native_id)) => {
                let file_content_start = file_path + 6; // Skip 'File:"'
                let file_end = self.title[file_content_start..].find('"');

                let native_content_start = native_id + 10; // Skip 'NativeID:"'
                let native_end = self.title[native_content_start..].find('"');

                match (file_end, native_end) {
                    (None, None) => SpectrumIds::None,
                    (None, Some(native_end)) => {
                        let native_id: SpectrumId = SpectrumId::Native(
                            self.title[native_content_start..native_content_start + native_end]
                                .to_string(),
                        );
                        SpectrumIds::FileNotKnown(vec![native_id])
                    }
                    (Some(file_end), None) => {
                        let raw_file = self.title
                            [file_content_start..file_content_start + file_end]
                            .to_string();

                        SpectrumIds::FileKnown(vec![(raw_file.into(), vec![])])
                    }
                    (Some(file_end), Some(native_end)) => {
                        let native_id: SpectrumId = SpectrumId::Native(
                            self.title[native_content_start..native_content_start + native_end]
                                .to_string(),
                        );
                        let raw_file = self.title
                            [file_content_start..file_content_start + file_end]
                            .to_string();
                        SpectrumIds::FileKnown(vec![(raw_file.into(), vec![native_id])])
                    }
                }
            }
        }
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
