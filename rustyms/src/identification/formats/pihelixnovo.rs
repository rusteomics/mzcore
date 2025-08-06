use std::{borrow::Cow, marker::PhantomData, ops::Range};

use crate::{
    error::CustomError,
    identification::{
        FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
        IdentifiedPeptidoformSource, KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId,
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
    "Invalid pi-HelixNovo line",
    "This column is not a number but it is required to be a number in this format",
);

format_family!(
    PiHelixNovo,
    SemiAmbiguous, PeptidoformPresent, [&PIHELIXNOVO_V1_1], b'\t', Some(vec!["title".to_string(),"peptide".to_string(),"score".to_string()]);
    required {
        title: String, |location: Location, _| Ok(location.get_string());
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters {
                allow_unwrapped_modifications: true,
                ..Default::default()
            },
        );
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
    }
    optional { }
);

/// The only known version of pi-HelixNovo
pub const PIHELIXNOVO_V1_1: PiHelixNovoFormat = PiHelixNovoFormat {
    version: PiHelixNovoVersion::V1_1,
    title: "title",
    peptide: "peptide",
    score: "score",
};

/// All possible pi-HelixNovo versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PiHelixNovoVersion {
    #[default]
    /// pi-HelixNovo version 1.1
    V1_1,
}

impl std::fmt::Display for PiHelixNovoVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PiHelixNovoFormat> for PiHelixNovoVersion {
    fn format(self) -> PiHelixNovoFormat {
        match self {
            Self::V1_1 => PIHELIXNOVO_V1_1,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V1_1 => "v1.1",
        }
    }
}

impl MetaData for PiHelixNovoData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PiHelixNovo(self.version)
    }

    fn id(&self) -> String {
        "-".to_string() // TODO: best would be to use the scan index in some way shape or form
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
        None
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    // TODO: maybe do this parsing as post processing right after parsing the line instead of every time this function is called.
    // Related: this looks like a mgf title line, so I like the setup here that is very resilient to lines that are not in the same format.
    // It might make sense though to have MGF title line processing somewhere more organised (maybe even call mzcore?)
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
