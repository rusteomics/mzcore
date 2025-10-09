use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
    IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
    KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn},
    csv::{CsvLine, parse_csv},
};
use mzcore::{
    ontology::CustomDatabase,
    sequence::{CompoundPeptidoformIon, Peptidoform, SemiAmbiguous, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Novor line",
    "This column is not a number but it is required to be a number in this Novor format",
);

format_family!(
    Novor,
    SemiAmbiguous, PeptidoformPresent, [&OLD_DENOVO, &OLD_PSM, &NEW_DENOVO, &NEW_PSM], b',', None;
    required {
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default(),
        ).map_err(BoxedError::to_owned);
    }
    optional {
        id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        spectra_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        fraction: usize, |location: Location, _| location.skip(1).parse::<usize>(NUMBER_ERROR); // Skip the F of the F{num} definition
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        peptide_no_ptm: String, |location: Location, _| Ok(Some(location.get_string()));
        protein: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        protein_start: u16, |location: Location, _| location.parse::<u16>(NUMBER_ERROR);
        protein_origin: String, |location: Location, _| Ok(Some(location.get_string()));
        protein_all: String, |location: Location, _| Ok(Some(location.get_string()));
        database_sequence: String, |location: Location, _| Ok(Some(location.get_string()));
        local_confidence: Vec<f64>, |location: Location, _| location.array('-')
                    .map(|l| l.parse::<f64>(NUMBER_ERROR))
                    .collect::<Result<Vec<_>, _>>();
    }
);

/// All available Novor versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum NovorVersion {
    /// An older version for the denovo file
    #[default]
    OldDenovo,
    /// An older version for the psms file
    OldPSM,
    /// Seen since v3.36.893 (not necessarily the time it was rolled out)
    NewDenovo,
    /// Seen since v3.36.893 (not necessarily the time it was rolled out)
    NewPSM,
}
impl std::fmt::Display for NovorVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<NovorFormat> for NovorVersion {
    fn format(self) -> NovorFormat {
        match self {
            Self::OldDenovo => OLD_DENOVO,
            Self::OldPSM => OLD_PSM,
            Self::NewDenovo => NEW_DENOVO,
            Self::NewPSM => NEW_PSM,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::OldDenovo => "Older Denovo",
            Self::OldPSM => "Older PSM",
            Self::NewDenovo => "New Denovo",
            Self::NewPSM => "New PSM",
        }
    }
}

/// The older supported format for denovo.csv files from Novor
///# -
/// 1       Fraction
/// 2       Scan #
/// 3       m/z
/// 4       z
/// 5       Score
/// 6       Peptide Mass
/// 7       Error (ppm)
/// 8       Length
/// 9       De Novo Peptide
/// 10      DB Sequence
/// <https://github.com/snijderlab/stitch/issues/156#issuecomment-1097862072>
pub const OLD_DENOVO: NovorFormat = NovorFormat {
    version: NovorVersion::OldDenovo,
    scan_number: "scan #",
    mz: "m/z",
    z: "z",
    mass: "peptide mass",
    score: "score",
    peptide: "de novo peptide",
    id: OptionalColumn::NotAvailable,
    spectra_id: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::Required("fraction"),
    rt: OptionalColumn::NotAvailable,
    peptide_no_ptm: OptionalColumn::NotAvailable,
    protein: OptionalColumn::NotAvailable,
    protein_start: OptionalColumn::NotAvailable,
    protein_origin: OptionalColumn::NotAvailable,
    protein_all: OptionalColumn::NotAvailable,
    database_sequence: OptionalColumn::Required("db sequence"),
    local_confidence: OptionalColumn::NotAvailable,
};

/// The older supported format for psms.csv files from Novor
///# -
/// 1       ID
/// 2       Fraction
/// 3       Scan
/// 4       m/z
/// 5       z
/// 6       Score
/// 7       Mass
/// 8       Error (ppm)
/// 9       # Proteins
/// 10      Sequence
/// <https://github.com/snijderlab/stitch/issues/156#issuecomment-1097862072>
pub const OLD_PSM: NovorFormat = NovorFormat {
    version: NovorVersion::OldPSM,
    scan_number: "scan",
    mz: "m/z",
    z: "z",
    mass: "mass",
    score: "score",
    peptide: "sequence",
    id: OptionalColumn::Required("id"),
    spectra_id: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::Required("fraction"),
    rt: OptionalColumn::NotAvailable,
    peptide_no_ptm: OptionalColumn::NotAvailable,
    protein: OptionalColumn::Required("# proteins"),
    protein_start: OptionalColumn::NotAvailable,
    protein_origin: OptionalColumn::NotAvailable,
    protein_all: OptionalColumn::NotAvailable,
    database_sequence: OptionalColumn::NotAvailable,
    local_confidence: OptionalColumn::NotAvailable,
};

/// denovo: `# id, scanNum, RT, mz(data), z, pepMass(denovo), err(data-denovo), ppm(1e6*err/(mz*z)), score, peptide, aaScore,`
pub const NEW_DENOVO: NovorFormat = NovorFormat {
    version: NovorVersion::NewDenovo,
    scan_number: "scannum",
    mz: "mz(data)",
    z: "z",
    mass: "pepmass(denovo)",
    score: "score",
    peptide: "peptide",
    id: OptionalColumn::Required("# id"),
    spectra_id: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::NotAvailable,
    rt: OptionalColumn::Required("rt"),
    peptide_no_ptm: OptionalColumn::NotAvailable,
    protein: OptionalColumn::NotAvailable,
    protein_start: OptionalColumn::NotAvailable,
    protein_origin: OptionalColumn::NotAvailable,
    protein_all: OptionalColumn::NotAvailable,
    database_sequence: OptionalColumn::NotAvailable,
    local_confidence: OptionalColumn::Required("aascore"),
};

/// PSM: `#id, spectraId, scanNum, RT, mz, z, pepMass, err, ppm, score, protein, start, length, origin, peptide, noPTMPeptide, aac, allProteins`
pub const NEW_PSM: NovorFormat = NovorFormat {
    version: NovorVersion::NewPSM,
    scan_number: "scannum",
    mz: "mz",
    z: "z",
    mass: "pepmass",
    score: "score",
    peptide: "peptide",
    id: OptionalColumn::Required("#id"),
    spectra_id: OptionalColumn::Required("spectraid"),
    fraction: OptionalColumn::NotAvailable,
    rt: OptionalColumn::Required("rt"),
    peptide_no_ptm: OptionalColumn::Required("noptmpeptide"),
    protein: OptionalColumn::Required("protein"),
    protein_start: OptionalColumn::Required("start"),
    protein_origin: OptionalColumn::Required("origin"),
    protein_all: OptionalColumn::Required("allproteins"),
    database_sequence: OptionalColumn::NotAvailable,
    local_confidence: OptionalColumn::Required("aac"),
};

impl MetaData for NovorData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Novor(self.version)
    }

    fn id(&self) -> String {
        self.id.unwrap_or(self.scan_number).to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some((self.score / 100.0).clamp(-1.0, 1.0))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| lc.iter().map(|v| *v / 100.0).collect())
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        self.local_confidence.as_deref()
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
        SpectrumIds::FileNotKnown(vec![SpectrumId::Number(self.scan_number)])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mass)
    }

    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        self.protein
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        self.protein_start.map(|s| s..s + self.peptide.len() as u16)
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }
}
