use std::{
    ops::Range,
    path::{Path, PathBuf},
};

use crate::identification::{
    common_parser::Location,
    csv::{parse_csv, CsvLine},
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, MetaData,
};
use crate::{
    error::{Context, CustomError},
    ontologies::CustomDatabase,
    peptidoform::{SemiAmbiguous, SloppyParsingParameters},
    system::{usize::Charge, Mass, MassOverCharge, Time},
    AminoAcid, Peptidoform,
};
use serde::{Deserialize, Serialize};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid OPair line",
    "This column is not a number but it is required to be a number in this OPair format",
);
format_family!(
    /// The format for OPair data
    OpairFormat,
    /// The data for OPair data
    OpairData,
    OpairVersion, [&O_PAIR], b'\t', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        accession: String, |location: Location, _| Ok(location.get_string());
        organism: String, |location: Location, _| Ok(location.get_string());
        protein_name: FastaIdentifier<String>, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_location: Range<usize>, |location: Location, _| location.parse_with(
            |loc| {
                if loc.location.len() < 3 {
                    return Err(CustomError::error(
                        "Invalid Opair line",
                        "The location is not defined, it should be defined like this [<start> to <end>]",
                        Context::line(
                            Some(loc.line.line_index()),
                            loc.line.line(),
                            loc.location.start,
                            loc.location.len(),
                        ),
                    ))
                }
                let bytes =
                    &loc.line.line().as_bytes()[loc.location.start + 1..loc.location.end-1];
                let start = bytes.iter().take_while(|c| c.is_ascii_digit()).count();
                let end = bytes
                    .iter()
                    .rev()
                    .take_while(|c| c.is_ascii_digit())
                    .count();
                Ok(
                    loc.line.line()[loc.location.start + 1..loc.location.start + 1 + start]
                        .parse()
                        .map_err(|_| {
                            CustomError::error(NUMBER_ERROR.0, NUMBER_ERROR.1, Context::line(
                                Some(loc.line.line_index()),
                                loc.line.line(),
                                loc.location.start + 1,
                                start,
                            ))
                        })?..
                    loc.line.line()[loc.location.end - 1 - end..loc.location.end - 1]
                        .parse()
                        .map_err(|_| {
                            CustomError::error(NUMBER_ERROR.0, NUMBER_ERROR.1, Context::line(
                                Some(loc.line.line_index()),
                                loc.line.line(),
                                loc.location.end - 1 - end,
                                end,
                            ))
                        })?
                )
            },
        );
        base_sequence: String, |location: Location, _| Ok(location.get_string());
        flanking_residues: (AminoAcid, AminoAcid),|location: Location, _| location.parse_with(
            |loc| {
                Ok((
                    AminoAcid::try_from(loc.line.line().as_bytes()[loc.location.start]).map_err(
                        |()| {
                            CustomError::error(
                                "Invalid Opair line",
                                "The flanking residues could not be parsed as amino acids",
                                Context::line(
                                    Some(loc.line.line_index()),
                                    loc.line.line(),
                                    loc.location.start,
                                    1,
                                ),
                            )
                        },
                    )?,
                    AminoAcid::try_from(loc.line.line().as_bytes()[loc.location.end - 1])
                        .map_err(|()| {
                            CustomError::error(
                                "Invalid Opair line",
                                "The flanking residues could not be parsed as amino acids",
                                Context::line(
                                    Some(loc.line.line_index()),
                                    loc.line.line(),
                                    loc.location.end - 1,
                                    1,
                                ),
                            )
                        })?
                ))
            },
        );
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default()
        );
        mod_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        rank: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        matched_ion_series: String, |location: Location, _| Ok(location.get_string());
        matched_ion_mz_ratios: String, |location: Location, _| Ok(location.get_string());
        matched_ion_intensities: String, |location: Location, _| Ok(location.get_string());
        matched_ion_mass_error: String, |location: Location, _| Ok(location.get_string());
        matched_ion_ppm: String, |location: Location, _| Ok(location.get_string());
        matched_ion_counts: String,|location: Location, _| Ok(location.get_string());
        kind: OpairMatchKind, |location: Location, _| location.parse_with(|loc| {
            match &loc.line.line()[loc.location.clone()] {
                "T" => Ok(OpairMatchKind::Target),
                "C" => Ok(OpairMatchKind::Contamination),
                "D" => Ok(OpairMatchKind::Decoy),
                _ => Err(CustomError::error(
                    "Invalid Opair line",
                    "The kind column does not contain a valid value (T/C/D)",
                    Context::line(
                        Some(loc.line.line_index()),
                        loc.line.line(),
                        loc.location.start,
                        loc.location.len(),
                    ),
                )),
            }
        });
        q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep_q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        localisation_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        yion_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        diagnostic_ion_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        plausible_glycan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        total_glycosylation_sites: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        glycan_mass:Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        plausible_glycan_composition: String, |location: Location, _| Ok(location.get_string());
        n_glycan_motif: bool, |location: Location, _| location.parse_with(|loc| {
            match &loc.line.line()[loc.location.clone()] {
                "TRUE" => Ok(true),
                "FALSE" => Ok(false),
                _ => Err(CustomError::error(
                    "Invalid Opair line",
                    "The N glycan motif check column does not contain a valid value (TRUE/FALSE)",
                    Context::line(
                        Some(loc.line.line_index()),
                        loc.line.line(),
                        loc.location.start,
                        loc.location.len(),
                    ),
                )),
            }
        });
        r138_144: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        plausible_glycan_structure: String, |location: Location, _| Ok(location.get_string());
        glycan_localisation_level: String, |location: Location, _| Ok(location
            .get_string()
            .trim_start_matches("Level")
            .to_string());
        glycan_peptide_site_specificity: String, |location: Location, _| Ok(location.get_string());
        glycan_protein_site_specificity:String, |location: Location, _| Ok(location.get_string());
        all_potential_glycan_localisations: String, |location: Location, _| Ok(location.get_string());
        all_site_specific_localisation_probabilities: String, |location: Location, _| Ok(location.get_string());
    }
    optional { }
);

impl From<OpairData> for IdentifiedPeptidoform {
    fn from(value: OpairData) -> Self {
        Self {
            score: Some(value.score / 100.0),
            local_confidence: None,
            metadata: MetaData::Opair(value),
        }
    }
}

/// All possible peaks versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum OpairVersion {
    /// The single known version
    #[default]
    Opair,
}

impl std::fmt::Display for OpairVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Opair => "",
            }
        )
    }
}

impl IdentifiedPeptidoformVersion<OpairFormat> for OpairVersion {
    fn format(self) -> OpairFormat {
        match self {
            Self::Opair => O_PAIR,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::Opair => "",
        }
    }
}

/// The only supported format for Opair data
pub const O_PAIR: OpairFormat = OpairFormat {
    version: OpairVersion::Opair,
    raw_file: "file name",
    scan_number: "scan number",
    rt: "scan retention time",
    precursor_scan_number: "precursor scan number",
    mz: "precursor mz",
    z: "precursor charge",
    mass: "precursor mass",
    accession: "protein accession",
    organism: "organism",
    protein_name: "protein name",
    protein_location: "start and end residues in protein",
    base_sequence: "base sequence",
    flanking_residues: "flankingresidues",
    peptide: "full sequence",
    mod_number: "number of mods",
    theoretical_mass: "peptide monoisotopic mass",
    score: "score",
    rank: "rank",
    matched_ion_series: "matched ion series",
    matched_ion_mz_ratios: "matched ion mass-to-charge ratios",
    matched_ion_mass_error: "matched ion mass diff (da)",
    matched_ion_ppm: "matched ion mass diff (ppm)",
    matched_ion_intensities: "matched ion intensities",
    matched_ion_counts: "matched ion counts",
    kind: "decoy/contaminant/target",
    q_value: "qvalue",
    pep: "pep",
    pep_q_value: "pep_qvalue",
    localisation_score: "localization score",
    yion_score: "yion score",
    diagnostic_ion_score: "diagonosticion score",
    plausible_glycan_number: "plausible number of glycans",
    total_glycosylation_sites: "total glycosylation sites",
    glycan_mass: "glycanmass",
    plausible_glycan_composition: "plausible glycancomposition",
    n_glycan_motif: "n-glycan motif check",
    r138_144: "r138/144",
    plausible_glycan_structure: "plausible glycanstructure",
    glycan_localisation_level: "glycanlocalizationlevel",
    glycan_peptide_site_specificity: "localized glycans with peptide site specific probability",
    glycan_protein_site_specificity: "localized glycans with protein site specific probability",
    all_potential_glycan_localisations: "all potential glycan localizations",
    all_site_specific_localisation_probabilities: "allsitespecificlocalizationprobability",
};

#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[expect(missing_docs)]
pub enum OpairMatchKind {
    #[default]
    Decoy,
    Contamination,
    Target,
}

impl std::fmt::Display for OpairMatchKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::Decoy => "Decoy",
                Self::Contamination => "Contamination",
                Self::Target => "Target",
            }
        )
    }
}
