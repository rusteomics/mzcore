use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use custom_error::*;
use serde::{Deserialize, Serialize};

use crate::{
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, Peptidoform, SemiAmbiguous, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid OPair line",
    "This column is not a number but it is required to be a number in this OPair format",
);
format_family!(
    Opair,
    SemiAmbiguous, PeptidoformPresent, [&O_PAIR], b'\t', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        accession: String, |location: Location, _| Ok(location.get_string());
        organism: String, |location: Location, _| Ok(location.get_string());
        protein_name: FastaIdentifier<String>, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_location: Option<Range<u16>>, |location: Location, _| location.or_empty().parse_with(
            |loc| {
                if loc.location.len() < 3 {
                    return Err(BoxedError::new(BasicKind::Error,
                        "Invalid Opair line",
                        "The location is not defined, it should be defined like this [<start> to <end>]",
                        Context::line(
                            Some(loc.line.line_index() as u32),
                            loc.line.line(),
                            loc.location.start,
                            loc.location.len(),
                        ).to_owned(),
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
                            BoxedError::new(BasicKind::Error,NUMBER_ERROR.0, NUMBER_ERROR.1, Context::line(
                                Some(loc.line.line_index() as u32),
                                loc.line.line(),
                                loc.location.start + 1,
                                start,
                            )).to_owned()
                        })?..
                    loc.line.line()[loc.location.end - 1 - end..loc.location.end - 1]
                        .parse()
                        .map_err(|_| {
                            BoxedError::new(BasicKind::Error,NUMBER_ERROR.0, NUMBER_ERROR.1, Context::line(
                                Some(loc.line.line_index() as u32),
                                loc.line.line(),
                                loc.location.end - 1 - end,
                                end,
                            ).to_owned())
                        })?
                )
            },
        );
        base_sequence: String, |location: Location, _| Ok(location.get_string());
        flanking_residues: (FlankingSequence, FlankingSequence), |location: Location, _| location.parse_with(
            |loc| {
                let n = loc.line.line().as_bytes()[loc.location.start];
                let c = loc.line.line().as_bytes()[loc.location.end - 1];
                Ok((
                    if n == b'-' {
                        FlankingSequence::Terminal
                    } else {
                        FlankingSequence::AminoAcid(AminoAcid::try_from(n).map_err(
                            |()| {
                                BoxedError::new(BasicKind::Error,
                                    "Invalid Opair line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        Some(loc.line.line_index() as u32),
                                        loc.line.line(),
                                        loc.location.start,
                                        1,
                                    ).to_owned(),
                                )
                            },
                        )?)
                    },
                    if c == b'-' {
                        FlankingSequence::Terminal
                    } else {
                        FlankingSequence::AminoAcid(AminoAcid::try_from(c).map_err(
                            |()| {
                                BoxedError::new(BasicKind::Error,
                                    "Invalid Opair line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        Some(loc.line.line_index() as u32),
                                        loc.line.line(),
                                        loc.location.end - 1,
                                        1,
                                    ).to_owned(),
                                )
                        })?)
                    }
                ))
            },
        );
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned);
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
                _ => Err(BoxedError::new(BasicKind::Error,
                    "Invalid Opair line",
                    "The kind column does not contain a valid value (T/C/D)",
                    Context::line(
                        Some(loc.line.line_index() as u32),
                        loc.line.line(),
                        loc.location.start,
                        loc.location.len(),
                    ).to_owned(),
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
                _ => Err(BoxedError::new(BasicKind::Error,
                    "Invalid Opair line",
                    "The N glycan motif check column does not contain a valid value (TRUE/FALSE)",
                    Context::line(
                        Some(loc.line.line_index() as u32),
                        loc.line.line(),
                        loc.location.start,
                        loc.location.len(),
                    ).to_owned(),
                )),
            }
        });
        r138_144: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        plausible_glycan_structure: String, |location: Location, _| Ok(location.get_string());
        glycan_localisation_level: String, |location: Location, _| Ok(location
            .trim_start_matches("Level")
            .get_string());
        glycan_peptide_site_specificity: String, |location: Location, _| Ok(location.get_string());
        glycan_protein_site_specificity:String, |location: Location, _| Ok(location.get_string());
        all_potential_glycan_localisations: String, |location: Location, _| Ok(location.get_string());
        all_site_specific_localisation_probabilities: String, |location: Location, _| Ok(location.get_string());
    }
    optional { }
);

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
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
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

#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[expect(missing_docs)]
pub enum OpairMatchKind {
    #[default]
    Decoy,
    Contamination,
    Target,
}

impl std::fmt::Display for OpairMatchKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
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

impl MetaData for OpairData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Opair(self.version)
    }

    fn id(&self) -> String {
        self.scan_number.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score / 100.0)
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
        Some(self.rt)
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileKnown(vec![(
            self.raw_file.clone(),
            vec![SpectrumId::Number(self.scan_number)],
        )])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.mass)
    }

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        Some(self.protein_name.clone())
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        self.protein_location.clone()
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&self.flanking_residues.0, &self.flanking_residues.1)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }
}
