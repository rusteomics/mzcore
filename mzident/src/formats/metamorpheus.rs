use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use context_error::*;
use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform, IdentifiedPeptidoformData,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, KnownFileFormat, MetaData,
    PeptidoformPresent, SpectrumId, SpectrumIds, common_parser::Location,
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid MetaMorpheus line",
    "This column is not a number but it is required to be a number in this MetaMorpheus format",
);
format_family!(
    MetaMorpheus,
    SemiAmbiguous, PeptidoformPresent, [&META_MORPHEUS], b'\t', None;
    required {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        /// If this is a contaminant, same number of entries as the vector of organisms.
        contaminant: Vec<bool>, |location: Location, _| Ok(location.array('|').map(|l| l.as_str() == "Y").collect());
        decoy: bool, |location: Location, _| Ok(location.as_str() == "Y");
        scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        missed_cleavages: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        precursor_scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        z: Charge, |location: Location, _| location.trim_end_matches(".00000").parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        protein_accession: Vec<String>, |location: Location, _| Ok(location.array('|').map(Location::get_string).collect());
        organism: Vec<String>, |location: Location, _| Ok(location.array('|').map(Location::get_string).collect());
        protein_name: FastaIdentifier<String>, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_location: Vec<Range<u16>>, |location: Location, _| location.array('|').map(|l| l.parse_with(
            |loc| {
                if loc.location.len() < 3 {
                    return Err(BoxedError::new(BasicKind::Error,
                        "Invalid MetaMorpheus line",
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
                        .parse::<u16>()
                        .map_err(|_| {
                            BoxedError::new(BasicKind::Error,NUMBER_ERROR.0, NUMBER_ERROR.1, Context::line(
                                Some(loc.line.line_index() as u32),
                                loc.line.line(),
                                loc.location.start + 1,
                                start,
                            )).to_owned()
                        })?..
                    loc.line.line()[loc.location.end - 1 - end..loc.location.end - 1]
                        .parse::<u16>()
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
        )).collect::<Result<Vec<_>,_>>();
        previous_residue: Vec<FlankingSequence>, |location: Location, _| location.array('|').map(|l| if l.as_str() == "-" {
                        Ok(FlankingSequence::Terminal)
                    } else {
                        Ok(FlankingSequence::AminoAcid(AminoAcid::try_from(l.as_str()).map_err(
                            |()| {
                                BoxedError::new(BasicKind::Error,
                                    "Invalid MetaMorpheus line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        Some(l.line.line_index() as u32),
                                        l.line.line(),
                                        l.location.start,
                                        1,
                                    ).to_owned(),
                                )
                            },
                        )?))
                    }).collect::<Result<Vec<_>,_>>();
        next_residue: Vec<FlankingSequence>, |location: Location, _| location.array('|').map(|l| if l.as_str() == "-" {
                        Ok(FlankingSequence::Terminal)
                    } else {
                        Ok(FlankingSequence::AminoAcid(AminoAcid::try_from(l.as_str()).map_err(
                            |()| {
                                BoxedError::new(BasicKind::Error,
                                    "Invalid MetaMorpheus line",
                                    "The flanking residues could not be parsed as amino acids",
                                    Context::line(
                                        Some(l.line.line_index() as u32),
                                        l.line.line(),
                                        l.location.start,
                                        1,
                                    ).to_owned(),
                                )
                            },
                        )?))
                    }).collect::<Result<Vec<_>,_>>();
        peptide: CompoundPeptidoformIon, |location: Location, ontologies: &Ontologies| location.array('|').map(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            ontologies,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned)).collect::<Result<CompoundPeptidoformIon,_>>();
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        delta_score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        notch: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        matched_ion_series: String, |location: Location, _| Ok(location.get_string());
        matched_ion_mz_ratios: String, |location: Location, _| Ok(location.get_string());
        matched_ion_intensities: String, |location: Location, _| Ok(location.get_string());
        matched_ion_mass_error: String, |location: Location, _| Ok(location.get_string());
        matched_ion_ppm: String, |location: Location, _| Ok(location.get_string());
        matched_ion_counts: String,|location: Location, _| Ok(location.get_string());
        kind: Vec<MetaMorpheusMatchKind>, |location: Location, _| location.array('|').map(|loc| {
            match &loc.line.line()[loc.location.clone()] {
                "T" => Ok(MetaMorpheusMatchKind::Target),
                "C" => Ok(MetaMorpheusMatchKind::Contamination),
                "D" => Ok(MetaMorpheusMatchKind::Decoy),
                _ => Err(BoxedError::new(BasicKind::Error,
                    "Invalid MetaMorpheus line",
                    "The kind column does not contain a valid value (T/C/D)",
                    Context::line(
                        Some(loc.line.line_index() as u32),
                        loc.line.line(),
                        loc.location.start,
                        loc.location.len(),
                    ).to_owned(),
                )),
            }
        }).collect::<Result<Vec<_>,_>>();
        q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep_q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional { }
);

/// All possible peaks versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum MetaMorpheusVersion {
    /// The single known version
    #[default]
    MetaMorpheus,
}

impl std::fmt::Display for MetaMorpheusVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::MetaMorpheus => "",
            }
        )
    }
}

impl IdentifiedPeptidoformVersion<MetaMorpheusFormat> for MetaMorpheusVersion {
    fn format(self) -> MetaMorpheusFormat {
        match self {
            Self::MetaMorpheus => META_MORPHEUS,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::MetaMorpheus => "",
        }
    }
}

/// The only supported format for MetaMorpheus data
pub const META_MORPHEUS: MetaMorpheusFormat = MetaMorpheusFormat {
    version: MetaMorpheusVersion::MetaMorpheus,
    contaminant: "contaminant",
    decoy: "decoy",
    delta_score: "delta score",
    previous_residue: "previous residue",
    next_residue: "next residue",
    kind: "decoy/contaminant/target",
    mass: "precursor mass",
    matched_ion_counts: "matched ion counts",
    matched_ion_intensities: "matched ion intensities",
    matched_ion_mass_error: "matched ion mass diff (da)",
    matched_ion_mz_ratios: "matched ion mass-to-charge ratios",
    matched_ion_ppm: "matched ion mass diff (ppm)",
    matched_ion_series: "matched ion series",
    missed_cleavages: "missed cleavages",
    mz: "precursor mz",
    notch: "notch",
    organism: "organism name",
    pep_q_value: "pep_qvalue",
    pep: "pep",
    peptide: "full sequence",
    precursor_intensity: "precursor intensity",
    precursor_scan_number: "precursor scan number",
    protein_accession: "accession",
    protein_location: "Start and End Residues In Full Sequence",
    protein_name: "name",
    q_value: "qvalue",
    raw_file: "file name",
    rt: "scan retention time",
    scan_number: "scan number",
    score: "score",
    z: "precursor charge",
};

#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[expect(missing_docs)]
pub enum MetaMorpheusMatchKind {
    #[default]
    Decoy,
    Contamination,
    Target,
}

impl std::fmt::Display for MetaMorpheusMatchKind {
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

impl MetaData for MetaMorpheusData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::MetaMorpheus(self.version)
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

    fn mode(&self) -> Option<Cow<'_, str>> {
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

    fn protein_names(&self) -> Option<Cow<'_, [FastaIdentifier<String>]>> {
        Some(Cow::Borrowed(std::slice::from_ref(&self.protein_name)))
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        self.protein_location.first().cloned()
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (
            self.previous_residue
                .first()
                .unwrap_or(&FlankingSequence::Unknown),
            self.next_residue
                .first()
                .unwrap_or(&FlankingSequence::Unknown),
        )
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }
}
