use std::{
    borrow::Cow,
    marker::PhantomData,
    ops::Range,
    path::{Path, PathBuf},
};

use context_error::*;
use serde::{Deserialize, Serialize};

use crate::{
    BoxedIdentifiedPeptideIter, CVTerm, FastaIdentifier, KnownFileFormat, MetaMorpheusMatchKind,
    PSM, PSMData, PSMFileFormatVersion, PSMMetaData, PSMSource, PeptidoformPresent,
    ProteinMetaData, Reliability, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalLocation},
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, FlankingSequence, Peptidoform, SemiAmbiguous,
        SequencePosition, SimpleModification, SloppyParsingParameters,
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
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        precursor_scan_number: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
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
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, ontologies: &Ontologies| Peptidoform::sloppy_pro_forma_inner(
            &location.base_context(),
            location.full_line(),
            location.location.clone(),
            ontologies,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned);
        mod_number: u8, |location: Location, _| location.parse(NUMBER_ERROR);
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        score: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        rank: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        matched_ion_series: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        matched_ion_mz_ratios: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        matched_ion_intensities: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        matched_ion_mass_error: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        matched_ion_ppm: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        matched_ion_counts: Box<str>,|location: Location, _| Ok(location.get_boxed_str());
        q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        pep_q_value: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        localisation_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        yion_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        diagnostic_ion_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        plausible_glycan_number: u8, |location: Location, _| location.parse(NUMBER_ERROR);
        total_glycosylation_sites: u8, |location: Location, _| location.parse(NUMBER_ERROR);
        glycan_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<mzcore::system::dalton>);
        plausible_glycan_composition: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
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
        plausible_glycan_structure: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        glycan_localisation_level: Box<str>, |location: Location, _| Ok(location
            .trim_start_matches("Level")
            .get_boxed_str());
        glycan_peptide_site_specificity: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        glycan_protein_site_specificity:Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        all_potential_glycan_localisations: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        all_site_specific_localisation_probabilities: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
    }
    optional { }

    fn post_process(_source: &CsvLine, parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        Ok(parsed)
    }

    protein {
        accession => (|location: Location, _| Ok(location.get_string()));

        required {
            organism: String, |location: Location, _| Ok(location.get_string());
            protein_name: FastaIdentifier<Box<str>>, |location: Location, _| location.parse(NUMBER_ERROR);
            kind: MetaMorpheusMatchKind, |location: Location, _| location.parse_with(|loc| {
            match &loc.line.line()[loc.location.clone()] {
                "T" => Ok(MetaMorpheusMatchKind::Target),
                "C" => Ok(MetaMorpheusMatchKind::Contamination),
                "D" => Ok(MetaMorpheusMatchKind::Decoy),
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
        }
        optional {}
    }
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

impl PSMFileFormatVersion<OpairFormat> for OpairVersion {
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

impl PSMMetaData for OpairPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Opair(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.scan_number)
    }

    fn id(&self) -> String {
        self.scan_number.to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        Some(mzcv::term!(MS:1002826|MetaMorpheus))
    }

    fn confidence(&self) -> Option<f64> {
        Some(self.score / 100.0)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        Some((self.score, mzcv::term!(MS:1002827|MetaMorpheus:score)))
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

    type Protein = OpairProtein;

    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(std::slice::from_ref(&self.accession))
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

    fn unique(&self) -> Option<bool> {
        None
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }
}

impl ProteinMetaData for OpairProtein {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<mzcore::sequence::Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.protein_name.as_str()
    }

    fn description(&self) -> Option<&str> {
        None
    }

    fn species(&self) -> Option<mzcv::Curie> {
        None
    }

    fn species_name(&self) -> Option<&str> {
        Some(&self.organism)
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &[] // TODO: PLGS does not have a search engine entry
    }

    fn ambiguity_members(&self) -> &[String] {
        &[]
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &[]
    }

    fn coverage(&self) -> Option<f64> {
        None
    }

    fn gene_ontology(&self) -> &[mzcv::Curie] {
        &[]
    }

    fn reliability(&self) -> Option<Reliability> {
        None
    }

    fn uri(&self) -> Option<&str> {
        None
    }
}
