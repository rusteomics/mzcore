use std::{
    borrow::Cow,
    marker::PhantomData,
    num::NonZeroU32,
    ops::Range,
    path::{Path, PathBuf},
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    helper_functions::end_of_enclosure,
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MaybePeptidoform, MetaData, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::{AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, SequenceElement},
    sequence::{Modification, Peptidoform, SemiAmbiguous, SimpleLinear, SloppyParsingParameters},
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid MaxQuant line",
    "This column is not a number but it is required to be a number in this MaxQuant format",
);
static BOOL_ERROR: (&str, &str) = (
    "Invalid MaxQuant line",
    "This column is not a boolean ('0' or '1') but it is required to be a boolean in this MaxQuant format",
);

format_family!(
    MaxQuant,
    SemiAmbiguous, MaybePeptidoform, [&MSMS, &NOVO_MSMS_SCANS, &MSMS_SCANS, &SILAC], b'\t', None;
    required {
        scan_number: ThinVec<usize>, |location: Location, _| location.or_empty().array(';').map(|s| s.parse(NUMBER_ERROR)).collect::<Result<ThinVec<usize>, BoxedError<'_, BasicKind>>>();
        modifications: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        proteins: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        peptide: Option<Peptidoform<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>| location.or_empty().parse_with(|location| Peptidoform::sloppy_pro_forma(
            location.full_line(),
            location.location.clone(),
            custom_database,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned));
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        ty: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        pep: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        all_modified_sequences: ThinVec<Peptidoform<SemiAmbiguous>>, |location: Location, custom_database: Option<&CustomDatabase>| location.array(';')
                .map(|s| Peptidoform::sloppy_pro_forma(s.line.line(), s.location, custom_database, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned))
                .collect::<Result<ThinVec<Peptidoform<SemiAmbiguous>>, BoxedError<'static, BasicKind>>>();
        base_peak_intensity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        carbamidomethyl_c_probabilities: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        carbamidomethyl_c_score_differences: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        collision_energy: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        delta_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        dn_c_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_combined_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        dn_complete: bool, |location: Location, _| Ok(location.as_str() == "+");
        dn_n_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        dn_sequence: Peptidoform<SimpleLinear>, |location: Location, custom_database: Option<&CustomDatabase>| parse_de_novo_sequence(location, custom_database);
        evidence_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        experiment: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        mode: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        genes: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        intensity_coverage: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        intensity_h: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        intensity_l: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        intensity: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        isotope_index: isize, |location: Location, _| location.or_empty().parse::<isize>(NUMBER_ERROR);
        labeling_state: bool, |location: Location, _| location.or_empty().ignore("-1").parse::<u8>(BOOL_ERROR).map(|n| n.map(|n| n != 0));
        localisation_probability: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        mass_analyser: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        missed_cleavages: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        modified_peptide_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        nem_probabilities: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        nem_score_differences: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        number_of_matches: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        oxidation_m_probabilities: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        oxidation_m_score_differences: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        peak_coverage: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        peptide_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        precursor: usize, |location: Location, _| location.ignore("-1").parse::<usize>(NUMBER_ERROR);
        precursor_intensity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_function: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_offset: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_offset_time: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        protein_group_ids: ThinVec<usize>, |location: Location, _| location.array(';').map(|p| p.parse::<usize>(NUMBER_ERROR)).collect::<Result<ThinVec<_>,_>>();
        ration_h_l_normalised: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        ration_h_l: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        scan_event_number: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        scan_index: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        score_diff: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        simple_mass_error_ppm: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        total_ion_current: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
    }
);

/// All possible MaxQuant versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum MaxQuantVersion {
    /// msms.txt
    #[default]
    #[expect(clippy::upper_case_acronyms)]
    MSMS,
    /// msmsScans.txt
    MSMSScans,
    /// MaxNovo msmsScans.txt
    NovoMSMSScans,
    /// MaxNovo SILAC evidence.txt
    Silac,
}

impl std::fmt::Display for MaxQuantVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<MaxQuantFormat> for MaxQuantVersion {
    fn format(self) -> MaxQuantFormat {
        match self {
            Self::MSMS => MSMS,
            Self::MSMSScans => MSMS_SCANS,
            Self::NovoMSMSScans => NOVO_MSMS_SCANS,
            Self::Silac => SILAC,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::MSMS => "msms",
            Self::MSMSScans => "msmsScans",
            Self::NovoMSMSScans => "de novo msmsScans",
            Self::Silac => "SILAC evidence",
        }
    }
}

/// msms.txt
pub const MSMS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMS,
    all_modified_sequences: OptionalColumn::Required("all modified sequences"),
    base_peak_intensity: OptionalColumn::NotAvailable,
    carbamidomethyl_c_probabilities: OptionalColumn::NotAvailable,
    carbamidomethyl_c_score_differences: OptionalColumn::NotAvailable,
    collision_energy: OptionalColumn::NotAvailable,
    delta_score: OptionalColumn::Required("delta score"),
    dn_c_mass: OptionalColumn::NotAvailable,
    dn_combined_score: OptionalColumn::NotAvailable,
    dn_complete: OptionalColumn::NotAvailable,
    dn_n_mass: OptionalColumn::NotAvailable,
    dn_sequence: OptionalColumn::NotAvailable,
    evidence_id: OptionalColumn::Required("evidence id"),
    experiment: OptionalColumn::NotAvailable,
    mode: OptionalColumn::Required("fragmentation"),
    genes: OptionalColumn::NotAvailable,
    id: OptionalColumn::Required("id"),
    intensity_coverage: OptionalColumn::Required("intensity coverage"),
    intensity_h: OptionalColumn::NotAvailable,
    intensity_l: OptionalColumn::NotAvailable,
    intensity: OptionalColumn::NotAvailable,
    isotope_index: OptionalColumn::Required("isotope index"),
    labeling_state: OptionalColumn::NotAvailable,
    localisation_probability: OptionalColumn::Required("localization prob"),
    mass_analyser: OptionalColumn::Required("mass analyzer"),
    mass: OptionalColumn::Required("mass"),
    missed_cleavages: OptionalColumn::Required("missed cleavages"),
    modifications: "modifications",
    modified_peptide_id: OptionalColumn::Required("mod. peptide id"),
    mz: OptionalColumn::Required("m/z"),
    nem_probabilities: OptionalColumn::NotAvailable,
    nem_score_differences: OptionalColumn::NotAvailable,
    number_of_matches: OptionalColumn::Required("number of matches"),
    oxidation_m_probabilities: OptionalColumn::NotAvailable,
    oxidation_m_score_differences: OptionalColumn::NotAvailable,
    peak_coverage: OptionalColumn::Required("peak coverage"),
    pep: "pep",
    peptide_id: OptionalColumn::Required("peptide id"),
    peptide: "modified sequence",
    precursor_apex_function: OptionalColumn::Required("precursor apex fraction"),
    precursor_apex_offset_time: OptionalColumn::Required("precursor apex offset time"),
    precursor_apex_offset: OptionalColumn::Required("precursor apex offset"),
    precursor_intensity: OptionalColumn::Required("precursor intensity"),
    precursor: OptionalColumn::Required("precursor full scan number"),
    protein_group_ids: OptionalColumn::Required("protein group ids"),
    proteins: "proteins",
    ration_h_l_normalised: OptionalColumn::NotAvailable,
    ration_h_l: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Optional("raw file"),
    rt: OptionalColumn::Required("retention time"),
    scan_event_number: OptionalColumn::Required("scan event number"),
    scan_index: OptionalColumn::Required("scan index"),
    scan_number: "scan number",
    score_diff: OptionalColumn::Required("score diff"),
    score: "score",
    simple_mass_error_ppm: OptionalColumn::Required("simple mass error [ppm]"),
    total_ion_current: OptionalColumn::NotAvailable,
    ty: "type",
    z: "charge",
};

/// msmsScans.txt
pub const MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMSScans,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::Required("base peak intensity"),
    carbamidomethyl_c_probabilities: OptionalColumn::NotAvailable,
    carbamidomethyl_c_score_differences: OptionalColumn::NotAvailable,
    collision_energy: OptionalColumn::Required("collision energy"),
    delta_score: OptionalColumn::NotAvailable,
    dn_c_mass: OptionalColumn::NotAvailable,
    dn_combined_score: OptionalColumn::NotAvailable,
    dn_complete: OptionalColumn::NotAvailable,
    dn_n_mass: OptionalColumn::NotAvailable,
    dn_sequence: OptionalColumn::NotAvailable,
    evidence_id: OptionalColumn::NotAvailable,
    experiment: OptionalColumn::NotAvailable,
    mode: OptionalColumn::Required("fragmentation"),
    genes: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    intensity_coverage: OptionalColumn::NotAvailable,
    intensity_h: OptionalColumn::NotAvailable,
    intensity_l: OptionalColumn::NotAvailable,
    intensity: OptionalColumn::NotAvailable,
    isotope_index: OptionalColumn::NotAvailable,
    labeling_state: OptionalColumn::NotAvailable,
    localisation_probability: OptionalColumn::NotAvailable,
    mass_analyser: OptionalColumn::Required("mass analyzer"),
    mass: OptionalColumn::Required("mass"),
    missed_cleavages: OptionalColumn::NotAvailable,
    modifications: "modifications",
    modified_peptide_id: OptionalColumn::NotAvailable,
    mz: OptionalColumn::Required("m/z"),
    nem_probabilities: OptionalColumn::NotAvailable,
    nem_score_differences: OptionalColumn::NotAvailable,
    number_of_matches: OptionalColumn::NotAvailable,
    oxidation_m_probabilities: OptionalColumn::NotAvailable,
    oxidation_m_score_differences: OptionalColumn::NotAvailable,
    peak_coverage: OptionalColumn::NotAvailable,
    pep: "pep",
    peptide_id: OptionalColumn::NotAvailable,
    peptide: "modified sequence",
    precursor_apex_function: OptionalColumn::Required("precursor apex fraction"),
    precursor_apex_offset_time: OptionalColumn::Required("precursor apex offset time"),
    precursor_apex_offset: OptionalColumn::Required("precursor apex offset"),
    precursor_intensity: OptionalColumn::Required("precursor intensity"),
    precursor: OptionalColumn::Required("precursor full scan number"),
    protein_group_ids: OptionalColumn::NotAvailable,
    proteins: "proteins",
    ration_h_l_normalised: OptionalColumn::NotAvailable,
    ration_h_l: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Optional("raw file"),
    rt: OptionalColumn::Required("retention time"),
    scan_event_number: OptionalColumn::Required("scan event number"),
    scan_index: OptionalColumn::Required("scan index"),
    scan_number: "scan number",
    score_diff: OptionalColumn::NotAvailable,
    score: "score",
    simple_mass_error_ppm: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::Required("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxNovo msmsScans.txt
pub const NOVO_MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::NovoMSMSScans,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::Required("base peak intensity"),
    carbamidomethyl_c_probabilities: OptionalColumn::NotAvailable,
    carbamidomethyl_c_score_differences: OptionalColumn::NotAvailable,
    collision_energy: OptionalColumn::Required("collision energy"),
    delta_score: OptionalColumn::NotAvailable,
    dn_c_mass: OptionalColumn::Required("dn cterm mass"),
    dn_combined_score: OptionalColumn::Required("dn combined score"),
    dn_complete: OptionalColumn::Required("dn complete"),
    dn_n_mass: OptionalColumn::Required("dn nterm mass"),
    dn_sequence: OptionalColumn::Required("dn sequence"),
    evidence_id: OptionalColumn::NotAvailable,
    experiment: OptionalColumn::Optional("experiment"),
    mode: OptionalColumn::Required("fragmentation"),
    genes: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    intensity_coverage: OptionalColumn::NotAvailable,
    intensity_h: OptionalColumn::NotAvailable,
    intensity_l: OptionalColumn::NotAvailable,
    intensity: OptionalColumn::NotAvailable,
    isotope_index: OptionalColumn::NotAvailable,
    labeling_state: OptionalColumn::NotAvailable,
    localisation_probability: OptionalColumn::NotAvailable,
    mass_analyser: OptionalColumn::Required("mass analyzer"),
    mass: OptionalColumn::Required("mass"),
    missed_cleavages: OptionalColumn::NotAvailable,
    modifications: "modifications",
    modified_peptide_id: OptionalColumn::NotAvailable,
    mz: OptionalColumn::Required("m/z"),
    nem_probabilities: OptionalColumn::NotAvailable,
    nem_score_differences: OptionalColumn::NotAvailable,
    number_of_matches: OptionalColumn::NotAvailable,
    oxidation_m_probabilities: OptionalColumn::NotAvailable,
    oxidation_m_score_differences: OptionalColumn::NotAvailable,
    peak_coverage: OptionalColumn::NotAvailable,
    pep: "pep",
    peptide_id: OptionalColumn::NotAvailable,
    peptide: "modified sequence",
    precursor_apex_function: OptionalColumn::Required("precursor apex fraction"),
    precursor_apex_offset_time: OptionalColumn::Required("precursor apex offset time"),
    precursor_apex_offset: OptionalColumn::Required("precursor apex offset"),
    precursor_intensity: OptionalColumn::Required("precursor intensity"),
    precursor: OptionalColumn::Required("precursor full scan number"),
    protein_group_ids: OptionalColumn::NotAvailable,
    proteins: "proteins",
    ration_h_l_normalised: OptionalColumn::NotAvailable,
    ration_h_l: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Optional("raw file"),
    rt: OptionalColumn::Required("retention time"),
    scan_event_number: OptionalColumn::Required("scan event number"),
    scan_index: OptionalColumn::Required("scan index"),
    scan_number: "scan number",
    score_diff: OptionalColumn::NotAvailable,
    score: "score",
    simple_mass_error_ppm: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::Required("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxQuant v2.4.14.0 SILAC evidence.txt
pub const SILAC: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::Silac,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::NotAvailable,
    carbamidomethyl_c_probabilities: OptionalColumn::Required("carbamidomethyl (c) probabilities"),
    carbamidomethyl_c_score_differences: OptionalColumn::Required(
        "carbamidomethyl (c) score diffs",
    ),
    collision_energy: OptionalColumn::NotAvailable,
    delta_score: OptionalColumn::Required("delta score"),
    dn_c_mass: OptionalColumn::NotAvailable,
    dn_combined_score: OptionalColumn::NotAvailable,
    dn_complete: OptionalColumn::NotAvailable,
    dn_n_mass: OptionalColumn::NotAvailable,
    dn_sequence: OptionalColumn::NotAvailable,
    evidence_id: OptionalColumn::NotAvailable,
    experiment: OptionalColumn::Required("experiment"),
    mode: OptionalColumn::NotAvailable,
    genes: OptionalColumn::Required("gene names"),
    id: OptionalColumn::Required("id"),
    intensity: OptionalColumn::Required("intensity"),
    intensity_coverage: OptionalColumn::NotAvailable,
    intensity_h: OptionalColumn::Required("intensity h"),
    intensity_l: OptionalColumn::Required("intensity l"),
    isotope_index: OptionalColumn::NotAvailable,
    labeling_state: OptionalColumn::Required("labeling state"),
    localisation_probability: OptionalColumn::NotAvailable,
    mass_analyser: OptionalColumn::NotAvailable,
    mass: OptionalColumn::Required("mass"),
    missed_cleavages: OptionalColumn::NotAvailable,
    modifications: "modifications",
    modified_peptide_id: OptionalColumn::Required("mod. peptide id"),
    mz: OptionalColumn::Required("m/z"),
    nem_probabilities: OptionalColumn::Required("nem probabilities"),
    nem_score_differences: OptionalColumn::Required("nem score diffs"),
    number_of_matches: OptionalColumn::NotAvailable,
    oxidation_m_probabilities: OptionalColumn::Required("oxidation (m) probabilities"),
    oxidation_m_score_differences: OptionalColumn::Required("oxidation (m) score diffs"),
    peak_coverage: OptionalColumn::NotAvailable,
    pep: "pep",
    peptide_id: OptionalColumn::Required("peptide id"),
    peptide: "modified sequence",
    precursor_apex_function: OptionalColumn::NotAvailable,
    precursor_apex_offset_time: OptionalColumn::NotAvailable,
    precursor_apex_offset: OptionalColumn::NotAvailable,
    precursor_intensity: OptionalColumn::NotAvailable,
    precursor: OptionalColumn::NotAvailable,
    protein_group_ids: OptionalColumn::Required("protein group ids"),
    proteins: "proteins",
    ration_h_l: OptionalColumn::Required("ratio h/l"),
    ration_h_l_normalised: OptionalColumn::Required("ratio h/l normalized"),
    raw_file: OptionalColumn::Optional("raw file"),
    rt: OptionalColumn::Required("retention time"),
    scan_event_number: OptionalColumn::NotAvailable,
    scan_index: OptionalColumn::NotAvailable,
    scan_number: "ms/ms scan numbers",
    score_diff: OptionalColumn::NotAvailable,
    score: "score",
    simple_mass_error_ppm: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::NotAvailable,
    ty: "type",
    z: "charge",
};

impl MetaData for MaxQuantData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::MaxQuant(self.version)
    }

    fn id(&self) -> String {
        self.id
            .map_or_else(|| self.scan_number.iter().join(";"), |id| id.to_string())
    }

    fn confidence(&self) -> Option<f64> {
        (!self.score.is_nan()).then(|| 2.0 * (1.0 / (1.0 + 1.01_f64.powf(-self.score)) - 0.5))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        (self.score.is_normal() || self.score.is_subnormal()).then_some(self.score)
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
        self.rt
    }

    fn scans(&self) -> SpectrumIds {
        self.raw_file.as_ref().map_or_else(
            || {
                SpectrumIds::FileNotKnown(
                    self.scan_number
                        .iter()
                        .copied()
                        .map(SpectrumId::Number)
                        .collect(),
                )
            },
            |raw_file| {
                SpectrumIds::FileKnown(vec![(
                    raw_file.clone(),
                    self.scan_number
                        .iter()
                        .copied()
                        .map(SpectrumId::Number)
                        .collect(),
                )])
            },
        )
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.mz
    }

    fn experimental_mass(&self) -> Option<Mass> {
        self.mass
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

fn parse_de_novo_sequence(
    location: Location,
    custom_database: Option<&CustomDatabase>,
) -> Result<Peptidoform<SimpleLinear>, BoxedError<'static, BasicKind>> {
    #[derive(Debug, PartialEq, Eq)]
    enum Element {
        Either(Vec<Vec<Element>>),
        UnknownOrder(Vec<Element>),
        Sequence(SequenceElement<SemiAmbiguous>),
    }

    impl Element {
        /// Parse the next element and return the left over location
        fn parse<'a>(
            location: Location<'a>,
            custom_database: Option<&CustomDatabase>,
        ) -> Result<(Self, Location<'a>), BoxedError<'static, BasicKind>> {
            match location.as_str().as_bytes()[0] {
                b'{' => {
                    if let Some(end) = end_of_enclosure(
                        location.full_line(),
                        location.location.start + 1,
                        b'{',
                        b'}',
                    ) {
                        let mut inner_location = Location {
                            location: location.location.start + 1..end,
                            ..location.clone()
                        };
                        let mut outer = Vec::new();
                        let mut inner = Vec::new();
                        while !inner_location.is_empty() {
                            let next = Self::parse(inner_location, custom_database)?;
                            inner.push(next.0);
                            inner_location = next.1;
                            if inner_location.as_str().as_bytes().first() == Some(&b'|') {
                                inner_location.location.start += 1;
                                outer.push(inner);
                                inner = Vec::new();
                            }
                        }
                        outer.push(inner);
                        if outer.is_empty() {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MaxNovo de novo sequence",
                                "No sequence within the curly braces",
                                location.context().to_owned(),
                            ))
                        } else {
                            Ok((
                                Self::Either(outer),
                                Location {
                                    location: end + 1..location.location.end,
                                    ..location.clone()
                                },
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MaxNovo de novo sequence",
                            "The curly bracket is never closed",
                            location.context().to_owned(),
                        ))
                    }
                }
                b'[' => {
                    if let Some(end) = end_of_enclosure(
                        location.full_line(),
                        location.location.start + 1,
                        b'[',
                        b']',
                    ) {
                        let mut inner_location = Location {
                            location: location.location.start + 1..end,
                            ..location.clone()
                        };
                        let mut inner = Vec::new();
                        while !inner_location.is_empty() {
                            let next = Self::parse(inner_location, custom_database)?;
                            inner.push(next.0);
                            inner_location = next.1;
                        }
                        if inner.is_empty() {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MaxNovo de novo sequence",
                                "No sequence within the square braces",
                                location.context().to_owned(),
                            ))
                        } else {
                            Ok((
                                Self::UnknownOrder(inner),
                                Location {
                                    location: end + 1..location.location.end,
                                    ..location.clone()
                                },
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MaxNovo de novo sequence",
                            "The square bracket is never closed",
                            location.context().to_owned(),
                        ))
                    }
                }
                aa => {
                    let mut aa = SequenceElement::new(
                        AminoAcid::try_from(aa)
                            .map_err(|()| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid MaxNovo de novo sequence",
                                    "Invalid amino acid",
                                    Context::line(
                                        Some(location.line.line_index() as u32),
                                        location.full_line(),
                                        location.location.start,
                                        1,
                                    )
                                    .to_owned(),
                                )
                            })?
                            .into(),
                        None,
                    );
                    let mut offset = location.location.start + 1;

                    if location.as_str().as_bytes().get(1) == Some(&b'(') {
                        if let Some(end) = end_of_enclosure(
                            location.full_line(),
                            location.location.start + 2,
                            b'(',
                            b')',
                        ) {
                            offset = end + 1;
                            let modification = Modification::sloppy_modification(
                                location.full_line(),
                                location.location.start + 2..end,
                                Some(&aa),
                                custom_database,
                            )
                            .map_err(BoxedError::to_owned)?;
                            aa.add_simple_modification(modification);
                        } else {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MaxNovo de novo sequence",
                                "The round bracket is never closed",
                                Context::line(
                                    Some(location.line.line_index() as u32),
                                    location.full_line(),
                                    location.location.start + 1,
                                    1,
                                )
                                .to_owned(),
                            ));
                        }
                    }

                    Ok((
                        Self::Sequence(aa),
                        Location {
                            location: offset..location.location.end,
                            ..location.clone()
                        },
                    ))
                }
            }
        }

        fn add_to_peptidoform(
            self,
            peptidoform: &mut Peptidoform<SimpleLinear>,
            ambiguous_group_id: &mut NonZeroU32,
            with_id: Option<NonZeroU32>,
        ) {
            match self {
                Self::Either(mut outer) => {
                    if outer.len() == 2
                        && outer[0].len() == 1
                        && outer[1].len() == 1
                        && ((outer[0][0]
                            == Self::Sequence(SequenceElement::new(
                                CheckedAminoAcid::Isoleucine.into(),
                                None,
                            ))
                            && outer[1][0]
                                == Self::Sequence(SequenceElement::new(
                                    CheckedAminoAcid::Leucine.into(),
                                    None,
                                )))
                            || (outer[0][0]
                                == Self::Sequence(SequenceElement::new(
                                    CheckedAminoAcid::Leucine.into(),
                                    None,
                                ))
                                && outer[1][0]
                                    == Self::Sequence(SequenceElement::new(
                                        CheckedAminoAcid::Isoleucine.into(),
                                        None,
                                    ))))
                    {
                        peptidoform.sequence_mut().push(SequenceElement::new(
                            CheckedAminoAcid::AmbiguousLeucine.into(),
                            with_id,
                        ));
                    } else {
                        // FOr now just add the last branch but mark them as ambiguous
                        let id = with_id.unwrap_or(*ambiguous_group_id);
                        *ambiguous_group_id = ambiguous_group_id.saturating_add(1);
                        for el in outer.pop().unwrap() {
                            // The unwrap is fine because it is guaranteed to have at least one element
                            el.add_to_peptidoform(peptidoform, ambiguous_group_id, Some(id));
                        }
                    }
                }
                Self::UnknownOrder(inner) => {
                    let id = with_id.unwrap_or(*ambiguous_group_id);
                    *ambiguous_group_id = ambiguous_group_id.saturating_add(1);
                    for el in inner {
                        el.add_to_peptidoform(peptidoform, ambiguous_group_id, Some(id));
                    }
                }
                Self::Sequence(mut seq) => {
                    seq.ambiguous = with_id;
                    peptidoform.sequence_mut().push(seq.into());
                }
            }
        }
    }

    if !location.as_str().is_ascii() {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid MaxNovo DN sequence",
            "A de novo sequence cannot contain non-ASCII characters",
            location.context().to_owned(),
        ));
    }

    let mut peptidoform = Peptidoform::default();
    let mut ambiguous_group_id = NonZeroU32::MIN;

    let mut inner_location = location.clone();
    while !inner_location.is_empty() {
        let next = Element::parse(inner_location, custom_database)?;
        next.0
            .add_to_peptidoform(&mut peptidoform, &mut ambiguous_group_id, None);
        inner_location = next.1;
    }

    Ok(peptidoform)
}

#[cfg(test)]
mod tests {
    use crate::identification::{
        common_parser::Location, csv::CsvLine, formats::maxquant::parse_de_novo_sequence,
    };

    #[test]
    fn de_novo_sequence() {
        let expected = vec![
            ("{I|L}{I|L}C", "JJC"),
            ("[AG][CR]HHAN", "(?AG)(?CR)HHAN"),
            (
                "{M(Oxidation (M))|F}G{[QN]|[{E|Q(Deamidation (NQ))}{I|L}]|[KN]}W[C{E|Q(Deamidation (NQ))}]{D|N(Deamidation (NQ))}YQ",
                "(?F)G(?KN)W(?CQ[U:Deamidated])(?N[U:Deamidated])YQ",
            ),
        ];

        for (test, expected) in expected {
            let test = parse_de_novo_sequence(
                Location {
                    line: &CsvLine {
                        line_index: 0,
                        line: test.to_string(),
                        fields: Vec::new(),
                    },
                    location: 0..test.len(),
                    column: None,
                },
                None,
            )
            .unwrap();
            assert_eq!(test.to_string(), expected);
        }
    }
}
