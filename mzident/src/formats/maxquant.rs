use std::{
    borrow::Cow, marker::PhantomData, num::NonZeroU32, ops::Range, path::PathBuf, str::FromStr,
};

use itertools::Itertools;
#[cfg(feature = "mzannotate")]
use mzannotate::prelude::AnnotatedPeak;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    BoxedIdentifiedPeptideIter, FastaIdentifier, KnownFileFormat, MaybePeptidoform, PSM, PSMData,
    PSMFileFormatVersion, PSMMetaData, PSMSource, SpectrumId, SpectrumIds,
    common_parser::{Location, OptionalColumn, OptionalLocation},
    helper_functions::end_of_enclosure,
};
use mzcore::{
    csv::{CsvLine, parse_csv},
    ontology::Ontologies,
    sequence::{
        AminoAcid, CheckedAminoAcid, CompoundPeptidoformIon, FlankingSequence, Modification,
        Peptidoform, SemiAmbiguous, SequenceElement, SimpleLinear, SimpleModificationInner,
        SloppyParsingParameters,
    },
    system::{MassOverCharge, Time, dalton, f32::Mass, isize::Charge},
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
    /// This can contain data from the database match and the _de novo_ match at the same time when
    /// run with MaxNovo. In that case the _de novo_ data will not be shown via the methods of the
    /// [`PSMMetaData`] trait. If access is needed solely to the _de novo_ data and not to the database
    /// data the easiest way is detecting this case and overwriting the data in place.
    /// ```rust
    /// # use mzcore::{prelude::*, sequence::Linked};
    /// # use mzident::*;
    /// # let mut psm: PSM<Linked, PeptidoformPresent> = BasicCSVPSM::default().into();
    /// if let PSMData::MaxQuant(ref mut mq) = psm.data
    ///     && mq.format().version() == Some(MaxQuantVersion::NovoMSMSScans.to_string())
    /// {
    ///     mq.peptide = mq.dn_sequence.clone();
    ///     psm.score = mq.dn_combined_score.map(|v| f64::from(v) / 100.0);
    ///     mq.score = mq.dn_combined_score.map_or(0.0, f64::from);
    /// }
    /// ```
    MaxQuant,
    SimpleLinear, MaybePeptidoform, [&MSMS, &NOVO_MSMS_SCANS, &MSMS_SCANS, &SILAC], b'\t', None;
    required {
        scan_number: ThinVec<usize>, |location: Location, _| location.or_empty().array(';').map(|s| s.parse(NUMBER_ERROR)).collect::<Result<ThinVec<usize>, BoxedError<'_, BasicKind>>>();
        proteins: ThinVec<FastaIdentifier<Box<str>>>, |location: Location, _| location.array(';').map(|v|
            FastaIdentifier::<Box<str>>::from_str(v.as_str())
            .map_err(|e| BoxedError::new(
                BasicKind::Error,
                "Could not parse MaxQuant line",
                format!("The protein identifier could not be parsed: {e}"),
                v.context().to_owned()
            ))).collect::<Result<ThinVec<FastaIdentifier<Box<str>>>, _>>();
        /// The database matched peptide, annotated as [`SimpleLinear`] to allow replacing it with the _de novo_ peptide, no features of the [`SimpleLinear`] complexity are used for the database peptides
        peptide: Option<Peptidoform<SimpleLinear>>, |location: Location, ontologies: &Ontologies| location.or_empty().parse_with(|location| Peptidoform::sloppy_pro_forma_inner(
            &location.base_context(),
            location.full_line(),
            location.location.clone(),
            ontologies,
            &SloppyParsingParameters::default()
        ).map_err(BoxedError::to_owned))
        .map(|p| p.map(Into::into));
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<mzcore::system::e>);
        ty: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        pep: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        score: f32, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        all_modified_sequences: ThinVec<Peptidoform<SemiAmbiguous>>, |location: Location, ontologies: &Ontologies| location.array(';')
                .map(|s| Peptidoform::sloppy_pro_forma_inner(&s.base_context(), s.line.line(), s.location.clone(), ontologies, &SloppyParsingParameters::default()).map_err(BoxedError::to_owned))
                .collect::<Result<ThinVec<Peptidoform<SemiAmbiguous>>, BoxedError<'static, BasicKind>>>();
        base_peak_intensity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        collision_energy: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        delta_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        dn_c_mass: Mass, |location: Location, _| location.parse::<f32>(NUMBER_ERROR).map(Mass::new::<dalton>);
        /// The combined score is a combination of the complete score and the gap score. Both of these scores are ranked, and then the sum of the two ranks is taken and normalized to lie between 0 and 100.
        dn_combined_score: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        dn_complete: bool, |location: Location, _| Ok(location.as_str() == "+");
        dn_n_mass: Mass, |location: Location, _| location.parse::<f32>(NUMBER_ERROR).map(Mass::new::<dalton>);
        dn_sequence: Peptidoform<SimpleLinear>, |location: Location, ontologies: &Ontologies| location.or_empty().map(|l| parse_de_novo_sequence(l, ontologies));
        evidence_id: u32, |location: Location, _| location.parse::<u32>(NUMBER_ERROR);
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
        mass: Mass, |location: Location, _| location.parse::<f32>(NUMBER_ERROR).map(Mass::new::<dalton>);
        missed_cleavages: u8, |location: Location, _| location.parse::<u8>(NUMBER_ERROR);
        modified_peptide_id: u32, |location: Location, _| location.parse::<u32>(NUMBER_ERROR);
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<mzcore::system::thomson>);
        number_of_matches: u16, |location: Location, _| location.parse::<u16>(NUMBER_ERROR);
        peak_coverage: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        peptide_id: u32, |location: Location, _| location.parse::<u32>(NUMBER_ERROR);
        precursor: u32, |location: Location, _| location.ignore("-1").parse::<u32>(NUMBER_ERROR);
        precursor_intensity: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_function: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_offset: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        precursor_apex_offset_time: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        protein_group_ids: ThinVec<usize>, |location: Location, _| location.array(';').map(|p| p.parse::<usize>(NUMBER_ERROR)).collect::<Result<ThinVec<_>,_>>();
        ration_h_l_normalised: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        ration_h_l: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| Ok(PathBuf::from(location.get_string()));
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<mzcore::system::time::min>);
        scan_event_number: u32, |location: Location, _| location.parse::<u32>(NUMBER_ERROR);
        scan_index: u32, |location: Location, _| location.parse::<u32>(NUMBER_ERROR);
        score_diff: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        simple_mass_error_ppm: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        total_ion_current: f32, |location: Location, _| location.parse::<f32>(NUMBER_ERROR);
        #[cfg(feature = "mzannotate")]
        spectrum: ThinVec<AnnotatedPeak<mzannotate::fragment::Fragment>>, |_, _| None;
    }

    fn post_process(source: &CsvLine, mut parsed: Self, _ontologies: &Ontologies) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Some(dn_sequence) = parsed.dn_sequence.as_mut() {
            if let Some(n_mass) = parsed.dn_n_mass && n_mass != Mass::default() && !n_mass.is_nan() {
                dn_sequence.add_simple_n_term(SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, n_mass.into(), None).into());
            }
            if let Some(c_mass) = parsed.dn_c_mass && c_mass != Mass::default() && !c_mass.is_nan() {
                dn_sequence.add_simple_c_term(SimpleModificationInner::Mass(mzcore::sequence::MassTag::None, c_mass.into(), None).into());
            }
        }

        #[cfg(feature = "mzannotate")]
        {
            if let Some(peptidoform) = &parsed.peptide && let Ok((annotation, match_range)) = source.index_column("matches") && let Ok((intensities, intensity_range)) = source.index_column("intensities") && let Ok((mzs, mz_range)) = source.index_column("masses") && !annotation.is_empty() && !intensities.is_empty() && !mzs.is_empty(){
                let mut peaks = ThinVec::with_capacity(annotation.chars().filter(|c| *c == ';').count());
                let mut match_offset = 0;
                let mut intensity_offset = 0;
                let mut mz_offset = 0;
                for (index, ((annotation, intensity), mz)) in annotation.split(';').zip(intensities.split(';')).zip(mzs.split(';')).enumerate() {
                    let annotation_parsed = mzannotate::fragment::Fragment::maxquant_inner(
                        &source.range_context(match_range.start+match_offset..match_range.start+match_offset+annotation.len(), None),
                        &source.line,
                        match_range.start+match_offset..match_range.start+match_offset+annotation.len(),
                        peptidoform,
                    ).map_err(BoxedError::to_owned)?;
                    let mz_num = mz.parse::<f64>().map_err(|e| BoxedError::new(BasicKind::Error, "Invalid m/z", format!("A peak m/z is expected to be a number but: {e}"), source.range_context(mz_range.start+mz_offset..mz_range.start+mz_offset+mz.len(), None).to_owned()))?;
                    let intensity_num = intensity.parse::<f32>().map_err(|e| BoxedError::new(BasicKind::Error, "Invalid intensity", format!("A peak intensity is expected to be a number but: {e}"), source.range_context(intensity_range.start+intensity_offset..intensity_range.start+intensity_offset+intensity.len(), None).to_owned()))?;
                    peaks.push(AnnotatedPeak {
                        mz: MassOverCharge::new::<mzcore::system::thomson>(mz_num),
                        intensity: intensity_num,
                        index: index as u32,
                        annotations: vec![annotation_parsed],
                        aggregations: Vec::new()});
                    match_offset += annotation.len() + 1;
                    intensity_offset += intensity.len() + 1;
                    mz_offset += mz.len() + 1;
                }
                parsed.spectrum = Some(peaks);
            }
        }

        Ok(parsed)
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

impl PSMFileFormatVersion<MaxQuantFormat> for MaxQuantVersion {
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
    missed_cleavages: OptionalColumn::Optional("missed cleavages"),
    modified_peptide_id: OptionalColumn::Required("mod. peptide id"),
    mz: OptionalColumn::Required("m/z"),
    number_of_matches: OptionalColumn::Required("number of matches"),
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
    #[cfg(feature = "mzannotate")]
    spectrum: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::NotAvailable,
    ty: "type",
    z: "charge",
};

/// msmsScans.txt
pub const MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::MSMSScans,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::Required("base peak intensity"),
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
    modified_peptide_id: OptionalColumn::NotAvailable,
    mz: OptionalColumn::Required("m/z"),
    number_of_matches: OptionalColumn::NotAvailable,
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
    #[cfg(feature = "mzannotate")]
    spectrum: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::Required("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxNovo msmsScans.txt
pub const NOVO_MSMS_SCANS: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::NovoMSMSScans,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::Required("base peak intensity"),
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
    modified_peptide_id: OptionalColumn::NotAvailable,
    mz: OptionalColumn::Required("m/z"),
    number_of_matches: OptionalColumn::NotAvailable,
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
    #[cfg(feature = "mzannotate")]
    spectrum: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::Required("total ion current"),
    ty: "type",
    z: "charge",
};

/// MaxQuant v2.4.14.0 SILAC evidence.txt
pub const SILAC: MaxQuantFormat = MaxQuantFormat {
    version: MaxQuantVersion::Silac,
    all_modified_sequences: OptionalColumn::NotAvailable,
    base_peak_intensity: OptionalColumn::NotAvailable,
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
    modified_peptide_id: OptionalColumn::Required("mod. peptide id"),
    mz: OptionalColumn::Required("m/z"),
    number_of_matches: OptionalColumn::NotAvailable,
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
    #[cfg(feature = "mzannotate")]
    spectrum: OptionalColumn::NotAvailable,
    total_ion_current: OptionalColumn::NotAvailable,
    ty: "type",
    z: "charge",
};

impl PSMMetaData for MaxQuantPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptide.as_ref().map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::MaxQuant(self.version)
    }

    fn numerical_id(&self) -> Option<usize> {
        self.id
    }

    fn id(&self) -> String {
        self.id
            .map_or_else(|| self.scan_number.iter().join(";"), |id| id.to_string())
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        Some(mzcv::term!(MS:1001583|MaxQuant))
    }

    fn confidence(&self) -> Option<f64> {
        (!self.score.is_nan())
            .then(|| 2.0 * (1.0 / (1.0 + 1.01_f64.powf(-f64::from(self.score))) - 0.5))
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        (self.score.is_normal() || self.score.is_subnormal())
            .then_some((f64::from(self.score), mzcv::term!(MS:1001171|Mascot:score)))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        self.mode.as_deref().map(Cow::Borrowed)
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

    fn experimental_mass(&self) -> Option<mzcore::system::f64::Mass> {
        self.mass
            .map(|v| mzcore::system::f64::Mass::new::<dalton>(v.value.into()))
    }

    type Protein = FastaIdentifier<Box<str>>;

    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(&self.proteins)
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

    fn unique(&self) -> Option<bool> {
        None
    }

    fn reliability(&self) -> Option<crate::Reliability> {
        None
    }

    fn uri(&self) -> Option<String> {
        None
    }

    #[cfg(feature = "mzannotate")]
    fn annotated_spectrum(&self) -> Option<Cow<'_, mzannotate::spectrum::AnnotatedSpectrum>> {
        // TODO: create a generic function to convert a 'MetaData' into a spectrum description + analytes etc
        self.spectrum.as_ref().map(|s| {
            use mzannotate::mzspeclib::Analyte;

            Cow::Owned(mzannotate::spectrum::AnnotatedSpectrum::new(
                0,
                mzannotate::mzdata::spectrum::SpectrumDescription::default(),
                vec![Vec::new(); 1],
                self.peptide
                    .as_ref()
                    .map(|pep| {
                        vec![Analyte {
                            id: 0,
                            target: mzannotate::mzspeclib::AnalyteTarget::PeptidoformIon(
                                pep.clone().into(),
                            ),
                            proteins: Vec::new(),
                            params: Vec::new(),
                        }]
                    })
                    .unwrap_or_default(),
                Vec::new(),
                Vec::from(s.clone()).into(),
            ))
        })
    }

    #[cfg(feature = "mzannotate")]
    fn has_annotated_spectrum(&self) -> bool {
        self.spectrum.is_some()
    }
}

/// Parse a MaxNovo sequence
/// # Errors
/// If the sequence does not follow the MaxNovo format
///
/// # Specification (provided by MaxQuant)
/// Three different kinds of brackets (round, square, and curly) can be found in a regex de-novo-
/// sequence. The round brackets “()” contain the modification of the amino acid that is located
/// on the left side of the brackets. For example M(Oxidation (M)). The square brackets “[]”
/// contain two amino acids that their order cannot be distinguished by MaxNovo due to the
/// absence of fragment ion peaks in the MS/MS scan. For example [WR] can be either WR or
/// RW. The curly brackets “{}” contain amino acids that are separated by the pipe character “|”.
/// The separated by “|” amino acids are equally possible to be at that position because they are
/// isobaric or almost isobaric (up to mass tolerance error). For example {I|L} can be either I or L.
#[expect(clippy::needless_pass_by_value)]
fn parse_de_novo_sequence(
    location: Location,
    ontologies: &Ontologies,
) -> Result<Peptidoform<SimpleLinear>, BoxedError<'static, BasicKind>> {
    #[derive(Debug, Eq, PartialEq)]
    enum Element {
        Either(Vec<Vec<Self>>),
        UnknownOrder(Vec<Self>),
        Sequence(SequenceElement<SemiAmbiguous>),
    }

    impl Element {
        /// Parse the next element and return the left over location
        /// # Errors
        /// If the next section could not be parsed as a single sequence element
        #[expect(clippy::needless_pass_by_value)]
        fn parse<'a>(
            location: Location<'a>,
            ontologies: &Ontologies,
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
                            let next = Self::parse(inner_location, ontologies)?;
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
                            let next = Self::parse(inner_location, ontologies)?;
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
                                ontologies,
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

        /// Add this sequence element to the given existing peptidoform
        #[allow(clippy::missing_panics_doc)] // Cannot panic
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
        let next = Element::parse(inner_location, ontologies)?;
        next.0
            .add_to_peptidoform(&mut peptidoform, &mut ambiguous_group_id, None);
        inner_location = next.1;
    }

    Ok(peptidoform)
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{common_parser::Location, formats::maxquant::parse_de_novo_sequence};

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
                    line: &mzcore::csv::CsvLine {
                        line_index: 0,
                        line: test.to_string(),
                        fields: Vec::new(),
                        file: None,
                    },
                    location: 0..test.len(),
                    column: None,
                },
                &mzcore::ontology::STATIC_ONTOLOGIES,
            )
            .unwrap();
            assert_eq!(test.to_string(), expected);
        }
    }
}
