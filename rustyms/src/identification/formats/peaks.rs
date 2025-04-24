use std::path::{Path, PathBuf};

use crate::{
    error::CustomError,
    identification::PeaksFamilyId,
    ontologies::CustomDatabase,
    peptidoform::{SemiAmbiguous, SloppyParsingParameters},
    system::{usize::Charge, Mass, MassOverCharge, Time},
    Peptidoform,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::identification::{
    common_parser::{Location, OptionalColumn, OptionalLocation},
    csv::{parse_csv, CsvLine},
    modification::SimpleModification,
    peptidoform::PeptideModificationSearch,
    BoxedIdentifiedPeptideIter, FastaIdentifier, IdentifiedPeptidoform,
    IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion, MetaData, Modification,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Peaks line",
    "This column is not a number but it is required to be a number in this peaks format",
);
static ID_ERROR: (&str, &str) =  (
    "Invalid Peaks line",
    "This column is not a valid peaks ID but it is required to be in this peaks format\nExamples of valid IDs: '1234', 'F2:1234', 'F2:1234 12345'"
);

format_family!(
    /// The format for any Peaks file
    PeaksFormat,
    /// The data from any peaks file
    PeaksData,
    PeaksVersion, [&V12, &V11, &V11_FEATURES, &XPLUS, &AB, &X_PATCHED, &X, &DB_PEPTIDE, &DB_PSM, &DB_PROTEIN_PEPTIDE], b',', None;
    required {
        peptide: (Option<crate::AminoAcid>, Vec<Peptidoform<SemiAmbiguous>>, Option<crate::AminoAcid>), |location: Location, custom_database: Option<&CustomDatabase>| {
            let n_flanking: Option<crate::AminoAcid> =
                (location.as_str().chars().nth(1) == Some('.'))
                .then(|| location.as_str().chars().next().unwrap().try_into().map_err(|()|
                    CustomError::error(
                        "Invalid amino acid",
                        "This flanking residue is not a valid amino acid",
                        crate::error::Context::line(Some(location.line.line_index()), location.full_line(), location.location.start, 1)))).transpose()?;

            let c_flanking: Option<crate::AminoAcid> =
            (location.as_str().chars().nth_back(1) == Some('.'))
            .then(|| location.as_str().chars().next_back().unwrap().try_into().map_err(|()|
                    CustomError::error(
                        "Invalid amino acid",
                        "This flanking residue is not a valid amino acid",
                        crate::error::Context::line(Some(location.line.line_index()), location.full_line(), location.location.end-1, location.location.end)))).transpose()?;
            if c_flanking.is_none() && n_flanking.is_none() {
                location.array(';').map(|l| Peptidoform::sloppy_pro_forma(
                    l.full_line(),
                    l.location.clone(),
                    custom_database,
                    &SloppyParsingParameters::default()
                )).unique()
                .collect::<Result<Vec<_>,_>>()
                .map(|sequences| (n_flanking, sequences, c_flanking))
            } else {
            Peptidoform::sloppy_pro_forma(
                location.full_line(),
                n_flanking.map_or(location.location.start, |_| location.location.start+2)..c_flanking.map_or(location.location.end, |_| location.location.end-2),
                custom_database,
                &SloppyParsingParameters::default()
            ).map(|p| (n_flanking, vec![p], c_flanking))
        }};
        mz: MassOverCharge, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mz>);
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        area: Option<f64>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
    }
    optional {
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        ptm: Vec<SimpleModification>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.or_empty().array(';').map(|v| {
                let v = v.trim();
                Modification::sloppy_modification(v.full_line(), v.location.clone(), None, custom_database)
            }).unique().collect::<Result<Vec<_>,_>>();
        scan_number: Vec<PeaksFamilyId>, |location: Location, _| location.or_empty()
                        .map_or(Ok(Vec::new()), |l| l.array(';').map(|v| v.parse(ID_ERROR)).collect::<Result<Vec<_>,_>>());
        z: Charge, |location: Location, _| location.parse::<usize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
        alc: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        local_confidence: Vec<f64>, |location: Location, _| location
            .array(' ')
            .map(|l| l.parse::<f64>(NUMBER_ERROR))
            .collect::<Result<Vec<_>, _>>();
        fraction: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        raw_file: PathBuf, |location: Location, _| Ok(Some(Path::new(&location.get_string()).to_owned()));
        feature: PeaksFamilyId, |location: Location, _| location.or_empty().parse(ID_ERROR);
        de_novo_score: f64, |location: Location, _| location
                .parse::<f64>(NUMBER_ERROR);
        predicted_rt: Time, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|o| o.map(Time::new::<crate::system::time::min>));
        accession: String, |location: Location, _|  Ok(Some(location.get_string()));
        tag: String, |location: Location, _| Ok(location.get_string());
        mode: String, |location: Location, _| Ok(location.get_string());
        logp: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        area_tryp_ead: f64, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR);
        ascore: String, |location: Location, _| Ok(location.get_string());
        found_by: String, |location: Location, _| Ok(location.get_string());
        feature_tryp_cid: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        feature_tryp_ead: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        id: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        from_chimera: bool, |location: Location, _| Ok(location.get_string().eq_ignore_ascii_case("yes"));
        unique: bool, |location: Location, _| Ok(location.get_string() == "Y");
        protein_group: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        protein_id: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        protein_accession: FastaIdentifier<String>, |location: Location, _|  location.parse(NUMBER_ERROR).map(Some);
        start: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        end: usize, |location: Location, _| location.parse(NUMBER_ERROR).map(Some);
        quality: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        rt_begin: Time, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|o| o.map(Time::new::<crate::system::time::s>));
        rt_end: Time, |location: Location, _| location.or_empty().parse::<f64>(NUMBER_ERROR).map(|o| o.map(Time::new::<crate::system::time::s>));
        precursor_id: isize, |location: Location, _| location.parse::<isize>(NUMBER_ERROR);
        k0_range: std::ops::RangeInclusive<f64>, |location: Location, _| location.split_once('-').map(|(start, end)| Ok(start.parse(NUMBER_ERROR)?..=end.parse(NUMBER_ERROR)?));
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        // Add the meaningful modifications to replace mass modifications
        if let Some(ptm) = parsed.ptm.clone() {
            for pep in &mut parsed.peptide.1 {
                *pep = PeptideModificationSearch::in_modifications(ptm.clone())
                    .tolerance(super::Tolerance::Absolute(super::system::da(0.05)))
                    .search(pep.clone());
            }
        }
        Ok(parsed)
    }
);

impl From<PeaksData> for IdentifiedPeptidoform {
    fn from(value: PeaksData) -> Self {
        Self {
            score: value
                .de_novo_score
                .or(value.alc)
                .map(|v| v / 100.0)
                .or_else(|| {
                    value
                        .logp
                        .map(|v| 2.0 * (1.0 / (1.0 + 1.025_f64.powf(-v)) - 0.5))
                }),
            local_confidence: value
                .local_confidence
                .as_ref()
                .map(|lc| lc.iter().map(|v| *v / 100.0).collect()),
            metadata: MetaData::Peaks(value),
        }
    }
}

/// An older version of a PEAKS export
pub const X: PeaksFormat = PeaksFormat {
    version: PeaksVersion::X,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag (>=0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::NotAvailable,
    feature: OptionalColumn::NotAvailable,
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version X of PEAKS export (made for build 31 January 2019)
pub const X_PATCHED: PeaksFormat = PeaksFormat {
    version: PeaksVersion::XPatched,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag (>=0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::Required("fraction"),
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("feature"),
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version X+ of PEAKS export (made for build 20 November 2019)
pub const XPLUS: PeaksFormat = PeaksFormat {
    version: PeaksVersion::XPlus,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag (>=0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::Required("fraction"),
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("feature"),
    de_novo_score: OptionalColumn::Required("denovo score"),
    predicted_rt: OptionalColumn::Required("predict rt"),
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version 11 of PEAKS export
pub const V11: PeaksFormat = PeaksFormat {
    version: PeaksVersion::V11,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag(>=0.0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("feature id"),
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::Optional("precursor id"),
    k0_range: OptionalColumn::Optional("1/k0 range"),
};
/// Version 11 of PEAKS export
pub const V11_FEATURES: PeaksFormat = PeaksFormat {
    version: PeaksVersion::V11Features,
    scan_number: OptionalColumn::NotAvailable,
    peptide: "denovo peptide",
    alc: OptionalColumn::Required("alc (%)"),
    quality: OptionalColumn::Required("quality"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::NotAvailable,
    rt: "rt",
    rt_begin: OptionalColumn::Required("rt begin"),
    rt_end: OptionalColumn::Required("rt end"),
    area: "area",
    ptm: OptionalColumn::NotAvailable,
    local_confidence: OptionalColumn::NotAvailable,
    tag: OptionalColumn::NotAvailable,
    mode: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("feature id"),
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version 12 of PEAKS export
pub const V12: PeaksFormat = PeaksFormat {
    version: PeaksVersion::V12,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag(>=0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("feature id"),
    de_novo_score: OptionalColumn::Required("deep novo score (%)"),
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version Ab of PEAKS export
pub const AB: PeaksFormat = PeaksFormat {
    version: PeaksVersion::Ab,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::Required("alc (%)"),
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::Required("local confidence (%)"),
    tag: OptionalColumn::Required("tag (>=0%)"),
    mode: OptionalColumn::Required("mode"),
    fraction: OptionalColumn::NotAvailable,
    raw_file: OptionalColumn::NotAvailable,
    feature: OptionalColumn::NotAvailable,
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::Required("accession"),
    ascore: OptionalColumn::NotAvailable,
    found_by: OptionalColumn::NotAvailable,
    logp: OptionalColumn::NotAvailable,
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version DB peptide of PEAKS export
pub const DB_PEPTIDE: PeaksFormat = PeaksFormat {
    version: PeaksVersion::DBPeptide,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::NotAvailable,
    mz: "m/z",
    z: OptionalColumn::NotAvailable,
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area tryp-cid",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::NotAvailable,
    tag: OptionalColumn::NotAvailable,
    mode: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::Required("fraction"),
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::Required("#feature"),
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::Required("accession"),
    ascore: OptionalColumn::Required("ascore"),
    found_by: OptionalColumn::Required("found by"),
    logp: OptionalColumn::Required("-10lgp"),
    feature_tryp_cid: OptionalColumn::Required("#feature tryp-cid"),
    feature_tryp_ead: OptionalColumn::Required("#feature tryp ead"),
    area_tryp_ead: OptionalColumn::Required("area tryp ead"),
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version DB psm of PEAKS export
pub const DB_PSM: PeaksFormat = PeaksFormat {
    version: PeaksVersion::DBPSM,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::NotAvailable,
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::NotAvailable,
    tag: OptionalColumn::NotAvailable,
    mode: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::Required("fraction"),
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::NotAvailable,
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::Required("accession"),
    ascore: OptionalColumn::Required("ascore"),
    found_by: OptionalColumn::Required("found by"),
    logp: OptionalColumn::Required("-10lgp"),
    feature_tryp_cid: OptionalColumn::NotAvailable,
    feature_tryp_ead: OptionalColumn::NotAvailable,
    area_tryp_ead: OptionalColumn::NotAvailable,
    id: OptionalColumn::Required("id"),
    from_chimera: OptionalColumn::Required("from chimera"),
    unique: OptionalColumn::NotAvailable,
    protein_group: OptionalColumn::NotAvailable,
    protein_id: OptionalColumn::NotAvailable,
    protein_accession: OptionalColumn::NotAvailable,
    start: OptionalColumn::NotAvailable,
    end: OptionalColumn::NotAvailable,
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};
/// Version DB protein peptide of PEAKS export
/// protein group, protein id, protein accession, unique, start, end,
pub const DB_PROTEIN_PEPTIDE: PeaksFormat = PeaksFormat {
    version: PeaksVersion::DBProteinPeptide,
    scan_number: OptionalColumn::Required("scan"),
    peptide: "peptide",
    alc: OptionalColumn::NotAvailable,
    mz: "m/z",
    z: OptionalColumn::Required("z"),
    mass: OptionalColumn::Required("mass"),
    rt: "rt",
    area: "area tryp-cid",
    ptm: OptionalColumn::Required("ptm"),
    local_confidence: OptionalColumn::NotAvailable,
    tag: OptionalColumn::NotAvailable,
    mode: OptionalColumn::NotAvailable,
    fraction: OptionalColumn::Required("fraction"),
    raw_file: OptionalColumn::Required("source file"),
    feature: OptionalColumn::NotAvailable,
    de_novo_score: OptionalColumn::NotAvailable,
    predicted_rt: OptionalColumn::NotAvailable,
    accession: OptionalColumn::NotAvailable,
    ascore: OptionalColumn::Required("ascore"),
    found_by: OptionalColumn::Required("found by"),
    logp: OptionalColumn::Required("-10lgp"),
    feature_tryp_cid: OptionalColumn::Required("#feature tryp-cid"),
    feature_tryp_ead: OptionalColumn::Required("#feature tryp ead"),
    area_tryp_ead: OptionalColumn::Required("area tryp ead"),
    id: OptionalColumn::NotAvailable,
    from_chimera: OptionalColumn::NotAvailable,
    unique: OptionalColumn::Required("unique"),
    protein_group: OptionalColumn::Required("protein group"),
    protein_id: OptionalColumn::Required("protein id"),
    protein_accession: OptionalColumn::Required("protein accession"),
    start: OptionalColumn::Required("start"),
    end: OptionalColumn::Required("end"),
    quality: OptionalColumn::NotAvailable,
    rt_begin: OptionalColumn::NotAvailable,
    rt_end: OptionalColumn::NotAvailable,
    precursor_id: OptionalColumn::NotAvailable,
    k0_range: OptionalColumn::NotAvailable,
};

/// All possible peaks versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
#[expect(clippy::upper_case_acronyms)]
pub enum PeaksVersion {
    /// An older version of a PEAKS export
    X,
    /// Version X of PEAKS export (made for build 31 January 2019)
    XPatched,
    /// Version X+ of PEAKS export (made for build 20 November 2019)
    XPlus,
    /// Version DB peptide of PEAKS export
    DBPeptide,
    /// Version DB PSM of PEAKS export
    DBPSM,
    /// Version DB Protein Peptide of PEAKS export
    DBProteinPeptide,
    /// Version Ab of PEAKS export
    Ab,
    /// Version 11 denovo file
    V11,
    /// Version 11 features file
    V11Features,
    /// Version 12
    #[default]
    V12,
}

impl std::fmt::Display for PeaksVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PeaksFormat> for PeaksVersion {
    fn format(self) -> PeaksFormat {
        match self {
            Self::X => X,
            Self::XPatched => X_PATCHED,
            Self::XPlus => XPLUS,
            Self::DBPeptide => DB_PEPTIDE,
            Self::DBPSM => DB_PSM,
            Self::DBProteinPeptide => DB_PROTEIN_PEPTIDE,
            Self::Ab => AB,
            Self::V11 => V11,
            Self::V11Features => V11_FEATURES,
            Self::V12 => V12,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::X => "X",
            Self::XPatched => "X patched",
            Self::XPlus => "X+",
            Self::DBPeptide => "DB peptide",
            Self::DBPSM => "DB PSM",
            Self::DBProteinPeptide => "DB protein peptide",
            Self::Ab => "Ab",
            Self::V11 => "11",
            Self::V11Features => "11 features",
            Self::V12 => "12",
        }
    }
}
