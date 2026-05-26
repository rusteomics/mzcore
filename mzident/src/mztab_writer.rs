//! Write mzTab files.
//!
//! The easiest method to use is [`MzTabWriter::write`] which just needs the lists of proteins and
//! PSMs and writes the file in one go. Additionally, it is possible to write the file piecemeal
//! which allows more complex patterns with less memory overhead but needs some more setup on the
//! library users side.

use crate::{CVTerm, PSMMetaData, ProteinMetaData, SpectrumId, SpectrumIds};
use itertools::Itertools;
#[cfg(feature = "mzdata")]
use mzannotate::mzdata;
use mzcore::{
    chemistry::Chemical,
    ontology::Ontology,
    prelude::{MolecularFormula, SequencePosition},
    sequence::{
        FlankingSequence, IsAminoAcid, Modification, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};
use mzcv::{CVIndex, CVSource, Term};
use serde::{Deserialize, Serialize};
use std::{
    borrow::Cow, fmt::Display, io::Write, marker::PhantomData, num::NonZeroUsize, path::PathBuf,
};

/// Write PSMs as an mzTab file. It will always output 'Identification' type files in 'Summary' mode but does a best effort to correctly store all info.
#[derive(Debug)]
pub struct MzTabWriter<Writer, State> {
    writer: Writer,
    prh: Option<String>,
    psh: Option<String>,
    metadata: MzTabMetadata,
    state: PhantomData<State>,
}

/// The metadata for an MS run for a mzTab file.
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabMSRun {
    /// The file format terms
    pub format: Option<CVTerm>,
    /// The id format term
    pub id_format: Option<CVTerm>,
    /// The file location
    pub location: PathBuf,
    /// The file hash value
    pub hash: Option<String>,
    /// The file hash method
    pub hash_method: Option<CVTerm>,
    /// The term for the fragmentation method
    pub fragmentation_method: Option<CVTerm>,
}

/// The mzTab file has been started but nothing has been written yet.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct Initial;
/// The mzTab file has already been written with a header.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct HeaderWritten;
/// The mzTab file has already been written with a header and proteins.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct ProteinsWritten;
/// The mzTab file has already been written with a header and proteins and PSMs.
#[allow(missing_debug_implementations, missing_copy_implementations)] // Marker ZST
pub struct PSMsWritten;

trait Sealed {}
impl Sealed for HeaderWritten {}
impl Sealed for ProteinsWritten {}
impl Sealed for PSMsWritten {}

/// A trait to help indicate when proteins can be written
pub trait CanWriteProteins: Sealed {
    /// A bool to indicate if the section header (PRH) is already written
    const SECTION_HEADER_PRESENT: bool;
}
impl CanWriteProteins for HeaderWritten {
    const SECTION_HEADER_PRESENT: bool = false;
}
impl CanWriteProteins for ProteinsWritten {
    const SECTION_HEADER_PRESENT: bool = true;
}

/// A trait to help indicate when PSMs can be written
pub trait CanWritePSMs: Sealed {
    /// A bool to indicate if the section header (PSH) is already written
    const SECTION_HEADER_PRESENT: bool;
}
impl CanWritePSMs for HeaderWritten {
    const SECTION_HEADER_PRESENT: bool = false;
}
impl CanWritePSMs for ProteinsWritten {
    const SECTION_HEADER_PRESENT: bool = false;
}
impl CanWritePSMs for PSMsWritten {
    const SECTION_HEADER_PRESENT: bool = true;
}

impl<W: Write> MzTabWriter<W, Initial> {
    /// Convenience function to easily write all information to an mzTab file.
    /// An [`MSRun`] has to be given for all PSMs that are [`SpectrumIds::FileNotKnown`] to still point to the correct location.
    /// # Errors
    /// If writing to the underlying writer failed.
    pub fn write<PSM: PSMMetaData>(
        writer: W,
        mut header: MzTabMetadata,
        psms: &[PSM],
        unknown_file_run: MzTabMSRun,
    ) -> Result<(), std::io::Error> {
        let mut unknown = false;
        for run in psms
            .iter()
            .flat_map(|p| match p.scans() {
                SpectrumIds::FileKnown(scans) => scans,
                SpectrumIds::FileNotKnown(_) => {
                    unknown = true;
                    Vec::new()
                }
                SpectrumIds::None => Vec::new(),
            })
            .map(|(p, _)| p)
            .unique()
        {
            if !header.ms_runs.iter().any(|r| r.location == run) {
                header.ms_runs.push(MzTabMSRun {
                    location: run,
                    format: None,
                    id_format: None,
                    hash: None,
                    hash_method: None,
                    fragmentation_method: None,
                });
            }
        }
        let proteins = psms
            .iter()
            .flat_map(|p| p.proteins().into_owned())
            .unique_by(|p| p.id().accession().to_string())
            .collect::<Vec<PSM::Protein>>();
        if unknown {
            header.ms_runs.push(unknown_file_run);
        }
        header.protein_search_engines.extend(
            proteins
                .iter()
                .flat_map(ProteinMetaData::search_engine)
                .filter_map(|p| p.1.as_ref().map(|(_, t)| t.term.clone().into()))
                .unique(),
        );
        header.psm_search_engines.extend(
            psms.iter()
                .filter_map(|p| p.original_confidence().map(|(_, t)| t.into()))
                .unique(),
        );
        let writer = Self::new(writer, header);
        let writer = writer.write_header()?;
        let highest_used_id = psms
            .iter()
            .filter_map(PSMMetaData::numerical_id)
            .max()
            .unwrap_or_default();
        if proteins.is_empty() {
            writer.write_psms(psms, highest_used_id, &[])?;
        } else {
            let writer = writer.write_proteins(&proteins, &[])?;
            writer.write_psms(psms, highest_used_id, &[])?;
        }
        Ok(())
    }

    /// Create a new mzTab file writer that will output to the given writer and with the given MS runs.
    pub const fn new(writer: W, header: MzTabMetadata) -> Self {
        Self {
            writer,
            prh: None,
            psh: None,
            metadata: header,
            state: PhantomData,
        }
    }

    /// Write the header. Adds the standard mzTab version, mode, and type headers and writes all
    /// keys for the [`MSRun`]s. All other keys can be added as a tuple of (key, value).
    /// # Errors
    /// If the underlying writer fails.
    pub fn write_header(mut self) -> Result<MzTabWriter<W, HeaderWritten>, std::io::Error> {
        writeln!(self.writer, "MTD\tmzTab-version\t1.0.0")?;
        writeln!(self.writer, "MTD\tmzTab-mode\tSummary")?;
        writeln!(self.writer, "MTD\tmzTab-type\tIdentification")?;

        if !self.metadata.id.is_empty() {
            writeln!(self.writer, "MTD\tmzTab-ID\t{}", self.metadata.id)?;
        }
        if !self.metadata.title.is_empty() {
            writeln!(self.writer, "MTD\ttitle\t{}", self.metadata.title)?;
        }
        if !self.metadata.description.is_empty() {
            writeln!(
                self.writer,
                "MTD\tdescription\t{}",
                self.metadata.description
            )?;
        }
        for (i, terms) in self.metadata.sample_processing.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tsample_processing[{i}]\t{}",
                terms.iter().join("|"),
            )?;
        }
        for (i, instrument) in self.metadata.instruments.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tinstrument[{i}]-name\t{}",
                instrument.name
            )?;
            writeln!(
                self.writer,
                "MTD\tinstrument[{i}]-source\t{}",
                instrument.source
            )?;
            for (j, analyser) in instrument.analyser.iter().enumerate() {
                let j = j + 1;
                writeln!(
                    self.writer,
                    "MTD\tinstrument[{i}]-analyzer[{j}]\t{analyser}",
                )?;
            }
            writeln!(
                self.writer,
                "MTD\tinstrument[{i}]-detector\t{}",
                instrument.detector
            )?;
        }
        for (i, software) in self.metadata.software.iter().enumerate() {
            let i = i + 1;
            writeln!(self.writer, "MTD\tsoftware[{i}]\t{}", software.name)?;
            for (j, setting) in software.setting.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsoftware[{i}]-setting[{j}]\t{setting}")?;
            }
        }
        if !self.metadata.false_discovery_rate.is_empty() {
            writeln!(
                self.writer,
                "MTD\tfalse_discovery_rate\t{}",
                self.metadata.false_discovery_rate.iter().join("|"),
            )?;
        }
        for (i, publication) in self.metadata.publication.iter().enumerate() {
            let i = i + 1;
            writeln!(self.writer, "MTD\tpublication[{i}]\t{publication}")?;
        }
        for (i, contact) in self.metadata.contact.iter().enumerate() {
            let i = i + 1;
            if !contact.name.is_empty() {
                writeln!(self.writer, "MTD\tcontact[{i}]-name\t{}", contact.name)?;
            }
            if !contact.affiliation.is_empty() {
                writeln!(
                    self.writer,
                    "MTD\tcontact[{i}]-affiliation\t{}",
                    contact.affiliation
                )?;
            }
            if !contact.email.is_empty() {
                writeln!(self.writer, "MTD\tcontact[{i}]-email\t{}", contact.email)?;
            }
        }
        for (i, uri) in self.metadata.uri.iter().enumerate() {
            writeln!(self.writer, "MTD\turi[{i}]\t{uri}")?;
        }
        for (i, m) in self.metadata.fixed_mods.iter().enumerate() {
            write_mod(&mut self.writer, &format!("MTD\tfixed_mod[{i}]"), m)?;
        }
        for (i, m) in self.metadata.variable_mods.iter().enumerate() {
            write_mod(&mut self.writer, &format!("MTD\tvariable_mod[{i}]"), m)?;
        }
        if let Some(term) = &self.metadata.quantification_method {
            writeln!(self.writer, "MTD\tquantification_method\t{term}")?;
        }
        if let Some(term) = &self.metadata.protein_quantification_unit {
            writeln!(self.writer, "MTD\tprotein_quantification_unit\t{term}")?;
        }
        for (i, run) in self.metadata.ms_runs.iter().enumerate() {
            let i = i + 1; // 1 based
            if let Some(format) = &run.format {
                writeln!(self.writer, "MTD\tms_run[{i}]-format\t{format}")?;
            }
            if let Some(id_format) = &run.id_format {
                writeln!(self.writer, "MTD\tms_run[{i}]-id_format\t{id_format}")?;
            }
            writeln!(
                self.writer,
                "MTD\tms_run[{i}]-location\t{}{}",
                if run.location.is_absolute() {
                    "file://"
                } else {
                    ""
                },
                run.location.display()
            )?;
            if let Some(term) = &run.hash_method {
                writeln!(self.writer, "MTD\tms_run[{i}]-hash_method\t{term}")?;
            }
            if let Some(value) = &run.hash {
                writeln!(self.writer, "MTD\tms_run[{i}]-hash\t{value}")?;
            }
            if let Some(term) = &run.fragmentation_method {
                writeln!(self.writer, "MTD\tms_run[{i}]-fragmentation_method\t{term}")?;
            }
        }
        for (i, search_engine) in self.metadata.protein_search_engines.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tprotein_search_engine_score[{i}]\t{search_engine}",
            )?;
        }
        for (i, search_engine) in self.metadata.psm_search_engines.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tpsm_search_engine_score[{i}]\t{search_engine}",
            )?;
        }
        for (i, sample) in self.metadata.sample.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tsample[{i}]-description\t{}",
                sample.description
            )?;
            for (j, species) in sample.species.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsample[{i}]-species[{j}]\t{species}")?;
            }
            for (j, tissue) in sample.tissue.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsample[{i}]-tissue[{j}]\t{tissue}")?;
            }
            for (j, cell_type) in sample.cell_type.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsample[{i}]-cell_type[{j}]\t{cell_type}")?;
            }
            for (j, disease) in sample.disease.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsample[{i}]-disease[{j}]\t{disease}")?;
            }
            for (j, custom) in sample.custom.iter().enumerate() {
                let j = j + 1;
                writeln!(self.writer, "MTD\tsample[{i}]-custom[{j}]\t{custom}")?;
            }
        }
        for (i, assay) in self.metadata.assay.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tassay[{i}]-quantification_reagent\t{}",
                assay.quantification_reagent
            )?;
            for (j, m) in assay.quantification_mod.iter().enumerate() {
                let j = j + 1;
                write_mod(
                    &mut self.writer,
                    &format!("MTD\tassay[{i}]-quantification_mod[{j}]"),
                    m,
                )?;
            }
            if let Some(r) = assay.sample_ref {
                writeln!(self.writer, "MTD\tassay[{i}]-sample_ref\tsample[{r}]")?;
            }
            if let Some(r) = assay.ms_run_ref {
                writeln!(self.writer, "MTD\tassay[{i}]-sample_ref\tms_run[{r}]")?;
            }
        }
        for (i, study_variable) in self.metadata.study_variable.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tstudy_variable[{i}]-description\t{}",
                study_variable.description
            )?;
            if !study_variable.assay_refs.is_empty() {
                writeln!(
                    self.writer,
                    "MTD\tstudy_variable[{i}]-assay_refs\t{}",
                    study_variable
                        .assay_refs
                        .iter()
                        .map(|a| format!("assay[{a}]"))
                        .join(",")
                )?;
            }
            if !study_variable.sample_refs.is_empty() {
                writeln!(
                    self.writer,
                    "MTD\tstudy_variable[{i}]-sample_refs\t{}",
                    study_variable
                        .sample_refs
                        .iter()
                        .map(|a| format!("sample[{a}]"))
                        .join(",")
                )?;
            }
        }
        for (i, cv) in self.metadata.cv.iter().enumerate() {
            let i = i + 1;
            if !cv.label.is_empty() {
                writeln!(self.writer, "MTD\tcv[{i}]-label\t{}", cv.label)?;
            }
            if !cv.full_name.is_empty() {
                writeln!(self.writer, "MTD\tcv[{i}]-full_name\t{}", cv.full_name)?;
            }
            if !cv.url.is_empty() {
                writeln!(self.writer, "MTD\tcv[{i}]-url\t{}", cv.url)?;
            }
            if !cv.version.is_empty() {
                writeln!(self.writer, "MTD\tcv[{i}]-version\t{}", cv.version)?;
            }
        }
        for (name, id, value) in &self.metadata.colunit_protein {
            writeln!(self.writer, "MTD\tcolunit-protein\topt_{id}_{name}={value}")?;
        }
        writeln!(
            self.writer,
            "MTD\tcolunit-psm\tretention_time=[UO,UO:0000010,second,]"
        )?;
        for (name, id, value) in &self.metadata.colunit_psm {
            writeln!(self.writer, "MTD\tcolunit-psm\topt_{id}_{name}={value}")?;
        }
        for (i, custom) in self.metadata.custom.iter().enumerate() {
            let i = i + 1;
            writeln!(self.writer, "MTD\tcustom[{i}]\t{custom}")?;
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            prh: self.prh,
            psh: self.psh,
            metadata: self.metadata,
            state: PhantomData,
        })
    }
}

/// Define the name for an optional column in an mzTab file
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabMetadata {
    /// An ID for the file
    pub id: String,
    /// A human readable title
    pub title: String,
    /// A human readable description
    pub description: String,
    /// The fixed modifications that where searched for
    pub fixed_mods: Vec<SimpleModification>,
    /// The variable modifications that where searched for
    pub variable_mods: Vec<SimpleModification>,
    /// Defines quantification method used in the file
    pub quantification_method: Option<CVTerm>,
    /// Defines the unit in the protein quantification field
    pub protein_quantification_unit: Option<CVTerm>,
    /// The software used to analyse the data
    pub software: Vec<MzTabSoftware>,
    /// The sample processing steps these should be defined in chronological order
    pub sample_processing: Vec<Vec<CVTerm>>,
    /// The instruments used in this file
    pub instruments: Vec<MzTabInstrument>,
    /// The false discovery rates used together with a term identifying the exact method.
    pub false_discovery_rate: Vec<CVTerm>,
    /// Publications associated with this file. PubMed ids must be prefixed with 'pubmed:', DOIs with 'doi:' and identifiers can be separated with '|'.
    pub publication: Vec<String>,
    /// Any contacts for the file
    pub contact: Vec<MzTabContact>,
    /// URIs to point to the file source data, eg from PRIDE or PeptideAtlas
    pub uri: Vec<String>,
    /// Any additional custom parameters
    pub custom: Vec<CVTerm>,
    /// Define the biological samples
    pub sample: Vec<MzTabSample>,
    /// Define the study variables
    pub study_variable: Vec<MzTabStudyVariable>,
    /// Define the assays
    pub assay: Vec<MzTabAssay>,
    /// Define which CVs are used in the file
    pub cv: Vec<MzTabCV>,
    /// Define the unit for a custom protein column
    pub colunit_protein: Vec<(MzTabOptionalColumnName, MzTabObjectIdentifier, CVTerm)>,
    /// Define the unit for a custom PSM column
    pub colunit_psm: Vec<(MzTabOptionalColumnName, MzTabObjectIdentifier, CVTerm)>,
    /// The MS runs
    pub ms_runs: Vec<MzTabMSRun>,
    /// The protein search engines
    pub protein_search_engines: Vec<CVTerm>,
    /// The PSM search engines
    pub psm_search_engines: Vec<CVTerm>,
}

/// Define an instrument for an mzTab file
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabCV {
    /// What is the label (start of a CURIE) for this CV
    pub label: String,
    /// What is the full name of the CV
    pub full_name: String,
    /// What is the version of the CV
    pub version: String,
    /// Where can the CV be found
    pub url: String,
}

impl<T: CVSource> From<&CVIndex<T>> for MzTabCV {
    fn from(value: &CVIndex<T>) -> Self {
        Self {
            label: T::cv_label().to_string(),
            full_name: T::cv_name().to_string(),
            version: value.version().version.clone().unwrap_or_default(),
            url: T::files()
                .iter()
                .find_map(|f| f.url)
                .unwrap_or_default()
                .to_string(),
        }
    }
}

/// Define an instrument for an mzTab file
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabContact {
    /// The name, has to be provided in first name + initials + last name eg: Joseph J. Thomson
    pub name: String,
    /// The affiliation
    pub affiliation: String,
    /// The email
    pub email: String,
}

/// Define an instrument for an mzTab file
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabSoftware {
    /// The software and its version
    pub name: CVTerm,
    /// Any settings for this software
    pub setting: Vec<String>,
}

/// Define an instrument for an mzTab file
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabInstrument {
    /// The instrument itself eg: `MS:1000449|LTQ Orbitrap`
    pub name: CVTerm,
    /// The source eg: `MS:1000073|ESI`
    pub source: CVTerm,
    /// The analyser(s) eg: `MS:1000291|linear ion trap`
    pub analyser: Vec<CVTerm>,
    /// The detector type eg `MS:1000253|electron multiplier`
    pub detector: CVTerm,
}

/// Define a sample for an mzTab file
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabSample {
    /// The species of the sample
    pub species: Vec<CVTerm>,
    /// The tissues of the sample
    pub tissue: Vec<CVTerm>,
    /// The cell types of the sample
    pub cell_type: Vec<CVTerm>,
    /// The diseases of the sample
    pub disease: Vec<CVTerm>,
    /// Human readable description
    pub description: String,
    /// Any additional properties
    pub custom: Vec<CVTerm>,
}

/// Define a study variable for an mzTab file
#[derive(Clone, Debug, Default, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabStudyVariable {
    /// Which samples make up this study variable
    pub sample_refs: Vec<NonZeroUsize>,
    /// Which assays make up this study variable
    pub assay_refs: Vec<NonZeroUsize>,
    /// Human textual description of the variable
    pub description: String,
}

/// Define an assay for an mzTab file
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub struct MzTabAssay {
    /// Which modifications where used for quantification
    pub quantification_mod: Vec<SimpleModification>,
    /// Which sample is associated with this assay
    pub sample_ref: Option<NonZeroUsize>,
    /// Which MS run is associated with this assay
    pub ms_run_ref: Option<NonZeroUsize>,
    /// What quantification reagent was used, if not labelled `MS:1002038|unlabeled sample` should be used
    pub quantification_reagent: CVTerm,
}

impl Default for MzTabAssay {
    fn default() -> Self {
        Self {
            quantification_reagent: mzcv::term!(MS:1002038|unlabeled sample).into(),
            quantification_mod: Vec::new(),
            sample_ref: None,
            ms_run_ref: None,
        }
    }
}

/// Define the name for an optional column in an mzTab file
#[derive(Clone, Debug, Eq, Hash, PartialEq, Serialize, Deserialize)]
pub enum MzTabOptionalColumnName {
    /// A term from a CV
    Term(Term),
    /// A free text name
    Name(Cow<'static, str>),
}

impl Display for MzTabOptionalColumnName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        /// Print a name while replacing any invalid character.
        /// # Errors
        /// If writing to the writer fails.
        fn print_safe(name: &str, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
            for c in name.chars() {
                if c.is_ascii_alphanumeric() || ['_', '-', '[', ']', ':'].contains(&c) {
                    write!(f, "{c}")?;
                } else {
                    f.write_str("_")?;
                }
            }

            Ok(())
        }
        match self {
            Self::Term(term) => {
                write!(f, "{}_{}_", term.accession.cv, term.accession.accession)?;
                print_safe(&term.name, f)
            }
            Self::Name(name) => print_safe(name, f),
        }
    }
}

/// An object identifier for an mzTab optional column
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum MzTabObjectIdentifier {
    /// An assay, with the assay id
    Assay(NonZeroUsize),
    /// A study variable, with the id
    StudyVariable(NonZeroUsize),
    /// An MS run, with the id
    MSRun(NonZeroUsize),
    /// A global optional column
    Global,
}

impl Display for MzTabObjectIdentifier {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Assay(i) => write!(f, "assay[{i}]"),
            Self::StudyVariable(i) => write!(f, "study_variable[{i}]"),
            Self::MSRun(i) => write!(f, "msrun[{i}]"),
            Self::Global => write!(f, "global"),
        }
    }
}

/// An optional column that contains the information needed to create the column name and the write the values.
pub type MzTabOptionalColumn<'a, P, W> = (
    MzTabOptionalColumnName,
    MzTabObjectIdentifier,
    Box<dyn Fn(&P, &mut W) -> Result<(), std::io::Error> + 'a>,
);

impl<W: Write, State: CanWriteProteins> MzTabWriter<W, State> {
    /// Write the proteins. This can only be done if the header is already written or if some
    /// proteins where written just before. If proteins where already written before and the custom
    /// columns differed this will result in an error.
    ///
    /// To write custom columns for the proteins these need to be defined with the name and
    /// identifier to create the correct column name and with a type erased closure. The closure
    /// gets direct access to the underlying writer to prevent unnecessary allocations. This means
    /// that the closure is expected to write a valid value. The value has to be the value for only
    /// this column excluding separators. If the value contains separators inside, it should be
    /// escaped or encased by the closure itself.
    ///
    /// # Errors
    /// If the underlying writer fails. Or if proteins where previously already written with
    /// different custom columns.
    pub fn write_proteins<Protein: ProteinMetaData>(
        mut self,
        proteins: impl IntoIterator<Item = Protein>,
        custom_columns: &[MzTabOptionalColumn<Protein, W>],
    ) -> Result<MzTabWriter<W, ProteinsWritten>, std::io::Error> {
        // TODO: allow custom additional columns (sequence for example, or regions/annotations)
        let prh = format!(
            "PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tambiguity_members\tmodifications\tprotein_coverage\tgo_terms\treliability\turi{}",
            custom_columns
                .iter()
                .map(|(term, id, _)| format!("\topt_{id}_{term}"))
                .join(""),
        );
        if self.prh.is_none() {
            writeln!(self.writer, "{prh}")?;
        } else if self.prh.as_ref().is_some_and(|written| *written != prh) {
            return Err(std::io::Error::other(
                "A different header is already written",
            ));
        }
        for protein in proteins {
            // TODO: think about how to handle dynamic numbers of search engine scores
            // TODO: think about how to handle mod locations
            write!(
                self.writer,
                "PRT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                protein.id().accession(),
                protein.description().map_or("null", |t| t),
                protein
                    .species()
                    .map_or_else(|| "null".to_string(), |t| t.accession.to_string()),
                protein.species_name().map_or("null", |t| t),
                protein.database().map_or(Cow::Borrowed("null"), |t| t.0),
                protein
                    .database()
                    .and_then(|d| d.1)
                    .map_or(Cow::Borrowed("null"), |t| t),
                protein
                    .search_engine()
                    .iter()
                    .map(|(term, _)| term)
                    .join("|"),
                protein.ambiguity_members().join(","),
                protein
                    .modifications()
                    .iter()
                    .map(|(locations, m)| if locations.is_empty() {
                        make_mztab_mod(m)
                    } else {
                        format!(
                            "{}-{}",
                            locations
                                .iter()
                                .map(|(l, score)| {
                                    let i = match l {
                                        SequencePosition::NTerm => 0,
                                        SequencePosition::Index(i) => 1 + i,
                                        SequencePosition::CTerm => usize::MAX,
                                    };
                                    score.as_ref().map_or_else(|| i.to_string(), |score| format!(
                                            "{i}[MS, MS:1001876, modification probability, {score}]"
                                        ))
                                })
                                .join("|"),
                            make_mztab_mod(m)
                        )
                    })
                    .join(","),
                protein
                    .coverage()
                    .map_or_else(|| "null".to_string(), |c| c.to_string()),
                protein.gene_ontology().iter().join("|"),
                match protein.reliability() {
                    Some(crate::Reliability::High) => "1",
                    Some(crate::Reliability::Medium) => "2",
                    Some(crate::Reliability::Poor) => "3",
                    None => "null",
                },
                protein.uri().map_or("null", |t| t),
            )?;
            for (_, _, col) in custom_columns {
                write!(self.writer, "\t")?;
                col(&protein, &mut self.writer)?;
            }
            writeln!(self.writer)?;
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            prh: Some(prh),
            psh: self.psh,
            metadata: self.metadata,
            state: PhantomData,
        })
    }
}

fn write_mod(
    mut w: impl Write,
    prefix: &str,
    m: &SimpleModificationInner,
) -> Result<(), std::io::Error> {
    writeln!(w, "{prefix}\t{}", make_mztab_mod_param(m))?;
    if let SimpleModificationInner::Database { specificities, .. } = m {
        let sites = specificities
            .iter()
            .flat_map(|r| &r.0)
            .flat_map(|r| match r {
                PlacementRule::AminoAcid(aa, _) => aa.iter().map(ToString::to_string).collect(),
                PlacementRule::PsiModification(i, _) => vec![format!("MOD:{i}")],
                PlacementRule::Terminal(term) => vec![term.to_string()],
                PlacementRule::Anywhere => vec!["Anywhere".to_string()],
            })
            .unique()
            .join(", ");
        if !sites.is_empty() {
            writeln!(w, "{prefix}-site\t{sites}")?;
        }
        let position = specificities
            .iter()
            .flat_map(|r| &r.0)
            .map(|r| match r {
                PlacementRule::AminoAcid(_, pos)
                | PlacementRule::PsiModification(_, pos)
                | PlacementRule::Terminal(pos) => pos.to_string(),
                PlacementRule::Anywhere => Position::Anywhere.to_string(),
            })
            .unique()
            .join(", ");
        if !position.is_empty() {
            writeln!(w, "{prefix}-position\t{position}")?;
        }
    }
    Ok(())
}

fn make_mztab_mod_param(m: &SimpleModificationInner) -> String {
    match m {
        SimpleModificationInner::Database { formula, id, .. } => match id.ontology {
            Ontology::Unimod => format!("[UNIMOD,UNIMOD:{},{},]", id.id(), id.name),
            Ontology::Psimod => format!("[MOD,MOD:{},{},]", id.id(), id.name),
            _ => format!(
                "[CHEMMOD,{},{},{}:{}]",
                mztab_chemmod(formula),
                id.name,
                id.ontology,
                id.id()
            ),
        },
        _ => format!("[CHEMMOD,{},,]", mztab_chemmod(&m.formula())),
    }
}

fn make_mztab_mod(m: &SimpleModificationInner) -> String {
    match m {
        SimpleModificationInner::Database { formula, id, .. } => match id.ontology {
            Ontology::Unimod => format!("UNIMOD:{}", id.id()),
            Ontology::Psimod => format!("MOD:{}", id.id()),
            _ => mztab_chemmod(formula),
        },
        _ => mztab_chemmod(&m.formula()),
    }
}

fn mztab_chemmod(f: &MolecularFormula) -> String {
    if f.additional_mass() == 0.0 {
        format!("CHEMMOD:+{}", f.hill_notation_core())
    } else {
        format!("CHEMMOD:{:+}", f.monoisotopic_mass().value)
    }
}

/// An error returned when writing a mzTab file
#[derive(Debug)]
pub enum MzTabWriteError {
    /// An IO error, meaning that writing to the underlying writer failed
    IO(std::io::Error),
    /// A formatting error, meaning that writing to a string did not work
    Fmt(std::fmt::Error),
    /// No [`MSRun`] is written in the header for this file, or if None there are no files and a [`SpectrumIds::FileNotKnown`] is given
    MissingMSRun(Option<PathBuf>),
    /// This PSM search engine term is not written in the header for this file
    MissingSearchEngine(Term),
    /// PSMs were already written before but the custom columns definition is different
    MultipleDifferentPSMHeaders,
}

impl From<std::io::Error> for MzTabWriteError {
    fn from(value: std::io::Error) -> Self {
        Self::IO(value)
    }
}

impl From<std::fmt::Error> for MzTabWriteError {
    fn from(value: std::fmt::Error) -> Self {
        Self::Fmt(value)
    }
}

impl From<MzTabWriteError> for std::io::Error {
    fn from(value: MzTabWriteError) -> Self {
        match value {
            MzTabWriteError::IO(err) => err,
            a => Self::other(a),
        }
    }
}

impl Display for MzTabWriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IO(err) => write!(f, "{err}"),
            Self::Fmt(err) => write!(f, "{err}"),
            Self::MissingMSRun(None) => {
                write!(f, "Missing MS run for raw file without a defined path")
            }
            Self::MissingMSRun(Some(run)) => write!(f, "Missing MS run: {}", run.display()),
            Self::MissingSearchEngine(engine) => write!(
                f,
                "Missing search engine: {}|{}",
                engine.accession, engine.name
            ),
            Self::MultipleDifferentPSMHeaders => {
                write!(f, "A different PSM header is already written")
            }
        }
    }
}

impl std::error::Error for MzTabWriteError {}

impl<W: Write, State: CanWritePSMs> MzTabWriter<W, State> {
    /// Write the given list of PSMs to this mzTab file. It will take as much information as
    /// possible via the [`MetaData`] trait. Any cross-linked peptides (inter and intra) are
    /// ignored as these cannot be written in mzTab. Any chimeric peptidoforms are written on
    /// separate lines. Any [`SpectrumId::RetentionTime`] ids are ignored as these cannot be stored
    /// in mzTab.
    ///
    /// This can only be done if the header is already written. It can be done with and without
    /// any proteins written, and if some PSMs are already written.
    ///
    /// # Errors
    /// * If writing to the underlying writer failed.
    /// * If writing to a string for formatting failed (not expected as formatting is seen as infallible see [`std::fmt::Error`]).
    /// * If a spectrum is referenced that is not defined as a [`MSRun`].
    pub fn write_psms<PSM: PSMMetaData>(
        mut self,
        psms: impl IntoIterator<Item = PSM>,
        mut highest_used_id: usize,
        custom_columns: &[MzTabOptionalColumn<PSM, W>],
    ) -> Result<MzTabWriter<W, PSMsWritten>, MzTabWriteError> {
        let psh = format!(
            "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{}\tmodifications\tspectra_ref\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tpre\tpost\tstart\tend\treliability\turi{}",
            (1..=self.metadata.psm_search_engines.len())
                .map(|i| format!("search_engine_score[{i}]"))
                .join("\t"),
            custom_columns
                .iter()
                .map(|(term, id, _)| format!("\topt_{id}_{term}"))
                .join(""),
        );
        if self.psh.is_none() {
            writeln!(self.writer, "{psh}")?;
        } else if self.psh.as_ref().is_some_and(|written| *written != psh) {
            return Err(MzTabWriteError::MultipleDifferentPSMHeaders);
        }
        for psm in psms {
            let mut first_peptidoform = true;
            for peptidoform_ion in psm
                .peptidoform_ion_set()
                .iter()
                .flat_map(|p| p.peptidoform_ions())
            {
                if let Some(peptidoform) =
                    peptidoform_ion.singular_ref().and_then(|p| p.as_linear())
                {
                    let mut mods = String::new();
                    let mut ambiguous =
                        vec![(None, Vec::new()); peptidoform.get_ambiguous_modifications().len()];
                    let mut first_mod = true;

                    for (p, s) in peptidoform.iter(..) {
                        let i = match p.sequence_index {
                            SequencePosition::NTerm => 0,
                            SequencePosition::Index(i) => 1 + i,
                            SequencePosition::CTerm => p.sequence_length + 1,
                        };
                        for m in &s.modifications {
                            match m {
                                Modification::CrossLink { .. } => (), //skip
                                Modification::Ambiguous {
                                    id,
                                    modification,
                                    localisation_score,
                                    ..
                                } => {
                                    if ambiguous[*id].0.is_none() {
                                        ambiguous[*id].0 = Some(make_mztab_mod(modification));
                                    }
                                    ambiguous[*id].1.push((i, localisation_score));
                                }
                                Modification::Simple(m) => {
                                    use std::fmt::Write;
                                    if first_mod {
                                        first_mod = false;
                                    } else {
                                        mods.push(',');
                                    }
                                    write!(mods, "{i}-{}", make_mztab_mod(m))?;
                                }
                            }
                        }
                    }

                    for (m, locations) in ambiguous {
                        if let Some(m) = m {
                            use std::fmt::Write;
                            if first_mod {
                                first_mod = false;
                            } else {
                                mods.push(',');
                            }
                            let mut first_loc = true;
                            for (i, score) in locations {
                                if first_loc {
                                    first_loc = false;
                                } else {
                                    mods.push('|');
                                }
                                if let Some(score) = score {
                                    write!(
                                        mods,
                                        "{i}[MS, MS:1001876, modification probability, {score}]"
                                    )?;
                                } else {
                                    write!(mods, "{i}")?;
                                }
                            }
                            write!(mods, "-{m}")?;
                        }
                    }

                    write!(
                        self.writer,
                        "PSM\t{sequence}\t{psm_id}\t{accession}\t{unique}\t{database}\t{database_version}\t{search_engine}\t{search_engine_score}\t{modifications}\t{spectra_ref}\t{retention_time}\t{charge}\t{exp_mass_to_charge}\t{calc_mass_to_charge}\t{pre}\t{post}\t{start}\t{end}\t{reliability}\t{uri}",
                        sequence = peptidoform
                            .sequence()
                            .iter()
                            .filter_map(|s| s.aminoacid.one_letter_code())
                            .collect::<String>(),
                        psm_id = if first_peptidoform {
                            psm.numerical_id().map_or_else(
                                || {
                                    highest_used_id += 1;
                                    highest_used_id.to_string()
                                },
                                |id| id.to_string(),
                            )
                        } else {
                            highest_used_id += 1;
                            highest_used_id.to_string()
                        },
                        accession = psm
                            .proteins()
                            .first()
                            .map_or(Cow::Borrowed("null"), |p| p.id().accession().clone()), // TODO: contains an I assume unnecessary .clone
                        unique = psm
                            .unique()
                            .map_or_else(|| "null".to_string(), |d| d.to_string()),
                        database = psm.database().map_or("null", |d| d.0),
                        database_version = psm.database().and_then(|d| d.1).unwrap_or("null"),
                        search_engine = psm.search_engine().map_or_else(
                            || "null".to_string(),
                            |d| format!(
                                "[{}, {0}:{}, {}, ]",
                                d.accession.cv, d.accession.accession, d.name
                            )
                        ),
                        search_engine_score = {
                            if let Some((v, term)) = psm.original_confidence() {
                                let Some(pos) = self
                                    .metadata
                                    .psm_search_engines
                                    .iter()
                                    .position(|s| s.term == term)
                                else {
                                    return Err(MzTabWriteError::MissingSearchEngine(term));
                                };

                                format!(
                                    "{}{}{v}{}{}",
                                    (1..pos).map(|_| "null").join("\t"),
                                    if pos > 1 { "\t" } else { "" },
                                    if pos < self.metadata.psm_search_engines.len() - 1 {
                                        "\t"
                                    } else {
                                        ""
                                    },
                                    (pos + 1..self.metadata.psm_search_engines.len())
                                        .map(|_| "null")
                                        .join("\t")
                                )
                            } else {
                                (0..self.metadata.psm_search_engines.len())
                                    .map(|_| "null")
                                    .join("\t")
                            }
                        },
                        modifications = if mods.is_empty() { "null" } else { &mods },
                        spectra_ref = match psm.scans() {
                            SpectrumIds::None => "null".to_string(),
                            SpectrumIds::FileNotKnown(ids) =>
                                if self.metadata.ms_runs.is_empty() {
                                    return Err(MzTabWriteError::MissingMSRun(None));
                                } else {
                                    let ids = ids
                                        .iter()
                                        .filter_map(|id| match id {
                                            SpectrumId::Index(i) => {
                                                Some(format!("ms_run[1]:index={i}"))
                                            }
                                            SpectrumId::Number(n) => {
                                                Some(format!("ms_run[1]:scan={n}"))
                                            }
                                            SpectrumId::Native(n) => Some(format!("ms_run[1]:{n}")),
                                            SpectrumId::RetentionTime(_) => None,
                                        })
                                        .join("|");
                                    if ids.is_empty() {
                                        "null".to_string()
                                    } else {
                                        ids
                                    }
                                },
                            SpectrumIds::FileKnown(ids) => {
                                let mut column = String::new();
                                for (file, ids) in ids {
                                    let Some(index) = self
                                        .metadata
                                        .ms_runs
                                        .iter()
                                        .position(|run| run.location == file)
                                        .or_else(|| {
                                            file.file_name().and_then(|f| {
                                                self.metadata.ms_runs.iter().position(|run| {
                                                    run.location
                                                        .file_name()
                                                        .is_some_and(|rf| f == rf)
                                                })
                                            })
                                        })
                                    else {
                                        return Err(MzTabWriteError::MissingMSRun(Some(file)));
                                    };
                                    let index = index + 1; // 1 based
                                    let ids = ids
                                        .iter()
                                        .filter_map(|id| match id {
                                            SpectrumId::Index(i) => {
                                                Some(format!("ms_run[{index}]:index={i}"))
                                            }
                                            SpectrumId::Number(n) => {
                                                Some(format!("ms_run[{index}]:scan={n}"))
                                            }
                                            SpectrumId::Native(n) => {
                                                Some(format!("ms_run[{index}]:{n}"))
                                            }
                                            SpectrumId::RetentionTime(_) => None,
                                        })
                                        .join("|");
                                    column = format!(
                                        "{column}{}{ids}",
                                        if column.is_empty() { "" } else { "|" }
                                    );
                                }
                                if column.is_empty() {
                                    "null".to_string()
                                } else {
                                    column
                                }
                            }
                        },
                        retention_time = psm.retention_time().map_or_else(
                            || "null".to_string(),
                            |d| d.get::<mzcore::system::time::s>().to_string()
                        ),
                        charge = psm
                            .charge()
                            .map_or_else(|| "null".to_string(), |d| d.value.to_string()),
                        exp_mass_to_charge = psm
                            .experimental_mz()
                            .map_or_else(|| "null".to_string(), |d| d.value.to_string()),
                        calc_mass_to_charge = psm
                            .charge()
                            .and_then(|c| peptidoform
                                .formulas()
                                .single()
                                .map(|f| f.monoisotopic_mass() / c.to_float()))
                            .map_or_else(|| "null".to_string(), |d| d.value.to_string()),
                        pre = match psm.flanking_sequences().0 {
                            FlankingSequence::Unknown => "null".to_string(),
                            FlankingSequence::Terminal => "-".to_string(),
                            FlankingSequence::AminoAcid(aa) => aa
                                .one_letter_code()
                                .map_or_else(|| "null".to_string(), |s| s.to_string()),
                            FlankingSequence::Sequence(seq) => seq
                                .sequence()
                                .last()
                                .and_then(|s| s.aminoacid.one_letter_code())
                                .map_or_else(|| "null".to_string(), |s| s.to_string()),
                        },
                        post = match psm.flanking_sequences().1 {
                            FlankingSequence::Unknown => "null".to_string(),
                            FlankingSequence::Terminal => "-".to_string(),
                            FlankingSequence::AminoAcid(aa) => aa
                                .one_letter_code()
                                .map_or_else(|| "null".to_string(), |s| s.to_string()),
                            FlankingSequence::Sequence(seq) => seq
                                .sequence()
                                .first()
                                .and_then(|s| s.aminoacid.one_letter_code())
                                .map_or_else(|| "null".to_string(), |s| s.to_string()),
                        },
                        start = psm
                            .protein_location()
                            .map_or_else(|| "null".to_string(), |r| r.start.to_string()),
                        end = psm
                            .protein_location()
                            .map_or_else(|| "null".to_string(), |r| r.end.to_string()),
                        reliability = match psm.reliability() {
                            Some(crate::Reliability::High) => "1",
                            Some(crate::Reliability::Medium) => "2",
                            Some(crate::Reliability::Poor) => "3",
                            None => "null",
                        },
                        uri = psm.uri().unwrap_or_else(|| "null".to_string()),
                    )?;
                    for (_, _, col) in custom_columns {
                        write!(self.writer, "\t")?;
                        col(&psm, &mut self.writer)?;
                    }
                    writeln!(self.writer)?;
                }
                first_peptidoform = false;
            }
        }
        Ok(MzTabWriter {
            writer: self.writer,
            prh: self.prh,
            psh: Some(psh),
            metadata: self.metadata,
            state: PhantomData,
        })
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use std::io::BufWriter;

    use crate::{
        mztab_writer::{MzTabMSRun, MzTabMetadata, MzTabWriter},
        open_psm_file,
    };

    #[test]
    fn convert() {
        if !std::fs::exists("src/test_files_out").unwrap() {
            std::fs::create_dir("src/test_files_out").unwrap();
        }
        for file in std::fs::read_dir("src/test_files").unwrap() {
            if let Ok(entry) = file
                && entry.file_type().is_ok_and(|t| t.is_file())
            {
                // Parse the file
                let psms = open_psm_file(entry.path(), &mzcore::ontology::STATIC_ONTOLOGIES, false)
                    .unwrap()
                    .collect::<Result<Vec<_>, _>>()
                    .unwrap();

                // Write a converted mzTab file
                let new_path = std::path::Path::new("src/test_files_out")
                    .join(entry.path().with_extension("mzTab").file_name().unwrap());

                MzTabWriter::write(
                    BufWriter::new(std::fs::File::create(&new_path).unwrap()),
                    MzTabMetadata::default(),
                    &psms,
                    MzTabMSRun::default(),
                )
                .unwrap();
                println!("Wrote: {}", new_path.display());

                // Check that the new file does not produce any errors
                for psm in
                    crate::MzTabPSM::parse_file(&new_path, &mzcore::ontology::STATIC_ONTOLOGIES)
                        .unwrap()
                {
                    psm.unwrap();
                }
            }
        }
    }
}
