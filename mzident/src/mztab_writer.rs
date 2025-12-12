//! Write mzTab files.
//!
//! The easiest method to use is [`MzTabWriter::write`] which just needs the lists of proteins and
//! PSMs and writes the file in one go. Additionally, it is possible to write the file piecemeal
//! which allows more complex patterns with less memory overhead but needs some more setup on the
//! library users side.

use crate::{PSMMetaData, ProteinMetaData, SpectrumId, SpectrumIds};
use itertools::Itertools;
use mzcore::{
    chemistry::Chemical,
    ontology::Ontology,
    prelude::{MolecularFormula, SequencePosition},
    sequence::{FlankingSequence, IsAminoAcid, Modification, SimpleModificationInner},
};
use mzcv::Term;
use std::{io::Write, marker::PhantomData, path::PathBuf};

/// Write identified peptidoform ions as an mzTab file. It will always output 'Identification' type files in 'Summary' mode but does a best effort to correctly store all info.
#[derive(Debug)]
pub struct MzTabWriter<Writer, State> {
    writer: Writer,
    ms_runs: Vec<MSRun>,
    protein_search_engines: Vec<Term>,
    psm_search_engines: Vec<Term>,
    state: PhantomData<State>,
}

/// The metadata for an MS run for a mzTab file.
#[derive(Clone, Debug, Default)]
pub struct MSRun {
    /// The file format and id format terms
    pub format: Option<(Term, Term)>,
    /// The file location
    pub location: PathBuf,
    /// The file hash (term, value)
    pub hash: Option<(Term, String)>,
    /// The term for the fragmentation method
    pub fragmentation_method: Option<Term>,
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
    pub fn write<PSM: PSMMetaData, Protein: ProteinMetaData>(
        writer: W,
        header: &[(String, String)],
        proteins: &[Protein],
        psms: &[PSM],
        unknown_file_run: MSRun,
    ) -> Result<(), std::io::Error> {
        let mut unknown = false;
        let mut ms_runs = psms
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
            .map(|p| MSRun {
                location: p,
                format: None,
                hash: None,
                fragmentation_method: None,
            })
            .collect::<Vec<_>>();
        if unknown {
            ms_runs.insert(0, unknown_file_run);
        }
        let protein_search_engines = proteins
            .iter()
            .flat_map(|p| p.search_engine())
            .filter_map(|p| p.1.as_ref().map(|(_, t)| t.term.clone()))
            .unique()
            .collect::<Vec<_>>();
        let psm_search_engines = psms
            .iter()
            .filter_map(|p| p.original_confidence().map(|(_, t)| t))
            .unique()
            .collect::<Vec<_>>();
        let writer = Self::new(writer, ms_runs, protein_search_engines, psm_search_engines);
        let writer = writer.write_header(header)?;
        match if proteins.is_empty() {
            writer.write_psms(psms)
        } else {
            let writer = writer.write_proteins(proteins)?;
            writer.write_psms(psms)
        } {
            Ok(_) => Ok(()),
            Err(MzTabWriteError::IO(err)) => Err(err),
            Err(MzTabWriteError::Fmt(err)) => Err(std::io::Error::other(err)),
            Err(MzTabWriteError::MissingMSRun(err)) => {
                unreachable!("A MSRun was not defined in the MzTabWriter::write for path {err:?}")
            }
            Err(MzTabWriteError::MissingSearchEngine(term)) => {
                unreachable!(
                    "A search engine term was not defined in the MzTabWriter::write for term {term:?}"
                )
            }
        }
    }

    /// Create a new mzTab file writer that will output to the given writer and with the given MS runs.
    pub const fn new(
        writer: W,
        ms_runs: Vec<MSRun>,
        protein_search_engines: Vec<Term>,
        psm_search_engines: Vec<Term>,
    ) -> Self {
        Self {
            writer,
            ms_runs,
            protein_search_engines,
            psm_search_engines,
            state: PhantomData,
        }
    }

    /// Write the header. Adds the standard mzTab version, mode, and type headers and writes all
    /// keys for the [`MSRun`]s. All other keys can be added as a tuple of (key, value).
    /// # Errors
    /// If the underlying writer fails.
    pub fn write_header(
        mut self,
        header: &[(String, String)],
    ) -> Result<MzTabWriter<W, HeaderWritten>, std::io::Error> {
        writeln!(self.writer, "MTD\tmzTab-version\t1.0.0")?;
        writeln!(self.writer, "MTD\tmzTab-mode\tSummary")?;
        writeln!(self.writer, "MTD\tmzTab-type\tIdentification")?;
        for (key, value) in header {
            writeln!(self.writer, "MTD\t{key}\t{value}")?;
        }
        for (i, run) in self.ms_runs.iter().enumerate() {
            let i = i + 1; // 1 based
            if let Some((format, id_format)) = &run.format {
                writeln!(
                    self.writer,
                    "MTD\tms_run[{i}]-format\t[{0}, {}:{}, {}, ]",
                    format.accession.cv, format.accession.accession, format.name
                )?;
                writeln!(
                    self.writer,
                    "MTD\tms_run[{i}]-id_format\t[{0}, {}:{}, {}, ]",
                    id_format.accession.cv, id_format.accession.accession, id_format.name
                )?;
            }
            writeln!(
                self.writer,
                "MTD\tms_run[{i}]-location\tfile://{}",
                run.location.display()
            )?;
            if let Some((term, value)) = &run.hash {
                writeln!(
                    self.writer,
                    "MTD\tms_run[{i}]-hash_method\t[{0}, {}:{}, {}, ]",
                    term.accession.cv, term.accession.accession, term.name
                )?;
                writeln!(self.writer, "MTD\tms_run[{i}]-hash\t{value}",)?;
            }
            if let Some(term) = &run.fragmentation_method {
                writeln!(
                    self.writer,
                    "MTD\tms_run[{i}]-fragmentation_method\t[{0}, {}:{}, {}, ]",
                    term.accession.cv, term.accession.accession, term.name
                )?;
            }
        }
        for (i, search_engine) in self.psm_search_engines.iter().enumerate() {
            let i = i + 1;
            writeln!(
                self.writer,
                "MTD\tpsm_search_engine_score[{i}]\t[{0}, {}:{}, {}, ]",
                search_engine.accession.cv, search_engine.accession.accession, search_engine.name
            )?;
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            protein_search_engines: self.protein_search_engines,
            psm_search_engines: self.psm_search_engines,
            state: PhantomData,
        })
    }
}

impl<W: Write, State: CanWriteProteins> MzTabWriter<W, State> {
    /// Write the proteins. This can only be done if the header is already written or if some proteins where written just before.
    /// # Errors
    /// If the underlying writer fails.
    pub fn write_proteins<Protein: ProteinMetaData>(
        mut self,
        proteins: &[Protein],
    ) -> Result<MzTabWriter<W, ProteinsWritten>, std::io::Error> {
        if !State::SECTION_HEADER_PRESENT {
            writeln!(
                self.writer,
                "PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tambiguity_members\tmodifications\tprotein_coverage\tgo_terms\treliability\turi"
            )?;
        }
        for protein in proteins {
            // TODO: think about how to handle dynamic numbers of search engine scores
            // TODO: think about how to handle mod locations
            writeln!(
                self.writer,
                "PRT\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                protein.id().accession(),
                protein.description().map_or("null", |t| t),
                protein
                    .species()
                    .map_or_else(|| "null".to_string(), |t| t.accession.to_string()),
                protein.species_name().map_or("null", |t| t),
                protein.database().map_or("null", |t| t.0),
                protein.database().and_then(|d| d.1).map_or("null", |t| t),
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
                                    if let Some(score) = score {
                                        format!(
                                            "{i}[MS, MS:1001876, modification probability, {score}]"
                                        )
                                    } else {
                                        i.to_string()
                                    }
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
                protein
                    .reliability()
                    .map_or_else(|| "null".to_string(), |c| c.to_string()),
                protein.uri().map_or("null", |t| t),
            )?;
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            protein_search_engines: self.protein_search_engines,
            psm_search_engines: self.psm_search_engines,
            state: PhantomData,
        })
    }
}

fn make_mztab_mod(m: &SimpleModificationInner) -> String {
    match m {
        SimpleModificationInner::Database { formula, id, .. } => id.id().map_or_else(
            || mztab_chemmod(formula),
            |index| match id.ontology {
                Ontology::Unimod => format!("UNIMOD:{index}"),
                Ontology::Psimod => format!("MOD:{index}"),
                _ => mztab_chemmod(formula),
            },
        ),
        _ => mztab_chemmod(&m.formula()),
    }
}

fn mztab_chemmod(f: &MolecularFormula) -> String {
    if f.additional_mass() == 0.0 {
        format!("CHEMMOD:{}", f.hill_notation_core())
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
        psms: &[PSM],
    ) -> Result<MzTabWriter<W, PSMsWritten>, MzTabWriteError> {
        if !State::SECTION_HEADER_PRESENT {
            // TODO: think about how to handle local confidences
            writeln!(
                self.writer,
                "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\t{}\tmodifications\tspectra_ref\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tpre\tpost\tstart\tend\treliability\turi",
                (1..=self.psm_search_engines.len())
                    .map(|i| format!("search_engine_score[{i}]"))
                    .join("\t"),
            )?;
        }
        let mut highest_id = psms
            .iter()
            .filter_map(PSMMetaData::numerical_id)
            .max()
            .unwrap_or_default();
        for psm in psms {
            let mut first_peptidoform = true;
            for peptidoform_ion in psm
                .compound_peptidoform_ion()
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

                    writeln!(
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
                                    highest_id += 1;
                                    highest_id.to_string()
                                },
                                |id| id.to_string(),
                            )
                        } else {
                            highest_id += 1;
                            highest_id.to_string()
                        },
                        accession = psm
                            .protein_names()
                            .as_ref()
                            .and_then(|n| n.first())
                            .map_or_else(|| "null".to_string(), |c| c.accession().clone()),
                        unique = psm
                            .unique()
                            .map_or_else(|| "null".to_string(), |d| d.to_string()),
                        database = psm.database().map_or_else(|| "null", |d| d.0),
                        database_version = psm
                            .database()
                            .and_then(|d| d.1)
                            .map_or_else(|| "null", |d| d),
                        search_engine = psm.search_engine().map_or_else(
                            || "null".to_string(),
                            |d| format!(
                                "[{}, {0}:{}, {}, ]",
                                d.accession.cv, d.accession.accession, d.name
                            )
                        ),
                        search_engine_score = {
                            if let Some((v, term)) = psm.original_confidence() {
                                let Some(pos) =
                                    self.psm_search_engines.iter().position(|s| *s == term)
                                else {
                                    return Err(MzTabWriteError::MissingSearchEngine(term));
                                };

                                format!(
                                    "{}{}{v}{}{}",
                                    (1..pos).map(|_| "null").join("\t"),
                                    if pos > 1 { "\t" } else { "" },
                                    if pos < self.psm_search_engines.len() - 1 {
                                        "\t"
                                    } else {
                                        ""
                                    },
                                    (pos + 1..self.psm_search_engines.len())
                                        .map(|_| "null")
                                        .join("\t")
                                )
                            } else {
                                (0..self.psm_search_engines.len())
                                    .map(|_| "null")
                                    .join("\t")
                            }
                        },
                        modifications = if mods.is_empty() { "null" } else { &mods },
                        spectra_ref = match psm.scans() {
                            SpectrumIds::None => "null".to_string(),
                            SpectrumIds::FileNotKnown(ids) =>
                                if self.ms_runs.is_empty() {
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
                                        .ms_runs
                                        .iter()
                                        .position(|run| run.location == file)
                                        .or_else(|| {
                                            file.file_name().and_then(|f| {
                                                self.ms_runs.iter().position(|run| {
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
                        uri = psm.uri().map_or_else(|| "null".to_string(), |r| r),
                    )?;
                }
                first_peptidoform = false;
            }
        }
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            protein_search_engines: self.protein_search_engines,
            psm_search_engines: self.psm_search_engines,
            state: PhantomData,
        })
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use std::io::BufWriter;

    use crate::{
        MzTabProtein,
        mztab_writer::{MSRun, MzTabWriter},
        open_identified_peptidoforms_file,
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
                let psms = open_identified_peptidoforms_file(
                    entry.path(),
                    &mzcore::ontology::STATIC_ONTOLOGIES,
                    false,
                )
                .unwrap()
                .collect::<Result<Vec<_>, _>>()
                .unwrap();

                // Write a converted mzTab file
                let new_path = std::path::Path::new("src/test_files_out")
                    .join(entry.path().with_extension("mzTab").file_name().unwrap());

                MzTabWriter::write::<_, MzTabProtein>(
                    BufWriter::new(std::fs::File::create(&new_path).unwrap()),
                    &[],
                    &[],
                    &psms,
                    MSRun::default(),
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
