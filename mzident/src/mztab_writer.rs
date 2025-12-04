use crate::{MetaData, Protein, SpectrumIds};
use itertools::Itertools;
use mzcore::{
    chemistry::Chemical,
    ontology::Ontology,
    prelude::{MolecularFormula, SequencePosition},
    sequence::{IsAminoAcid, Modification, SimpleModification, SimpleModificationInner},
};
use mzcv::Term;
use std::{io::Write, marker::PhantomData, path::PathBuf, sync::Arc};

/// Write identified peptidoform ions as an mzTab file. It will always output 'Identification' type files in 'Summary' mode but does a best effort to correctly store all info.
#[derive(Debug)]
pub struct MzTabWriter<Writer, State> {
    writer: Writer,
    ms_runs: Vec<MSRun>,
    state: PhantomData<State>,
}

#[derive(Debug, Clone)]
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

// #[sealed]
pub trait CanWriteProteins {
    const SECTION_HEADER_PRESENT: bool;
}
impl CanWriteProteins for HeaderWritten {
    const SECTION_HEADER_PRESENT: bool = false;
}
impl CanWriteProteins for ProteinsWritten {
    const SECTION_HEADER_PRESENT: bool = true;
}

// #[sealed]
pub trait CanWritePSMs {
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
    // TODO: detect headers from PSMs
    pub fn write_header(
        mut self,
        header: &[(String, String)],
        ms_runs: Vec<MSRun>,
    ) -> Result<MzTabWriter<W, HeaderWritten>, std::io::Error> {
        self.ms_runs = ms_runs;
        writeln!(self.writer, "MTD\tmzTab-version\t1.0.0")?;
        writeln!(self.writer, "MTD\tmzTab-mode\tSummary")?;
        writeln!(self.writer, "MTD\tmzTab-type\tIdentification")?;
        for (key, value) in header {
            writeln!(self.writer, "MTD\t{key}\t{value}")?;
        }
        for (i, run) in self.ms_runs.iter().enumerate() {
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
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            state: PhantomData,
        })
    }
}

impl<W: Write, State: CanWriteProteins> MzTabWriter<W, State> {
    pub fn write_proteins(
        mut self,
        proteins: &[(String, Arc<Protein>)],
    ) -> Result<MzTabWriter<W, ProteinsWritten>, std::io::Error> {
        if !State::SECTION_HEADER_PRESENT {
            writeln!(
                self.writer,
                "PRH\taccession\tdescription\ttaxid\tspecies\tdatabase\tdatabase_version\tsearch_engine\tambiguity_members\tmodifications\tprotein_coverage\tgo_terms\treliability\turi"
            )?;
        }
        for (accession, protein) in proteins {
            // TODO: think about how to handle dynamic numbers of search engine scores
            writeln!(
                self.writer,
                "PRT\t{accession}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                protein.description,
                protein.taxid,
                protein.species,
                protein.database,
                protein.database_version,
                protein.search_engine.iter().map(|(term, _)| term).join("|"),
                protein.ambiguity_members.join(","),
                protein.modifications,
                protein
                    .coverage
                    .map_or_else(|| "null".to_string(), |c| c.to_string()),
                protein
                    .go_terms
                    .iter()
                    .map(|id| format!("GO:{id}"))
                    .join("|"),
                protein
                    .reliability
                    .map_or_else(|| "null".to_string(), |c| c.to_string()),
                protein
                    .uri
                    .as_ref()
                    .map_or_else(|| "null".to_string(), ToString::to_string),
            )?;
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            state: PhantomData,
        })
    }
}

fn make_mztab_mod(m: SimpleModification) -> String {
    match &*m {
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

#[derive(Debug)]
pub enum MzTabPSMWriteError {
    IO(std::io::Error),
    FMT(std::fmt::Error),
    /// No MSRun is written in the header for this file, or if none there are no files and a [`SpectrumIds::FileNotKnown`] is given
    MissingMSRun(Option<PathBuf>),
}

impl From<std::io::Error> for MzTabPSMWriteError {
    fn from(value: std::io::Error) -> Self {
        Self::IO(value)
    }
}
impl From<std::fmt::Error> for MzTabPSMWriteError {
    fn from(value: std::fmt::Error) -> Self {
        Self::FMT(value)
    }
}

impl<W: Write, State: CanWritePSMs> MzTabWriter<W, State> {
    /// Write the given list of PSMs to this mzTab file. It will take as much information as possible via the [`MetaData`] trait. Any cross-linked peptides (inter and intra) are ignored as these cannot be written in mzTab. Any chimeric peptidoforms are written on separate lines.
    pub fn write_psms<PSM: MetaData>(
        mut self,
        psms: impl Iterator<Item = PSM>,
    ) -> Result<MzTabWriter<W, PSMsWritten>, MzTabPSMWriteError> {
        if !State::SECTION_HEADER_PRESENT {
            // TODO: think about how to handle dynamic numbers of search engine scores
            // TODO: think about how to handle local confidences
            writeln!(
                self.writer,
                "PSH\tsequence\tPSM_ID\taccession\tunique\tdatabase\tdatabase_version\tsearch_engine\tsearch_engine_score[1]\tmodifications\tspectra_ref\tretention_time\tcharge\texp_mass_to_charge\tcalc_mass_to_charge\tpre\tpost\tstart\tend\reliability\turi"
            )?;
        }
        for psm in psms {
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
                                        ambiguous[*id].0 =
                                            Some(make_mztab_mod(modification.clone()));
                                    }
                                    ambiguous[*id].1.push((i, localisation_score));
                                }
                                Modification::Simple(m) => {
                                    use std::fmt::Write;
                                    if !first_mod {
                                        mods.push(',');
                                        first_mod = false;
                                    }
                                    write!(mods, "{i}-{}", make_mztab_mod(m.clone()))?;
                                }
                            }
                        }
                    }

                    for (m, locations) in ambiguous {
                        if let Some(m) = m {
                            use std::fmt::Write;
                            if !first_mod {
                                mods.push(',');
                                first_mod = false;
                            }
                            let mut first_loc = true;
                            for (i, score) in locations {
                                if !first_loc {
                                    mods.push('|');
                                    first_loc = false;
                                }
                                if let Some(score) = score {
                                    write!(
                                        mods,
                                        "{i}[MS,MS:1001876, modification probability, {score}]"
                                    )?;
                                } else {
                                    write!(mods, "{i}")?;
                                }
                            }
                            write!(mods, "-{m}")?;
                        }
                    }

                    // TODO: unique and search_engine are missing
                    writeln!(
                        self.writer,
                        "PSM\t{}\t{}\t{}\tnull\t{}\t{}\tnull\t{}\t{mods}\t{}\t{}\t{}\t{}\tcalc_mass_to_charge\tpre\tpost\tstart\tend\reliability\turi",
                        peptidoform
                            .sequence()
                            .iter()
                            .filter_map(|s| s.aminoacid.one_letter_code())
                            .collect::<String>(),
                        psm.id(), // TODO: this could be duplicated when multiple chimeric peptides are given
                        psm.protein_names()
                            .as_ref()
                            .and_then(|n| n.first())
                            .map_or_else(|| "null".to_string(), |c| c.accession().clone()),
                        psm.database().map_or_else(|| "null", |d| d.0),
                        psm.original_confidence()
                            .map_or_else(|| "null".to_string(), |d| d.to_string()),
                        psm.database()
                            .and_then(|d| d.1)
                            .map_or_else(|| "null", |d| d),
                        match psm.scans() {
                            SpectrumIds::None => "null",
                            SpectrumIds::FileNotKnown(ids) => "null", // TODO: fix
                            SpectrumIds::FileKnown(ids) => "null",
                        },
                        psm.retention_time().map_or_else(
                            || "null".to_string(),
                            |d| d.get::<mzcore::system::time::s>().to_string()
                        ),
                        psm.charge()
                            .map_or_else(|| "null".to_string(), |d| d.value.to_string()),
                        psm.experimental_mz()
                            .map_or_else(|| "null".to_string(), |d| d.value.to_string()),
                        // TODO: add rest peptidoform.formulas()
                    )?;
                }
            }
        }
        writeln!(self.writer)?;
        Ok(MzTabWriter {
            writer: self.writer,
            ms_runs: self.ms_runs,
            state: PhantomData,
        })
    }
}
