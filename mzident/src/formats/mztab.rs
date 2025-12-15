use std::{
    borrow::Cow,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    ops::Range,
    str::FromStr,
    sync::Arc,
};

use context_error::*;
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use mzcv::{ControlledVocabulary, Term};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    FastaIdentifier, KnownFileFormat, MaybePeptidoform, PSM, PSMData, PSMMetaData, ProteinMetaData,
    SpectrumId, SpectrumIds,
    helper_functions::{
        check_extension, explain_number_error, float_digits, next_number, split_with_brackets,
    },
};
use mzcore::{
    chemistry::MolecularFormula,
    quantities::Tolerance,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, FlankingSequence, MUPSettings,
        PeptideModificationSearch, Peptidoform, SequencePosition, SimpleLinear, SimpleModification,
        SimpleModificationInner, SloppyParsingParameters,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge, usize},
};
use mzcv::Curie;

/// Peptidoform data from a mzTab file
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct MzTabPSM {
    /// The sequence corresponding to the PSM
    pub peptidoform: Option<Peptidoform<SimpleLinear>>,
    /// A unique identifier for a PSM within the file. If a PSM can be matched to
    /// multiple proteins, the same PSM should be represented on multiple rows with
    /// different accessions and the same PSM_ID.
    pub id: usize,
    /// The protein's accession the corresponding peptide sequence (coming from the
    /// PSM) is associated with.
    pub protein: Option<(String, Option<Arc<MzTabProtein>>)>,
    /// Indicates whether the peptide sequence (coming from the PSM) is unique for
    /// this protein in respect to the searched database.
    pub unique: Option<bool>,
    /// The protein database used for the search together with the version of the
    /// database.
    pub database: Option<(String, Option<String>)>,
    /// The search engines that identified this peptide, alongside their identified score and the CV term describing the score
    pub search_engine: Vec<(CVTerm, Option<f64>, CVTerm)>,
    /// The estimated reliability of the PSM. (Optional parameter)
    pub reliability: Option<Reliability>,
    /// The retention time for this peptide.
    pub rt: Option<Time>,
    /// The charge for this peptide.
    pub z: Charge,
    /// The experimental mz
    pub mz: Option<MassOverCharge>,
    /// A URI pointing to the PSMs entry in the experiment it was identified in (e.g. the peptidoform PRIDE entry).
    pub uri: Option<String>,
    /// The spectra references grouped by raw file
    pub spectra_ref: SpectrumIds,
    /// The flanking sequences around this peptide
    pub flanking_sequence: (FlankingSequence, FlankingSequence),
    /// The start of this peptide in the containing protein (0-based)
    pub protein_location: Option<Range<u16>>,
    /// Casanovo specific additional metadata with the amino acid confidence
    pub local_confidence: Option<Vec<f64>>,
    /// Any additional metadata
    pub additional: HashMap<String, String>,
}

impl mzcore::space::Space for MzTabPSM {
    fn space(&self) -> mzcore::space::UsedSpace {
        (self.peptidoform.space()
            + self.id.space()
            + self.protein.space()
            + self.unique.space()
            + self.database.space()
            + self.search_engine.space()
            + self.reliability.space()
            + self.rt.space()
            + self.z.space()
            + self.mz.space()
            + self.uri.space()
            + self.spectra_ref.space()
            + self.flanking_sequence.space()
            + self.protein_location.space()
            + self.local_confidence.space()
            + self.additional.space())
        .set_total::<Self>()
    }
}

impl MzTabPSM {
    /// Parse a mzTab file.
    /// # Errors
    /// If the file is not in the correct format
    pub fn parse_file(
        path: impl AsRef<std::path::Path>,
        ontologies: &mzcore::ontology::Ontologies,
    ) -> Result<
        Box<dyn Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + '_>,
        BoxedError<'static, BasicKind>,
    > {
        let file = File::open(path.as_ref()).map_err(|e| {
            BoxedError::new(
                BasicKind::Error,
                "Could not open file",
                e.to_string(),
                Context::default()
                    .source(path.as_ref().to_string_lossy())
                    .to_owned(),
            )
        })?;
        if check_extension(path, "gz") {
            Ok(Box::new(Self::parse_reader(
                BufReader::new(GzDecoder::new(BufReader::new(file))),
                ontologies,
            )))
        } else {
            Ok(Box::new(Self::parse_reader(
                BufReader::new(file),
                ontologies,
            )))
        }
    }

    /// Parse a mzTab file directly from a buffered reader
    pub fn parse_reader<'a, T: BufRead + 'a>(
        reader: T,
        ontologies: &'a mzcore::ontology::Ontologies,
    ) -> impl Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + 'a {
        let mut search_engine_score_type: Vec<CVTerm> = Vec::new();
        let mut modifications: Vec<SimpleModification> = Vec::new();
        let mut raw_files: Vec<(Option<String>, Option<CVTerm>, Option<CVTerm>)> = Vec::new(); //path, file format, identifier type
        let mut peptide_header: Option<Vec<String>> = None;
        let mut protein_header: Option<Vec<String>> = None;
        let mut proteins: HashMap<String, Arc<MzTabProtein>> = HashMap::new();

        parse_mztab_reader(reader).filter_map(move |item| {
            item.transpose().and_then(|item| match item {
                Ok(MzTabLine::MTD(line_index, line, fields)) => {
                    if fields.len() == 3 {
                        match line[fields[1].clone()].to_ascii_lowercase().as_str() {
                            m if (m.starts_with("variable_mod[")
                                || m.starts_with("fixed_mod["))
                                && m.ends_with(']') =>
                            {
                                if line[fields[2].clone()].starts_with('[')
                                    && line[fields[2].clone()].ends_with(']')
                                {
                                    let text_range = fields[2].start + 1..fields[2].end - 1;
                                    let mut split = line[text_range.clone()].split(',');
                                    let first = split.next();
                                    let offset = split.next();
                                    let m_range = text_range.start
                                        + first.map_or(0, |f| f.len() + 1)
                                        + offset.map_or(0, |s| {
                                            s.bytes().take_while(u8::is_ascii_whitespace).count()
                                        })
                                        ..text_range.start
                                            + first.map_or(0, |f| f.len() + 1)
                                            + offset.map_or(0, str::len);
                                    if &line[m_range.clone()] != "MS:1002453"
                                        && &line[m_range.clone()] != "MS:1002454"
                                    {
                                        match parse_single_modification(
                                            &line,
                                            m_range,
                                            ontologies,
                                            line_index,
                                        ) {
                                            Ok(m) => {
                                                if !modifications.contains(&m) {
                                                    modifications.push(m);
                                                }
                                            }
                                            Err(e) => return Some(Err(e.to_owned())),
                                        }
                                    }
                                } else {
                                    return Some(Err(BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid mzTab modification",
                                        "A modification has to be enclosed by square brackets",
                                        Context::line_range(
                                            Some(line_index),
                                            &line,
                                            fields[2].clone(),
                                        )
                                        .to_owned(),
                                    )));
                                }
                            }
                            m if m.starts_with("psm_search_engine_score[") && m.ends_with(']') => {
                                match CVTerm::from_str(&line[fields[2].clone()]) {
                                    Ok(term) => search_engine_score_type.push(term),
                                    Err(err) => return Some(Err(err)),
                                }
                            }
                            m if m.starts_with("ms_run[") && m.ends_with("]-location") => {
                                let index = match m
                                    .trim_start_matches("ms_run[")
                                    .trim_end_matches("]-location")
                                    .parse::<usize>()
                                    .map_err(|err| {
                                        BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid mzTab ms_run identifier",
                                            format!(
                                                "The ms_run identifier {}",
                                                explain_number_error(&err)
                                            ),
                                            Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned(),
                                        )
                                    }) {
                                        Ok(0) => return Some(Err(BoxedError::new(BasicKind::Error, "Invalid mzTab mz_run identifier", "A ms_run identifier uses 1 based indexing and so cannot be 0.", Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned()))),
                                        Ok(i) => i - 1,
                                        Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].0 = Some(line[fields[2].clone()].to_string());
                            }
                            m if m.starts_with("ms_run[") && m.ends_with("]-format") => {
                                let index = match m
                                    .trim_start_matches("ms_run[")
                                    .trim_end_matches("]-format")
                                    .parse::<usize>()
                                    .map_err(|err| {
                                        BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid mzTab ms_run identifier",
                                            format!(
                                                "The ms_run identifier {}",
                                                explain_number_error(&err)
                                            ),
                                            Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned(),
                                        )
                                    }) {
                                        Ok(0) => return Some(Err(BoxedError::new(BasicKind::Error, "Invalid mzTab mz_run identifier", "A ms_run identifier uses 1 based indexing and so cannot be 0.", Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned()))),
                                        Ok(i) => i - 1,
                                        Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].1 =
                                    Some(match CVTerm::from_str(&line[fields[2].clone()]) {
                                        Ok(i) => i,
                                        Err(err) => return Some(Err(err)),
                                    });
                            }
                            m if m.starts_with("ms_run[") && m.ends_with("]-id_format") => {
                                let index = match m
                                    .trim_start_matches("ms_run[")
                                    .trim_end_matches("]-id_format")
                                    .parse::<usize>()
                                    .map_err(|err| {
                                        BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid mzTab ms_run identifier",
                                            format!(
                                                "The ms_run identifier {}",
                                                explain_number_error(&err)
                                            ),
                                            Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned(),
                                        )
                                    }) {
                                        Ok(0) => return Some(Err(BoxedError::new(BasicKind::Error, "Invalid mzTab mz_run identifier", "A ms_run identifier uses 1 based indexing and so cannot be 0.", Context::line_range(
                                                Some(line_index),
                                                &line,
                                                fields[1].clone(),
                                            )
                                            .to_owned()))),
                                        Ok(i) => i - 1,
                                        Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].2 =
                                    Some(match CVTerm::from_str(&line[fields[2].clone()]) {
                                        Ok(i) => i,
                                        Err(err) => return Some(Err(err)),
                                    });
                            }
                            _ => (),
                        }
                        None
                    } else {
                        Some(Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MTD line",
                            "MTD lines should contain three columns (the tag, key, and value)",
                            Context::full_line(line_index, line),
                        )))
                    }
                }
                Ok(MzTabLine::PRH(line_index, line, fields)) => {
                    let header = fields
                        .into_iter()
                        .map(|field| line[field].to_ascii_lowercase())
                        .collect_vec();
                    // not checked: best_search_engine_score[n]
                    for required in [
                        "accession",
                        "description",
                        "taxid",
                        "species",
                        "database",
                        "database_version",
                        "search_engine",
                        "ambiguity_members",
                        "modifications",
                    ] {
                        if !header.contains(&required.to_string()) {
                            return Some(Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid protein table",
                                format!("The required column '{required}' is not present"),
                                Context::full_line(line_index, line),
                            )));
                        }
                    }
                    protein_header = Some(header);
                    None
                }
                Ok(MzTabLine::PRT(line_index, line, fields)) => {
                    match PSMLine::new(line_index, protein_header.as_deref(), &line, &fields)
                        .and_then(|line| MzTabProtein::from_line(line))
                    {
                        Ok(protein) => {
                            for name in &protein.ambiguity_members {
                                proteins.insert(name.clone(), protein.clone());
                            }
                            proteins.insert(protein.accession.clone(), protein);
                            None
                        }
                        Err(err) => Some(Err(err.to_owned())),
                    }
                }

                Ok(MzTabLine::PSH(line_index, line, fields)) => {
                    let header = fields
                        .into_iter()
                        .map(|field| line[field].to_ascii_lowercase())
                        .collect_vec();
                    // optional: opt_*, reliability, uri,
                    // not checked: search_engine_score[n]
                    for required in [
                        "sequence",
                        "psm_id",
                        "accession",
                        "unique",
                        "database",
                        "database_version",
                        "search_engine",
                        "modifications",
                        "retention_time",
                        "charge",
                        "exp_mass_to_charge",
                        "spectra_ref",
                        "pre",
                        "post",
                        "start",
                        "end",
                    ] {
                        if !header.contains(&required.to_string()) {
                            return Some(Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid peptide table",
                                format!("The required column '{required}' is not present"),
                                Context::full_line(line_index, line),
                            )));
                        }
                    }
                    peptide_header = Some(header);
                    None
                }
                Ok(MzTabLine::PSM(line_index, line, fields)) => Some(
                    PSMLine::new(line_index, peptide_header.as_deref(), &line, &fields)
                        .and_then(|line| {
                            Self::from_line(
                                line,
                                &modifications,
                                &search_engine_score_type,
                                &raw_files,
                                ontologies,
                                &proteins,
                            )
                        })
                        .map_err(BoxedError::to_owned),
                ),
                Err(e) => Some(Err(e)),
            })
        })
    }

    /// Parse a single PSM line
    /// # Errors
    /// When not in the correct format
    #[expect(clippy::missing_panics_doc)]
    fn from_line<'a>(
        line: PSMLine<'a>,
        global_modifications: &[SimpleModification],
        search_engine_score_types: &[CVTerm],
        raw_files: &[(Option<String>, Option<CVTerm>, Option<CVTerm>)],
        ontologies: &mzcore::ontology::Ontologies,
        proteins: &HashMap<String, Arc<MzTabProtein>>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let (mod_column, mod_range) = line
            .required_column("modifications")
            .map_err(BoxedError::to_owned)?;
        let modifications: Vec<MzTabReturnModification> =
            if mod_column.eq_ignore_ascii_case("null") || mod_column == "0" {
                Vec::new()
            } else {
                split_with_brackets(line.line, mod_range, b',', b'[', b']')
                    .into_iter()
                    .map(|field| parse_modification(line.line, field, ontologies, line.line_index))
                    .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()?
            };

        let mut result = Self {
            peptidoform: {
                let range = line
                    .required_column("sequence")
                    .map_err(BoxedError::to_owned)?
                    .1;

                if range.is_empty() {
                    None
                } else {
                    let mut peptide: Peptidoform<SimpleLinear> = Peptidoform::pro_forma_or_sloppy(
                        line.line,
                        range.clone(),
                        ontologies,
                        &SloppyParsingParameters {
                            allow_unwrapped_modifications: true,
                            ..Default::default()
                        },
                    )?
                    .into();
                    for modification in modifications {
                        match modification {
                            MzTabReturnModification::Defined(location, modification) => {
                                peptide.add_simple_modification(
                                    SequencePosition::from_index(location, peptide.len()).ok_or_else(||
                                        BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid modification position",
                                            format!("Index {location} falls outside the peptide of length {}", peptide.len()),
                                            Context::default().line_index(line.line_index).lines(0, line.line).add_highlight((0, range.clone()))))?,
                                    modification,
                                );
                            }
                            MzTabReturnModification::GlobalAmbiguous(modification) => {
                                let _possible = peptide.add_unknown_position_modification(
                                    modification,
                                    0..peptide.len(),
                                    &MUPSettings::default(),
                                );
                            }
                            MzTabReturnModification::Ambiguous(pos, modification) => {
                                let locations = pos
                                    .into_iter()
                                    .map(|(index, score)| {
                                        Ok((SequencePosition::from_index(index, peptide.len()).ok_or_else(||
                                        BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid modification position",
                                            format!("Index {index} falls outside the peptide of length {}", peptide.len()),
                                            Context::default().line_index(line.line_index).lines(0, line.line).add_highlight((0, range.clone()))))?, score))
                                    })
                                    .collect::<Result<Vec<_>, BoxedError<'a, BasicKind>>>()?;
                                let _possible = peptide.add_ambiguous_modification(
                                    modification,
                                    None,
                                    &locations,
                                    None,
                                    None,
                                    true,
                                );
                            }
                            MzTabReturnModification::NeutralLoss(_, _) => {
                                // TODO: handle neutral losses
                            }
                        }
                    }
                    Some(
                        PeptideModificationSearch::in_modifications(global_modifications.to_vec())
                            .tolerance(Tolerance::new_ppm(20.0))
                            .search(peptide),
                    )
                }
            },
            id: line
                .required_column("psm_id")
                .map_err(BoxedError::to_owned)?
                .0
                .parse()
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab PSM_ID",
                        format!("The PSM_ID {}", explain_number_error(&err)),
                        Context::line_range(
                            Some(line.line_index),
                            line.line,
                            line.optional_column("psm_id").unwrap().1,
                        ),
                    )
                })?,
            protein: line.optional_column("accession").and_then(|(v, _)| {
                (!v.eq_ignore_ascii_case("null"))
                    .then(|| (v.to_string(), proteins.get(v).map(Clone::clone)))
            }), // TODO: when warnings are hooked up add a warning about the missing reference protein
            unique: line
                .optional_column("unique")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v == "1")),
            database: line
                .optional_column("database")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v.to_string()))
                .map(|db| {
                    (
                        db,
                        line.optional_column("database_version").and_then(|(v, _)| {
                            (!v.eq_ignore_ascii_case("null")).then(|| v.to_string())
                        }),
                    )
                }),
            search_engine: {
                let (value, range) = line
                    .required_column("search_engine")
                    .map_err(BoxedError::to_owned)?;

                if value.trim().eq_ignore_ascii_case("null") {
                    Vec::new()
                } else {
                    value
                        .split('|')
                        .enumerate()
                        .map(|(i, s)| {
                            line.optional_column(&format!("search_engine_score[{}]", i + 1))
                                .and_then(|(v, _)| {
                                    (!v.eq_ignore_ascii_case("null")).then(|| {
                                        v.parse::<f64>().map_err(|err| {
                                            BoxedError::new(BasicKind::Error,
                                                "Invalid mzTab search engine score",
                                                format!(
                                        "The search engine score can not be parsed as f64: {err}"
                                    ),
                                                Context::line_range(
                                                    Some(line.line_index),
                                                    line.line,
                                                    line.optional_column(&format!(
                                                        "search_engine_score[{}]",
                                                        i + 1
                                                    ))
                                                    .unwrap()
                                                    .1,
                                                ),
                                            )
                                        })
                                    })
                                })
                                .transpose()
                                .and_then(|score| {
                                    CVTerm::from_str(s)
                                        .map_err(|e| {
                                            e.replace_context(Context::line_range(
                                                Some(line.line_index),
                                                line.line,
                                                range.clone(),
                                            ))
                                        })
                                        .and_then(|engine| Ok((engine, score, search_engine_score_types.get(i).ok_or_else(|| BoxedError::new(BasicKind::Error,"Missing search engine score type", "All search engines require a defined search type", Context::line_range(
                                            Some(line.line_index),
                                            line.line,
                                            range.clone(),
                                        )))?.clone())))
                                })
                        })
                        .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()?
                }
            },
            reliability: line
                .optional_column("reliability")
                .map(|(v, range)| match v {
                    "1" => Ok(Some(Reliability::High)),
                    "2" => Ok(Some(Reliability::Medium)),
                    "3" => Ok(Some(Reliability::Poor)),
                    s if s.eq_ignore_ascii_case("null") => Ok(None),
                    _ => Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid PSM reliability",
                        format!("A reliability should be 1, 2, 3, or null, '{v}' is invalid"),
                        Context::line_range(Some(line.line_index), line.line, range),
                    )),
                })
                .transpose()?
                .flatten(),
            rt: line
                .optional_column("retention_time")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid mzTab retention time",
                                    format!("The retention time can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| Time::new::<mzcore::system::s>(v))
                    })
                })
                .transpose()?,
            z: {
                let (value, range) = line
                    .required_column("charge")
                    .map_err(BoxedError::to_owned)?;

                if value.trim().eq_ignore_ascii_case("null") {
                    Charge::new::<mzcore::system::e>(1)
                } else {
                    value
                        .trim_end_matches(".0")
                        .parse::<isize>()
                        .map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid mzTab charge",
                                format!("The charge {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, range),
                            )
                        })
                        .map(|v| Charge::new::<mzcore::system::e>(v))?
                }
            },
            mz: line
                .optional_column("exp_mass_to_charge")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid mzTab experimental mz",
                                    format!("The experimental mz can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| MassOverCharge::new::<mzcore::system::thomson>(v))
                    })
                })
                .transpose()?,
            uri: line.optional_column("uri").map(|(v, _)| v.to_string()),
            spectra_ref: {
                let (value, range) = line
                    .required_column("spectra_ref")
                    .map_err(BoxedError::to_owned)?;
                if value.eq_ignore_ascii_case("null") {
                    SpectrumIds::None
                } else {
                    let grouped = value
                .split('|')
                .map(|value|
                value.split_once(':')
                .ok_or_else(|| {
                    BoxedError::new(BasicKind::Error,
                        "Invalid mzTab spectra_ref",
                        "The spectra_ref should be 'ms_run[x]:id'",
                        Context::line_range(Some(line.line_index), line.line, range.clone()),
                    )
                })
                .and_then(|(run, scan_id)| {
                    let index = run
                        .trim_start_matches("ms_run[")
                        .trim_end_matches(']')
                        .parse::<usize>()
                        .map_err(|err| {
                            BoxedError::new(BasicKind::Error,
                                "Invalid mzTab ms_run",
                                format!("The ms_run identifier {}", explain_number_error(&err)),
                                Context::line_range(
                                    Some(line.line_index),
                                    line.line,
                                    range.clone(),
                                ),
                            )
                        })
                        .and_then(|v| if v == 0 {
                            Err(BoxedError::new(BasicKind::Error,
                                "Invalid mzTab ms_run",
                                "A ms_run identifier uses 1 based indexing and so cannot be 0.",
                                Context::line_range(
                                    Some(line.line_index),
                                    line.line,
                                    range.clone(),
                                ),
                            ))
                        } else {Ok(v-1)})?;
                    let path = raw_files.get(index).ok_or_else(|| BoxedError::new(BasicKind::Error,"Missing raw file definition", "All raw files should be defined in the MTD section before being used in the PSM Section", Context::line_range(
                        Some(line.line_index),
                        line.line,
                        range.clone(),
                    )))?.0.as_ref().ok_or_else(|| BoxedError::new(BasicKind::Error,"Missing raw file path definition", "The path is not defined for this raw file", Context::line_range(
                        Some(line.line_index),
                        line.line,
                        range.clone(),
                    )))?;

                    let id = match scan_id.split_once('=') {
                        Some(("scan", num)) if num.chars().all(|c| c.is_ascii_digit()) => SpectrumId::Number(num.parse().map_err(|err| {
                            BoxedError::new(BasicKind::Error,
                                "Invalid mzTab spectra_ref scan number",
                                format!("The spectra_ref scan number {}", explain_number_error(&err)),
                                Context::line_range(
                                    Some(line.line_index),
                                    line.line,
                                    range.clone(),
                                ),
                            )
                        })?),
                        Some(("index", index)) if index.chars().all(|c| c.is_ascii_digit()) => SpectrumId::Index(index.parse().map_err(|err| {
                            BoxedError::new(BasicKind::Error,
                                "Invalid mzTab spectra_ref index",
                                format!("The spectra_ref index {}", explain_number_error(&err)),
                                Context::line_range(
                                    Some(line.line_index),
                                    line.line,
                                    range.clone(),
                                ),
                            )
                        })?),
                        _ =>  SpectrumId::Native(scan_id.to_owned()),
                    };

                    Ok((std::path::PathBuf::from(path), id))
                })).collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()?
                .into_iter()
                .sorted_by(|(a, _), (b,_)| a.cmp(b))
                .into_group_map_by(|(path, _)| path.clone());

                    if grouped.is_empty() {
                        SpectrumIds::None
                    } else if grouped.len() == 1
                        && let Some((path, ids)) = grouped.iter().next()
                        && *path == std::path::PathBuf::default()
                    {
                        SpectrumIds::FileNotKnown(ids.iter().map(|(_, i)| i.clone()).collect())
                    } else {
                        SpectrumIds::FileKnown(
                            grouped
                                .into_iter()
                                .map(|(path, ids)| {
                                    (path, ids.into_iter().map(|(_, i)| i).collect())
                                })
                                .collect(),
                        )
                    }
                }
            },
            flanking_sequence: {
                let pre = line.required_column("pre").map_err(BoxedError::to_owned)?;
                let post = line.required_column("post").map_err(BoxedError::to_owned)?;

                let pre = match pre.0 {
                    "null" => FlankingSequence::Unknown,
                    "-" => FlankingSequence::Terminal,
                    aa if aa.len() == 1 => {
                        FlankingSequence::AminoAcid(AminoAcid::from_str(aa).map_err(|()| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid preceding amino acid",
                                "This is not a valid amino acid code",
                                Context::line_range(Some(line.line_index), line.line, pre.1),
                            )
                        })?)
                    }
                    _ => FlankingSequence::Sequence(Box::new(Peptidoform::sloppy_pro_forma(
                        line.line,
                        pre.1,
                        ontologies,
                        &SloppyParsingParameters::default(),
                    )?)),
                };

                let post = match post.0 {
                    "null" => FlankingSequence::Unknown,
                    "-" => FlankingSequence::Terminal,
                    aa if aa.len() == 1 => {
                        FlankingSequence::AminoAcid(AminoAcid::from_str(aa).map_err(|()| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid following amino acid",
                                "This is not a valid amino acid code",
                                Context::line_range(Some(line.line_index), line.line, post.1),
                            )
                        })?)
                    }
                    _ => FlankingSequence::Sequence(Box::new(Peptidoform::sloppy_pro_forma(
                        line.line,
                        post.1,
                        ontologies,
                        &SloppyParsingParameters::default(),
                    )?)),
                };
                (pre, post)
            },
            protein_location: line
                .optional_column("start")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<u16>().map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid mzTab start",
                                format!("The start {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, r),
                            )
                        })
                    })
                })
                .transpose()
                .and_then(|start| {
                    let end = line
                        .optional_column("end")
                        .and_then(|(v, r)| {
                            (!v.eq_ignore_ascii_case("null")).then(|| {
                                v.parse::<u16>().map_err(|err| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid mzTab end",
                                        format!("The end {}", explain_number_error(&err)),
                                        Context::line_range(Some(line.line_index), line.line, r),
                                    )
                                })
                            })
                        })
                        .transpose()?;
                    Ok(start.and_then(|s| end.map(|e| s..e)))
                })?,
            local_confidence: line
                .optional_column("opt_ms_run[1]_aa_scores")
                .filter(|(lc, _)| !lc.trim().is_empty())
                .map(|(v, r)| {
                    v.split(',')
                        .map(|score| {
                            score.parse::<f64>().map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid mzTab local confidence",
                                    format!("The local confidence can not be parsed: {err}"),
                                    Context::line_range(
                                        Some(line.line_index),
                                        line.line,
                                        r.clone(),
                                    ),
                                )
                            })
                        })
                        .collect()
                })
                .transpose()?,
            additional: line
                .header
                .iter()
                .enumerate()
                .filter(|(_, column)| {
                    column.starts_with("opt") && *column != "opt_ms_run[1]_aa_scores"
                })
                .map(|(index, column)| {
                    (
                        column.clone(),
                        line.line[line.fields[index].clone()].to_string(),
                    )
                })
                .collect(),
        };

        result.local_confidence = result.local_confidence.as_ref().map(|lc| {
            let pep_len = result.peptidoform.as_ref().map_or(0, Peptidoform::len);
            if lc.len() > pep_len {
                // Casanovo stores the confidence for N and C terminal modifications.
                // As Casanovo has a double N terminal modification (+43.006-17.027) which could also
                // exist as two separate modifications the number of N terminal modifications is not a
                // reliable measure to determine how many local confidence scores to ignore.
                let c = result
                    .peptidoform
                    .as_ref()
                    .map_or(0, |p| p.get_c_term().len());
                let n = lc.len() - c - pep_len;
                let range = n..lc.len() - c;
                if range.len() == pep_len {
                    Ok(lc[n..lc.len() - c].to_vec())
                } else {
                    Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid local confidence", 
                        format!("The number of elements in the local confidence ({}) does not match the number of amino acids ({pep_len}).", lc.len()), 
                        Context::default().line_index(line.line_index).lines(0, line.line)))
                }
            } else if lc.len() == pep_len {
                Ok(lc.to_vec())
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid local confidence", 
                    format!("The number of elements in the local confidence ({}) does not match the number of amino acids ({pep_len}).", lc.len()), 
                    Context::default().line_index(line.line_index).lines(0, line.line)))
            }
        }).transpose()?;

        Ok(result)
    }
}

/// A protein definition from mzTab
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct MzTabProtein {
    /// The accession number, like 'Q340U4'
    pub accession: String,
    /// The protein’s name and or description line.
    pub description: String,
    /// The NCBI/NEWT taxonomy id for the species the protein was identified in.
    pub taxid: u32,
    /// The human readable species the protein was identified in - this SHOULD be the NCBI entry’s name.
    pub species: String,
    /// The protein database used for the search (could theoretically come from a different species).
    pub database: String,
    /// The version of the database used, if there is no version the data.
    pub database_version: String,
    /// The search engines that identified this protein along with their scores
    pub search_engine: Vec<(CVTerm, Option<(f64, CVTerm)>)>,
    /// A list of all proteins that cannot be separated based on peptide evidence from this main protein.
    pub ambiguity_members: Vec<String>,
    /// The reported modifications on this protein.
    pub modifications: Vec<(Vec<(SequencePosition, Option<f64>)>, SimpleModification)>,
    /// The coverage of this protein based on the peptidoforms in this file.
    pub coverage: Option<f64>,
    /// The GO terms for this protein.
    pub go_terms: Vec<Curie>,
    /// The reliability of this protein
    pub reliability: Option<Reliability>,
    /// A URI pointing to the protein's source entry in the unit it was identified in (e.g., the PRIDE database or a local database / file identifier).
    pub uri: Option<String>,
}

impl From<MzTabProtein> for crate::ProteinData {
    fn from(value: MzTabProtein) -> Self {
        Self::MzTab(value)
    }
}

impl mzcore::space::Space for MzTabProtein {
    fn space(&self) -> mzcore::space::UsedSpace {
        (self.accession.space()
            + self.description.space()
            + self.taxid.space()
            + self.species.space()
            + self.database.space()
            + self.database_version.space()
            + self.search_engine.space()
            + self.ambiguity_members.space()
            + self.modifications.space()
            + self.coverage.space()
            + self.go_terms.space()
            + self.reliability.space()
            + self.uri.space())
        .set_total::<Self>()
    }
}

impl MzTabProtein {
    /// Parse a single PRT line
    /// # Errors
    /// When not in the correct format
    fn from_line(line: PSMLine<'_>) -> Result<Arc<Self>, BoxedError<'static, BasicKind>> {
        Ok(Arc::new(Self {
            accession: line
                .required_column("accession")
                .map(|(v, _)| v.to_string())
                .map_err(BoxedError::to_owned)?,
            description: line
                .required_column("description")
                .map(|(v, _)| v.to_string())
                .map_err(BoxedError::to_owned)?,
            taxid: line
                .required_column("taxid")
                .and_then(|(v, location)| {
                    v.parse::<u32>().map_err(|e| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid PRT Line",
                            format!("The taxid number {}", explain_number_error(&e)),
                            Context::line_range(Some(line.line_index), line.line, location),
                        )
                    })
                })
                .map_err(BoxedError::to_owned)?,
            species: line
                .required_column("species")
                .map(|(v, _)| v.to_string())
                .map_err(BoxedError::to_owned)?,
            database: line
                .required_column("database")
                .map(|(v, _)| v.to_string())
                .map_err(BoxedError::to_owned)?,
            database_version: line
                .required_column("database_version")
                .map(|(v, _)| v.to_string())
                .map_err(BoxedError::to_owned)?,
            search_engine: Vec::new(), // TODO: actually parse
            reliability: line
                .optional_column("reliability")
                .map(|(v, range)| match v {
                    "1" => Ok(Some(Reliability::High)),
                    "2" => Ok(Some(Reliability::Medium)),
                    "3" => Ok(Some(Reliability::Poor)),
                    s if s.eq_ignore_ascii_case("null") => Ok(None),
                    _ => Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid PRT reliability",
                        format!("A reliability should be 1, 2, 3, or null, '{v}' is invalid"),
                        Context::line_range(Some(line.line_index), line.line, range).to_owned(),
                    )),
                })
                .transpose()?
                .flatten(),
            ambiguity_members: line
                .required_column("ambiguity_members")
                .map(|(v, _)| v.split(',').map(|s| s.trim().to_string()).collect())
                .map_err(BoxedError::to_owned)?,
            modifications: Vec::new(), // TODO: actually parse
            uri: line
                .optional_column("uri")
                .filter(|(v, _)| !v.eq_ignore_ascii_case("null"))
                .map(|(v, _)| v.to_string()),
            go_terms: line
                .optional_column("go_terms")
                .map(|(v, l)| {
                    let mut offset = 0;
                    v.split('|')
                        .map(|s| {
                            let r = s.parse::<Curie>().map_err(|e| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid CV accession",
                                    e.description(),
                                    Context::none()
                                        .line_index(line.line_index)
                                        .lines(0, line.line)
                                        .add_highlight((0, l.start + offset, s.len()))
                                        .to_owned(),
                                )
                            });
                            offset += s.len() + 1;
                            r
                        })
                        .collect::<Result<Vec<_>, _>>()
                })
                .transpose()?
                .unwrap_or_default(),
            coverage: line
                .optional_column("protein_coverage")
                .map(|(v, location)| {
                    v.parse::<f64>().map_err(|e| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid PRT Line",
                            format!("The protein coverage {e}"),
                            Context::line_range(Some(line.line_index), line.line, location)
                                .to_owned(),
                        )
                    })
                })
                .transpose()?,
        }))
    }
}

#[allow(dead_code)] // Neutral losses not yet handled
enum MzTabReturnModification {
    GlobalAmbiguous(SimpleModification),
    Ambiguous(Vec<(usize, Option<OrderedFloat<f64>>)>, SimpleModification),
    Defined(usize, SimpleModification),
    NeutralLoss(Option<usize>, Mass),
}

/// Parse a single modification definition. These are quite complex in mzTab. Rough schema:
/// ```ebnfish
/// Modification: (<position>|'null')'-'<definition>
/// Position: (<location>'|')*<location>
/// Location: <int>('[' <CV term> ']')?
/// Definition: 'unimod:'<int>|'mod:'<int>|'chemmod:'(<formula>|<float>)
/// Neutral loss: (<int>'-')? <CV term>
/// ```
///
/// # Errors
/// If this is not a valid modification definition
fn parse_modification<'a>(
    line: &'a str,
    range: Range<usize>,
    ontologies: &mzcore::ontology::Ontologies,
    line_index: u32,
) -> Result<MzTabReturnModification, BoxedError<'a, BasicKind>> {
    if let Some((pos, modification)) = line[range.clone()].split_once('-') {
        let position: Option<Vec<(usize, Option<OrderedFloat<f64>>)>> = if pos
            .eq_ignore_ascii_case("null")
        {
            Ok::<Option<Vec<(usize, Option<OrderedFloat<f64>>)>>, BoxedError<'a, BasicKind>>(None)
        } else {
            let mut index = range.start;
            Ok(Some(
                pos.split('|')
                    .map(|single| {
                        if let Some((offset, _, number)) =
                            next_number::<false, false, usize>(line, index..index + single.len())
                        {
                            let res = number.map_err(|err| {
                                BoxedError::new(BasicKind::Error,
                                    "Invalid modification position",
                                    format!("The position {}", explain_number_error(&err)),
                                    Context::line_range(
                                        Some(line_index),
                                        line,
                                        index..index + offset,
                                    ),
                                )
                            })?;

                            let parameter = &single[offset..];
                            let score = if !parameter.is_empty()
                                && single[offset..].starts_with('[')
                                && single.ends_with(']')
                            {
                                let value_range = CVTerm::parse_and_identify(line, index + offset..index+single.len(), "MS:1001876", "modification probability")?;
                                Some(line[value_range.clone()].parse::<f64>().map(OrderedFloat).map_err(|err| BoxedError::new(BasicKind::Error,
                                                "Invalid modification position",
                                                format!("The modification probability is not a valid number: {err}"),
                                                Context::line_range(
                                        Some(line_index),
                                        line,
                                        value_range
                                    )
                                )))
                            } else if parameter.is_empty() {
                                None
                            } else {
                                Some(Err(BoxedError::new(BasicKind::Error,
                                    "Invalid modification position",
                                    "A modification position parameter should be enclosed in square brackets '[]'",
                                    Context::line_range(
                                        Some(line_index),
                                        line,
                                        index + offset..index + single.len(),
                                    ),
                                )))
                            }.transpose()?;
                            index += single.len() + 1;

                            Ok((res, score))
                        } else {
                            Err(BoxedError::new(BasicKind::Error,
                                "Invalid modification position",
                                "A modification position should start with a number",
                                Context::line_range(
                                    Some(line_index),
                                    line,
                                    index..index + single.len(),
                                ),
                            ))
                        }
                    })
                    .collect::<Result<_, _>>()?,
            ))
        }?;

        if modification.starts_with('[') && modification.ends_with(']') {
            let value_range = CVTerm::parse_and_identify(
                line,
                range.start + pos.len() + 1..range.end,
                "MS:1001524",
                "fragment neutral loss",
            )?;
            let loss = line[value_range.clone()]
                .parse::<f64>()
                .map(|v| Mass::new::<mzcore::system::dalton>(v))
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid neutral loss",
                        format!("The fragment neutral loss is not a valid number: {err}"),
                        Context::line_range(Some(line_index), line, value_range),
                    )
                })?;

            position.map_or_else(
                || Ok(MzTabReturnModification::NeutralLoss(None, loss)),
                |pos| {
                    if pos.len() == 1 {
                        Ok(MzTabReturnModification::NeutralLoss(Some(pos[0].0), loss))
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid neutral loss",
                            "A neutral loss cannot be placed on multiple positions",
                            Context::line_range(
                                Some(line_index),
                                line,
                                range.start..range.start + pos.len(),
                            ),
                        ))
                    }
                },
            )
        } else {
            let modification = parse_single_modification(
                line,
                range.start + pos.len() + 1..range.end,
                ontologies,
                line_index,
            )?;

            Ok(match position {
                None => MzTabReturnModification::GlobalAmbiguous(modification),
                Some(pos) => {
                    if pos.len() == 1 {
                        MzTabReturnModification::Defined(pos[0].0, modification)
                    } else {
                        MzTabReturnModification::Ambiguous(pos, modification)
                    }
                }
            })
        }
    } else if line[range.clone()].starts_with('[') && line[range.clone()].ends_with(']') {
        let value_range =
            CVTerm::parse_and_identify(line, range, "MS:1001524", "fragment neutral loss")?;
        line[value_range.clone()]
            .parse::<f64>()
            .map(|v| Mass::new::<mzcore::system::dalton>(v))
            .map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid neutral loss",
                    format!("The fragment neutral loss is not a valid number: {err}"),
                    Context::line_range(Some(line_index), line, value_range),
                )
            })
            .map(|loss| MzTabReturnModification::NeutralLoss(None, loss))
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid modification",
            "A modification should be the position followed by a hyphen ('-') followed by the modification",
            Context::line_range(Some(line_index), line, range),
        ))
    }
}

/// Parse a single mzTab modification.
/// # Errors
/// if it does not follow the UNIMOD/MOD/CHEMMOD rules.
fn parse_single_modification<'a>(
    line: &'a str,
    range: Range<usize>,
    ontologies: &mzcore::ontology::Ontologies,
    line_index: u32,
) -> Result<SimpleModification, BoxedError<'a, BasicKind>> {
    if let Some((tag, value)) = line[range.clone()].split_once(':') {
        let value_context = Context::line_range(Some(line_index), line, range.clone());
        let modification = if tag.eq_ignore_ascii_case("unimod") {
            ontologies
                .unimod()
                .get_by_index(&value.parse::<u32>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid unimod code",
                        format!("The unimod modification {}", explain_number_error(&err)),
                        value_context.clone(),
                    )
                })?)
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid unimod code",
                        "The given unimod modification does not exist",
                        value_context.clone(),
                    )
                })?
        } else if tag.eq_ignore_ascii_case("mod") {
            ontologies
                .psimod()
                .get_by_index(&value.parse::<u32>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid PSI-MOD code",
                        format!("The PSI-MOD modification {}", explain_number_error(&err)),
                        value_context.clone(),
                    )
                })?)
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid PSI-MOD code",
                        "The given PSI-MOD modification does not exist",
                        value_context.clone(),
                    )
                })?
        } else if tag.eq_ignore_ascii_case("custom") {
            ontologies
                .custom()
                .get_by_index(&value.parse::<u32>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid custom code",
                        format!("The custom modification {}", explain_number_error(&err)),
                        value_context.clone(),
                    )
                })?)
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid custom code",
                        "The given custom modification does not exist",
                        value_context.clone(),
                    )
                })?
        } else if tag.eq_ignore_ascii_case("chemmod") {
            if let Ok(mass) = value.parse::<f64>() {
                SimpleModificationInner::Mass(
                    mzcore::sequence::MassTag::None,
                    Mass::new::<mzcore::system::dalton>(mass).into(),
                    float_digits(value),
                )
            } else {
                let factor = match line.as_bytes()[range.start + tag.len() + 1] {
                    b'-' => -1,
                    b'+' => 1,
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzTab modification",
                            "A chemmod formula modification should be prepended by a sign",
                            Context::line(Some(line_index), line, range.start + tag.len() + 1, 1),
                        ));
                    }
                };
                MolecularFormula::pro_forma_inner::<false, false>(
                    &Context::none().lines(0, line),
                    line,
                    range.start + tag.len() + 2..range.end,
                )
                .map(|f| SimpleModificationInner::Formula(f * factor))?
            }
            .into()
        } else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid mzTab modification",
                "The modification should be prepended by a tag describing the kind of modification, the possible tags are: 'unimod', 'mod', and 'chemmod'",
                Context::line_range(Some(line_index), line, range),
            ));
        };

        Ok(modification)
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid mzTab modification",
            "An mzTab modification should be in format 'tag:value' but the colon (':') is missing",
            Context::line_range(Some(line_index), line, range),
        ))
    }
}

#[derive(Clone, Copy, Debug)]
struct PSMLine<'a> {
    line_index: u32,
    header: &'a [String],
    pub line: &'a str,
    fields: &'a [Range<usize>],
}

impl<'a> PSMLine<'a> {
    /// Form an indexable line out of a set of fields
    /// # Errors
    /// When there is no header or the line has a different number of columns
    fn new(
        line_index: u32,
        header: Option<&'a [String]>,
        line: &'a str,
        fields: &'a [Range<usize>],
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let header = header.ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Missing PSH line",
                "The PSH peptide header line should precede any PSM line",
                Context::full_line(line_index, line),
            )
        })?;
        if header.len() == fields.len() {
            Ok(Self {
                line_index,
                header,
                line,
                fields,
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid PSM line",
                format!(
                    "This PSM line does not have the same number of columns as the PSH line (PSH: {}, PSM: {})",
                    header.len(),
                    fields.len(),
                ),
                Context::full_line(line_index, line),
            ))
        }
    }

    fn optional_column(&self, column: &str) -> Option<(&str, Range<usize>)> {
        self.header
            .iter()
            .position(|h| h == column)
            .map(|i| (&self.line[self.fields[i].clone()], self.fields[i].clone()))
    }

    /// Get a required columns
    /// # Errors
    /// If the column is not available
    fn required_column<'b>(
        &'b self,
        column: &str,
    ) -> Result<(&'b str, Range<usize>), BoxedError<'b, BasicKind>> {
        self.optional_column(column).ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Missing column",
                format!("The column '{column}' is required but not present"),
                Context::full_line(self.line_index, self.line),
            )
        })
    }
}

impl std::fmt::Display for PSMLine<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            Context::default()
                .line_index(self.line_index)
                .lines(0, self.line)
                .add_highlights(self.fields.iter().map(|r| (0, r.clone())))
        )
    }
}

impl From<MzTabPSM> for PSM<SimpleLinear, MaybePeptidoform> {
    fn from(value: MzTabPSM) -> Self {
        Self {
            score: (!value.search_engine.is_empty())
                .then(|| {
                    (value
                        .search_engine
                        .iter()
                        .filter_map(|(_, s, _)| *s)
                        .sum::<f64>()
                        / value.search_engine.len() as f64)
                        .clamp(-1.0, 1.0)
                })
                .filter(|v| !v.is_nan()),
            local_confidence: value.local_confidence.clone(),
            data: PSMData::MzTab(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

/// A CV term
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub struct CVTerm {
    /// The term
    pub term: Term,
    /// Additional comments on the term, eg additional specification or version
    pub comment: String,
}

impl std::fmt::Display for CVTerm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{}, {}, {}, {}]",
            self.term.accession.cv, self.term.accession, self.term.name, self.comment
        )
    }
}

impl mzcore::space::Space for CVTerm {
    fn space(&self) -> mzcore::space::UsedSpace {
        (self.term.space() + self.comment.space()).set_total::<Self>()
    }
}

impl CVTerm {
    /// Get a line parse the range as a cv term, check the id, and give back the comment/value (if any).
    /// All done without any allocations.
    /// # Errors
    /// If the selected section is not enclosed in square brackets. Or if the Line has a different term id.
    fn parse_and_identify<'a>(
        line: &'a str,
        range: Range<usize>,
        required_id: &str,
        required_term: &str,
    ) -> Result<Range<usize>, BoxedError<'a, BasicKind>> {
        let value = &line[range.clone()];
        if value.starts_with('[') && value.ends_with(']') {
            let value = &value[1..value.len() - 1];
            let mut field = 0;
            let mut field_index: Option<Range<usize>> = None;
            let mut id = range.start..range.start;
            for (index, c) in value.char_indices() {
                let index = index + range.start + 1;
                if c == ',' {
                    if field == 1 {
                        id = field_index.as_ref().map_or(index..index, Clone::clone);
                    }
                    if field == 3 {
                        // Just continue
                    } else {
                        field += 1;
                        field_index = None;
                    }
                } else if !c.is_ascii_whitespace() {
                    field_index = field_index.map_or_else(
                        || Some(index..index + c.len_utf8()),
                        |range| Some(range.start..index + c.len_utf8()),
                    );
                }
            }
            let comment = field_index.map_or(range.end..range.end, |range| range);

            if line[id.clone()].eq_ignore_ascii_case(required_id) {
                Ok(comment)
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid CV term",
                    format!(
                        "A CV term with id {required_id} \"{required_term}\" was expected but id {} was found",
                        &line[id]
                    ),
                    Context::line_range(None, line, range),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid CV term",
                "A CV term should be enclosed by '[]'",
                Context::line_range(None, line, range),
            ))
        }
    }
}

impl FromStr for CVTerm {
    type Err = BoxedError<'static, BasicKind>;
    fn from_str(value: &str) -> Result<Self, BoxedError<'static, BasicKind>> {
        let value = value.trim();
        if value.starts_with('[') && value.ends_with(']') {
            let value = &value[1..value.len() - 1];
            let mut split = value.splitn(4, ',');
            let _ontology = split.next().unwrap_or_default().trim().to_string();
            let accession = split
                .next()
                .unwrap_or_default()
                .trim()
                .to_string()
                .parse::<Curie>()
                .map_err(|e| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid CV accession",
                        e.description(),
                        Context::none().lines(0, value).to_owned(),
                    )
                })?;
            let name = split.next().unwrap_or_default().trim().to_string();
            Ok(Self {
                term: Term {
                    accession,
                    name: name.into(),
                },
                comment: split.next().unwrap_or_default().trim().to_string(),
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid CV term",
                "A CV term should be enclosed by '[]'",
                Context::show(value).to_owned(),
            ))
        }
    }
}

/// The reliability of a PSM
#[expect(missing_docs)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub enum Reliability {
    High,
    Medium,
    #[default]
    Poor,
}

impl std::fmt::Display for Reliability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::High => write!(f, "High"),
            Self::Medium => write!(f, "Medium"),
            Self::Poor => write!(f, "Poor"),
        }
    }
}

impl mzcore::space::Space for Reliability {
    fn space(&self) -> mzcore::space::UsedSpace {
        mzcore::space::UsedSpace::stack(size_of::<Self>())
    }
}

/// A basic structure for a mzTab file line
#[expect(clippy::upper_case_acronyms)]
enum MzTabLine {
    /// Metadata line
    MTD(u32, String, Vec<Range<usize>>),
    /// Protein header line
    PRH(u32, String, Vec<Range<usize>>),
    /// Protein line
    PRT(u32, String, Vec<Range<usize>>),
    /// Peptide header line
    PSH(u32, String, Vec<Range<usize>>),
    /// Peptide line, stored as hash map with the columns names from PSH
    PSM(u32, String, Vec<Range<usize>>),
}

/// Parse a mzTab file
/// # Errors
/// If the file is not a valid mzTab file
fn parse_mztab_reader<T: BufRead>(
    reader: T,
) -> impl Iterator<Item = Result<Option<MzTabLine>, BoxedError<'static, BasicKind>>> {
    reader.lines().enumerate().map(move |(line_index, line)| {
        let line_index = line_index as u32;
        line.map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Could not read line",
                err.to_string(),
                Context::default().line_index(line_index),
            )
        })
        .and_then(move |line| {
            if line.trim().is_empty() {
                Ok(None)
            } else {
                mzcore::csv::csv_separate(&line, b'\t')
                    .map_err(BoxedError::to_owned)
                    .map(|fields| match &line[fields[0].clone()] {
                        "MTD" => Some(MzTabLine::MTD(line_index, line, fields)),
                        "PRH" => Some(MzTabLine::PRH(line_index, line, fields)),
                        "PRT" => Some(MzTabLine::PRT(line_index, line, fields)),
                        "PSH" => Some(MzTabLine::PSH(line_index, line, fields)),
                        "PSM" => Some(MzTabLine::PSM(line_index, line, fields)),
                        _ => None,
                    })
            }
        })
    })
}

impl PSMMetaData for MzTabPSM {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        self.peptidoform
            .as_ref()
            .map(|p| Cow::Owned(p.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::MzTab
    }

    fn numerical_id(&self) -> Option<usize> {
        Some(self.id)
    }

    fn id(&self) -> String {
        self.id.to_string()
    }

    fn search_engine(&self) -> Option<Term> {
        self.search_engine.first().map(|(t, _, _)| t.term.clone())
    }

    fn confidence(&self) -> Option<f64> {
        (!self.search_engine.is_empty())
            .then(|| {
                (self
                    .search_engine
                    .iter()
                    .filter_map(|(_, s, _)| *s)
                    .sum::<f64>()
                    / self.search_engine.len() as f64)
                    .clamp(-1.0, 1.0)
            })
            .filter(|v| !v.is_nan())
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        self.local_confidence
            .as_ref()
            .map(|lc| Cow::Borrowed(lc.as_slice()))
    }

    fn original_confidence(&self) -> Option<(f64, Term)> {
        self.search_engine
            .first()
            .and_then(|(_, s, t)| s.map(|s| (s, t.term.clone())))
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        self.local_confidence.as_deref()
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.z)
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        self.rt
    }

    fn scans(&self) -> SpectrumIds {
        self.spectra_ref.clone()
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        self.mz
    }

    fn experimental_mass(&self) -> Option<Mass> {
        self.mz.map(|mz| mz * self.z.to_float())
    }

    type Protein = MzTabProtein;
    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(
            self.protein
                .as_ref()
                .and_then(|p| p.1.as_ref().map(|a| std::slice::from_ref(&**a)))
                .unwrap_or_default(),
        )
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        self.protein_location.clone()
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&self.flanking_sequence.0, &self.flanking_sequence.1)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        self.database
            .as_ref()
            .map(|(db, version)| (db.as_str(), version.as_deref()))
    }

    fn unique(&self) -> Option<bool> {
        self.unique
    }

    fn reliability(&self) -> Option<Reliability> {
        self.reliability
    }

    fn uri(&self) -> Option<String> {
        self.uri.clone()
    }
}

impl ProteinMetaData for MzTabProtein {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<mzcore::sequence::Linear>>> {
        None
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        FastaIdentifier::Undefined(false, self.accession.as_str())
    }

    fn description(&self) -> Option<&str> {
        Some(&self.description)
    }

    fn species(&self) -> Option<Curie> {
        Some(Curie {
            cv: ControlledVocabulary::NCBITaxon,
            accession: mzcv::AccessionCode::Numeric(self.taxid),
        })
    }

    fn species_name(&self) -> Option<&str> {
        Some(&self.species)
    }

    fn search_engine(&self) -> &[(CVTerm, Option<(f64, CVTerm)>)] {
        &self.search_engine
    }

    fn ambiguity_members(&self) -> &[String] {
        &self.ambiguity_members
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        Some((&self.database, Some(&self.database_version)))
    }

    fn modifications(&self) -> &[(Vec<(SequencePosition, Option<f64>)>, SimpleModification)] {
        &self.modifications
    }

    fn coverage(&self) -> Option<f64> {
        self.coverage
    }

    fn gene_ontology(&self) -> &[Curie] {
        &self.go_terms
    }

    fn reliability(&self) -> Option<Reliability> {
        self.reliability
    }

    fn uri(&self) -> Option<&str> {
        self.uri.as_deref()
    }
}
