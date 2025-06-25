use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    ops::Range,
    str::FromStr,
};

use flate2::bufread::GzDecoder;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::{check_extension, explain_number_error},
    identification::{IdentifiedPeptidoform, MaybePeptidoform, IdentifiedPeptidoformData, SpectrumId, SpectrumIds},
    ontology::CustomDatabase,
    quantities::Tolerance,
    sequence::{
        AminoAcid, PeptideModificationSearch, Peptidoform, ReturnModification, SemiAmbiguous,
        SimpleModification, SimpleModificationInner, SloppyParsingParameters,
    },
    system::{MassOverCharge, Time, isize::Charge},
};

/// Peptide data from a mzTab file
#[derive(Clone, Debug, Default, Deserialize, PartialEq, Serialize)]
pub struct MZTabData {
    /// The peptide's sequence corresponding to the PSM
    pub peptide: Option<Peptidoform<SemiAmbiguous>>,
    /// A unique identifier for a PSM within the file. If a PSM can be matched to
    /// multiple proteins, the same PSM should be represented on multiple rows with
    /// different accessions and the same PSM_ID.
    pub id: usize,
    /// The protein's accession the corresponding peptide sequence (coming from the
    /// PSM) is associated with.
    pub accession: Option<String>,
    /// Indicates whether the peptide sequence (coming from the PSM) is unique for
    /// this protein in respect to the searched database.
    pub unique: Option<bool>,
    /// The protein database used for the search (could theoretically come from a
    /// different species) and the peptide sequence comes from.
    pub database: Option<String>,
    /// The protein database's version
    pub database_version: Option<String>,
    /// The search engines that identified this peptide, alongside their identified score and the CV term describing the score
    pub search_engine: Vec<(CVTerm, Option<f64>, CVTerm)>,
    /// If available the estaimated reliability of the PSM.
    pub reliability: Option<PSMReliability>,
    /// The retention time for this peptide.
    pub rt: Option<Time>,
    /// The charge for this peptide.
    pub z: Charge,
    /// The experimental mz
    pub mz: Option<MassOverCharge>,
    /// A URI pointing to the PSM's entry in the experiment it was identified in (e.g. the peptide’s PRIDE entry).
    pub uri: Option<String>,
    /// The spectra references grouped by raw file
    pub spectra_ref: SpectrumIds,
    /// The amino acide before this peptide
    pub preceding_aa: FlankingResidue,
    /// The amino acide after this peptide
    pub following_aa: FlankingResidue,
    /// The start of this peptide in the containing protein (0-based)
    pub start: Option<usize>,
    /// The end of this peptide in the containing protein (0-based)
    pub end: Option<usize>,
    /// Casanovo specific additional metadata with the amino acid confidence
    pub local_confidence: Option<Vec<f64>>,
    /// Any additional metadata
    pub additional: HashMap<String, String>,
}

impl MZTabData {
    /// Parse a mzTab file.
    /// # Errors
    /// If the file is not in the correct format
    pub fn parse_file(
        path: impl AsRef<std::path::Path>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Box<dyn Iterator<Item = Result<Self, CustomError>> + '_>, CustomError> {
        let file = File::open(path.as_ref()).map_err(|e| {
            CustomError::error(
                "Could not open file",
                e,
                Context::Show {
                    line: path.as_ref().to_string_lossy().to_string(),
                },
            )
        })?;
        if check_extension(path, "gz") {
            Ok(Box::new(Self::parse_reader(
                BufReader::new(GzDecoder::new(BufReader::new(file))),
                custom_database,
            )))
        } else {
            Ok(Box::new(Self::parse_reader(
                BufReader::new(file),
                custom_database,
            )))
        }
    }

    /// Parse a mzTab file directly from a buffered reader
    pub fn parse_reader<'a, T: BufRead + 'a>(
        reader: T,
        custom_database: Option<&'a CustomDatabase>,
    ) -> impl Iterator<Item = Result<Self, CustomError>> + 'a {
        let mut search_engine_score_type: Vec<CVTerm> = Vec::new();
        let mut modifications: Vec<SimpleModification> = Vec::new();
        let mut raw_files: Vec<(Option<String>, Option<CVTerm>, Option<CVTerm>)> = Vec::new(); //path, file format, identifier type
        let mut peptide_header: Option<Vec<String>> = None;

        parse_mztab_reader(reader).filter_map(move |item| {
            item.transpose().and_then(|item| match item {
                Ok(MZTabLine::MTD(line_index, line, fields)) => {
                    if fields.len() == 3 {
                        match line[fields[1].clone()].to_ascii_lowercase().as_str() {
                            m if (m.starts_with("variable_mod[") || m.starts_with("fixed_mod[")) && m.ends_with(']') => {
                                match CVTerm::from_str(&line[fields[2].clone()]).and_then(|term|
                                        (term.id.trim() != "MS:1002453" && term.id.trim()  != "MS:1002454").then(||
                                            SimpleModificationInner::parse_pro_forma(term.id.trim(), 0..term.id.trim().len(), &mut Vec::new(), &mut Vec::new(), custom_database)).transpose()) {
                                    Ok(Some((ReturnModification::Defined(modification), _))) => if !modifications.contains(&modification) { modifications.push(modification)},
                                    Ok(Some(_)) => return Some(Err(CustomError::error("Invalid modification in mzTab", "Modifications in mzTab have to be defeined, not ambiguous or cross-linkers", Context::line_range(Some(line_index), line, fields[2].clone())))),
                                    Err(err) => return Some(Err(err)),
                                    Ok(None) => (),
                                }
                            },
                            m if m.starts_with("psm_search_engine_score[") && m.ends_with(']') => {
                                match CVTerm::from_str(&line[fields[2].clone()]) {
                                    Ok(term) => search_engine_score_type.push(term),
                                    Err(err) => return Some(Err(err)),
                                }
                            }
                            m if m.starts_with("ms_run[") && m.ends_with("]-location") => {
                                let index = match m.trim_start_matches("ms_run[").trim_end_matches("]-location").parse::<usize>().map_err(|err| {
                                    CustomError::error(
                                        "Invalid mzTab ms_run identifier",
                                        format!("The ms_run identifier {}", explain_number_error(&err)),
                                        Context::line_range(
                                            Some(line_index),
                                            &line,
                                            fields[1].clone(),
                                        ),
                                    )
                                }) {
                                    Ok(i) => i - 1,
                                    Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].0 = Some(line[fields[2].clone()].to_string());
                            },
                            m if m.starts_with("ms_run[") && m.ends_with("]-format") => {
                                let index = match m.trim_start_matches("ms_run[").trim_end_matches("]-format").parse::<usize>().map_err(|err| {
                                    CustomError::error(
                                        "Invalid mzTab ms_run identifier",
                                        format!("The ms_run identifier {}", explain_number_error(&err)),
                                        Context::line_range(
                                            Some(line_index),
                                            &line,
                                            fields[1].clone(),
                                        ),
                                    )
                                }) {
                                    Ok(i) => i - 1,
                                    Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].1 = Some(match CVTerm::from_str(&line[fields[2].clone()]) {
                                        Ok(i) => i,
                                        Err(err) => return Some(Err(err)),
                                });
                            },
                            m if m.starts_with("ms_run[") && m.ends_with("]-id_format") => {
                                let index = match m.trim_start_matches("ms_run[").trim_end_matches("]-id_format").parse::<usize>().map_err(|err| {
                                    CustomError::error(
                                        "Invalid mzTab ms_run identifier",
                                        format!("The ms_run identifier {}", explain_number_error(&err)),
                                        Context::line_range(
                                            Some(line_index),
                                            &line,
                                            fields[1].clone(),
                                        ),
                                    )
                                }) {
                                    Ok(i) => i - 1,
                                    Err(err) => return Some(Err(err)),
                                };

                                while raw_files.len() <= index {
                                    raw_files.push((None, None, None));
                                }

                                raw_files[index].2 = Some(match CVTerm::from_str(&line[fields[2].clone()]) {
                                        Ok(i) => i,
                                        Err(err) => return Some(Err(err)),
                                });
                            },
                            _ => (),
                        }
                        None
                    } else {
                        Some(Err(CustomError::error(
                            "Invalid MTD line",
                            "MTD lines should contain three columns (the tag, key, and value)",
                            Context::full_line(line_index, line),
                        )))
                    }
                }
                Ok(MZTabLine::PSH(line_index, line, fields)) => {
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
                            return Some(Err(CustomError::error(
                                "Invalid peptide table",
                                format!("The required column '{required}' is not present"),
                                Context::full_line(line_index, line),
                            )));
                        }
                    }
                    peptide_header = Some(header);
                    None
                }
                Ok(MZTabLine::PSM(line_index, line, fields)) => Some(
                    PSMLine::new(line_index, peptide_header.as_deref(), &line, &fields)
                        .and_then(|line| Self::from_line(line, &modifications, &search_engine_score_type, &raw_files, custom_database)),
                ),
                Err(e) => Some(Err(e)),
            })
        })
    }

    /// Parse a single PSM line
    /// # Errors
    /// When not in the correct format
    #[expect(clippy::missing_panics_doc)]
    fn from_line(
        line: PSMLine<'_>,
        global_modifications: &[SimpleModification],
        search_engine_score_types: &[CVTerm],
        raw_files: &[(Option<String>, Option<CVTerm>, Option<CVTerm>)],
        custom_database: Option<&CustomDatabase>,
    ) -> Result<Self, CustomError> {
        let (mod_column, mod_range) = line.required_column("modifications")?;
        let mut mod_index = mod_range.start;
        let modifications: Vec<(usize, SimpleModification)> = mod_column
            .split(',')
            .flat_map(|definition| {
                let pair = definition
                    .split_once('-')
                    .map(|(pos, _)| Ok((
                        pos.parse::<usize>().map_err(|err| CustomError::error(
                        "Invalid modification position",
                        format!("The position {}", explain_number_error(&err)),
                        Context::line_range(Some(line.line_index), line.line, mod_range.clone()),
                    ))?,
                    SimpleModificationInner::parse_pro_forma(
                        line.line,
                        mod_index+1+pos.len()..mod_index+definition.len(),
                        &mut Vec::new(),
                        &mut Vec::new(),
                        custom_database)?.0.defined()
                        .ok_or_else(
                            || CustomError::error(
                                "Invalid modification",
                                "A modification should be a fully defined modification, no cross-link or ambiguous modification",
                                Context::line_range(Some(line.line_index), line.line, mod_range.clone())))?)))
                    .ok_or_else(
                        || CustomError::error(
                            "Invalid modification",
                            "A modification should be the position followed by a hyphen ('-') followed by the modification",
                            Context::line_range(Some(line.line_index), line.line, mod_range.clone())));
                            mod_index += definition.len() + 1;
                            pair
            })
            .collect::<Result<Vec<_>, CustomError>>()?;

        let mut result = Self {
            peptide: {
                let range = line.required_column("sequence")?.1;

                if range.is_empty() {
                    None
                } else {
                    let mut peptide = Peptidoform::sloppy_pro_forma(
                        line.line,
                        range,
                        custom_database,
                        &SloppyParsingParameters {
                            allow_unwrapped_modifications: true,
                            ..Default::default()
                        },
                    )?;
                    for (location, modification) in modifications {
                        match location {
                            0 => peptide.add_simple_n_term(modification),
                            c if c == peptide.len() + 1 => {
                                peptide.add_simple_c_term(modification);
                            }
                            i => {
                                peptide.sequence_mut()[i - 1].add_simple_modification(modification);
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
            id: line.required_column("psm_id")?.0.parse().map_err(|err| {
                CustomError::error(
                    "Invalid mzTab PSM_ID",
                    format!("The PSM_ID {}", explain_number_error(&err)),
                    Context::line_range(
                        Some(line.line_index),
                        line.line,
                        line.optional_column("psm_id").unwrap().1,
                    ),
                )
            })?,
            accession: line
                .optional_column("accession")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v.to_string())),
            unique: line
                .optional_column("unique")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v == "1")),
            database: line
                .optional_column("database")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v.to_string())),
            database_version: line
                .optional_column("database_version")
                .and_then(|(v, _)| (!v.eq_ignore_ascii_case("null")).then(|| v.to_string())),
            search_engine: {
                let (value, range) = line.required_column("search_engine")?;

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
                                            CustomError::error(
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
                                            e.with_context(Context::line_range(
                                                Some(line.line_index),
                                                line.line,
                                                range.clone(),
                                            ))
                                        })
                                        .and_then(|engine| Ok((engine, score, search_engine_score_types.get(i).ok_or_else(|| CustomError::error("Missing search engine score type", "All search engines require a defined search type", Context::line_range(
                                            Some(line.line_index),
                                            line.line,
                                            range.clone(),
                                        )))?.clone())))
                                })
                        })
                        .collect::<Result<Vec<_>, CustomError>>()?
                }
            },
            reliability: line
                .optional_column("reliability")
                .map(|(v, range)| match v {
                    "1" => Ok(PSMReliability::High),
                    "2" => Ok(PSMReliability::Medium),
                    "3" => Ok(PSMReliability::Poor),
                    _ => Err(CustomError::error(
                        "Invalid PSM reliability",
                        format!("A reliability should be 1, 2, or 3, '{v}' is invalid"),
                        Context::line_range(Some(line.line_index), line.line, range),
                    )),
                })
                .transpose()?,
            rt: line
                .optional_column("retention_time")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid mzTab retention time",
                                    format!("The retention time can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| Time::new::<crate::system::s>(v))
                    })
                })
                .transpose()?,
            z: {
                let (value, range) = line.required_column("charge")?;

                if value.trim().eq_ignore_ascii_case("null") {
                    Charge::new::<crate::system::e>(1)
                } else {
                    value
                        .trim_end_matches(".0")
                        .parse::<isize>()
                        .map_err(|err| {
                            CustomError::error(
                                "Invalid mzTab charge",
                                format!("The charge {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, range),
                            )
                        })
                        .map(|v| Charge::new::<crate::system::e>(v))?
                }
            },
            mz: line
                .optional_column("exp_mass_to_charge")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<f64>()
                            .map_err(|err| {
                                CustomError::error(
                                    "Invalid mzTab experimental mz",
                                    format!("The experimental mz can not be parsed as f64: {err}"),
                                    Context::line_range(Some(line.line_index), line.line, r),
                                )
                            })
                            .map(|v| MassOverCharge::new::<crate::system::mz>(v))
                    })
                })
                .transpose()?,
            uri: line.optional_column("uri").map(|(v, _)| v.to_string()),
            spectra_ref: {
                let (value, range) = line.required_column("spectra_ref")?;
                let grouped = value
                .split('|')
                .map(|value|
                value.split_once(':')
                .ok_or_else(|| {
                    CustomError::error(
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
                            CustomError::error(
                                "Invalid mzTab ms_run",
                                format!("The ms_run identifier {}", explain_number_error(&err)),
                                Context::line_range(
                                    Some(line.line_index),
                                    line.line,
                                    range.clone(),
                                ),
                            )
                        })? - 1;
                    let path = raw_files.get(index).ok_or_else(|| CustomError::error("Missing raw file definition", "All raw files should be defined in the MTD section before being used in the PSM Section", Context::line_range(
                        Some(line.line_index),
                        line.line,
                        range.clone(),
                    )))?.0.as_ref().ok_or_else(|| CustomError::error("Missing raw file path definition", "The path is not defined for this raw file", Context::line_range(
                        Some(line.line_index),
                        line.line,
                        range.clone(),
                    )))?;

                    let id = match scan_id.split_once('=') {
                        Some(("scan", num)) if num.chars().all(|c| c.is_ascii_digit()) => SpectrumId::Number(num.parse().map_err(|err| {
                            CustomError::error(
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
                            CustomError::error(
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
                })).collect::<Result<Vec<_>, CustomError>>()?
                .into_iter()
                .sorted_by(|(a, _), (b,_)| a.cmp(b))
                .chunk_by(|(path, _)| path.clone());

                SpectrumIds::FileKnown(
                    grouped
                        .into_iter()
                        .map(|(path, ids)| (path, ids.into_iter().map(|(_, i)| i).collect()))
                        .collect(),
                )
            },
            preceding_aa: line.required_column("pre")?.0.parse().map_err(|()| {
                CustomError::error(
                    "Invalid preceding amino acid",
                    "The pre column should contain null, -, or an aminoacid",
                    Context::line_range(
                        Some(line.line_index),
                        line.line,
                        line.optional_column("pre").unwrap().1,
                    ),
                )
            })?,
            following_aa: line.required_column("post")?.0.parse().map_err(|()| {
                CustomError::error(
                    "Invalid following amino acid",
                    "The post column should contain null, -, or an aminoacid",
                    Context::line_range(
                        Some(line.line_index),
                        line.line,
                        line.optional_column("post").unwrap().1,
                    ),
                )
            })?,
            start: line
                .optional_column("start")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<usize>().map_err(|err| {
                            CustomError::error(
                                "Invalid mzTab start",
                                format!("The start {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, r),
                            )
                        })
                    })
                })
                .transpose()?,
            end: line
                .optional_column("end")
                .and_then(|(v, r)| {
                    (!v.eq_ignore_ascii_case("null")).then(|| {
                        v.parse::<usize>().map_err(|err| {
                            CustomError::error(
                                "Invalid mzTab end",
                                format!("The end {}", explain_number_error(&err)),
                                Context::line_range(Some(line.line_index), line.line, r),
                            )
                        })
                    })
                })
                .transpose()?,
            local_confidence: line
                .optional_column("opt_ms_run[1]_aa_scores")
                .filter(|(lc, _)| !lc.trim().is_empty())
                .map(|(v, r)| {
                    v.split(',')
                        .map(|score| {
                            score.parse::<f64>().map_err(|err| {
                                CustomError::error(
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
                        column.to_string(),
                        line.line[line.fields[index].clone()].to_string(),
                    )
                })
                .collect(),
        };

        result.local_confidence = result.local_confidence.as_ref().map(|lc| {
            // Casanovo stores the confidence for N and C terminal modifications
            // As Casanovo has a double N terminal modification (+43.006-17.027) which could also
            // exist as two separate modifications the number of N terminal modifications is not a
            // reliable measure to detrmine how many local confidence scores to ignore.
            let c = result.peptide.as_ref().map_or(0, |p| p.get_c_term().len());
            let n = lc.len() - c - result.peptide.as_ref().map_or(0, Peptidoform::len);
            lc[n..lc.len() - c].to_vec()
        });

        Ok(result)
    }
}

#[derive(Clone, Copy, Debug)]
struct PSMLine<'a> {
    line_index: usize,
    header: &'a [String],
    pub line: &'a str,
    fields: &'a [Range<usize>],
}

impl<'a> PSMLine<'a> {
    /// Form a indexable line out of a set of fields
    /// # Errors
    /// When there is no header or the line has a different number of columns
    fn new(
        line_index: usize,
        header: Option<&'a [String]>,
        line: &'a str,
        fields: &'a [Range<usize>],
    ) -> Result<Self, CustomError> {
        let header = header.ok_or_else(|| {
            CustomError::error(
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
            Err(CustomError::error(
                "Invalid PSM line",
                "This PSM line does not have the same number of columns as the PSH line",
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
    fn required_column(&self, column: &str) -> Result<(&str, Range<usize>), CustomError> {
        self.optional_column(column).ok_or_else(|| {
            CustomError::error(
                "Missing column",
                format!("The column '{column}' is required but not present"),
                Context::full_line(self.line_index, self.line),
            )
        })
    }
}

impl From<MZTabData> for IdentifiedPeptidoform<SemiAmbiguous, MaybePeptidoform> {
    fn from(value: MZTabData) -> Self {
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
            metadata: IdentifiedPeptidoformData::MZTab(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

/// A flanking residue for a sequence, N or C terminal agnostic
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub enum FlankingResidue {
    /// The flanking residue is unknown (for example in de novo data)
    #[default]
    Unknown,
    /// The residue is terminal
    Terminal,
    /// The flanking residue
    AminoAcid(AminoAcid),
}

impl FromStr for FlankingResidue {
    type Err = <AminoAcid as FromStr>::Err;
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim() {
            "null" => Ok(Self::Unknown),
            "-" => Ok(Self::Terminal),
            _ => AminoAcid::from_str(value).map(Self::AminoAcid),
        }
    }
}

impl std::fmt::Display for FlankingResidue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unknown => write!(f, "Unknown"),
            Self::Terminal => write!(f, "Terminal"),
            Self::AminoAcid(a) => write!(f, "{a}"),
        }
    }
}

/// A CV term
#[derive(Clone, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub struct CVTerm {
    /// The ontology
    pub ontology: String,
    /// The id within the ontology
    pub id: String,
    /// The human name for the term
    pub term: String,
    /// Additional comments on the term, eg additional specification
    pub comment: String,
}

impl FromStr for CVTerm {
    type Err = CustomError;
    fn from_str(value: &str) -> Result<Self, CustomError> {
        let value = value.trim();
        if value.starts_with('[') && value.ends_with(']') {
            let value = &value[1..value.len() - 1];
            let mut split = value.splitn(4, ',');
            Ok(Self {
                ontology: split.next().unwrap_or_default().trim().to_string(),
                id: split.next().unwrap_or_default().trim().to_string(),
                term: split.next().unwrap_or_default().trim().to_string(),
                comment: split.next().unwrap_or_default().trim().to_string(),
            })
        } else {
            Err(CustomError::error(
                "Invalid CV term",
                "A CV term should be encolsed by '[]'",
                Context::None,
            ))
        }
    }
}

/// The reliability of a PSM
#[expect(missing_docs)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, PartialEq, Serialize)]
pub enum PSMReliability {
    High,
    Medium,
    #[default]
    Poor,
}

impl std::fmt::Display for PSMReliability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::High => write!(f, "High"),
            Self::Medium => write!(f, "Medium"),
            Self::Poor => write!(f, "Poor"),
        }
    }
}

/// A basic structure for a mzTab file line
#[expect(clippy::upper_case_acronyms)]
enum MZTabLine {
    /// Metadata line
    MTD(usize, String, Vec<Range<usize>>),
    /// Peptide header line
    PSH(usize, String, Vec<Range<usize>>),
    /// Peptide line, stored as hash map with the columns names from PSH
    PSM(usize, String, Vec<Range<usize>>),
}

/// Parse a mzTab file
/// # Errors
/// If the file is not a valid mzTab file
fn parse_mztab_reader<T: BufRead>(
    reader: T,
) -> impl Iterator<Item = Result<Option<MZTabLine>, CustomError>> {
    reader.lines().enumerate().map(move |(line_index, line)| {
        line.map_err(|err| {
            CustomError::error(
                "Could not read line",
                err,
                Context::full_line(line_index, "(failed)"),
            )
        })
        .and_then(|line| {
            if line.trim().is_empty() {
                Ok(None)
            } else {
                crate::identification::csv::csv_separate(&line, b'\t').map(|fields| {
                    match &line[fields[0].clone()] {
                        "MTD" => Some(MZTabLine::MTD(line_index, line, fields)),
                        "PSH" => Some(MZTabLine::PSH(line_index, line, fields)),
                        "PSM" => Some(MZTabLine::PSM(line_index, line, fields)),
                        _ => None,
                    }
                })
            }
        })
    })
}
