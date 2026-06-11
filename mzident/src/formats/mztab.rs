use std::{
    borrow::Cow,
    cmp::Ordering,
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    marker::PhantomData,
    num::NonZeroUsize,
    ops::Range,
    str::FromStr,
    sync::Arc,
};

use context_error::*;
use flate2::bufread::GzDecoder;
use itertools::Itertools;
use mzcv::{ControlledVocabulary, Term, term};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    CVTerm, FastaIdentifier, KnownFileFormat, MaybePeptidoform, PSM, PSMData, PSMMetaData,
    ProteinMetaData, SpectrumId, SpectrumIds,
    helper_functions::{
        check_extension, explain_number_error, float_digits, next_number, split_with_brackets,
    },
    mztab_writer::{
        MzTabAssay, MzTabCV, MzTabContact, MzTabInstrument, MzTabKind, MzTabMSRun, MzTabMetadata,
        MzTabMode, MzTabSample, MzTabSoftware, MzTabStudyVariable,
    },
};
use mzcore::{
    chemistry::MolecularFormula,
    ontology::Ontologies,
    quantities::Tolerance,
    sequence::{
        AminoAcid, FlankingSequence, MUPSettings, PeptideModificationSearch, Peptidoform,
        PeptidoformIonSet, SequencePosition, SimpleLinear, SimpleModification,
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
    /// The metadata of the file
    pub metadata: Arc<MzTabMetadata>,
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
            + self.additional.space()
            + self.metadata.space())
        .set_total::<Self>()
    }
}

impl MzTabPSM {
    /// Parse a mzTab file which can be gzipped which is detected from the path. Returns the metadata, all proteins, and an iterator over the PSMs. The file is consumed until the first PSM line is detected and is consumed lazily after that.
    /// # Errors
    /// If the file is not in the correct format
    pub fn parse_file(
        path: impl AsRef<std::path::Path>,
        ontologies: &Ontologies,
    ) -> Result<
        (
            Arc<MzTabMetadata>,
            Arc<HashMap<String, Arc<MzTabProtein>>>,
            Box<dyn Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + '_>,
        ),
        BoxedError<'static, BasicKind>,
    > {
        let path = path.as_ref();
        let file = File::open(path).map_err(|e| {
            BoxedError::new(
                BasicKind::Error,
                "Could not open file",
                e.to_string(),
                Context::default().source(path.to_string_lossy()).to_owned(),
            )
        })?;
        if check_extension(path, "gz") {
            Self::parse_reader(
                BufReader::new(GzDecoder::new(BufReader::new(file))),
                ontologies,
                Context::default().source(path.to_string_lossy()).to_owned(),
            )
            .map(|(m, p, i)| {
                let i: Box<dyn Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + '_> =
                    Box::new(i);
                (m, p, i)
            })
        } else {
            Self::parse_reader(
                BufReader::new(file),
                ontologies,
                Context::default().source(path.to_string_lossy()).to_owned(),
            )
            .map(|(m, p, i)| {
                let i: Box<dyn Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + '_> =
                    Box::new(i);
                (m, p, i)
            })
        }
    }

    /// Parse a mzTab file directly from a buffered reader (needed for the `.lines()` function). Returns the metadata, all proteins, and an iterator over the PSMs. The file is consumed until the first PSM line is detected and is consumed lazily after that.
    pub fn parse_reader<'a, T: BufRead + 'a>(
        reader: T,
        ontologies: &'a Ontologies,
        base: Context<'a>,
    ) -> Result<
        (
            Arc<MzTabMetadata>,
            Arc<HashMap<String, Arc<MzTabProtein>>>,
            impl Iterator<Item = Result<Self, BoxedError<'static, BasicKind>>> + 'a,
        ),
        BoxedError<'static, BasicKind>,
    > {
        let mut peptide_header: Option<Vec<String>> = None;
        let mut protein_header: Option<Vec<String>> = None;
        let mut proteins: HashMap<String, Arc<MzTabProtein>> = HashMap::new();
        let mut metadata = Arc::new(MzTabMetadata::default());
        let mut line_iter = parse_mztab_reader(reader).peekable();

        while let Some(Ok(item)) = line_iter.peek()
            && item
                .as_ref()
                .is_none_or(|i| !matches!(i, MzTabLine::PSM(..)))
        {
            if let Some(Ok(Some(item))) = line_iter.next() {
                match item {
                    MzTabLine::MTD(line_index, line, fields) => {
                        let Some(metadata) = Arc::get_mut(&mut metadata) else {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MTD line",
                                "Metadata lines need to be defined before all other line types",
                                base.clone()
                                    .lines(0, line)
                                    .line_index(line_index)
                                    .to_owned(),
                            ));
                        };
                        match parse_metadata(
                            metadata,
                            &line,
                            &fields,
                            &base.clone().lines(0, line.clone()).line_index(line_index),
                            ontologies,
                        ) {
                            Ok(()) => (),
                            Err(err) => return Err(err.to_owned()),
                        }
                    }
                    MzTabLine::PRH(line_index, line, fields) => {
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
                                return Err(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid protein table",
                                    format!("The required column '{required}' is not present"),
                                    base.clone()
                                        .lines(0, line)
                                        .line_index(line_index)
                                        .to_owned(),
                                ));
                            }
                        }
                        protein_header = Some(header);
                    }
                    MzTabLine::PRT(line_index, line, fields) => {
                        match PSMLine::new(
                            base.clone()
                                .lines(0, line.clone())
                                .line_index(line_index)
                                .to_owned(),
                            protein_header.as_deref(),
                            &line,
                            &fields,
                        )
                        .and_then(|line| MzTabProtein::from_line(&line, metadata.clone()))
                        {
                            Ok(protein) => {
                                for name in &protein.ambiguity_members {
                                    proteins.insert(name.clone(), protein.clone());
                                }
                                proteins.insert(protein.accession.clone(), protein);
                            }
                            Err(err) => return Err(err.to_owned()),
                        }
                    }

                    MzTabLine::PSH(line_index, line, fields) => {
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
                                return Err(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid peptide table",
                                    format!("The required column '{required}' is not present"),
                                    base.clone()
                                        .lines(0, line)
                                        .line_index(line_index)
                                        .to_owned(),
                                ));
                            }
                        }
                        peptide_header = Some(header);
                    }
                    MzTabLine::PSM(..) => unreachable!(),
                }
            }
        }
        let proteins = Arc::new(proteins);

        Ok((
            metadata.clone(),
            proteins.clone(),
            line_iter.filter_map(move |item| {
                item.transpose().map(|item| match item {
                    Ok(MzTabLine::PSM(line_index, line, fields)) => PSMLine::new(
                            base.clone()
                                .lines(0, line.clone())
                                .line_index(line_index)
                                .to_owned(),
                            peptide_header.as_deref(),
                            &line,
                            &fields,
                        )
                        .and_then(|line| {
                            Self::from_line(&line, ontologies, &proteins, metadata.clone())
                                .map_err(BoxedError::to_owned)
                        }),
                    Ok(MzTabLine::MTD(line_index, ..) | MzTabLine::PRH(line_index, ..) | MzTabLine::PRT(line_index,.. ) | MzTabLine::PSH(line_index,.. )) =>
                        Err(BoxedError::new(BasicKind::Warning, "Invalid line type", "After the first PSM line no MTD, PRH, PRT, or PSH lines are allowed anymore", base.clone().line_index(line_index).to_owned())),
                    Err(e) => Err(e),
                })
            }),
        ))
    }

    /// Parse a single PSM line
    /// # Errors
    /// When not in the correct format
    fn from_line<'a>(
        line: &PSMLine<'a>,
        ontologies: &Ontologies,
        proteins: &HashMap<String, Arc<MzTabProtein>>,
        metadata: Arc<MzTabMetadata>,
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
                    .map(|field| parse_modification(line.line, field, ontologies, &line.context))
                    .collect::<Result<Vec<_>, BoxedError<'_, BasicKind>>>()?
            };

        let mut result = Self {
            metadata: metadata.clone(),
            peptidoform: {
                let range = line
                    .required_column("sequence")
                    .map_err(BoxedError::to_owned)?
                    .1;

                if range.is_empty() {
                    None
                } else {
                    let mut peptide: Peptidoform<SimpleLinear> = Peptidoform::pro_forma_or_sloppy(
                        &line.context,
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
                                            line.context.clone().add_highlight((0, range.clone()))))?,
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
                                            line.context.clone().add_highlight((0, range.clone()))))?, score))
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
                        PeptideModificationSearch::in_modifications(
                            metadata
                                .fixed_mods
                                .iter()
                                .chain(metadata.variable_mods.iter())
                                .cloned()
                                .collect(),
                        )
                        .tolerance(Tolerance::new_ppm(20.0))
                        .search(peptide),
                    )
                }
            },
            id: {
                let col = line
                    .required_column("psm_id")
                    .map_err(BoxedError::to_owned)?;
                col.0.parse().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab PSM_ID",
                        format!("The PSM_ID {}", explain_number_error(&err)),
                        line.context.clone().add_highlight((0, col.1)),
                    )
                })?
            },
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
                                .and_then(|(v, inner_range)| {
                                    (!v.eq_ignore_ascii_case("null")).then(|| {
                                        v.parse::<f64>().map_err(|err| {
                                            BoxedError::new(BasicKind::Error,
                                                "Invalid mzTab search engine score",
                                                format!(
                                        "The search engine score can not be parsed as f64: {err}"
                                    ),
                                    line.context.clone().add_highlight((0, inner_range))
                                            )
                                        })
                                    })
                                })
                                .transpose()
                                .and_then(|score| {
                                    CVTerm::from_str(s)
                                        .map_err(|e| {
                                            e.replace_context(line.context.clone().add_highlight((0, range.clone())))
                                        })
                                        .and_then(|engine| Ok((engine, score, metadata.psm_search_engines.get(i).ok_or_else(|| BoxedError::new(BasicKind::Error,"Missing search engine score type", "All search engines require a defined search type", 
                                        line.context.clone().add_highlight((0, range.clone()))))?.clone())))
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
                        line.context.clone().add_highlight((0, range)),
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
                                    line.context.clone().add_highlight((0, r)),
                                )
                            })
                            .map(|v| Time::new::<mzcore::system::s>(v)) // TODO: use unit from metadata if present
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
                                line.context.clone().add_highlight((0, range)),
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
                                    line.context.clone().add_highlight((0, r)),
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
                        line.context.clone().add_highlight((0, range.clone()))
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
                                line.context.clone().add_highlight((0, range.clone()))
                            )
                        })
                        .and_then(|v| if v == 0 {
                            Err(BoxedError::new(BasicKind::Error,
                                "Invalid mzTab ms_run",
                                "A ms_run identifier uses 1 based indexing and so cannot be 0.",
                                line.context.clone().add_highlight((0, range.clone()))
                            ))
                        } else {Ok(v-1)})?;
                    let path = metadata.ms_runs.get(index).ok_or_else(|| BoxedError::new(BasicKind::Error,"Missing raw file definition", "All raw files should be defined in the MTD section before being used in the PSM Section", 
                    line.context.clone().add_highlight((0, range.clone()))))?.location.clone();

                    let id = match scan_id.split_once('=') {
                        Some(("scan", num)) if num.chars().all(|c| c.is_ascii_digit()) => SpectrumId::Number(num.parse().map_err(|err| {
                            BoxedError::new(BasicKind::Error,
                                "Invalid mzTab spectra_ref scan number",
                                format!("The spectra_ref scan number {}", explain_number_error(&err)),
                                line.context.clone().add_highlight((0, range.clone()))
                            )
                        })?),
                        Some(("index", index)) if index.chars().all(|c| c.is_ascii_digit()) => SpectrumId::Index(index.parse().map_err(|err| {
                            BoxedError::new(BasicKind::Error,
                                "Invalid mzTab spectra_ref index",
                                format!("The spectra_ref index {}", explain_number_error(&err)),
                                line.context.clone().add_highlight((0, range.clone()))
                            )
                        })?),
                        _ =>  SpectrumId::Native(scan_id.to_owned()),
                    };

                    Ok((path, id))
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
                                line.context.clone().add_highlight((0, pre.1)),
                            )
                        })?)
                    }
                    _ => FlankingSequence::Sequence(Box::new(Peptidoform::sloppy_pro_forma_inner(
                        &line.context,
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
                                line.context.clone().add_highlight((0, post.1)),
                            )
                        })?)
                    }
                    _ => FlankingSequence::Sequence(Box::new(Peptidoform::sloppy_pro_forma_inner(
                        &line.context,
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
                                line.context.clone().add_highlight((0, r)),
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
                                        line.context.clone().add_highlight((0, r)),
                                    )
                                })
                            })
                        })
                        .transpose()?;
                    Ok(start.and_then(|s| end.map(|e| s..e)))
                })?,
            local_confidence: line
                .optional_column("opt_global_ms_1003984_amino_acid_confidence_level") // TODO: update if the final accepted name changes
                .or_else(|| line.optional_column("opt_ms_run[1]_aa_scores"))
                .filter(|(lc, _)| !lc.trim().is_empty() && !lc.trim().eq_ignore_ascii_case("null"))
                .map(|(v, r)| {
                    v.split(',')
                        .map(|score| {
                            score.parse::<f64>().map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid mzTab local confidence",
                                    format!("The local confidence can not be parsed: {err}"),
                                    line.context.clone().add_highlight((0, r.clone())),
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
            match lc.len().cmp(&pep_len) {
                Ordering::Greater => {
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
                            line.context.clone()))
                    }
                }
                Ordering::Equal => {
                    Ok(lc.clone())
                }
                Ordering::Less => {
                    Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid local confidence", 
                    format!("The number of elements in the local confidence ({}) does not match the number of amino acids ({pep_len}).", lc.len()), 
                    line.context.clone()))
                }
            }
        }).transpose()?;

        Ok(result)
    }
}

fn parse_metadata<'a>(
    metadata: &mut MzTabMetadata,
    line: &'a str,
    fields: &'a [Range<usize>],
    context: &Context<'a>,
    ontologies: &Ontologies,
) -> Result<(), BoxedError<'a, BasicKind>> {
    if fields.len() == 3 {
        match line[fields[1].clone()].to_ascii_lowercase().as_str() {
            m if (m.starts_with("variable_mod[") || m.starts_with("fixed_mod["))
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
                        + offset
                            .map_or(0, |s| s.bytes().take_while(u8::is_ascii_whitespace).count())
                        ..text_range.start
                            + first.map_or(0, |f| f.len() + 1)
                            + offset.map_or(0, str::len);
                    if &line[m_range.clone()] != "MS:1002453"
                        && &line[m_range.clone()] != "MS:1002454"
                    {
                        match parse_single_modification(line, m_range, ontologies, context) {
                            Ok(modification) => {
                                if m.starts_with("variable") {
                                    metadata.variable_mods.push(modification);
                                } else {
                                    metadata.fixed_mods.push(modification);
                                }
                            }
                            Err(e) => return Err(e.to_owned()),
                        }
                    }
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab modification",
                        "A modification has to be enclosed by square brackets",
                        context.clone().add_highlight((0, fields[2].clone())),
                    ));
                }
            }
            m if m.starts_with("psm_search_engine_score[") && m.ends_with(']') => {
                match CVTerm::from_str(&line[fields[2].clone()]) {
                    Ok(term) => metadata.psm_search_engines.push(term),
                    Err(err) => return Err(err),
                }
            }
            m if m.starts_with("protein_search_engine_score[") && m.ends_with(']') => {
                match CVTerm::from_str(&line[fields[2].clone()]) {
                    Ok(term) => metadata.protein_search_engines.push(term),
                    Err(err) => return Err(err),
                }
            }
            "description" => metadata.description = line[fields[2].clone()].into(),
            "title" => metadata.title = line[fields[2].clone()].into(),
            "mztab-id" => metadata.id = line[fields[2].clone()].into(),
            "mztab-mode" => {
                metadata.mode = match &line[fields[2].clone()] {
                    "Summary" => MzTabMode::Summary,
                    "Complete" => MzTabMode::Complete,
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzTab mode",
                            "Only Summary and Complete modes are allowed",
                            context.clone().add_highlight((0, fields[2].clone())),
                        ));
                    }
                }
            }
            "mztab-type" => {
                metadata.kind = match &line[fields[2].clone()] {
                    "Identification" => MzTabKind::Identification,
                    "Quantification" => MzTabKind::Quantification,
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzTab type",
                            "Only Identification and Quantification types are allowed",
                            context.clone().add_highlight((0, fields[2].clone())),
                        ));
                    }
                }
            }
            m if m.starts_with("publication[") && m.ends_with(']') => {
                metadata.publication.push(line[fields[2].clone()].into());
            }
            m if m.starts_with("uri[") && m.ends_with(']') => {
                metadata.uri.push(line[fields[2].clone()].into());
            }
            "protein-quantification-unit" => {
                metadata.protein_quantification_unit =
                    match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(term) => Some(term),
                        Err(err) => return Err(err),
                    }
            }
            "custom" => {
                metadata
                    .custom
                    .push(CVTerm::from_str(&line[fields[2].clone()])?);
            }
            "false_discovery_rate" => {
                metadata.false_discovery_rate = line[fields[2].clone()]
                    .split('|')
                    .map(CVTerm::from_str)
                    .collect::<Result<ThinVec<_>, _>>()
                    .map_err(|err| {
                        err.replace_context(context.clone().add_highlight((0, fields[2].clone())))
                    })?;
            }
            "quantification_method" => {
                metadata.quantification_method = match CVTerm::from_str(&line[fields[2].clone()]) {
                    Ok(term) => Some(term),
                    Err(err) => return Err(err),
                }
            }
            m if m.starts_with("software[") => {
                let Some((index, tag)) = m.trim_start_matches("software[").split_once(']') else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab software parameter",
                        "An software parameter needs to be 'software[n]' or 'software[n]-setting[n]'.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab software identifier",
                        format!("The software identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.software.len() <= index {
                    metadata.software.push(MzTabSoftware {
                        name: CVTerm {
                            term: term!(MS:1000531|software),
                            value: Box::default(),
                        },
                        settings: ThinVec::new(),
                    });
                }

                if tag.is_empty() {
                    metadata.software[index].name = CVTerm::from_str(&line[fields[2].clone()])?;
                } else if tag.starts_with("-setting[") {
                    metadata.software[index]
                        .settings
                        .push(line[fields[2].clone()].to_string());
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab software parameter",
                        "An software parameter needs to be 'software[n]' or 'software[n]-setting[n]'.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                }
            }
            m if m.starts_with("ms_run[") => {
                let Some((index, tag)) = m.trim_start_matches("ms_run[").split_once("]-") else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab ms_run parameter",
                        "An ms_run parameter needs to be 'ms_run[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab ms_run identifier",
                        format!("The ms_run identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.ms_runs.len() <= index {
                    metadata.ms_runs.push(MzTabMSRun::default());
                }

                let elem = &mut metadata.ms_runs[index];

                match tag {
                    "location" => elem.location = line[fields[2].clone()].into(),
                    "format" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.format = Some(i),
                        Err(err) => return Err(err),
                    },
                    "id_format" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.id_format = Some(i),
                        Err(err) => return Err(err),
                    },
                    "fragmentation_method" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.fragmentation_method = Some(i),
                        Err(err) => return Err(err),
                    },
                    "hash_method" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.hash_method = Some(i),
                        Err(err) => return Err(err),
                    },
                    "hash" => elem.hash = Some(line[fields[2].clone()].to_string()),
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MTD line",
                            "This parameter does not exist for an ms_run",
                            context.clone().add_highlight((0, fields[0].clone())),
                        ));
                    }
                }
            }
            m if m.starts_with("assay[") => {
                let Some((index, tag)) = m.trim_start_matches("assay[").split_once("]-") else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab assay parameter",
                        "An assay parameter needs to be 'assay[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab assay identifier",
                        format!("The assay identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.assay.len() <= index {
                    metadata.assay.push(MzTabAssay::default());
                }

                let elem = &mut metadata.assay[index];

                match tag {
                    "ms_run_ref" => {
                        elem.ms_run_ref = Some(parse_ref(
                            "ms_run",
                            line[fields[2].clone()].trim(),
                            &context.clone().add_highlight((0, fields[2].clone())),
                        )?);
                    }
                    "quantification_reagent" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.quantification_reagent = i,
                        Err(err) => return Err(err),
                    },
                    "sample_ref" => {
                        elem.sample_ref = Some(parse_ref(
                            "sample",
                            line[fields[2].clone()].trim(),
                            &context.clone().add_highlight((0, fields[2].clone())),
                        )?);
                    }
                    _ => {
                        if let Some(m) = tag.strip_prefix("quantification_mod[") {
                            if m.ends_with(']') {
                                elem.quantification_mod.push(parse_single_modification(
                                    line,
                                    fields[2].clone(),
                                    ontologies,
                                    context,
                                )?);
                            }
                        } else {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MTD line",
                                "This parameter does not exist for an ms_run",
                                context.clone().add_highlight((0, fields[0].clone())),
                            ));
                        }
                    }
                }
            }
            m if m.starts_with("instrument[") => {
                let Some((index, tag)) = m.trim_start_matches("instrument[").split_once("]-")
                else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab instrument parameter",
                        "An instrument parameter needs to be 'instrument[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab instrument identifier",
                        format!("The instrument identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.instruments.len() <= index {
                    metadata.instruments.push(MzTabInstrument {
                        name: term!(MS:1000463|instrument).into(),
                        source: term!(MS:1000458|source).into(),
                        detector: term!(MS:1000453|detector).into(),
                        analyser: Vec::new(),
                    });
                }

                let elem = &mut metadata.instruments[index];

                match tag {
                    "name" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.name = i,
                        Err(err) => return Err(err),
                    },
                    "source" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.source = i,
                        Err(err) => return Err(err),
                    },
                    "detector" => match CVTerm::from_str(&line[fields[2].clone()]) {
                        Ok(i) => elem.detector = i,
                        Err(err) => return Err(err),
                    },
                    _ => {
                        if let Some(m) = tag.strip_prefix("analyzer[") {
                            if m.ends_with(']') {
                                elem.analyser
                                    .push(CVTerm::from_str(&line[fields[2].clone()])?);
                            }
                        } else {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MTD line",
                                "This parameter does not exist for an instrument",
                                context.clone().add_highlight((0, fields[0].clone())),
                            ));
                        }
                    }
                }
            }
            m if m.starts_with("contact[") => {
                let Some((index, tag)) = m.trim_start_matches("contact[").split_once("]-") else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab contact parameter",
                        "An contact parameter needs to be 'contact[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab contact identifier",
                        format!("The contact identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.contact.len() <= index {
                    metadata.contact.push(MzTabContact::default());
                }

                let elem = &mut metadata.contact[index];

                match tag {
                    "name" => elem.name = line[fields[2].clone()].into(),
                    "affiliation" => elem.affiliation = line[fields[2].clone()].into(),
                    "email" => elem.email = line[fields[2].clone()].into(),
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MTD line",
                            "This parameter does not exist for a contact",
                            context.clone().add_highlight((0, fields[0].clone())),
                        ));
                    }
                }
            }
            m if m.starts_with("sample[") => {
                let Some((index, tag)) = m.trim_start_matches("sample[").split_once("]-") else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab sample parameter",
                        "An sample parameter needs to be 'sample[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab sample identifier",
                        format!("The sample identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.sample.len() <= index {
                    metadata.sample.push(MzTabSample::default());
                }

                let elem = &mut metadata.sample[index];

                match tag {
                    "description" => elem.description = line[fields[2].clone()].into(),
                    _ => match tag.split_once('[') {
                        Some(("species", _)) => elem
                            .species
                            .push(CVTerm::from_str(&line[fields[2].clone()])?),
                        Some(("tissue", _)) => elem
                            .tissue
                            .push(CVTerm::from_str(&line[fields[2].clone()])?),
                        Some(("cell_type", _)) => elem
                            .cell_type
                            .push(CVTerm::from_str(&line[fields[2].clone()])?),
                        Some(("disease", _)) => elem
                            .disease
                            .push(CVTerm::from_str(&line[fields[2].clone()])?),
                        Some(("custom", _)) => elem
                            .custom
                            .push(CVTerm::from_str(&line[fields[2].clone()])?),
                        _ => {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MTD line",
                                "This parameter does not exist for a sample",
                                context.clone().add_highlight((0, fields[0].clone())),
                            ));
                        }
                    },
                }
            }
            m if m.starts_with("study_variable[") => {
                let Some((index, tag)) = m.trim_start_matches("study_variable[").split_once("]-")
                else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab study_variable parameter",
                        "An study_variable parameter needs to be 'study_variable[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab study_variable identifier",
                        format!(
                            "The study_variable identifier {}",
                            explain_number_error(&err)
                        ),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.study_variable.len() <= index {
                    metadata.study_variable.push(MzTabStudyVariable::default());
                }

                let elem = &mut metadata.study_variable[index];

                match tag {
                    "description" => elem.description = line[fields[2].clone()].into(),
                    "sample_refs" => {
                        elem.sample_refs = parse_refs(
                            "sample",
                            line[fields[2].clone()].trim(),
                            &context.clone().add_highlight((0, fields[2].clone())),
                        )?;
                    }
                    "assay_refs" => {
                        elem.assay_refs = parse_refs(
                            "assay",
                            line[fields[2].clone()].trim(),
                            &context.clone().add_highlight((0, fields[2].clone())),
                        )?;
                    }
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MTD line",
                            "This parameter does not exist for a study_variable",
                            context.clone().add_highlight((0, fields[0].clone())),
                        ));
                    }
                }
            }
            m if m.starts_with("cv[") => {
                let Some((index, tag)) = m.trim_start_matches("cv[").split_once("]-") else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab cv parameter",
                        "An cv parameter needs to be 'cv[0]-<name>' where the 0 is the index and the name can be any of the known parameters.",
                        context.clone().add_highlight((0, fields[1].clone())),
                    ));
                };
                let index = match index.parse::<NonZeroUsize>().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzTab cv identifier",
                        format!("The cv identifier {}", explain_number_error(&err)),
                        context.clone().add_highlight((0, fields[1].clone())),
                    )
                }) {
                    Ok(i) => i.get() - 1,
                    Err(err) => return Err(err),
                };

                while metadata.cv.len() <= index {
                    metadata.cv.push(MzTabCV::default());
                }

                let elem = &mut metadata.cv[index];

                match tag {
                    "label" => elem.label = line[fields[2].clone()].into(),
                    "full_name" => elem.full_name = line[fields[2].clone()].into(),
                    "version" => elem.version = line[fields[2].clone()].into(),
                    "url" => elem.url = line[fields[2].clone()].into(),
                    _ => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid MTD line",
                            "This parameter does not exist for a cv",
                            context.clone().add_highlight((0, fields[0].clone())),
                        ));
                    }
                }
            }
            m if m.starts_with("sample_processing[") && m.ends_with(']') => {
                metadata.sample_processing.push(
                    line[fields[2].clone()]
                        .split('|')
                        .map(CVTerm::from_str)
                        .collect::<Result<Vec<_>, _>>()
                        .map_err(|err| {
                            err.replace_context(
                                context.clone().add_highlight((0, fields[2].clone())),
                            )
                        })?,
                );
            }
            "colunit-protein" => {
                let Some((column, param)) = line[fields[2].clone()].split_once('=') else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid colunit-protein line",
                        "A colunit line must contain '{column}={param}' but the '=' is missing",
                        context.clone().add_highlight((0, fields[0].clone())),
                    ));
                };
                metadata.colunit_protein.push((
                    column.parse().map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid column name",
                            "A column name has to be an exsting column or an optional column",
                            context.clone().add_highlight((0, fields[0].clone())),
                        )
                    })?,
                    param.parse()?,
                ));
            }
            "colunit-psm" => {
                let Some((column, param)) = line[fields[2].clone()].split_once('=') else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid colunit-psm line",
                        "A colunit line must contain '{column}={param}' but the '=' is missing",
                        context.clone().add_highlight((0, fields[0].clone())),
                    ));
                };
                metadata.colunit_psm.push((
                    column.parse().map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid column name",
                            "A column name has to be an exsting column or an optional column",
                            context.clone().add_highlight((0, fields[0].clone())),
                        )
                    })?,
                    param.parse()?,
                ));
            }
            _ => (),
        }
        Ok(())
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid MTD line",
            "MTD lines should contain three columns (the tag, key, and value)",
            context.clone(),
        ))
    }
}

fn parse_refs<'a>(
    ty: &'static str,
    value: &str,
    context: &Context<'a>,
) -> Result<ThinVec<NonZeroUsize>, BoxedError<'a, BasicKind>> {
    value
        .split(',')
        .map(|s| parse_ref(ty, s, context))
        .collect()
}

fn parse_ref<'a>(
    ty: &'static str,
    value: &str,
    context: &Context<'a>,
) -> Result<NonZeroUsize, BoxedError<'a, BasicKind>> {
    if let Some(boxed) = value.strip_prefix(ty)
        && let Some(tail) = boxed.strip_prefix('[')
        && let Some(num) = tail.strip_suffix(']')
    {
        num.parse().map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                format!("Invalid {ty} reference"),
                format!("A {ty} reference {}", explain_number_error(&err)),
                context.clone(),
            )
        })
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            format!("Invalid {ty} reference"),
            format!("A {ty} reference should look like '{ty}[n]'"),
            context.clone(),
        ))
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
    pub taxid: Option<u32>,
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
    /// Any additional metadata
    pub additional: HashMap<String, String>,
    /// The metadata of the file
    pub metadata: Arc<MzTabMetadata>,
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
            + self.uri.space()
            + self.additional.space()
            + self.metadata.space())
        .set_total::<Self>()
    }
}

impl MzTabProtein {
    /// Parse a single PRT line
    /// # Errors
    /// When not in the correct format
    fn from_line(
        line: &PSMLine<'_>,
        metadata: Arc<MzTabMetadata>,
    ) -> Result<Arc<Self>, BoxedError<'static, BasicKind>> {
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
                    if v.eq_ignore_ascii_case("null") {
                        Ok(None)
                    } else {
                        v.parse::<u32>().map(Some).map_err(|e| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid PRT Line",
                                format!("The taxid number {}", explain_number_error(&e)),
                                line.context.clone().add_highlight((0, location)),
                            )
                        })
                    }
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
                        line.context.clone().add_highlight((0, range)).to_owned(),
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
                .filter(|(v, _)| !v.eq_ignore_ascii_case("null") && !v.is_empty())
                .map(|(v, _)| v.to_string()),
            go_terms: line
                .optional_column("go_terms")
                .filter(|(v, _)| !v.eq_ignore_ascii_case("null") && !v.is_empty())
                .map(|(v, l)| {
                    let mut offset = 0;
                    v.split('|')
                        .map(|s| {
                            let r = s.parse::<Curie>().map_err(|e| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid CV accession",
                                    e.description(),
                                    line.context
                                        .clone()
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
                .filter(|(v, _)| !v.eq_ignore_ascii_case("null") && !v.is_empty())
                .map(|(v, location)| {
                    v.parse::<f64>().map_err(|e| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid PRT Line",
                            format!("The protein coverage {e}"),
                            line.context.clone().add_highlight((0, location)).to_owned(),
                        )
                    })
                })
                .transpose()?,
            additional: line
                .header
                .iter()
                .enumerate()
                .filter(|(_, column)| column.starts_with("opt"))
                .map(|(index, column)| {
                    (
                        column.clone(),
                        line.line[line.fields[index].clone()].to_string(),
                    )
                })
                .collect(),
            metadata,
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
    ontologies: &Ontologies,
    context: &Context<'a>,
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
                                    context.clone().add_highlight((0, index..index + offset))
                                )
                            })?;

                            let parameter = &single[offset..];
                            let score = if !parameter.is_empty()
                                && single[offset..].starts_with('[')
                                && single.ends_with(']')
                            {
                                let value_range = CVTerm::parse_and_identify(line, index + offset..index+single.len(), "MS:1001876", "modification probability", context)?;
                                Some(line[value_range.clone()].parse::<f64>().map(OrderedFloat).map_err(|err| BoxedError::new(BasicKind::Error,
                                                "Invalid modification position",
                                                format!("The modification probability is not a valid number: {err}"),
                                                context.clone().add_highlight((0, value_range))
                                )))
                            } else if parameter.is_empty() {
                                None
                            } else {
                                Some(Err(BoxedError::new(BasicKind::Error,
                                    "Invalid modification position",
                                    "A modification position parameter should be enclosed in square brackets '[]'",
                                    context.clone().add_highlight((0, index + offset..index + single.len())),
                                )))
                            }.transpose()?;
                            index += single.len() + 1;

                            Ok((res, score))
                        } else {
                            Err(BoxedError::new(BasicKind::Error,
                                "Invalid modification position",
                                "A modification position should start with a number",
                                context.clone().add_highlight((0, index..index + single.len()))
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
                context,
            )?;
            let loss = line[value_range.clone()]
                .parse::<f64>()
                .map(|v| Mass::new::<mzcore::system::dalton>(v))
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid neutral loss",
                        format!("The fragment neutral loss is not a valid number: {err}"),
                        context.clone().add_highlight((0, value_range)),
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
                            context
                                .clone()
                                .add_highlight((0, range.start..range.start + pos.len())),
                        ))
                    }
                },
            )
        } else {
            let modification = parse_single_modification(
                line,
                range.start + pos.len() + 1..range.end,
                ontologies,
                context,
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
        let value_range = CVTerm::parse_and_identify(
            line,
            range,
            "MS:1001524",
            "fragment neutral loss",
            context,
        )?;
        line[value_range.clone()]
            .parse::<f64>()
            .map(|v| Mass::new::<mzcore::system::dalton>(v))
            .map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid neutral loss",
                    format!("The fragment neutral loss is not a valid number: {err}"),
                    context.clone().add_highlight((0, value_range)),
                )
            })
            .map(|loss| MzTabReturnModification::NeutralLoss(None, loss))
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid modification",
            "A modification should be the position followed by a hyphen ('-') followed by the modification",
            context.clone().add_highlight((0, range)),
        ))
    }
}

/// Parse a single mzTab modification.
/// # Errors
/// if it does not follow the UNIMOD/MOD/CHEMMOD rules.
fn parse_single_modification<'a>(
    line: &'a str,
    range: Range<usize>,
    ontologies: &Ontologies,
    context: &Context<'a>,
) -> Result<SimpleModification, BoxedError<'a, BasicKind>> {
    if let Some((tag, value)) = line[range.clone()].split_once(':') {
        let value_context = context.clone().add_highlight((0, range.clone()));
        let modification = if tag.eq_ignore_ascii_case("unimod") {
            ontologies
                .unimod()
                .get_by_index(&mzcv::AccessionCode::Numeric(
                    value.parse::<u32>().map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid unimod code",
                            format!("The unimod modification {}", explain_number_error(&err)),
                            value_context.clone(),
                        )
                    })?,
                ))
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
                .get_by_index(&mzcv::AccessionCode::Numeric(
                    value.parse::<u32>().map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid PSI-MOD code",
                            format!("The PSI-MOD modification {}", explain_number_error(&err)),
                            value_context.clone(),
                        )
                    })?,
                ))
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
                .get_by_index(&mzcv::AccessionCode::Numeric(
                    value.parse::<u32>().map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid custom code",
                            format!("The custom modification {}", explain_number_error(&err)),
                            value_context.clone(),
                        )
                    })?,
                ))
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
                            context
                                .clone()
                                .add_highlight((0, range.start + tag.len() + 1, 1)),
                        ));
                    }
                };
                MolecularFormula::pro_forma_inner::<false, false>(
                    &Context::default().lines(0, line),
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
                context.clone().add_highlight((0, range)),
            ));
        };

        Ok(modification)
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid mzTab modification",
            "An mzTab modification should be in format 'tag:value' but the colon (':') is missing",
            context.clone().add_highlight((0, range)),
        ))
    }
}

#[derive(Clone, Debug)]
struct PSMLine<'a> {
    context: Context<'a>,
    header: &'a [String],
    pub line: &'a str,
    fields: &'a [Range<usize>],
}

impl<'a> PSMLine<'a> {
    /// Form an indexable line out of a set of fields
    /// # Errors
    /// When there is no header or the line has a different number of columns
    fn new(
        context: Context<'static>,
        header: Option<&'a [String]>,
        line: &'a str,
        fields: &'a [Range<usize>],
    ) -> Result<Self, BoxedError<'static, BasicKind>> {
        let header = header.ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Missing PSH line",
                "The PSH peptide header line should precede any PSM line",
                context.clone().to_owned(),
            )
        })?;
        if header.len() == fields.len() {
            Ok(Self {
                context,
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
                context.to_owned(),
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
                self.context.clone(),
            )
        })
    }
}

impl std::fmt::Display for PSMLine<'_> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.context
                .clone()
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
        context: &Context<'a>,
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
                    context.clone().add_highlight((0, range)),
                ))
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid CV term",
                "A CV term should be enclosed by '[]'",
                context.clone().add_highlight((0, range)),
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
                        Context::default().lines(0, value).to_owned(),
                    )
                })?;
            let name = split.next().unwrap_or_default().trim().to_string();
            Ok(Self {
                term: Term {
                    accession,
                    name: name.into(),
                },
                value: split.next().unwrap_or_default().trim().into(),
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid CV term",
                "A CV term should be enclosed by '[]'",
                Context::default().lines(0, value).to_owned(),
            ))
        }
    }
}

/// The reliability of a PSM
#[expect(missing_docs)]
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
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
    fn peptidoform_ion_set(&self) -> Option<Cow<'_, PeptidoformIonSet>> {
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

    fn id(&self) -> FastaIdentifier<Cow<'_, str>> {
        FastaIdentifier::Undefined(false, Cow::Borrowed(&self.accession))
    }

    fn description(&self) -> Option<&str> {
        Some(&self.description)
    }

    fn species(&self) -> Option<Curie> {
        self.taxid.map(|taxid| Curie {
            cv: ControlledVocabulary::NCBITaxon,
            accession: mzcv::AccessionCode::Numeric(taxid),
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

    fn database(&self) -> Option<(Cow<'_, str>, Option<Cow<'_, str>>)> {
        Some((
            Cow::Borrowed(&self.database),
            Some(Cow::Borrowed(&self.database_version)),
        ))
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
