use std::{
    borrow::Cow,
    io::{BufRead, BufReader},
    marker::PhantomData,
    num::ParseIntError,
    ops::Range,
    path::Path,
    str::FromStr,
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    helper_functions::explain_number_error,
    identification::{
        IdentifiedPeptidoform, IdentifiedPeptidoformData, KnownFileFormat, MetaData,
        PeptidoformPresent, SpectrumIds,
    },
    prelude::{CompoundPeptidoformIon, HasPeptidoformImpl},
    sequence::{
        AminoAcid, AnnotatedPeptide, Annotation, Peptidoform, Region, SemiAmbiguous,
        SequenceElement,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

/// A single parsed line of a fasta file
#[allow(missing_docs)]
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct FastaData {
    identifier: FastaIdentifier<Range<usize>>,
    description: Range<usize>,
    tags: Vec<(Range<usize>, Range<usize>)>,
    line_index: usize,
    full_header: String,
    peptide: Peptidoform<SemiAmbiguous>,
    regions: Vec<(Region, usize)>,
    annotations: Vec<(Annotation, usize)>,
}

impl Ord for FastaData {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.line_index.cmp(&other.line_index)
    }
}

impl PartialOrd for FastaData {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl HasPeptidoformImpl for FastaData {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
        self.peptide()
    }
}

impl HasPeptidoformImpl for &FastaData {
    type Complexity = SemiAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
        self.peptide()
    }
}

impl AnnotatedPeptide for FastaData {
    fn regions(&self) -> &[(Region, usize)] {
        &self.regions
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        &self.annotations
    }
}

/// A fasta identifier following the NCBI identifier definition
#[expect(missing_docs)]
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[expect(clippy::upper_case_acronyms)]
pub enum FastaIdentifier<T> {
    Undefined(T),
    Local(T),
    GenInfoBackboneSeqID(T),
    GenInfoBackboneMolType(T),
    GenInfoImportID(T),
    GenBank(T, T),
    EMBL(T, T),
    PIR(T, T),
    SwissProt(T, T),
    Patent(T, T, T),
    PrePatent(T, T, T),
    RefSeq(T, T),
    GeneralDatabase(T, T),
    GenInfoIntegratedDatabase(T),
    DDBJ(T, T),
    PRF(T, T),
    PDB(T, T),
    ThirdPartyGenBank(T, T),
    ThirdPartyEMBL(T, T),
    ThirdPartyDDJ(T, T),
    TrEMBL(T, T),
}

impl<T: Default> Default for FastaIdentifier<T> {
    fn default() -> Self {
        Self::Undefined(T::default())
    }
}

impl FastaIdentifier<Range<usize>> {
    fn as_str<'a>(&'a self, header: &'a str) -> FastaIdentifier<&'a str> {
        match self {
            Self::GenInfoBackboneSeqID(a) => {
                FastaIdentifier::GenInfoBackboneSeqID(&header[a.clone()])
            }
            Self::GenInfoBackboneMolType(a) => {
                FastaIdentifier::GenInfoBackboneMolType(&header[a.clone()])
            }
            Self::GenInfoImportID(a) => FastaIdentifier::GenInfoImportID(&header[a.clone()]),
            Self::GenInfoIntegratedDatabase(a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(&header[a.clone()])
            }
            Self::Undefined(a) => FastaIdentifier::Undefined(&header[a.clone()]),
            Self::Local(a) => FastaIdentifier::Local(&header[a.clone()]),
            Self::GenBank(a, b) => FastaIdentifier::GenBank(&header[a.clone()], &header[b.clone()]),
            Self::EMBL(a, b) => FastaIdentifier::EMBL(&header[a.clone()], &header[b.clone()]),
            Self::PIR(a, b) => FastaIdentifier::PIR(&header[a.clone()], &header[b.clone()]),
            Self::SwissProt(a, b) => {
                FastaIdentifier::SwissProt(&header[a.clone()], &header[b.clone()])
            }
            Self::RefSeq(a, b) => FastaIdentifier::RefSeq(&header[a.clone()], &header[b.clone()]),
            Self::GeneralDatabase(a, b) => {
                FastaIdentifier::GeneralDatabase(&header[a.clone()], &header[b.clone()])
            }
            Self::DDBJ(a, b) => FastaIdentifier::DDBJ(&header[a.clone()], &header[b.clone()]),
            Self::PRF(a, b) => FastaIdentifier::PRF(&header[a.clone()], &header[b.clone()]),
            Self::ThirdPartyGenBank(a, b) => {
                FastaIdentifier::ThirdPartyGenBank(&header[a.clone()], &header[b.clone()])
            }
            Self::ThirdPartyEMBL(a, b) => {
                FastaIdentifier::ThirdPartyEMBL(&header[a.clone()], &header[b.clone()])
            }
            Self::ThirdPartyDDJ(a, b) => {
                FastaIdentifier::ThirdPartyDDJ(&header[a.clone()], &header[b.clone()])
            }
            Self::TrEMBL(a, b) => FastaIdentifier::TrEMBL(&header[a.clone()], &header[b.clone()]),
            Self::PDB(a, b) => FastaIdentifier::PDB(&header[a.clone()], &header[b.clone()]),
            Self::Patent(a, b, c) => {
                FastaIdentifier::Patent(&header[a.clone()], &header[b.clone()], &header[c.clone()])
            }
            Self::PrePatent(a, b, c) => FastaIdentifier::PrePatent(
                &header[a.clone()],
                &header[b.clone()],
                &header[c.clone()],
            ),
        }
    }

    fn as_string(&self, header: &str) -> FastaIdentifier<String> {
        match self {
            Self::GenInfoBackboneSeqID(a) => {
                FastaIdentifier::GenInfoBackboneSeqID(header[a.clone()].to_string())
            }
            Self::GenInfoBackboneMolType(a) => {
                FastaIdentifier::GenInfoBackboneMolType(header[a.clone()].to_string())
            }
            Self::GenInfoImportID(a) => {
                FastaIdentifier::GenInfoImportID(header[a.clone()].to_string())
            }
            Self::GenInfoIntegratedDatabase(a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(header[a.clone()].to_string())
            }
            Self::Undefined(a) => FastaIdentifier::Undefined(
                header
                    .get(a.clone())
                    .map_or(String::new(), ToString::to_string),
            ),
            Self::Local(a) => FastaIdentifier::Local(header[a.clone()].to_string()),
            Self::GenBank(a, b) => FastaIdentifier::GenBank(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::EMBL(a, b) => {
                FastaIdentifier::EMBL(header[a.clone()].to_string(), header[b.clone()].to_string())
            }
            Self::PIR(a, b) => {
                FastaIdentifier::PIR(header[a.clone()].to_string(), header[b.clone()].to_string())
            }
            Self::SwissProt(a, b) => FastaIdentifier::SwissProt(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::RefSeq(a, b) => FastaIdentifier::RefSeq(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::GeneralDatabase(a, b) => FastaIdentifier::GeneralDatabase(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::DDBJ(a, b) => {
                FastaIdentifier::DDBJ(header[a.clone()].to_string(), header[b.clone()].to_string())
            }
            Self::PRF(a, b) => {
                FastaIdentifier::PRF(header[a.clone()].to_string(), header[b.clone()].to_string())
            }
            Self::ThirdPartyGenBank(a, b) => FastaIdentifier::ThirdPartyGenBank(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::ThirdPartyEMBL(a, b) => FastaIdentifier::ThirdPartyEMBL(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::ThirdPartyDDJ(a, b) => FastaIdentifier::ThirdPartyDDJ(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::TrEMBL(a, b) => FastaIdentifier::TrEMBL(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::PDB(a, b) => {
                FastaIdentifier::PDB(header[a.clone()].to_string(), header[b.clone()].to_string())
            }
            Self::Patent(a, b, c) => FastaIdentifier::Patent(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
                header[c.clone()].to_string(),
            ),
            Self::PrePatent(a, b, c) => FastaIdentifier::PrePatent(
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
                header[c.clone()].to_string(),
            ),
        }
    }
}

impl std::fmt::Display for FastaIdentifier<&'_ str> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GenInfoBackboneSeqID(a) => write!(f, "bbs|{a}"),
            Self::GenInfoBackboneMolType(a) => write!(f, "bbm|{a}"),
            Self::GenInfoImportID(a) => write!(f, "gim|{a}"),
            Self::GenInfoIntegratedDatabase(a) => write!(f, "gi|{a}"),
            Self::GenBank(a, b) => write!(f, "gb|{a}|{b}"),
            Self::EMBL(a, b) => write!(f, "emb|{a}|{b}"),
            Self::PIR(a, b) => write!(f, "pir|{a}|{b}"),
            Self::SwissProt(a, b) => write!(f, "sp|{a}|{b}"),
            Self::Patent(a, b, c) => write!(f, "pat|{a}|{b}|{c}"),
            Self::PrePatent(a, b, c) => write!(f, "pgp|{a}|{b}|{c}"),
            Self::RefSeq(a, b) => write!(f, "ref|{a}|{b}"),
            Self::GeneralDatabase(b, a) => write!(f, "gnl|{a}|{b}"),
            Self::DDBJ(a, b) => write!(f, "dbj|{a}|{b}"),
            Self::PRF(a, b) => write!(f, "prf|{a}|{b}"),
            Self::ThirdPartyGenBank(a, b) => write!(f, "tpg|{a}|{b}"),
            Self::ThirdPartyEMBL(a, b) => write!(f, "tpe|{a}|{b}"),
            Self::ThirdPartyDDJ(a, b) => write!(f, "tpd|{a}|{b}"),
            Self::TrEMBL(a, b) => write!(f, "tr|{a}|{b}"),
            Self::Undefined(a) => write!(f, "{a}"),
            Self::Local(a) => write!(f, "lcl|{a}"),
            Self::PDB(a, b) => write!(f, "pdb|{a}|{b}"),
        }
    }
}

impl std::fmt::Display for FastaIdentifier<String> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::GenInfoBackboneSeqID(a) => write!(f, "bbs|{a}"),
            Self::GenInfoBackboneMolType(a) => write!(f, "bbm|{a}"),
            Self::GenInfoImportID(a) => write!(f, "gim|{a}"),
            Self::GenInfoIntegratedDatabase(a) => write!(f, "gi|{a}"),
            Self::GenBank(a, b) => write!(f, "gb|{a}|{b}"),
            Self::EMBL(a, b) => write!(f, "emb|{a}|{b}"),
            Self::PIR(a, b) => write!(f, "pir|{a}|{b}"),
            Self::SwissProt(a, b) => write!(f, "sp|{a}|{b}"),
            Self::Patent(a, b, c) => write!(f, "pat|{a}|{b}|{c}"),
            Self::PrePatent(a, b, c) => write!(f, "pgp|{a}|{b}|{c}"),
            Self::RefSeq(a, b) => write!(f, "ref|{a}|{b}"),
            Self::GeneralDatabase(b, a) => write!(f, "gnl|{a}|{b}"),
            Self::DDBJ(a, b) => write!(f, "dbj|{a}|{b}"),
            Self::PRF(a, b) => write!(f, "prf|{a}|{b}"),
            Self::ThirdPartyGenBank(a, b) => write!(f, "tpg|{a}|{b}"),
            Self::ThirdPartyEMBL(a, b) => write!(f, "tpe|{a}|{b}"),
            Self::ThirdPartyDDJ(a, b) => write!(f, "tpd|{a}|{b}"),
            Self::TrEMBL(a, b) => write!(f, "tr|{a}|{b}"),
            Self::Undefined(a) => write!(f, "{a}"),
            Self::Local(a) => write!(f, "lcl|{a}"),
            Self::PDB(a, b) => write!(f, "pdb|{a}|{b}"),
        }
    }
}

impl<T> FastaIdentifier<T> {
    /// Get the accession or ID for this sequence
    pub const fn accession(&self) -> &T {
        match self {
            Self::GenInfoBackboneSeqID(a)
            | Self::GenInfoBackboneMolType(a)
            | Self::GenInfoImportID(a)
            | Self::GenInfoIntegratedDatabase(a)
            | Self::GenBank(a, _)
            | Self::EMBL(a, _)
            | Self::PIR(a, _)
            | Self::SwissProt(a, _)
            | Self::Patent(_, _, a)
            | Self::PrePatent(_, _, a)
            | Self::RefSeq(a, _)
            | Self::GeneralDatabase(_, a)
            | Self::DDBJ(a, _)
            | Self::PRF(a, _)
            | Self::ThirdPartyGenBank(a, _)
            | Self::ThirdPartyEMBL(a, _)
            | Self::ThirdPartyDDJ(a, _)
            | Self::TrEMBL(a, _)
            | Self::Undefined(a)
            | Self::Local(a)
            | Self::PDB(_, a) => a,
        }
    }

    /// Get the name, if no name is defined in this schema take the accession
    pub const fn name(&self) -> &T {
        match self {
            Self::GenInfoBackboneSeqID(n)
            | Self::GenInfoBackboneMolType(n)
            | Self::GenInfoImportID(n)
            | Self::GenInfoIntegratedDatabase(n)
            | Self::GenBank(n, _)
            | Self::EMBL(n, _)
            | Self::PIR(_, n)
            | Self::SwissProt(_, n)
            | Self::Patent(_, _, n)
            | Self::PrePatent(_, _, n)
            | Self::RefSeq(_, n)
            | Self::GeneralDatabase(_, n)
            | Self::DDBJ(n, _)
            | Self::PRF(_, n)
            | Self::ThirdPartyGenBank(_, n)
            | Self::ThirdPartyEMBL(_, n)
            | Self::ThirdPartyDDJ(_, n)
            | Self::TrEMBL(_, n)
            | Self::Undefined(n)
            | Self::Local(n)
            | Self::PDB(_, n) => n,
        }
    }
}

impl FromStr for FastaIdentifier<String> {
    type Err = ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(FastaIdentifier::<Range<usize>>::from_str(s)?.as_string(s))
    }
}

impl FromStr for FastaIdentifier<Range<usize>> {
    type Err = ParseIntError;
    /// Get the header string as ">header|stuff", so including the '>' until the first space
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let pipes = s
            .char_indices()
            .filter_map(|(i, c)| (c == '|').then_some(i))
            .collect_vec();
        let len = s.len();
        if pipes.is_empty() {
            Ok(Self::Undefined(usize::from(s.starts_with('>'))..len))
        } else {
            match s[usize::from(s.starts_with('>'))..pipes[0]]
                .to_ascii_lowercase()
                .as_str()
            {
                "lcl" => Ok(Self::Local(pipes[0] + 1..len)),
                "bbs" => Ok(Self::GenInfoBackboneSeqID(pipes[0] + 1..len)),
                "bbm" => Ok(Self::GenInfoBackboneMolType(pipes[0] + 1..len)),
                "gim" => Ok(Self::GenInfoImportID(pipes[0] + 1..len)),
                "gi" => Ok(Self::GenInfoIntegratedDatabase(pipes[0] + 1..len)),
                "gb" => Ok(Self::GenBank(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "emb" => Ok(Self::EMBL(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pir" => Ok(Self::PIR(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "sp" => Ok(Self::SwissProt(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "ref" => Ok(Self::RefSeq(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "gnl" => Ok(Self::GeneralDatabase(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "dbj" => Ok(Self::DDBJ(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "prf" => Ok(Self::PRF(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pdb" => Ok(Self::PDB(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpg" => Ok(Self::ThirdPartyGenBank(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpe" => Ok(Self::ThirdPartyEMBL(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpd" => Ok(Self::ThirdPartyDDJ(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tr" => Ok(Self::TrEMBL(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pat" => Ok(Self::Patent(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..pipes.get(2).map_or(len, |s| *s),
                    pipes.get(2).map_or(len, |s| *s + 1)..len,
                )),
                "pgp" => Ok(Self::PrePatent(
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..pipes.get(2).map_or(len, |s| *s),
                    pipes.get(2).map_or(len, |s| *s + 1)..len,
                )),
                _ => Ok(Self::Undefined(1..len)),
            }
        }
    }
}

impl FastaData {
    /// The identifier
    pub fn identifier(&self) -> FastaIdentifier<&str> {
        self.identifier.as_str(&self.full_header)
    }

    /// The description
    pub fn description(&self) -> &str {
        &self.full_header[self.description.clone()]
    }

    /// Get the tags, as key/value pairs, the keys are defined to be in uppercase
    pub fn tags(&self) -> impl DoubleEndedIterator<Item = (&str, &str)> + '_ {
        self.tags
            .iter()
            .map(|(k, v)| (&self.full_header[k.clone()], &self.full_header[v.clone()]))
    }

    /// Get the full header line
    pub fn header(&self) -> &str {
        &self.full_header
    }

    /// Get the sequence
    pub const fn peptide(&self) -> &Peptidoform<SemiAmbiguous> {
        &self.peptide
    }

    /// Parse a single fasta file
    /// # Errors
    /// A custom error when it is not a valid fasta file
    pub fn parse_file(path: impl AsRef<Path>) -> Result<Vec<Self>, CustomError> {
        let path = path.as_ref();
        let file = std::fs::File::open(path).map_err(|_| {
            CustomError::error(
                "Failed reading fasta file",
                "Error occurred while opening the file",
                Context::show(path.to_string_lossy()),
            )
        })?;
        let reader = BufReader::new(file);
        Self::parse_reader(reader, Some(path))
    }
    /// Parse a single fasta file from a reader
    /// # Errors
    /// A custom error when it is not a valid fasta file
    pub fn parse_reader(
        reader: impl BufRead,
        path: Option<&Path>,
    ) -> Result<Vec<Self>, CustomError> {
        let mut sequences = Vec::new();
        let mut last_header = None;
        let mut last_sequence: Vec<SequenceElement<SemiAmbiguous>> = Vec::new();

        for (line_index, line) in reader.lines().enumerate() {
            let line = line.map_err(|_| {
                CustomError::error(
                    "Failed reading fasta file",
                    format!("Error occurred while reading line {}", line_index + 1),
                    path.map_or(Context::None, |p| Context::show(p.to_string_lossy())),
                )
            })?;
            if line.starts_with('>') {
                if let Some(last_header) = last_header {
                    sequences.push(
                        Self {
                            peptide: last_sequence.into(),
                            ..last_header
                        }
                        .validate()?,
                    );
                }
                last_header = Some(Self::parse_header(line_index, line)?);
                last_sequence = Vec::new();
            } else {
                last_sequence.extend(
                    line.char_indices()
                        .filter(|(_, c)| !c.is_ascii_whitespace())
                        .map(|(i, c)| {
                            c.try_into()
                                .map(|aa: AminoAcid| SequenceElement::new(aa.into(), None))
                                .map_err(|()| {
                                    CustomError::error(
                                        "Failed reading fasta file",
                                        "Character is not an amino acid",
                                        Context::line(Some(line_index), &line, i, 1),
                                    )
                                })
                        })
                        .collect::<Result<Vec<SequenceElement<_>>, _>>()?,
                );
            }
        }
        if let Some(last_header) = last_header {
            sequences.push(
                Self {
                    peptide: last_sequence.into(),
                    ..last_header
                }
                .validate()?,
            );
        }

        Ok(sequences)
    }

    /// # Errors
    /// If the total length of the regions is not identical to the length of the peptide, or if any of the annotations is outside of the peptide
    fn validate(self) -> Result<Self, CustomError> {
        let total_regions_len: usize = self.regions.iter().map(|(_, l)| *l).sum();
        if total_regions_len > 0 && total_regions_len != self.peptide.len() {
            Err(CustomError::error(
                "Invalid regions definition",
                format!(
                    "The 'REGIONS' definition is invalid, the total length of the regions ({}) has to be identical to the length of the peptide ({})",
                    total_regions_len,
                    self.peptide.len()
                ),
                Context::full_line(self.line_index, &self.full_header),
            ))
        } else if self
            .annotations
            .iter()
            .any(|(_, p)| *p >= self.peptide.len())
        {
            Err(CustomError::error(
                "Invalid annotations definition",
                format!(
                    "The 'ANNOTATIONS' definition is invalid, on of the annotations is out of range of the peptide (length {})",
                    self.peptide.len()
                ),
                Context::full_line(self.line_index, &self.full_header),
            ))
        } else if total_regions_len > 0 {
            Ok(self)
        } else {
            // Add unannotated region annotation
            Ok(Self {
                regions: vec![(Region::None, self.peptide.len())],
                ..self
            })
        }
    }

    /// # Errors
    /// When the parsing of the fasta identifier is not succesful
    #[expect(clippy::missing_panics_doc)] // Regions and annotation parse cannot fail
    fn parse_header(line_index: usize, full_header: String) -> Result<Self, CustomError> {
        // thread 'main' panicked at C:\Users\5803969\src\rustyms\rustyms\src\identification\fasta.rs:301:26:
        // begin <= end (16 <= 15) when slicing `>Trastuzumab_HC REGIONS=FR1:25;CDR1:8;FR2:17;CDR2:8;FR3:38;CDR3:13;FR4:11;CH1:98;H:15;CH2:110;CH3:105;CHS:2 ANNOTATIONS=C:21;C:5;C:80;C:95;C:109;C:110;C:112;C:146;C:160;C:202;C:217;C:263;C:279;N:299;C:323;C:338;C:369;C:383;C:412;C:427;C:443`
        // note: run with `RUST_BACKTRACE=1` environment variable to display a backtrace
        // thread 'main' panicked at core\src\panicking.rs:221:5:
        // panic in a function that cannot unwind
        let first_space = full_header.find(' ').unwrap_or(full_header.len());
        let mut description = 0..0;
        let mut last_equals = None;
        let mut tags = Vec::new();
        let mut last_tag = None;

        loop {
            let start = last_equals.unwrap_or(first_space);
            let slice = &full_header[start..];
            if let Some(equals_position) = slice.find('=') {
                let tag_end = slice[..equals_position]
                    .char_indices()
                    .rev()
                    .take_while(|(_, c)| c.is_ascii_uppercase())
                    .last()
                    .map(|(i, _)| i)
                    .unwrap_or_default();
                if let Some(last_tag) = last_tag.take() {
                    tags.push((
                        last_tag,
                        trim_whitespace(&full_header, start..start + tag_end),
                    ));
                } else {
                    description = trim_whitespace(&full_header, start..start + tag_end);
                }
                last_tag = Some(start + tag_end..start + equals_position);
                last_equals = Some(start + equals_position + 1);
            } else {
                if let Some(last_tag) = last_tag.take() {
                    tags.push((
                        last_tag,
                        trim_whitespace(&full_header, start..full_header.len()),
                    ));
                } else {
                    description = trim_whitespace(&full_header, start..full_header.len());
                }
                break;
            }
        }

        let mut regions = Vec::new();
        let mut annotations = Vec::new();
        for tag in &tags {
            match &full_header[tag.0.clone()] {
                "REGIONS" => {
                    let mut index = 0;
                    regions =full_header[tag.1.clone()].split(';').map(|region| {
                        let last = index;
                        index += region.len() + usize::from(index != 0);
                    if let Some((region, n)) = region.split_once(':') {
                        Ok((
                            region.parse::<Region>().unwrap(),
                            n.parse::<usize>().map_err(|err| CustomError::error(
                            "Invalid regions definition", 
                            format!("The fasta header 'REGIONS' key, should contain regions followed by a colon, e.g. 'CDR3:6', but the number is {}", explain_number_error(&err)), 
                            Context::line(Some(line_index), &full_header, tag.1.start + last, tag.1.start+index)))?
                        ))
                    } else {
                        Err(CustomError::error(
                            "Invalid regions definition", 
                            "The fasta header 'REGIONS' key, should contain regions followed by a colon, e.g. 'CDR3:6'", 
                            Context::line(Some(line_index), &full_header, tag.1.start + last, tag.1.start+index)))
                    }
                }).collect::<Result<Vec<_>,_>>()?;
                }
                "ANNOTATIONS" => {
                    let mut index = 0;
                    annotations =full_header[tag.1.clone()].split(';').map(|region| {
                        let last = index;
                        index += region.len() + usize::from(index != 0);
                    if let Some((region, n)) = region.split_once(':') {
                        Ok((
                            region.parse::<Annotation>().unwrap(),
                            n.parse::<usize>().map_err(|err| CustomError::error(
                            "Invalid annotations definition", 
                            format!("The fasta header 'ANNOTATIONS' key, should contain annotations followed by a colon, e.g. 'Conserved:6', but the number is {}", explain_number_error(&err)), 
                            Context::line(Some(line_index), &full_header, tag.1.start + last, tag.1.start+index)))?
                        ))
                    } else {
                        Err(CustomError::error(
                            "Invalid annotations definition", 
                            "The fasta header 'ANNOTATIONS' key, should contain annotations followed by a colon, e.g. 'Conserved:6'", 
                            Context::line(Some(line_index), &full_header, tag.1.start + last, tag.1.start+index)))
                    }
                }).collect::<Result<Vec<_>,_>>()?;
                }
                _ => (),
            }
        }

        Ok(Self {
            identifier: full_header[0..first_space]
                .parse::<FastaIdentifier<Range<usize>>>()
                .map_err(|err| {
                    CustomError::error(
                        "Failed reading fasta file",
                        format!(
                            "Error occurred parsing NCBI identifier: number {}",
                            explain_number_error(&err)
                        ),
                        Context::line(Some(line_index), &full_header, 1, first_space - 1),
                    )
                })?,
            description,
            tags,
            regions,
            annotations,
            full_header,
            line_index,
            peptide: Peptidoform::default(),
        })
    }
}

fn trim_whitespace(line: &str, range: Range<usize>) -> Range<usize> {
    let start = range.len() - line[range.clone()].trim_start().len();
    let end = range.len() - line[range.clone()].trim_end().len();
    range.start + start..range.end - end
}

impl From<FastaData> for IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent> {
    fn from(value: FastaData) -> Self {
        Self {
            score: None,
            local_confidence: None,
            data: IdentifiedPeptidoformData::Fasta(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn empty_lines() {
    let file = ">A\naaa\n\naaa";
    let fasta = FastaData::parse_reader(BufReader::new(file.as_bytes()), None).unwrap();
    assert_eq!(fasta.len(), 1);
    assert_eq!(
        fasta[0].peptide,
        Peptidoform::pro_forma("AAAAAA", None).unwrap()
    );
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn parse_header() {
    let header = ">sp|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier PE=ProteinExistence SV=SequenceVersion REGIONS=FR1:12;CDR1:6;FR2:13 ANNOTATIONS=C:12;Conserved:25";
    let header = FastaData::parse_header(0, header.to_string()).unwrap();
    let identifier = header.identifier();
    assert_eq!(*identifier.name(), "EntryName");
    assert_eq!(*identifier.accession(), "UniqueIdentifier");
    assert_eq!(header.description(), "ProteinName");
    assert!(
        header
            .tags()
            .any(|(k, v)| k == "PE" && v == "ProteinExistence")
    );
    assert_eq!(header.regions().len(), 3);
    assert_eq!(header.regions()[0], (Region::Framework(1), 12));
    assert_eq!(header.annotations().len(), 2);
    assert_eq!(header.annotations()[0], (Annotation::Conserved, 12));
}

impl MetaData for FastaData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide().clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Fasta
    }

    fn id(&self) -> String {
        (*self.identifier().accession()).to_string()
    }

    fn confidence(&self) -> Option<f64> {
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        None
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        None
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        None
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::None
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        None
    }

    fn experimental_mass(&self) -> Option<Mass> {
        None
    }

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        None
    }

    fn protein_id(&self) -> Option<usize> {
        None
    }

    fn protein_location(&self) -> Option<Range<usize>> {
        None
    }
}
