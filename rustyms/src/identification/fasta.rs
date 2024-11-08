use crate::{
    error::{Context, CustomError},
    peptide::SemiAmbiguous,
    LinearPeptide, SequenceElement,
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::{
    io::{BufRead, BufReader},
    num::ParseIntError,
    ops::Range,
    path::Path,
    str::FromStr,
};

use super::{helper_functions::explain_number_error, AminoAcid, IdentifiedPeptide, MetaData};

/// A single parsed line of a fasta file
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize, Hash)]
pub struct FastaData {
    identifier: FastaIdentifier<Range<usize>>,
    description: Range<usize>,
    tags: Vec<(Range<usize>, Range<usize>)>,
    full_header: String,
    peptide: LinearPeptide<SemiAmbiguous>,
}

/// A fasta identifier following the NCBI identifier definition
#[allow(missing_docs)]
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
#[allow(clippy::upper_case_acronyms)]
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
}

impl<'a> std::fmt::Display for FastaIdentifier<&'a str> {
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

impl<T: Copy> FastaIdentifier<T> {
    /// Get the accession or ID for this sequence
    pub const fn accession(&self) -> T {
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
            | Self::PDB(_, a) => *a,
        }
    }

    /// Get the name, if no name is defined in this schema take the accession
    pub const fn name(&self) -> T {
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
            | Self::PDB(_, n) => *n,
        }
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
            Ok(Self::Undefined(1..len))
        } else {
            match s[1..pipes[0]].to_ascii_lowercase().as_str() {
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
    pub const fn peptide(&self) -> &LinearPeptide<SemiAmbiguous> {
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
            #[allow(clippy::manual_strip)]
            if line.starts_with('>') {
                if let Some(((identifier, description, tags), full_header)) = last_header {
                    sequences.push(Self {
                        identifier,
                        description,
                        tags,
                        full_header,
                        peptide: last_sequence.into(),
                    });
                }
                last_header = Some((
                    Self::parse_header(line_index, &line)?,
                    line[1..].to_string(),
                ));
                last_sequence = Vec::new();
            } else {
                last_sequence.extend(
                    line.char_indices()
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
        if let Some(((identifier, description, tags), full_header)) = last_header {
            sequences.push(Self {
                identifier,
                description,
                tags,
                full_header,
                peptide: last_sequence.into(),
            });
        }

        Ok(sequences)
    }

    /// # Errors
    /// When the parsing of the fasta identifier is not succesful
    fn parse_header(
        line_index: usize,
        header: &str,
    ) -> Result<
        (
            FastaIdentifier<Range<usize>>,
            Range<usize>,
            Vec<(Range<usize>, Range<usize>)>,
        ),
        CustomError,
    > {
        let first_space = header.find(' ').unwrap_or(header.len());
        let mut description = 0..0;
        let mut last_equals = None;
        let mut tags = Vec::new();
        let mut last_tag = None;

        loop {
            let start = last_equals.unwrap_or(first_space);
            let slice = &header[start..];
            if let Some(equals_position) = slice.find('=') {
                let tag_end = slice[..equals_position]
                    .char_indices()
                    .rev()
                    .take_while(|(_, c)| c.is_ascii_uppercase())
                    .last()
                    .map(|(i, _)| i)
                    .unwrap_or_default();
                if let Some(last_tag) = last_tag.take() {
                    tags.push((last_tag, trim_whitespace(header, start..start + tag_end)));
                } else {
                    description = trim_whitespace(header, start..start + tag_end);
                }
                last_tag = Some(start + tag_end..start + equals_position);
                last_equals = Some(start + equals_position + 1);
            } else {
                if let Some(last_tag) = last_tag.take() {
                    tags.push((last_tag, trim_whitespace(header, start..header.len())));
                } else {
                    description = trim_whitespace(header, start..header.len());
                }
                break;
            }
        }

        Ok((
            header[0..first_space]
                .parse::<FastaIdentifier<Range<usize>>>()
                .map_err(|err| {
                    CustomError::error(
                        "Failed reading fasta file",
                        format!(
                            "Error occurred parsing NCBI identifier: number {}",
                            explain_number_error(&err)
                        ),
                        Context::line(Some(line_index), header, 1, first_space - 1),
                    )
                })?,
            description,
            tags,
        ))
    }
}

fn trim_whitespace(line: &str, range: Range<usize>) -> Range<usize> {
    let start = range.len() - line[range.clone()].trim_ascii_start().len();
    let end = range.len() - line[range.clone()].trim_ascii_end().len();
    range.start + start..range.end - end
}

impl From<FastaData> for IdentifiedPeptide {
    fn from(value: FastaData) -> Self {
        Self {
            score: None,
            metadata: MetaData::Fasta(value),
        }
    }
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn empty_lines() {
    let file = ">A\naaa\n\naaa";
    let fasta = FastaData::parse_reader(BufReader::new(file.as_bytes()), None).unwrap();
    assert_eq!(fasta.len(), 1);
    assert_eq!(
        fasta[0].peptide,
        LinearPeptide::pro_forma("AAAAAA", None).unwrap()
    );
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn parse_header() {
    let header = ">sp|UniqueIdentifier|EntryName ProteinName OS=OrganismName OX=OrganismIdentifier PE=ProteinExistence SV=SequenceVersion";
    let (identifier, description, tags) = FastaData::parse_header(0, header).unwrap();
    let identifier = identifier.as_str(header);
    assert_eq!(identifier.name(), "EntryName");
    assert_eq!(identifier.accession(), "UniqueIdentifier");
    assert_eq!(&header[description], "ProteinName");
    assert!(tags
        .iter()
        .any(|(k, v)| &header[k.clone()] == "PE" && &header[v.clone()] == "ProteinExistence"));
}
