use std::{
    borrow::Cow,
    io::{BufRead, BufReader},
    marker::PhantomData,
    num::ParseIntError,
    ops::Range,
    path::Path,
    str::FromStr,
};

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    KnownFileFormat, PSM, PSMData, PSMMetaData, PeptidoformPresent, ProteinMetaData, SpectrumIds,
    helper_functions::explain_number_error,
};
use mzcore::{
    sequence::{
        AminoAcid, AnnotatedPeptidoform, Annotation, CompoundPeptidoformIon, FlankingSequence,
        HasPeptidoformImpl, Peptidoform, Region, SemiAmbiguous, SequenceElement,
    },
    system::{Mass, MassOverCharge, Time, isize::Charge},
};

/// A single parsed line of a fasta file
#[allow(missing_docs)]
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
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

impl AnnotatedPeptidoform for FastaData {
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
    Undefined(bool, T),
    Local(bool, T),
    GenInfoBackboneSeqID(bool, T),
    GenInfoBackboneMolType(bool, T),
    GenInfoImportID(bool, T),
    GenBank(bool, T, T),
    EMBL(bool, T, T),
    PIR(bool, T, T),
    SwissProt(bool, T, T),
    Patent(bool, T, T, T),
    PrePatent(bool, T, T, T),
    RefSeq(bool, T, T),
    GeneralDatabase(bool, T, T),
    GenInfoIntegratedDatabase(bool, T),
    DDBJ(bool, T, T),
    PRF(bool, T, T),
    PDB(bool, T, T),
    ThirdPartyGenBank(bool, T, T),
    ThirdPartyEMBL(bool, T, T),
    ThirdPartyDDJ(bool, T, T),
    TrEMBL(bool, T, T),
}

pub(crate) const EMPTY_FASTA_IDENTIFIER: FastaIdentifier<&str> =
    FastaIdentifier::Undefined(false, "");

impl<T: Default> Default for FastaIdentifier<T> {
    fn default() -> Self {
        Self::Undefined(false, T::default())
    }
}

impl FastaIdentifier<Range<usize>> {
    fn as_str<'a>(&'a self, header: &'a str) -> FastaIdentifier<&'a str> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, &header[a.clone()])
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, &header[a.clone()])
            }
            Self::GenInfoImportID(decoy, a) => {
                FastaIdentifier::GenInfoImportID(*decoy, &header[a.clone()])
            }
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, &header[a.clone()])
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(*decoy, &header[a.clone()]),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, &header[a.clone()]),
            Self::GenBank(decoy, a, b) => {
                FastaIdentifier::GenBank(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::EMBL(decoy, a, b) => {
                FastaIdentifier::EMBL(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::PIR(decoy, a, b) => {
                FastaIdentifier::PIR(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::SwissProt(decoy, a, b) => {
                FastaIdentifier::SwissProt(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::RefSeq(decoy, a, b) => {
                FastaIdentifier::RefSeq(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::GeneralDatabase(decoy, a, b) => {
                FastaIdentifier::GeneralDatabase(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::DDBJ(decoy, a, b) => {
                FastaIdentifier::DDBJ(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::PRF(decoy, a, b) => {
                FastaIdentifier::PRF(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::ThirdPartyGenBank(decoy, a, b) => {
                FastaIdentifier::ThirdPartyGenBank(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::ThirdPartyEMBL(decoy, a, b) => {
                FastaIdentifier::ThirdPartyEMBL(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::ThirdPartyDDJ(decoy, a, b) => {
                FastaIdentifier::ThirdPartyDDJ(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::TrEMBL(decoy, a, b) => {
                FastaIdentifier::TrEMBL(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::PDB(decoy, a, b) => {
                FastaIdentifier::PDB(*decoy, &header[a.clone()], &header[b.clone()])
            }
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(
                *decoy,
                &header[a.clone()],
                &header[b.clone()],
                &header[c.clone()],
            ),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(
                *decoy,
                &header[a.clone()],
                &header[b.clone()],
                &header[c.clone()],
            ),
        }
    }

    fn as_boxed_str(&self, header: &str) -> FastaIdentifier<Box<str>> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, header[a.clone()].into())
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, header[a.clone()].into())
            }
            Self::GenInfoImportID(decoy, a) => {
                FastaIdentifier::GenInfoImportID(*decoy, header[a.clone()].into())
            }
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, header[a.clone()].into())
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(
                *decoy,
                header
                    .get(a.clone())
                    .map_or(String::new(), ToString::to_string)
                    .into_boxed_str(),
            ),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, header[a.clone()].into()),
            Self::GenBank(decoy, a, b) => {
                FastaIdentifier::GenBank(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::EMBL(decoy, a, b) => {
                FastaIdentifier::EMBL(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::PIR(decoy, a, b) => {
                FastaIdentifier::PIR(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::SwissProt(decoy, a, b) => FastaIdentifier::SwissProt(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
            ),
            Self::RefSeq(decoy, a, b) => {
                FastaIdentifier::RefSeq(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::GeneralDatabase(decoy, a, b) => FastaIdentifier::GeneralDatabase(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
            ),
            Self::DDBJ(decoy, a, b) => {
                FastaIdentifier::DDBJ(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::PRF(decoy, a, b) => {
                FastaIdentifier::PRF(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::ThirdPartyGenBank(decoy, a, b) => FastaIdentifier::ThirdPartyGenBank(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
            ),
            Self::ThirdPartyEMBL(decoy, a, b) => FastaIdentifier::ThirdPartyEMBL(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
            ),
            Self::ThirdPartyDDJ(decoy, a, b) => FastaIdentifier::ThirdPartyDDJ(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
            ),
            Self::TrEMBL(decoy, a, b) => {
                FastaIdentifier::TrEMBL(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::PDB(decoy, a, b) => {
                FastaIdentifier::PDB(*decoy, header[a.clone()].into(), header[b.clone()].into())
            }
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
                header[c.clone()].into(),
            ),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(
                *decoy,
                header[a.clone()].into(),
                header[b.clone()].into(),
                header[c.clone()].into(),
            ),
        }
    }

    fn as_string(&self, header: &str) -> FastaIdentifier<String> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, header[a.clone()].to_string())
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, header[a.clone()].to_string())
            }
            Self::GenInfoImportID(decoy, a) => {
                FastaIdentifier::GenInfoImportID(*decoy, header[a.clone()].to_string())
            }
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, header[a.clone()].to_string())
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(
                *decoy,
                header
                    .get(a.clone())
                    .map_or(String::new(), ToString::to_string),
            ),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, header[a.clone()].to_string()),
            Self::GenBank(decoy, a, b) => FastaIdentifier::GenBank(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::EMBL(decoy, a, b) => FastaIdentifier::EMBL(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::PIR(decoy, a, b) => FastaIdentifier::PIR(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::SwissProt(decoy, a, b) => FastaIdentifier::SwissProt(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::RefSeq(decoy, a, b) => FastaIdentifier::RefSeq(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::GeneralDatabase(decoy, a, b) => FastaIdentifier::GeneralDatabase(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::DDBJ(decoy, a, b) => FastaIdentifier::DDBJ(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::PRF(decoy, a, b) => FastaIdentifier::PRF(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::ThirdPartyGenBank(decoy, a, b) => FastaIdentifier::ThirdPartyGenBank(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::ThirdPartyEMBL(decoy, a, b) => FastaIdentifier::ThirdPartyEMBL(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::ThirdPartyDDJ(decoy, a, b) => FastaIdentifier::ThirdPartyDDJ(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::TrEMBL(decoy, a, b) => FastaIdentifier::TrEMBL(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::PDB(decoy, a, b) => FastaIdentifier::PDB(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
            ),
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
                header[c.clone()].to_string(),
            ),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(
                *decoy,
                header[a.clone()].to_string(),
                header[b.clone()].to_string(),
                header[c.clone()].to_string(),
            ),
        }
    }
}

impl FastaIdentifier<String> {
    /// Create a string slice fasta identifier from an owned version
    pub fn as_str(&self) -> FastaIdentifier<&str> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, a)
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, a)
            }
            Self::GenInfoImportID(decoy, a) => FastaIdentifier::GenInfoImportID(*decoy, a),
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, a)
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(*decoy, a),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, a),
            Self::GenBank(decoy, a, b) => FastaIdentifier::GenBank(*decoy, a, b),
            Self::EMBL(decoy, a, b) => FastaIdentifier::EMBL(*decoy, a, b),
            Self::PIR(decoy, a, b) => FastaIdentifier::PIR(*decoy, a, b),
            Self::SwissProt(decoy, a, b) => FastaIdentifier::SwissProt(*decoy, a, b),
            Self::RefSeq(decoy, a, b) => FastaIdentifier::RefSeq(*decoy, a, b),
            Self::GeneralDatabase(decoy, a, b) => FastaIdentifier::GeneralDatabase(*decoy, a, b),
            Self::DDBJ(decoy, a, b) => FastaIdentifier::DDBJ(*decoy, a, b),
            Self::PRF(decoy, a, b) => FastaIdentifier::PRF(*decoy, a, b),
            Self::ThirdPartyGenBank(decoy, a, b) => {
                FastaIdentifier::ThirdPartyGenBank(*decoy, a, b)
            }
            Self::ThirdPartyEMBL(decoy, a, b) => FastaIdentifier::ThirdPartyEMBL(*decoy, a, b),
            Self::ThirdPartyDDJ(decoy, a, b) => FastaIdentifier::ThirdPartyDDJ(*decoy, a, b),
            Self::TrEMBL(decoy, a, b) => FastaIdentifier::TrEMBL(*decoy, a, b),
            Self::PDB(decoy, a, b) => FastaIdentifier::PDB(*decoy, a, b),
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(*decoy, a, b, c),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(*decoy, a, b, c),
        }
    }

    /// Create a string slice fasta identifier from an owned version
    pub fn as_boxed_str(&self) -> FastaIdentifier<Box<str>> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, a.clone().into_boxed_str())
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, a.clone().into_boxed_str())
            }
            Self::GenInfoImportID(decoy, a) => {
                FastaIdentifier::GenInfoImportID(*decoy, a.clone().into_boxed_str())
            }
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, a.clone().into_boxed_str())
            }
            Self::Undefined(decoy, a) => {
                FastaIdentifier::Undefined(*decoy, a.clone().into_boxed_str())
            }
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, a.clone().into_boxed_str()),
            Self::GenBank(decoy, a, b) => FastaIdentifier::GenBank(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::EMBL(decoy, a, b) => FastaIdentifier::EMBL(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::PIR(decoy, a, b) => FastaIdentifier::PIR(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::SwissProt(decoy, a, b) => FastaIdentifier::SwissProt(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::RefSeq(decoy, a, b) => FastaIdentifier::RefSeq(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::GeneralDatabase(decoy, a, b) => FastaIdentifier::GeneralDatabase(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::DDBJ(decoy, a, b) => FastaIdentifier::DDBJ(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::PRF(decoy, a, b) => FastaIdentifier::PRF(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::ThirdPartyGenBank(decoy, a, b) => FastaIdentifier::ThirdPartyGenBank(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::ThirdPartyEMBL(decoy, a, b) => FastaIdentifier::ThirdPartyEMBL(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::ThirdPartyDDJ(decoy, a, b) => FastaIdentifier::ThirdPartyDDJ(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::TrEMBL(decoy, a, b) => FastaIdentifier::TrEMBL(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::PDB(decoy, a, b) => FastaIdentifier::PDB(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
            ),
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
                c.clone().into_boxed_str(),
            ),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(
                *decoy,
                a.clone().into_boxed_str(),
                b.clone().into_boxed_str(),
                c.clone().into_boxed_str(),
            ),
        }
    }
}

impl FastaIdentifier<Box<str>> {
    /// Create a string slice fasta identifier from an owned version
    pub fn as_str(&self) -> FastaIdentifier<&str> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, a)
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, a)
            }
            Self::GenInfoImportID(decoy, a) => FastaIdentifier::GenInfoImportID(*decoy, a),
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, a)
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(*decoy, a),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, a),
            Self::GenBank(decoy, a, b) => FastaIdentifier::GenBank(*decoy, a, b),
            Self::EMBL(decoy, a, b) => FastaIdentifier::EMBL(*decoy, a, b),
            Self::PIR(decoy, a, b) => FastaIdentifier::PIR(*decoy, a, b),
            Self::SwissProt(decoy, a, b) => FastaIdentifier::SwissProt(*decoy, a, b),
            Self::RefSeq(decoy, a, b) => FastaIdentifier::RefSeq(*decoy, a, b),
            Self::GeneralDatabase(decoy, a, b) => FastaIdentifier::GeneralDatabase(*decoy, a, b),
            Self::DDBJ(decoy, a, b) => FastaIdentifier::DDBJ(*decoy, a, b),
            Self::PRF(decoy, a, b) => FastaIdentifier::PRF(*decoy, a, b),
            Self::ThirdPartyGenBank(decoy, a, b) => {
                FastaIdentifier::ThirdPartyGenBank(*decoy, a, b)
            }
            Self::ThirdPartyEMBL(decoy, a, b) => FastaIdentifier::ThirdPartyEMBL(*decoy, a, b),
            Self::ThirdPartyDDJ(decoy, a, b) => FastaIdentifier::ThirdPartyDDJ(*decoy, a, b),
            Self::TrEMBL(decoy, a, b) => FastaIdentifier::TrEMBL(*decoy, a, b),
            Self::PDB(decoy, a, b) => FastaIdentifier::PDB(*decoy, a, b),
            Self::Patent(decoy, a, b, c) => FastaIdentifier::Patent(*decoy, a, b, c),
            Self::PrePatent(decoy, a, b, c) => FastaIdentifier::PrePatent(*decoy, a, b, c),
        }
    }

    /// Create a string slice fasta identifier from an owned version
    pub fn as_string(&self) -> FastaIdentifier<String> {
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => {
                FastaIdentifier::GenInfoBackboneSeqID(*decoy, a.to_string())
            }
            Self::GenInfoBackboneMolType(decoy, a) => {
                FastaIdentifier::GenInfoBackboneMolType(*decoy, a.to_string())
            }
            Self::GenInfoImportID(decoy, a) => {
                FastaIdentifier::GenInfoImportID(*decoy, a.to_string())
            }
            Self::GenInfoIntegratedDatabase(decoy, a) => {
                FastaIdentifier::GenInfoIntegratedDatabase(*decoy, a.to_string())
            }
            Self::Undefined(decoy, a) => FastaIdentifier::Undefined(*decoy, a.to_string()),
            Self::Local(decoy, a) => FastaIdentifier::Local(*decoy, a.to_string()),
            Self::GenBank(decoy, a, b) => {
                FastaIdentifier::GenBank(*decoy, a.to_string(), b.to_string())
            }
            Self::EMBL(decoy, a, b) => FastaIdentifier::EMBL(*decoy, a.to_string(), b.to_string()),
            Self::PIR(decoy, a, b) => FastaIdentifier::PIR(*decoy, a.to_string(), b.to_string()),
            Self::SwissProt(decoy, a, b) => {
                FastaIdentifier::SwissProt(*decoy, a.to_string(), b.to_string())
            }
            Self::RefSeq(decoy, a, b) => {
                FastaIdentifier::RefSeq(*decoy, a.to_string(), b.to_string())
            }
            Self::GeneralDatabase(decoy, a, b) => {
                FastaIdentifier::GeneralDatabase(*decoy, a.to_string(), b.to_string())
            }
            Self::DDBJ(decoy, a, b) => FastaIdentifier::DDBJ(*decoy, a.to_string(), b.to_string()),
            Self::PRF(decoy, a, b) => FastaIdentifier::PRF(*decoy, a.to_string(), b.to_string()),
            Self::ThirdPartyGenBank(decoy, a, b) => {
                FastaIdentifier::ThirdPartyGenBank(*decoy, a.to_string(), b.to_string())
            }
            Self::ThirdPartyEMBL(decoy, a, b) => {
                FastaIdentifier::ThirdPartyEMBL(*decoy, a.to_string(), b.to_string())
            }
            Self::ThirdPartyDDJ(decoy, a, b) => {
                FastaIdentifier::ThirdPartyDDJ(*decoy, a.to_string(), b.to_string())
            }
            Self::TrEMBL(decoy, a, b) => {
                FastaIdentifier::TrEMBL(*decoy, a.to_string(), b.to_string())
            }
            Self::PDB(decoy, a, b) => FastaIdentifier::PDB(*decoy, a.to_string(), b.to_string()),
            Self::Patent(decoy, a, b, c) => {
                FastaIdentifier::Patent(*decoy, a.to_string(), b.to_string(), c.to_string())
            }
            Self::PrePatent(decoy, a, b, c) => {
                FastaIdentifier::PrePatent(*decoy, a.to_string(), b.to_string(), c.to_string())
            }
        }
    }
}

impl std::fmt::Display for FastaIdentifier<&'_ str> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = |decoy: bool| if decoy { "rev_" } else { "" };
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => write!(f, "{}bbs|{a}", prefix(*decoy)),
            Self::GenInfoBackboneMolType(decoy, a) => write!(f, "{}bbm|{a}", prefix(*decoy)),
            Self::GenInfoImportID(decoy, a) => write!(f, "{}gim|{a}", prefix(*decoy)),
            Self::GenInfoIntegratedDatabase(decoy, a) => write!(f, "{}gi|{a}", prefix(*decoy)),
            Self::GenBank(decoy, a, b) => write!(f, "{}gb|{a}|{b}", prefix(*decoy)),
            Self::EMBL(decoy, a, b) => write!(f, "{}emb|{a}|{b}", prefix(*decoy)),
            Self::PIR(decoy, a, b) => write!(f, "{}pir|{a}|{b}", prefix(*decoy)),
            Self::SwissProt(decoy, a, b) => write!(f, "{}sp|{a}|{b}", prefix(*decoy)),
            Self::Patent(decoy, a, b, c) => write!(f, "{}pat|{a}|{b}|{c}", prefix(*decoy)),
            Self::PrePatent(decoy, a, b, c) => write!(f, "{}pgp|{a}|{b}|{c}", prefix(*decoy)),
            Self::RefSeq(decoy, a, b) => write!(f, "{}ref|{a}|{b}", prefix(*decoy)),
            Self::GeneralDatabase(decoy, b, a) => write!(f, "{}gnl|{a}|{b}", prefix(*decoy)),
            Self::DDBJ(decoy, a, b) => write!(f, "{}dbj|{a}|{b}", prefix(*decoy)),
            Self::PRF(decoy, a, b) => write!(f, "{}prf|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyGenBank(decoy, a, b) => write!(f, "{}tpg|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyEMBL(decoy, a, b) => write!(f, "{}tpe|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyDDJ(decoy, a, b) => write!(f, "{}tpd|{a}|{b}", prefix(*decoy)),
            Self::TrEMBL(decoy, a, b) => write!(f, "{}tr|{a}|{b}", prefix(*decoy)),
            Self::Undefined(decoy, a) => write!(f, "{}{a}", prefix(*decoy)),
            Self::Local(decoy, a) => write!(f, "{}lcl|{a}", prefix(*decoy)),
            Self::PDB(decoy, a, b) => write!(f, "{}pdb|{a}|{b}", prefix(*decoy)),
        }
    }
}

impl std::fmt::Display for FastaIdentifier<Box<str>> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = |decoy: bool| if decoy { "rev_" } else { "" };
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => write!(f, "{}bbs|{a}", prefix(*decoy)),
            Self::GenInfoBackboneMolType(decoy, a) => write!(f, "{}bbm|{a}", prefix(*decoy)),
            Self::GenInfoImportID(decoy, a) => write!(f, "{}gim|{a}", prefix(*decoy)),
            Self::GenInfoIntegratedDatabase(decoy, a) => write!(f, "{}gi|{a}", prefix(*decoy)),
            Self::GenBank(decoy, a, b) => write!(f, "{}gb|{a}|{b}", prefix(*decoy)),
            Self::EMBL(decoy, a, b) => write!(f, "{}emb|{a}|{b}", prefix(*decoy)),
            Self::PIR(decoy, a, b) => write!(f, "{}pir|{a}|{b}", prefix(*decoy)),
            Self::SwissProt(decoy, a, b) => write!(f, "{}sp|{a}|{b}", prefix(*decoy)),
            Self::Patent(decoy, a, b, c) => write!(f, "{}pat|{a}|{b}|{c}", prefix(*decoy)),
            Self::PrePatent(decoy, a, b, c) => write!(f, "{}pgp|{a}|{b}|{c}", prefix(*decoy)),
            Self::RefSeq(decoy, a, b) => write!(f, "{}ref|{a}|{b}", prefix(*decoy)),
            Self::GeneralDatabase(decoy, b, a) => write!(f, "{}gnl|{a}|{b}", prefix(*decoy)),
            Self::DDBJ(decoy, a, b) => write!(f, "{}dbj|{a}|{b}", prefix(*decoy)),
            Self::PRF(decoy, a, b) => write!(f, "{}prf|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyGenBank(decoy, a, b) => write!(f, "{}tpg|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyEMBL(decoy, a, b) => write!(f, "{}tpe|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyDDJ(decoy, a, b) => write!(f, "{}tpd|{a}|{b}", prefix(*decoy)),
            Self::TrEMBL(decoy, a, b) => write!(f, "{}tr|{a}|{b}", prefix(*decoy)),
            Self::Undefined(decoy, a) => write!(f, "{}{a}", prefix(*decoy)),
            Self::Local(decoy, a) => write!(f, "{}lcl|{a}", prefix(*decoy)),
            Self::PDB(decoy, a, b) => write!(f, "{}pdb|{a}|{b}", prefix(*decoy)),
        }
    }
}

impl std::fmt::Display for FastaIdentifier<String> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let prefix = |decoy: bool| if decoy { "rev_" } else { "" };
        match self {
            Self::GenInfoBackboneSeqID(decoy, a) => write!(f, "{}bbs|{a}", prefix(*decoy)),
            Self::GenInfoBackboneMolType(decoy, a) => write!(f, "{}bbm|{a}", prefix(*decoy)),
            Self::GenInfoImportID(decoy, a) => write!(f, "{}gim|{a}", prefix(*decoy)),
            Self::GenInfoIntegratedDatabase(decoy, a) => write!(f, "{}gi|{a}", prefix(*decoy)),
            Self::GenBank(decoy, a, b) => write!(f, "{}gb|{a}|{b}", prefix(*decoy)),
            Self::EMBL(decoy, a, b) => write!(f, "{}emb|{a}|{b}", prefix(*decoy)),
            Self::PIR(decoy, a, b) => write!(f, "{}pir|{a}|{b}", prefix(*decoy)),
            Self::SwissProt(decoy, a, b) => write!(f, "{}sp|{a}|{b}", prefix(*decoy)),
            Self::Patent(decoy, a, b, c) => write!(f, "{}pat|{a}|{b}|{c}", prefix(*decoy)),
            Self::PrePatent(decoy, a, b, c) => write!(f, "{}pgp|{a}|{b}|{c}", prefix(*decoy)),
            Self::RefSeq(decoy, a, b) => write!(f, "{}ref|{a}|{b}", prefix(*decoy)),
            Self::GeneralDatabase(decoy, b, a) => write!(f, "{}gnl|{a}|{b}", prefix(*decoy)),
            Self::DDBJ(decoy, a, b) => write!(f, "{}dbj|{a}|{b}", prefix(*decoy)),
            Self::PRF(decoy, a, b) => write!(f, "{}prf|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyGenBank(decoy, a, b) => write!(f, "{}tpg|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyEMBL(decoy, a, b) => write!(f, "{}tpe|{a}|{b}", prefix(*decoy)),
            Self::ThirdPartyDDJ(decoy, a, b) => write!(f, "{}tpd|{a}|{b}", prefix(*decoy)),
            Self::TrEMBL(decoy, a, b) => write!(f, "{}tr|{a}|{b}", prefix(*decoy)),
            Self::Undefined(decoy, a) => write!(f, "{}{a}", prefix(*decoy)),
            Self::Local(decoy, a) => write!(f, "{}lcl|{a}", prefix(*decoy)),
            Self::PDB(decoy, a, b) => write!(f, "{}pdb|{a}|{b}", prefix(*decoy)),
        }
    }
}

impl<T> FastaIdentifier<T> {
    /// Get the accession or ID for this sequence
    pub const fn accession(&self) -> &T {
        match self {
            Self::GenInfoBackboneSeqID(_, a)
            | Self::GenInfoBackboneMolType(_, a)
            | Self::GenInfoImportID(_, a)
            | Self::GenInfoIntegratedDatabase(_, a)
            | Self::GenBank(_, a, _)
            | Self::EMBL(_, a, _)
            | Self::PIR(_, a, _)
            | Self::SwissProt(_, a, _)
            | Self::Patent(_, _, _, a)
            | Self::PrePatent(_, _, _, a)
            | Self::RefSeq(_, a, _)
            | Self::GeneralDatabase(_, _, a)
            | Self::DDBJ(_, a, _)
            | Self::PRF(_, a, _)
            | Self::ThirdPartyGenBank(_, a, _)
            | Self::ThirdPartyEMBL(_, a, _)
            | Self::ThirdPartyDDJ(_, a, _)
            | Self::TrEMBL(_, a, _)
            | Self::Undefined(_, a)
            | Self::Local(_, a)
            | Self::PDB(_, _, a) => a,
        }
    }

    /// Get the name, if no name is defined in this schema take the accession
    pub const fn name(&self) -> &T {
        match self {
            Self::GenInfoBackboneSeqID(_, n)
            | Self::GenInfoBackboneMolType(_, n)
            | Self::GenInfoImportID(_, n)
            | Self::GenInfoIntegratedDatabase(_, n)
            | Self::GenBank(_, n, _)
            | Self::EMBL(_, n, _)
            | Self::PIR(_, _, n)
            | Self::SwissProt(_, _, n)
            | Self::Patent(_, _, _, n)
            | Self::PrePatent(_, _, _, n)
            | Self::RefSeq(_, _, n)
            | Self::GeneralDatabase(_, _, n)
            | Self::DDBJ(_, n, _)
            | Self::PRF(_, _, n)
            | Self::ThirdPartyGenBank(_, _, n)
            | Self::ThirdPartyEMBL(_, _, n)
            | Self::ThirdPartyDDJ(_, _, n)
            | Self::TrEMBL(_, _, n)
            | Self::Undefined(_, n)
            | Self::Local(_, n)
            | Self::PDB(_, _, n) => n,
        }
    }

    /// See if this name is a decoy protein
    pub const fn decoy(&self) -> bool {
        match self {
            Self::GenInfoBackboneSeqID(decoy, _)
            | Self::GenInfoBackboneMolType(decoy, _)
            | Self::GenInfoImportID(decoy, _)
            | Self::GenInfoIntegratedDatabase(decoy, _)
            | Self::GenBank(decoy, _, _)
            | Self::EMBL(decoy, _, _)
            | Self::PIR(decoy, _, _)
            | Self::SwissProt(decoy, _, _)
            | Self::Patent(decoy, _, _, _)
            | Self::PrePatent(decoy, _, _, _)
            | Self::RefSeq(decoy, _, _)
            | Self::GeneralDatabase(decoy, _, _)
            | Self::DDBJ(decoy, _, _)
            | Self::PRF(decoy, _, _)
            | Self::ThirdPartyGenBank(decoy, _, _)
            | Self::ThirdPartyEMBL(decoy, _, _)
            | Self::ThirdPartyDDJ(decoy, _, _)
            | Self::TrEMBL(decoy, _, _)
            | Self::Undefined(decoy, _)
            | Self::Local(decoy, _)
            | Self::PDB(decoy, _, _) => *decoy,
        }
    }
}

impl FromStr for FastaIdentifier<String> {
    type Err = ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(FastaIdentifier::<Range<usize>>::from_str(s)?.as_string(s))
    }
}

impl FromStr for FastaIdentifier<Box<str>> {
    type Err = ParseIntError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(FastaIdentifier::<Range<usize>>::from_str(s)?.as_boxed_str(s))
    }
}

impl FromStr for FastaIdentifier<Range<usize>> {
    type Err = ParseIntError;
    /// Get the header string as ">header|stuff", so including the '>' until the first space
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut offset = usize::from(s.starts_with('>'));
        let decoy = s[offset..].len() > 4 && s[offset..offset + 4].eq_ignore_ascii_case("rev_");
        if decoy {
            offset += 4;
        }
        let pipes = s
            .char_indices()
            .filter_map(|(i, c)| (c == '|').then_some(i))
            .collect_vec();
        let len = s.len();
        if pipes.is_empty() {
            Ok(Self::Undefined(decoy, offset..len))
        } else {
            match s[offset..pipes[0]].to_ascii_lowercase().as_str() {
                "lcl" => Ok(Self::Local(decoy, pipes[0] + 1..len)),
                "bbs" => Ok(Self::GenInfoBackboneSeqID(decoy, pipes[0] + 1..len)),
                "bbm" => Ok(Self::GenInfoBackboneMolType(decoy, pipes[0] + 1..len)),
                "gim" => Ok(Self::GenInfoImportID(decoy, pipes[0] + 1..len)),
                "gi" => Ok(Self::GenInfoIntegratedDatabase(decoy, pipes[0] + 1..len)),
                "gb" => Ok(Self::GenBank(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "emb" => Ok(Self::EMBL(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pir" => Ok(Self::PIR(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "sp" => Ok(Self::SwissProt(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "ref" => Ok(Self::RefSeq(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "gnl" => Ok(Self::GeneralDatabase(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "dbj" => Ok(Self::DDBJ(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "prf" => Ok(Self::PRF(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pdb" => Ok(Self::PDB(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpg" => Ok(Self::ThirdPartyGenBank(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpe" => Ok(Self::ThirdPartyEMBL(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tpd" => Ok(Self::ThirdPartyDDJ(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "tr" => Ok(Self::TrEMBL(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..len,
                )),
                "pat" => Ok(Self::Patent(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..pipes.get(2).map_or(len, |s| *s),
                    pipes.get(2).map_or(len, |s| *s + 1)..len,
                )),
                "pgp" => Ok(Self::PrePatent(
                    decoy,
                    pipes[0] + 1..pipes.get(1).map_or(len, |s| *s),
                    pipes.get(1).map_or(len, |s| *s + 1)..pipes.get(2).map_or(len, |s| *s),
                    pipes.get(2).map_or(len, |s| *s + 1)..len,
                )),
                _ => Ok(Self::Undefined(decoy, 0..len)),
            }
        }
    }
}

impl<T: mzcore::space::Space> mzcore::space::Space for FastaIdentifier<T> {
    fn space(&self) -> mzcore::space::UsedSpace {
        match self {
            Self::Patent(decoy, a, b, c) | Self::PrePatent(decoy, a, b, c) => {
                decoy.space() + a.space() + b.space() + c.space()
            }
            Self::GenInfoBackboneSeqID(decoy, a)
            | Self::GenInfoBackboneMolType(decoy, a)
            | Self::GenInfoImportID(decoy, a)
            | Self::GenInfoIntegratedDatabase(decoy, a)
            | Self::Undefined(decoy, a)
            | Self::Local(decoy, a) => decoy.space() + a.space(),
            Self::GenBank(decoy, a, b)
            | Self::EMBL(decoy, a, b)
            | Self::PIR(decoy, a, b)
            | Self::SwissProt(decoy, a, b)
            | Self::RefSeq(decoy, a, b)
            | Self::GeneralDatabase(decoy, b, a)
            | Self::DDBJ(decoy, a, b)
            | Self::PRF(decoy, a, b)
            | Self::ThirdPartyGenBank(decoy, a, b)
            | Self::ThirdPartyEMBL(decoy, a, b)
            | Self::ThirdPartyDDJ(decoy, a, b)
            | Self::TrEMBL(decoy, a, b)
            | Self::PDB(decoy, a, b) => decoy.space() + a.space() + b.space(),
        }
        .set_total::<Self>()
    }
}

impl FastaData {
    /// The identifier
    pub fn identifier(&self) -> FastaIdentifier<&str> {
        self.identifier.as_str(&self.full_header[1..])
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
    pub fn parse_file(path: impl AsRef<Path>) -> Result<Vec<Self>, BoxedError<'static, BasicKind>> {
        let path = path.as_ref();
        let file = std::fs::File::open(path).map_err(|_| {
            BoxedError::new(
                BasicKind::Error,
                "Failed reading fasta file",
                "Error occurred while opening the file",
                Context::default().source(path.to_string_lossy()).to_owned(),
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
    ) -> Result<Vec<Self>, BoxedError<'static, BasicKind>> {
        let mut sequences = Vec::new();
        let mut last_header = None;
        let mut last_sequence: Vec<SequenceElement<SemiAmbiguous>> = Vec::new();

        for (line_index, line) in reader.lines().enumerate() {
            let line = line.map_err(|_| {
                BoxedError::new(
                    BasicKind::Error,
                    "Failed reading fasta file",
                    format!("Error occurred while reading line {}", line_index + 1),
                    path.map_or_else(Context::none, |p| {
                        Context::default().source(p.to_string_lossy()).to_owned()
                    }),
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
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Failed reading fasta file",
                                        "Character is not an amino acid",
                                        Context::line(Some(line_index as u32), &line, i, 1)
                                            .to_owned(),
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
    fn validate(self) -> Result<Self, BoxedError<'static, BasicKind>> {
        let total_regions_len: usize = self.regions.iter().map(|(_, l)| *l).sum();
        if total_regions_len > 0 && total_regions_len != self.peptide.len() {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid regions definition",
                format!(
                    "The 'REGIONS' definition is invalid, the total length of the regions ({}) has to be identical to the length of the peptide ({})",
                    total_regions_len,
                    self.peptide.len()
                ),
                Context::full_line(self.line_index as u32, &self.full_header).to_owned(),
            ))
        } else if self
            .annotations
            .iter()
            .any(|(_, p)| *p >= self.peptide.len())
        {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid annotations definition",
                format!(
                    "The 'ANNOTATIONS' definition is invalid, on of the annotations is out of range of the peptide (length {})",
                    self.peptide.len()
                ),
                Context::full_line(self.line_index as u32, &self.full_header).to_owned(),
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
    fn parse_header(
        line_index: usize,
        full_header: String,
    ) -> Result<Self, BoxedError<'static, BasicKind>> {
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
                            n.parse::<usize>().map_err(|err| BoxedError::new(BasicKind::Error,
                            "Invalid regions definition",
                            format!("The fasta header 'REGIONS' key, should contain regions followed by a colon, e.g. 'CDR3:6', but the number is {}", explain_number_error(&err)),
                            Context::line(Some(line_index as u32), &full_header, tag.1.start + last, tag.1.start+index).to_owned()))?
                        ))
                    } else {
                        Err(BoxedError::new(BasicKind::Error,
                            "Invalid regions definition",
                            "The fasta header 'REGIONS' key, should contain regions followed by a colon, e.g. 'CDR3:6'",
                            Context::line(Some(line_index as u32), &full_header, tag.1.start + last, tag.1.start+index).to_owned()))
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
                            n.parse::<usize>().map_err(|err| BoxedError::new(BasicKind::Error,
                            "Invalid annotations definition",
                            format!("The fasta header 'ANNOTATIONS' key, should contain annotations followed by a colon, e.g. 'Conserved:6', but the number is {}", explain_number_error(&err)),
                            Context::line(Some(line_index as u32), &full_header, tag.1.start + last, tag.1.start+index).to_owned()))?
                        ))
                    } else {
                        Err(BoxedError::new(BasicKind::Error,
                            "Invalid annotations definition",
                            "The fasta header 'ANNOTATIONS' key, should contain annotations followed by a colon, e.g. 'Conserved:6'",
                            Context::line(Some(line_index as u32), &full_header, tag.1.start + last, tag.1.start+index).to_owned()))
                    }
                }).collect::<Result<Vec<_>,_>>()?;
                }
                _ => (),
            }
        }

        Ok(Self {
            identifier: full_header[1..first_space]
                .parse::<FastaIdentifier<Range<usize>>>()
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Failed reading fasta file",
                        format!(
                            "Error occurred parsing NCBI identifier: number {}",
                            explain_number_error(&err)
                        ),
                        Context::line(Some(line_index as u32), &full_header, 1, first_space - 1)
                            .to_owned(),
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

impl From<FastaData> for PSM<SemiAmbiguous, PeptidoformPresent> {
    fn from(value: FastaData) -> Self {
        Self {
            score: None,
            local_confidence: None,
            data: PSMData::Fasta(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

impl From<FastaIdentifier<String>> for crate::ProteinData {
    fn from(value: FastaIdentifier<String>) -> Self {
        Self::FastaId(value.as_boxed_str())
    }
}

impl From<FastaIdentifier<Box<str>>> for crate::ProteinData {
    fn from(value: FastaIdentifier<Box<str>>) -> Self {
        Self::FastaId(value)
    }
}

impl From<FastaData> for crate::ProteinData {
    fn from(value: FastaData) -> Self {
        Self::Fasta(value)
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
        Peptidoform::pro_forma("AAAAAA", &mzcore::ontology::STATIC_ONTOLOGIES)
            .unwrap()
            .0
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

impl PSMMetaData for FastaData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide().clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::Fasta
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> String {
        (*self.identifier().accession()).to_string()
    }

    fn search_engine(&self) -> Option<mzcv::Term> {
        None
    }

    fn confidence(&self) -> Option<f64> {
        None
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<(f64, mzcv::Term)> {
        None
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        None
    }

    fn mode(&self) -> Option<Cow<'_, str>> {
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

    type Protein = Self;
    fn proteins(&self) -> Cow<'_, [Self::Protein]> {
        Cow::Borrowed(std::slice::from_ref(self))
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
}

impl ProteinMetaData for FastaData {
    fn sequence(&self) -> Option<Cow<'_, Peptidoform<mzcore::sequence::Linear>>> {
        Some(Cow::Borrowed(self.peptide.as_ref()))
    }

    fn numerical_id(&self) -> Option<usize> {
        None
    }

    fn id(&self) -> FastaIdentifier<&str> {
        self.identifier.as_str(&self.full_header)
    }

    fn description(&self) -> Option<&str> {
        Some(&self.full_header[self.description.clone()])
    }

    fn species(&self) -> Option<mzcv::Curie> {
        todo!()
    }

    fn species_name(&self) -> Option<&str> {
        todo!()
    }

    fn search_engine(&self) -> &[(super::CVTerm, Option<(f64, super::CVTerm)>)] {
        todo!()
    }

    fn ambiguity_members(&self) -> &[String] {
        todo!()
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        todo!()
    }

    fn modifications(
        &self,
    ) -> &[(
        Vec<(mzcore::prelude::SequencePosition, Option<f64>)>,
        mzcore::sequence::SimpleModification,
    )] {
        todo!()
    }

    fn coverage(&self) -> Option<f64> {
        todo!()
    }

    fn gene_ontology(&self) -> &[mzcv::Curie] {
        todo!()
    }

    fn reliability(&self) -> Option<super::Reliability> {
        todo!()
    }

    fn uri(&self) -> Option<&str> {
        todo!()
    }
}

impl mzcore::space::Space for FastaData {
    fn space(&self) -> mzcore::space::UsedSpace {
        (self.identifier.space()
            + self.description.space()
            + self.tags.space()
            + self.line_index.space()
            + self.full_header.space()
            + self.peptide.space()
            + self.regions.space()
            + self.annotations.space())
        .set_total::<Self>()
    }
}

#[allow(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use std::str::FromStr;

    use crate::FastaIdentifier;

    #[test]
    fn identifiers() {
        let a = "rev_sp|P0C1U8|SSPA_STAAU";
        let parsed = FastaIdentifier::<String>::from_str(a).unwrap();
        assert!(parsed.decoy());
        assert_eq!(parsed.accession(), "P0C1U8");
        assert_eq!(parsed.name(), "SSPA_STAAU");
        let a = "sp|P00778|PRLA_LYSEN";
        let parsed = FastaIdentifier::<String>::from_str(a).unwrap();
        assert!(!parsed.decoy());
        assert_eq!(parsed.accession(), "P00778");
        assert_eq!(parsed.name(), "PRLA_LYSEN");
    }
}
