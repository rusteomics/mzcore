//! Handle modification related issues, access provided if you want to dive deeply into
//! modifications in your own code.

use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fmt::{Display, Write},
    path::Path,
    str::FromStr,
    sync::Arc,
};

use context_error::*;
use itertools::Itertools;
use mzcv::{AccessionCode, SynonymScope};
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use thin_vec::ThinVec;

use crate::{
    chemistry::{AmbiguousLabel, Chemical, DiagnosticIon, MolecularFormula, NeutralLoss},
    glycan::{BackboneFragmentKind, GlycanAttachement, GlycanStructure, MonoSaccharide},
    helper_functions::{explain_number_error, merge_hashmap},
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    parse_json::{ParseJson, use_serde},
    quantities::Multi,
    sequence::{
        AminoAcid, CrossLinkName, CrossLinkSide, Linked, LinkerSpecificity, MUPSettings,
        Peptidoform, PlacementRule, Position, SequenceElement, SequencePosition,
        SimpleModification, SimpleModificationInner,
    },
    space::{Space, UsedSpace},
    system::OrderedMass,
};

/// A modification on an amino acid
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum Modification {
    /// Any of the simple modifications
    Simple(SimpleModification),
    /// A cross link to another (or the same) peptide, a branch is also seen as a cross-link but
    /// then the name is None.
    CrossLink {
        /// The index of the peptide this cross-link is bound to (can be the index for this peptide
        /// if it is an intra link)
        peptide: usize,
        /// The sequence index where this cross-link is bound to
        sequence_index: SequencePosition,
        /// The linker that defines the chemical structure that is the actual linker
        linker: SimpleModification,
        /// The name of the cross-linker, if [`CrossLinkName::Branch`] it is a branch instead of
        /// cross-link
        name: CrossLinkName,
        /// To determine if the cross-link is placed symmetrically or if asymmetrically if this is
        /// the left or right side
        side: CrossLinkSide,
    },
    /// An ambiguous modification, that can be placed at multiple locations
    Ambiguous {
        /// The name of the group
        group: String,
        /// The id to compare be able to find the other locations where this modifications can be
        /// placed
        id: usize,
        /// The modification itself
        modification: SimpleModification,
        /// If present the localisation score, meaning the chance/ratio for this modification to
        /// show up on this exact spot
        localisation_score: Option<OrderedFloat<f64>>,
        /// If this is the preferred location or not
        preferred: bool,
    },
}

impl Space for Modification {
    fn space(&self) -> UsedSpace {
        match self {
            Self::Simple(d) => d.space(),
            Self::CrossLink {
                peptide,
                sequence_index,
                linker,
                name,
                side,
            } => {
                peptide.space()
                    + sequence_index.space()
                    + linker.space()
                    + name.space()
                    + side.space()
            }
            Self::Ambiguous {
                group,
                id,
                modification,
                localisation_score,
                preferred,
            } => {
                group.space()
                    + id.space()
                    + modification.space()
                    + localisation_score.space()
                    + preferred.space()
            }
        }
        .set_total::<Self>()
    }
}

/// A modification id/name
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct ModificationId {
    /// The ontology where this linker is defined
    pub ontology: Ontology,
    /// The name
    pub name: Box<str>,
    /// The ID, `u32::MAX` is defined as not having an ID
    id: AccessionCode,
    /// The description or definition
    pub description: Box<str>,
    /// Any synonyms
    pub synonyms: ThinVec<(SynonymScope, Box<str>)>,
    /// Cross reference IDs
    pub cross_ids: ThinVec<CrossId>,
    /// Indicate if this modification is marked as obsolete
    pub obsolete: bool,
    /// Parent terms following 'is_a' relationships
    pub parents: ThinVec<AccessionCode>,
    /// Child terms following 'is_a' relationships
    pub children: ThinVec<AccessionCode>,
}

impl ModificationId {
    /// Create a new Modification ID
    pub fn new(
        ontology: Ontology,
        name: Box<str>,
        id: AccessionCode,
        description: Box<str>,
        synonyms: ThinVec<(SynonymScope, Box<str>)>,
        cross_ids: ThinVec<CrossId>,
        obsolete: bool,
    ) -> Self {
        Self {
            ontology,
            name,
            id,
            description,
            synonyms,
            cross_ids,
            obsolete,
            parents: ThinVec::default(),
            children: ThinVec::default(),
        }
    }

    /// Get the ID of this modification
    pub const fn id(&self) -> AccessionCode {
        self.id
    }

    /// Set the ID of this modification
    pub const fn set_id(&mut self, id: AccessionCode) {
        self.id = id;
    }
}

impl Space for ModificationId {
    fn space(&self) -> UsedSpace {
        (self.ontology.space()
            + self.name.space()
            + self.id.space()
            + self.description.space()
            + self.synonyms.space()
            + self.cross_ids.space()
            + self.obsolete.space())
        .set_total::<Self>()
    }
}

impl ParseJson for ModificationId {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::default().lines(
                    0,
                    map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                )
            };
            Ok(Self {
                ontology: Ontology::from_json_value(map.remove("ontology").ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ModificationID",
                        "The required property 'ontology' is missing",
                        context(&map),
                    )
                })?)?,
                name: Box::from_json_value(map.remove("name").ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ModificationID",
                        "The required property 'name' is missing",
                        context(&map),
                    )
                })?)?,
                id: match map.remove("id") {
                    None => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ModificationID",
                            "The required property 'id' is missing",
                            context(&map),
                        ));
                    }
                    Some(v) => AccessionCode::from_json_value(v)?,
                },
                description: Box::from_json_value(map.remove("description").ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ModificationID",
                        "The required property 'description' is missing",
                        context(&map),
                    )
                })?)?,
                synonyms: {
                    let element = map.remove("synonyms").ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid ModificationID",
                            "The required property 'synonyms' is missing",
                            context(&map),
                        )
                    })?;
                    if let Value::Array(arr) = element {
                        let mut output = ThinVec::with_capacity(arr.len());
                        for el in arr {
                            match el {
                                Value::String(str) => {
                                    output.push((SynonymScope::Exact, str.into_boxed_str()));
                                }
                                Value::Array(mut tup) if tup.len() == 2 => {
                                    let b = Box::from_json_value(tup.pop().unwrap())?;
                                    let a = SynonymScope::from_json_value(tup.pop().unwrap())?;
                                    output.push((a, b));
                                }
                                el => {
                                    return Err(BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid ModificationID",
                                        "A synonym could not be parsed, it has to be either a string or a tuple with two elements",
                                        Context::default().lines(0, el.to_string()),
                                    ));
                                }
                            }
                        }
                        output
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid ModificationID",
                            "Synonyms should be an array",
                            Context::default().lines(0, element.to_string()),
                        ));
                    }
                },
                cross_ids: ThinVec::from_json_value(map.remove("cross_ids").ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ModificationID",
                        "The required property 'cross_ids' is missing",
                        context(&map),
                    )
                })?)?,
                obsolete: map.remove("obsolete").map_or(Ok(false), |v| bool::from_json_value(v))?,
                parents: map
                    .remove("parents")
                    .map(|v| ThinVec::from_json_value(v))
                    .transpose()?
                    .unwrap_or_default(),
                children: map
                    .remove("children")
                    .map(|v| ThinVec::from_json_value(v))
                    .transpose()?
                    .unwrap_or_default(),
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid ModificationID",
                "The JSON value has to be a map",
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

/// A cross-identifier for a modification
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[allow(clippy::upper_case_acronyms, clippy::doc_markdown)]
pub enum CrossId {
    /// A reference to an article
    Article(Box<str>),
    /// A reference to the Beilstein organic chemistry database
    Beilstein(usize),
    /// A reference to a book
    Book(Box<str>),
    /// A CAS registry number
    CAS(u32, u8, u8),
    /// A reference to the chemical entities of biological interest database
    ChEBI(usize),
    /// A reference to the chemicalbook database
    ChemicalBook(Box<str>),
    /// A reference to the chemspider database
    ChemSpider(usize),
    /// Reference a COMe entry, a database of metalloproteins (10.1186/1472-6807-4-3)
    COMe(Box<str>),
    /// Enzyme Commission number
    EC(u16, u16, u16, u16, bool),
    /// A gene ontology term
    GO(AccessionCode),
    /// A deltamass ID
    Deltamass(usize),
    /// A digital object identifier
    DOI(Box<str>),
    /// A findmod ID
    Findmod(Box<str>),
    /// An MDL chemical database ID
    MDL(Box<str>),
    /// A modification in another modification ontology
    Mod(Ontology, AccessionCode, PlacementRule),
    /// A reference to the open mass spectrometry search algorithm
    OMSSA(usize),
    /// Any other reference, note that this can still contain a key in '{key}:{value}' format
    Other(Box<str>),
    /// A patent reference
    Patent(Box<str>),
    /// Reference a PDB entry
    PDB(Box<str>),
    /// Reference a PDBHet entry
    PDBHet(Box<str>),
    /// Reference a PubMed entry
    PubMed(usize),
    /// Reference a PubChem Compound entry
    PubChemCompound(usize),
    /// Reference a PubChem Substance entry
    PubChemSubstance(usize),
    /// Reference a uniprot PTMList entry by name
    Uniprot(Box<str>),
    /// Reference a URL with possibly a name for the link
    URL(Option<Box<str>>, Box<str>),
}

impl TryFrom<(Option<Box<str>>, Box<str>)> for CrossId {
    type Error = BoxedError<'static, BasicKind>;

    fn try_from(value: (Option<Box<str>>, Box<str>)) -> Result<Self, Self::Error> {
        let parse_ref = |r: &str| {
            if let Some((loc, rule)) = r.split_once('#') {
                let rule = rule.trim();
                loc.trim()
                    .parse()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid modification location",
                            format!("Contains invalid accession code {err:?}"),
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })
                    .and_then(|v| {
                        let r = if rule == "CYS2" || rule == "CYS1" {
                            PlacementRule::AminoAcid(
                                vec![AminoAcid::Cysteine].into(),
                                Position::Anywhere,
                            )
                        } else if rule == "N-term" {
                            PlacementRule::Position(Position::AnyNTerm)
                        } else if rule == "C-term" {
                            PlacementRule::Position(Position::AnyCTerm)
                        } else {
                            let loc = rule
                                .chars()
                                .map(TryInto::try_into)
                                .collect::<Result<ThinVec<AminoAcid>, _>>()
                                .map_err(|()| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid modification location",
                                        "Contains invalid amino acid",
                                        Context::default().lines(0, value.1.to_string()),
                                    )
                                })?;
                            PlacementRule::AminoAcid(loc, Position::Anywhere)
                        };

                        Ok((v, r))
                    })
            } else {
                r.trim()
                    .parse()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid modification location",
                            format!("Contains invalid accession code {err:?}"),
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })
                    .map(|v| (v, PlacementRule::Position(Position::Anywhere)))
            }
        };

        let t = value.0.map(|t| t.to_ascii_lowercase());
        Ok(match t.as_deref() {
            Some("article" | "journal") => Self::Article(value.1),
            Some("book") => Self::Book(value.1),
            Some("chemicalbookno") => Self::ChemicalBook(value.1),
            Some("come") => Self::COMe(value.1),
            Some("doi") => Self::DOI(value.1.trim_start_matches("https://doi.org/").into()),
            Some("epo" | "wipo" | "uspto") => Self::Patent(value.1),
            Some("findmod") => Self::Findmod(value.1),
            Some("http") => Self::URL(None, format!("http:{}", value.1).into_boxed_str()),
            Some("https") => Self::URL(None, format!("https:{}", value.1).into_boxed_str()),
            Some("mdl") => Self::MDL(value.1),
            Some("pdb") => Self::PDB(value.1),
            Some("pdbhet") => Self::PDBHet(value.1),
            Some("uniprot") => Self::Uniprot(value.1),
            Some("url" | "misc. url") => Self::URL(None, value.1),
            Some("cas" | "cas registry" | "cas registry number" | "cs") => {
                let mut split = value.1.trim().split('-');
                let a = split
                    .next()
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            "Missing first number",
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?
                    .parse()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            format!("The first number {}", explain_number_error(&err)),
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?;
                let b = split
                    .next()
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            "Missing second number",
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?
                    .parse()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            format!("The second number {}", explain_number_error(&err)),
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?;
                let c = split
                    .next()
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            "Missing third number",
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?
                    .parse()
                    .map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid CAS number",
                            format!("The third number {}", explain_number_error(&err)),
                            Context::default().lines(0, value.1.to_string()),
                        )
                    })?;
                if split.next().is_some() {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid CAS number",
                        "Too many numbers",
                        Context::default().lines(0, value.1.to_string()),
                    ));
                }

                // TODO: check the last number, this is a single digit that is a checksum of the
                // entire number
                Self::CAS(a, b, c)
            } // TODO: fix the one XLMOD that uses cs
            Some("chebi") => Self::ChEBI(value.1.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid ChEBI identifier",
                    format!("The number {}", explain_number_error(&err)),
                    Context::default().lines(0, value.1.to_string()),
                )
            })?),
            Some("deltamass") => Self::Deltamass(value.1.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid deltamass identifier",
                    format!("The number {}", explain_number_error(&err)),
                    Context::default().lines(0, value.1.to_string()),
                )
            })?),
            Some("omssa") => Self::OMSSA(value.1.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid OMSSA identifier",
                    format!("The number {}", explain_number_error(&err)),
                    Context::default().lines(0, value.1.to_string()),
                )
            })?),
            Some("beilstein") => Self::Beilstein(value.1.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid Beilstein identifier",
                    format!("The number {}", explain_number_error(&err)),
                    Context::default().lines(0, value.1.to_string()),
                )
            })?),
            Some("chemspider" | "chemspiderid" | "chemspider id") => {
                Self::ChemSpider(value.1.parse().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid ChemSpider identifier",
                        format!("The number {}", explain_number_error(&err)),
                        Context::default().lines(0, value.1.to_string()),
                    )
                })?)
            }
            Some("pubchemid" | "pubchem cid" | "pubchem" | "pubchem_compound") => {
                Self::PubChemCompound(value.1.parse().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid PubChem identifier",
                        format!("The number {}", explain_number_error(&err)),
                        Context::default().lines(0, value.1.to_string()),
                    )
                })?)
            }
            Some("pubchem_substance") => {
                Self::PubChemSubstance(value.1.parse().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid PubChem identifier",
                        format!("The number {}", explain_number_error(&err)),
                        Context::default().lines(0, value.1.to_string()),
                    )
                })?)
            }
            Some("pmid" | "pubmed" | "pubmed pmid") => {
                Self::PubMed(value.1.parse().map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid PubMed identifier",
                        format!("The number {}", explain_number_error(&err)),
                        Context::default().lines(0, value.1.to_string()),
                    )
                })?)
            }
            Some("go") => Self::GO(value.1.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid GO term",
                    format!("Contains invalid accession code {err:?}"),
                    Context::default().lines(0, value.1.to_string()),
                )
            })?),
            Some("unimod") => {
                let (id, rule) = parse_ref(&value.1)?;
                Self::Mod(Ontology::Unimod, id, rule)
            }
            Some("resid") => {
                let (id, rule) = parse_ref(&value.1)?;
                Self::Mod(Ontology::Resid, id, rule)
            }
            Some("psi-mod") => {
                let (id, rule) = parse_ref(&value.1)?;
                Self::Mod(Ontology::Psimod, id, rule)
            }
            Some("ec" | "brenda" | "expasy") => {
                let split = value.1.split('.').collect::<Vec<_>>();

                if split.len() == 4 {
                    let mut nums = Vec::with_capacity(4);
                    let mut preliminary = false;
                    let mut offset = 0;

                    for (n, element) in split.iter().enumerate() {
                        let num = if n == 3 && element.starts_with('n') {
                            preliminary = true;
                            element[1..].trim()
                        } else {
                            element.trim()
                        };
                        let a =
                            num.parse::<u16>().map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid EC number",
                                    format!("The EC number {}", explain_number_error(&err)),
                                    Context::default()
                                        .lines(0, value.1.to_string())
                                        .add_highlight((0, offset, element.len())),
                                )
                            })?;
                        nums.push(a);
                        offset += element.len() + 1;
                    }
                    Self::EC(nums[0], nums[1], nums[2], nums[3], preliminary)
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid EC number",
                        "An EC number consists of 4 parts separated by '.'",
                        Context::default().lines(0, value.1.to_string()),
                    ));
                }
            }
            Some(name) if value.1.starts_with("http://") | value.1.starts_with("https://") => {
                Self::URL(Some(name.into()), value.1)
            }
            Some(t) => {
                // println!("Unknown cross-id tag: {t}:{}", value.1);
                Self::Other(format!("{t}:{}", value.1).into_boxed_str())
            }
            None => Self::Other(value.1),
        })
    }
}

impl FromStr for CrossId {
    type Err = BoxedError<'static, BasicKind>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        s.split_once(':')
            .map_or_else(|| (None, s.into()), |(t, v)| (Some(t.into()), v.into()))
            .try_into()
    }
}

impl Display for CrossId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Article(a) => write!(f, "Article:{a}"),
            Self::Beilstein(a) => write!(f, "Beilstein:{a}"),
            Self::Book(a) => write!(f, "Book:{a}"),
            Self::CAS(a, b, c) => write!(f, "CAS:{a:02}-{b:02}-{c}"),
            Self::ChEBI(a) => write!(f, "CHEBI:{a}"),
            Self::ChemicalBook(a) => write!(f, "ChemicalBook:{a}"),
            Self::ChemSpider(a) => write!(f, "ChemSpider:{a}"),
            Self::COMe(a) => write!(f, "COMe:{a}"),
            Self::Deltamass(a) => write!(f, "DeltaMass:{a}"),
            Self::DOI(a) => write!(f, "DOI:{a}"),
            Self::EC(a, b, c, d, p) => write!(f, "EC:{a}.{b}.{c}.{}{d}", if *p { "n" } else { "" }),
            Self::Findmod(a) => write!(f, "FindMod:{a}"),
            Self::GO(a) => write!(f, "GO:{a}"),
            Self::MDL(a) => write!(f, "MDL:{a}"),
            Self::Mod(a, b, _) => write!(f, "{}:{b}", a.name()),
            Self::OMSSA(a) => write!(f, "OMSSA:{a}"),
            Self::Other(a) => write!(f, "{a}"),
            Self::Patent(a) => write!(f, "Patent:{a}"),
            Self::PDB(a) => write!(f, "PDB:{a}"),
            Self::PDBHet(a) => write!(f, "PDBHet:{a}"),
            Self::PubChemCompound(a) => write!(f, "PubChem_Compound:{a}"),
            Self::PubChemSubstance(a) => write!(f, "PubChem_Substance:{a}"),
            Self::PubMed(a) => write!(f, "PMID:{a}"),
            Self::Uniprot(a) => write!(f, "Uniprot:{a}"),
            Self::URL(name, b) => write!(f, "{}:{b}", name.as_deref().unwrap_or("URL")),
        }
    }
}

impl Space for CrossId {
    fn space(&self) -> UsedSpace {
        match self {
            Self::CAS(a, b, c) => a.space() + b.space() + c.space(),
            Self::EC(a, b, c, d, p) => a.space() + b.space() + c.space() + d.space() + p.space(),
            Self::GO(a) => a.space(),
            Self::Mod(a, b, c) => a.space() + b.space() + c.space(),
            Self::Beilstein(a)
            | Self::ChEBI(a)
            | Self::ChemSpider(a)
            | Self::Deltamass(a)
            | Self::OMSSA(a)
            | Self::PubChemCompound(a)
            | Self::PubChemSubstance(a)
            | Self::PubMed(a) => a.space(),
            Self::Article(a)
            | Self::Book(a)
            | Self::ChemicalBook(a)
            | Self::COMe(a)
            | Self::DOI(a)
            | Self::Findmod(a)
            | Self::MDL(a)
            | Self::Other(a)
            | Self::Patent(a)
            | Self::PDB(a)
            | Self::PDBHet(a)
            | Self::Uniprot(a) => a.space(),
            Self::URL(a, b) => a.space() + b.space(),
        }
        .set_total::<Self>()
    }
}

impl CrossId {
    /// Get a URL to resolve this cross-id
    pub fn url(&self) -> Option<String> {
        match self {
            Self::CAS(a, b, c) => Some(format!(
                "https://commonchemistry.cas.org/detail?cas_rn={a}-{b}-{c}"
            )),
            Self::ChEBI(a) => Some(format!("https://www.ebi.ac.uk/chebi/CHEBI:{a}")),
            Self::ChemicalBook(a) => Some(format!(
                "https://www.chemicalbook.com/ProductChemicalProperties{a}_EN.htm"
            )),
            Self::ChemSpider(a) => Some(format!(
                "https://www.chemspider.com/Chemical-Structure.{a}.html"
            )),
            Self::DOI(a) => Some(format!("https://doi.org/{a}")),
            Self::EC(a, b, c, d, _p) => Some(format!(
                "https://iubmb.qmul.ac.uk/enzyme/EC{a}/{b}/{c}/{d}.html"
            )), //TODO: not sure how the preliminary works here
            Self::Findmod(a) => Some(format!("https://web.expasy.org/findmod/{a}.html")),
            Self::GO(a) => Some(format!("https://amigo.geneontology.org/amigo/term/GO:{a}")), /* TODO: test if this works (has enough 0s) */
            Self::PDB(a) => Some(format!("https://www.rcsb.org/structure/{a}")),
            Self::PDBHet(a) => Some(format!("https://www.rcsb.org/ligand/{a}")),
            Self::PubChemCompound(a) => {
                Some(format!("https://pubchem.ncbi.nlm.nih.gov/compound/{a}"))
            }
            Self::PubChemSubstance(a) => {
                Some(format!("https://pubchem.ncbi.nlm.nih.gov/substance/{a}"))
            }
            Self::PubMed(a) => Some(format!("https://pubmed.ncbi.nlm.nih.gov/{a}/")),
            Self::URL(_, b) => Some(b.to_string()),
            Self::Article(_)
            | Self::Beilstein(_)
            | Self::Book(_)
            | Self::COMe(_)
            | Self::Deltamass(_)
            | Self::MDL(_)
            | Self::Mod(..)
            | Self::OMSSA(_)
            | Self::Other(_)
            | Self::Patent(_)
            | Self::Uniprot(_) => None,
        }
    }
}

impl ParseJson for CrossId {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        match value {
            value @ Value::Object(_) => use_serde(value),
            Value::Array(value) => {
                if value.len() == 2 {
                    Self::try_from((
                        (value[0] != Value::Null).then(|| value[0].to_string().into()),
                        value[1].to_string().into(),
                    ))
                } else {
                    todo!()
                }
            }
            Value::String(value) => Self::from_str(&value),
            value => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid CrossId",
                "A CrossId has to be an object",
                Context::default().lines(0, value.to_string()),
            )),
        }
    }
}

/// All possible compositions in the GNO ontology
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum GnoComposition {
    /// Only the mass is known
    Weight(OrderedMass),
    /// The composition,
    Composition(Vec<(MonoSaccharide, isize)>),
    /// The (full) structure is known
    Topology(GlycanStructure),
}

impl Space for GnoComposition {
    fn space(&self) -> UsedSpace {
        (UsedSpace::stack(1)
            + match self {
                Self::Weight(w) => w.space(),
                Self::Composition(c) => c.space(),
                Self::Topology(t) => t.space(),
            })
        .set_total::<Self>()
    }
}

/// All possible subsumption levels in the GNOme database indicating different levels of description
/// for a glycan species
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash, Serialize, Deserialize,
)]
pub enum GnoSubsumption {
    /// Indicates only the average weight is defined
    #[default]
    AverageWeight,
    /// Indicates the basic composition, without isomeric information
    BaseComposition,
    /// Indicates the composition, with isomeric information
    Composition,
    /// Indicates the topology, without linkage and anomeric information
    Topology,
    /// Indicates the topology, without reducing end ring and anomeric information
    Saccharide,
}

impl Display for GnoSubsumption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AverageWeight => write!(f, "Average weight"),
            Self::BaseComposition => write!(f, "Base composition (no isomeric information)"),
            Self::Composition => write!(f, "Composition"),
            Self::Topology => write!(f, "Topology (no linkage)"),
            Self::Saccharide => write!(f, "Saccharide"),
        }
    }
}

impl ModificationId {
    /// Get the accession number name for the ontology
    pub fn url(&self) -> Option<String> {
        match self.ontology {
            Ontology::Unimod => Some(format!(
                "https://www.unimod.org/modifications_view.php?editid1={}",
                self.id
            )),
            Ontology::Psimod => Some(format!(
                "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{}",
                self.id
            )),
            Ontology::Gnome => Some(format!(
                "http://glytoucan.org/Structures/Glycans/{}",
                self.name.to_ascii_uppercase()
            )),
            Ontology::Resid => Some(format!(
                "https://proteininformationresource.org/cgi-bin/resid?id={}",
                self.id
            )),
            Ontology::Xlmod | Ontology::Custom => None,
        }
    }
}

impl Display for ModificationId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.ontology == Ontology::Gnome {
            write!(
                f,
                "{}:{}",
                self.ontology.char(),
                self.name.to_ascii_uppercase()
            )
        } else if self.ontology == Ontology::Xlmod
            && (self.id() == AccessionCode::Numeric(1711)
                || self.id() == AccessionCode::Numeric(9070))
        {
            write!(f, "{}:{}", self.ontology.name(), self.id())
        } else {
            write!(f, "{}:{}", self.ontology.char(), self.name)
        }
    }
}

/// The result of checking if a modification can be placed somewhere.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum RulePossible {
    /// This modification cannot be placed
    No,
    /// This modification can be placed and if it is a cross-link it can be placed on both ends
    Symmetric(BTreeSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the
    /// 'left' side of the cross-link
    AsymmetricLeft(BTreeSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the
    /// 'right' side of the cross-link
    AsymmetricRight(BTreeSet<usize>),
}

impl RulePossible {
    /// Flatten this into a bool, check if the rule is not [`Self::No`]
    pub fn any_possible(self) -> bool {
        self != Self::No
    }
}

impl std::ops::Add for RulePossible {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Symmetric(a), _) | (_, Self::Symmetric(a)) => Self::Symmetric(a),
            (Self::AsymmetricLeft(l), Self::AsymmetricRight(r))
            | (Self::AsymmetricRight(l), Self::AsymmetricLeft(r)) => {
                let overlap: BTreeSet<usize> = l.intersection(&r).copied().collect();
                if overlap.is_empty() {
                    Self::No
                } else {
                    Self::Symmetric(overlap)
                }
            }
            (Self::AsymmetricLeft(l), _) | (_, Self::AsymmetricLeft(l)) => Self::AsymmetricLeft(l),
            (Self::AsymmetricRight(r), _) | (_, Self::AsymmetricRight(r)) => {
                Self::AsymmetricRight(r)
            }
            _ => Self::No,
        }
    }
}

impl std::iter::Sum for RulePossible {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::No, |acc, i| acc + i)
    }
}

impl CrossLinkSide {
    /// Get all allowed placement rules with all applicable neutral losses, stubs, and diagnostic
    /// ions.
    pub(crate) fn allowed_rules(
        &self,
        linker: &SimpleModification,
    ) -> (
        Vec<NeutralLoss>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ) {
        let selected_rules = match self {
            Self::Left(r) | Self::Right(r) | Self::Symmetric(r) => r,
        };
        let mut all_neutral = Vec::new();
        let mut all_stubs = Vec::new();
        let mut all_diagnostic = Vec::new();

        match &**linker {
            SimpleModificationInner::Linker { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    match rule {
                        LinkerSpecificity::Asymmetric {
                            diagnostic,
                            stubs,
                            neutral_losses,
                            ..
                        } => {
                            all_neutral.extend_from_slice(neutral_losses);
                            all_diagnostic.extend_from_slice(diagnostic);
                            match self {
                                Self::Left(_) => all_stubs.extend(stubs.iter().cloned()),
                                Self::Right(_) => {
                                    all_stubs
                                        .extend(stubs.iter().map(|(l, r)| (r.clone(), l.clone())));
                                }
                                Self::Symmetric(_) => {
                                    all_stubs.extend(stubs.iter().flat_map(|(l, r)| {
                                        vec![(l.clone(), r.clone()), (r.clone(), l.clone())]
                                    }));
                                }
                            }
                        }
                        LinkerSpecificity::Symmetric {
                            stubs,
                            neutral_losses,
                            diagnostic,
                            ..
                        } => {
                            all_stubs.extend_from_slice(stubs);
                            all_neutral.extend_from_slice(neutral_losses);
                            all_diagnostic.extend_from_slice(diagnostic);
                        }
                    }
                }
            }
            SimpleModificationInner::Database { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    all_neutral.extend_from_slice(&rule.1);
                    all_diagnostic.extend_from_slice(&rule.2);
                }
            }
            _ => (),
        }
        (all_neutral, all_stubs, all_diagnostic)
    }
}

impl Modification {
    /// Check if this modification is a simple modification.
    pub const fn is_simple(&self) -> bool {
        matches!(self, Self::Simple(_))
    }

    /// Check if this modification is a cross-link.
    pub const fn is_cross_link(&self) -> bool {
        matches!(self, Self::CrossLink { .. })
    }

    /// Check if this modification is an ambiguous modification.
    pub const fn is_ambiguous(&self) -> bool {
        matches!(self, Self::Ambiguous { .. })
    }

    /// Get the formula for the whole addition (or subtraction) for this modification
    #[doc(hidden)]
    pub fn formula_inner(
        &self,
        all_peptides: &[Peptidoform<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
        glycan_model: &impl GlycanAttachement,
        attachment: Option<AminoAcid>,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        match self {
            Self::Simple(modification) | Self::Ambiguous { modification, .. } => {
                match &**modification {
                    // A linker that is not cross-linked is hydrolysed
                    SimpleModificationInner::Linker { formula, .. } => (
                        (formula.clone() + molecular_formula!(H 2 O 1)).into(),
                        HashMap::new(),
                        HashSet::new(),
                    ),

                    s => {
                        let default_rules = glycan_model.get_default_fragments(attachment);
                        let specific_rules = glycan_model.get_specific_fragments(attachment);

                        let f = s.formula_inner(
                            sequence_index,
                            peptidoform_index,
                            default_rules,
                            attachment,
                        );
                        let specific = specific_rules
                            .into_iter()
                            .map(|(k, settings)| {
                                (
                                    k,
                                    s.formula_inner(
                                        sequence_index,
                                        peptidoform_index,
                                        settings,
                                        attachment,
                                    ),
                                )
                            })
                            .collect();
                        (f, specific, HashSet::new())
                    }
                }
            }
            Self::CrossLink {
                peptide: other_peptide,
                linker,
                name,
                side,
                ..
            } => {
                if applied_cross_links.contains(name) {
                    (Multi::default(), HashMap::new(), HashSet::default())
                } else if visited_peptides.contains(other_peptide) {
                    applied_cross_links.push(name.clone());
                    (
                        linker
                            .formula_inner(
                                sequence_index,
                                peptidoform_index,
                                glycan_model.get_default_fragments(attachment),
                                attachment,
                            )
                            .with_label(&AmbiguousLabel::CrossLinkBound(name.clone())),
                        HashMap::new(),
                        HashSet::from([name.clone()]),
                    )
                } else {
                    applied_cross_links.push(name.clone());
                    let link = linker.formula_inner(
                        sequence_index,
                        peptidoform_index,
                        glycan_model.get_default_fragments(attachment),
                        attachment,
                    );
                    let (_, stubs, _) = side.allowed_rules(linker);

                    if allow_ms_cleavable && !stubs.is_empty() {
                        let mut options: Vec<MolecularFormula> = stubs
                            .iter()
                            .map(|s| {
                                s.0.clone().with_label(AmbiguousLabel::CrossLinkBroken(
                                    name.clone(),
                                    s.0.clone(),
                                ))
                            })
                            .unique()
                            .collect();
                        let mut seen_peptides = HashSet::from([name.clone()]);
                        let mut specific = HashMap::new();

                        options.extend_from_slice(&{
                            let (f, f_specific, seen) = all_peptides[*other_peptide]
                                .formulas_inner(
                                    *other_peptide,
                                    peptidoform_ion_index,
                                    all_peptides,
                                    visited_peptides,
                                    applied_cross_links,
                                    false,
                                    glycan_model,
                                );
                            seen_peptides.extend(seen);
                            specific = merge_hashmap(specific, &f_specific, &link, &f);
                            (f * link)
                                .with_label(&AmbiguousLabel::CrossLinkBound(name.clone()))
                                .to_vec()
                        });

                        (options.into(), specific, seen_peptides)
                    } else {
                        let (f, specific, mut seen) = all_peptides[*other_peptide].formulas_inner(
                            *other_peptide,
                            peptidoform_ion_index,
                            all_peptides,
                            visited_peptides,
                            applied_cross_links,
                            false,
                            glycan_model,
                        );
                        seen.insert(name.clone());
                        (
                            (f * link).with_label(&AmbiguousLabel::CrossLinkBound(name.clone())),
                            specific,
                            seen,
                        )
                    }
                }
            }
        }
    }

    /// Get the formula for a modification, if it is a cross linked modification only get the cross
    /// link
    pub fn formula(&self) -> MolecularFormula {
        match self {
            Self::Simple(s) => s.formula(),
            Self::CrossLink { linker, .. } => linker.formula(),
            Self::Ambiguous { modification, .. } => modification.formula(),
        }
    }

    /// Check if this is a simple modification
    pub const fn simple(&self) -> Option<&SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } | Self::Ambiguous { .. } => None,
        }
    }

    /// Check if this is a simple modification
    pub fn into_simple(self) -> Option<SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } | Self::Ambiguous { .. } => None,
        }
    }

    /// Get the underlying simple mod, even if this is a cross-link or ambiguous
    pub fn get_simple(&self) -> SimpleModification {
        match self {
            Self::Simple(modification)
            | Self::CrossLink {
                linker: modification,
                ..
            }
            | Self::Ambiguous { modification, .. } => modification.clone(),
        }
    }

    /// Get a url for more information on this modification. Only defined for modifications from
    /// ontologies.
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Simple(modification)
            | Self::Ambiguous { modification, .. }
            | Self::CrossLink {
                linker: modification,
                ..
            } => modification.ontology_url(),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible<T>(
        &self,
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> RulePossible {
        self.simple().map_or(RulePossible::Symmetric(BTreeSet::new()), |s| {
            s.is_possible(seq, position)
        })
    }

    /// Get the name if this is an Unimod modification (for use in mzPAF)
    pub fn unimod_name(&self) -> Option<&str> {
        match self {
            Self::Simple(s)
            | Self::CrossLink { linker: s, .. }
            | Self::Ambiguous {
                modification: s, ..
            } => match &**s {
                SimpleModificationInner::Database {
                    id:
                        ModificationId {
                            ontology: Ontology::Unimod,
                            name,
                            ..
                        },
                    ..
                } => Some(name),
                _ => None,
            },
        }
    }
}

/// The structure to lookup ambiguous modifications, with a list of all modifications (the order is
/// fixed) with for each modification their name and the actual modification itself (if already
/// defined)
pub type AmbiguousLookup = Vec<AmbiguousLookupEntry>;

/// An entry in the ambiguous lookup
#[derive(Clone, Debug)]
pub struct AmbiguousLookupEntry {
    /// The name of the modification
    pub name: String,
    /// The group of the modification
    pub group: Option<usize>,
    /// The modification itself
    pub modification: Option<SimpleModification>,
    /// The allowed locations, the actual allowed locations is the intersection of this set with
    /// the ruleset from the modification
    position: Option<Vec<PlacementRule>>,
    /// The maximal number of this modification on one place
    limit: Option<usize>,
    /// Determines if this modification can colocalise with placed modifications eg if the
    /// modification of unknown position is allowed at the second M '[Oxidation]?MM[Dioxidation]M'
    comkp: bool,
    /// Determines if this modification can colocalise with other modifications of unknown position
    comup: bool,
}

impl AmbiguousLookupEntry {
    /// Create a new ambiguous lookup entry
    pub const fn new(name: String, modification: Option<SimpleModification>) -> Self {
        Self {
            name,
            group: None,
            modification,
            limit: None,
            position: None,
            comkp: true,
            comup: true,
        }
    }

    /// Copy settings into this lookup entry
    pub fn copy_settings(&mut self, settings: &MUPSettings) {
        self.position.clone_from(&settings.position);
        self.limit = settings.limit;
        self.comkp = settings.comkp;
        self.comup = settings.comup;
    }

    /// Get the settings for this modification of unknown position
    pub fn as_settings(&self) -> MUPSettings {
        MUPSettings {
            position: self.position.clone(),
            limit: self.limit,
            comkp: self.comkp,
            comup: self.comup,
        }
    }
}

/// The structure to lookup cross-links, with a list of all linkers (the order is fixed) with for
/// each linker their name or None if it is a branch and the actual linker itself (if already
/// defined)
pub type CrossLinkLookup = Vec<(CrossLinkName, Option<SimpleModification>)>;

impl Modification {
    /// Display a modification either normalised to the internal representation or as fully valid
    /// ProForma (no glycan structure or custom modifications). `display_ambiguous` shows or
    /// hides the modification definition of any ambiguous modifications (eg true results in
    /// '+1#1' false in '#1').
    /// # Errors
    /// When the given writer errors.
    pub fn display(
        &self,
        f: &mut impl Write,
        specification_compliant: bool,
        display_ambiguous: bool,
    ) -> std::fmt::Result {
        match self {
            Self::Simple(sim) => sim.display(f, specification_compliant),
            Self::CrossLink { name, linker, .. } => {
                linker.display(f, specification_compliant)?;
                write!(f, "{name}")?;
                Ok(())
            }
            Self::Ambiguous {
                group,
                modification,
                localisation_score,
                ..
            } => {
                if display_ambiguous {
                    modification.display(f, specification_compliant)?;
                }
                write!(
                    f,
                    "\x23{group}{}",
                    localisation_score.map(|v| format!("({v})")).unwrap_or_default()
                )?;
                Ok(())
            }
        }
    }
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true, true)
    }
}

impl From<SimpleModification> for Modification {
    fn from(value: SimpleModification) -> Self {
        Self::Simple(value)
    }
}

impl From<SimpleModificationInner> for Modification {
    fn from(value: SimpleModificationInner) -> Self {
        Self::Simple(Arc::new(value))
    }
}

impl Display for CrossLinkName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Branch => write!(f, "#BRANCH"),
            Self::Name(n) => write!(f, "#XL{n}"),
        }
    }
}

/// Parse a custom modifications JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom database in previous version of
/// the library.
/// # Errors
/// If the string could not be parsed.
pub fn parse_custom_modifications(
    path: &Path,
) -> Result<CustomDatabase, BoxedError<'static, BasicKind>> {
    let string = std::fs::read_to_string(path).map_err(|err| {
        BoxedError::new(
            BasicKind::Error,
            "Could not parse custom modifications file",
            err.to_string(),
            Context::default().source(path.to_string_lossy()).to_owned(),
        )
    })?;
    CustomDatabase::from_json(&string)
}

/// Parse a custom modifications JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom database in previous version of
/// the library.
/// # Errors
/// If the string could not be parsed.
pub fn parse_custom_modifications_str(
    value: &str,
) -> Result<CustomDatabase, BoxedError<'static, BasicKind>> {
    CustomDatabase::from_json(value)
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use super::*;
    use crate::ontology::STATIC_ONTOLOGIES;

    #[test]
    fn test_reading_custom_modifications_json_2024() {
        let data = include_str!("custom_modifications_2024.json");
        let mods = parse_custom_modifications_str(data).unwrap();
        assert!(mods.len() > 1);
    }

    #[test]
    fn test_reading_custom_modifications_json_2025() {
        let data = include_str!("custom_modifications_20250207.json");
        let mods = parse_custom_modifications_str(data).unwrap();
        assert!(mods.len() > 1);
    }

    #[test]
    fn parse_mod_json() {
        let modification: SimpleModificationInner = dbg!(SimpleModificationInner::from_json_value(serde_json::from_str(
            r#"{"Database":{"specificities":[[[{"AminoAcid":[["Arginine","Asparagine","AsparticAcid","Cysteine","Glutamine","GlutamicAcid","Histidine","Isoleucine","Leucine","Lysine","Phenylalanine","Proline","Serine","Threonine","Tryptophan","Tyrosine","Valine","Selenocysteine"],"Anywhere"]},{"AminoAcid":[["Glycine"],"AnyCTerm"]}],[],[]],[[{"AminoAcid":[["Methionine"],"Anywhere"]}],[{"Loss":[1,{"elements":[["H",null,4],["C",null,1],["O",null,1],["S",null,1]],"additional_mass":0.0,"labels":[]}]}],[]]],"formula":{"elements":[["O",null,1]],"additional_mass":0.0,"labels":[]},"id":{"ontology":"Unimod","name":"Oxidation","id":{"Numeric":35},"description":"Oxidation or Hydroxylation","synonyms":[["Exact","Oxidation or Hydroxylation"]],"cross_ids":[{"Mod":["Resid",{"Alphanumeric":[65,[65,48,48,50,55,32,32]]},{"Position":"Anywhere"}]},{"PubMed":11461766},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,48,50,57,32,32]]},{"Position":"Anywhere"}]},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,48,50,56,32,32]]},{"Position":"Anywhere"}]},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,48,51,48,32,32]]},{"Position":"Anywhere"}]},{"PubMed":9004526},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,50,48,53,32,32]]},{"Position":"Anywhere"}]},{"Article":"Lagerwerf FM, van de Weert M, Heerma W, Haverkamp J, Rapid Commun Mass Spectrom. 1996;10(15):1905-10"},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,49,52,54,32,32]]},{"Position":"Anywhere"}]},{"Findmod":"DOPA"},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,50,49,53,32,32]]},{"Position":"Anywhere"}]},{"Findmod":"CSEA"},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,48,50,54,32,32]]},{"Position":"Anywhere"}]},{"PubMed":14661084},{"Article":"Berlett, Barbara S.; Stadtman, Earl R. Journal of Biological Chemistry (1997), 272(33), 20313-20316."},{"PubMed":15569593},{"PubMed":11120890},{"PubMed":11212008},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,51,50,50,32,32]]},{"Position":"Anywhere"}]},{"Mod":["Resid",{"Alphanumeric":[65,[65,48,50,51,53,32,32]]},{"Position":"Anywhere"}]},{"Findmod":"HYDR"},{"PubMed":14661085},{"PubMed":12781462},{"PubMed":2057999},{"Mod":["Psimod",{"Numeric":35},{"Position":"Anywhere"}]}],"obsolete":false,"parents":[],"children":[]}}}"#,
        ).unwrap()).unwrap());
        assert_eq!(
            modification,
            *STATIC_ONTOLOGIES
                .get_by_index(Ontology::Unimod, &AccessionCode::Numeric(35))
                .unwrap()
        );
        let modification: SimpleModificationInner = dbg!(SimpleModificationInner::from_json_value(serde_json::from_str(
            r#"{"Database":{"specificities":[[[{"AminoAcid":[["Phenylalanine"],"Anywhere"]}],[],[]]],"formula":{"elements":[["O",null,1]],"additional_mass":0.0,"labels":[]},"id":{"ontology":"Psimod","name":"monohydroxylated phenylalanine","id":{"Numeric":2116},"description":"A protein modification that effectively converts an L-phenylalanine residue to a monohydroxylated phenylalanine.","synonyms":[["Related","Oxidation"]],"cross_ids":[{"Mod":["Unimod",{"Numeric":35},{"AminoAcid":[["Phenylalanine"],"Anywhere"]}]}],"obsolete":false,"parents":[{"Numeric":425},{"Numeric":914}],"children":[{"Numeric":1385},{"Numeric":2102},{"Numeric":2104}]}}}"#,
        ).unwrap()).unwrap());
        assert_eq!(
            modification,
            *STATIC_ONTOLOGIES
                .get_by_index(Ontology::Psimod, &AccessionCode::Numeric(2116))
                .unwrap()
        );
    }
}
