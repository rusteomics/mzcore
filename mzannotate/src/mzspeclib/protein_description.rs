use context_error::{BoxedError, Context, CreateError};
use mzcore::{
    prelude::{AminoAcid, IsAminoAcid},
    sequence::FlankingSequence,
};
use mzdata::{
    curie,
    params::{ParamValue, Value},
};

use crate::{
    mzspeclib::{Attribute, AttributeValue, MzSpecLibErrorKind},
    term,
};

/// All protein level metadata
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct ProteinDescription {
    /// Protein accession (MS:1000885) `sp|P12955|PEPD_HUMAN`
    pub accession: Option<Box<str>>,
    /// Protein name (MS:1000886) `Alpha-enolase`
    pub name: Option<Box<str>>,
    /// Database name (MS:1001013)
    pub database_name: Option<Box<str>>,
    /// Database version (MS:1001016)
    pub database_version: Option<Box<str>>,
    /// The cleavage agent or protease (MS:1001045)
    pub cleavage_agent: CleaveAgent,
    /// Protein description (MS:1001088)
    pub description: Option<Box<str>>,
    /// N (MS:1001112) and C (MS:1001113) terminal flanking sequences
    pub flanking_sequences: (FlankingSequence, FlankingSequence),
    /// NCBI species accession (MS:1001467)
    pub species_accession: Option<u32>,
    /// The source species common name (MS:1001468)
    pub species_common_name: Option<Box<str>>,
    /// The source species scientific name (MS:1001469)
    pub species_scientific_name: Option<Box<str>>,
    /// The number of missed cleavages in the peptide (MS:1003044)
    pub missed_cleavages: Option<u16>,
    /// The number of enzymatic termini (0, 1, or 2) (MS:1003048)
    pub enzymatic_termini: Option<u8>,
    /// All set names for this protein (MS:1003212)
    pub set_names: Vec<Box<str>>,
    /// Other attributes
    pub attributes: Vec<Attribute>,
}

impl ProteinDescription {
    /// Check if this description is any if any information is present
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }

    /// Generate the attributes needed to describe to contents of this protein
    pub fn attributes(&self) -> Vec<Attribute> {
        let mut attributes = Vec::new();
        if let Some(acc) = &self.accession {
            attributes.push(Attribute {
                name: term!(MS:1000885|protein accession),
                value: Value::String(acc.to_string()).into(),
            });
        }
        if let Some(name) = &self.name {
            attributes.push(Attribute {
                name: term!(MS:1000886|protein name),
                value: Value::String(name.to_string()).into(),
            });
        }
        if let Some(name) = &self.database_name {
            attributes.push(Attribute {
                name: term!(MS:1001013|database name),
                value: Value::String(name.to_string()).into(),
            });
        }
        if let Some(version) = &self.database_version {
            attributes.push(Attribute {
                name: term!(MS:1001016|database version),
                value: Value::String(version.to_string()).into(),
            });
        }
        match &self.cleavage_agent {
            CleaveAgent::Unknown => (),
            CleaveAgent::Name(name) => attributes.push(Attribute {
                name: term!(MS:1001045|cleavage agent name),
                value: Value::String(name.to_string()).into(),
            }),
            CleaveAgent::Term(term) => attributes.push(Attribute {
                name: term!(MS:1001045|cleavage agent name),
                value: term.clone().into(),
            }),
        }
        if let Some(desc) = &self.description {
            attributes.push(Attribute {
                name: term!(MS:1001088|protein description),
                value: Value::String(desc.to_string()).into(),
            });
        }
        if let Some(string) = match &self.flanking_sequences.0 {
            FlankingSequence::Unknown => None,
            FlankingSequence::Terminal => Some("N-term".to_string()),
            FlankingSequence::AminoAcid(aa) => Some(aa.to_string()),
            FlankingSequence::Sequence(pep) => pep
                .sequence()
                .first()
                .and_then(|s| s.aminoacid.one_letter_code())
                .map(|c| c.to_string()),
        } {
            attributes.push(Attribute {
                name: term!(MS:1001112|n-terminal flanking residue),
                value: Value::String(string).into(),
            });
        }
        if let Some(string) = match &self.flanking_sequences.1 {
            FlankingSequence::Unknown => None,
            FlankingSequence::Terminal => Some("C-term".to_string()),
            FlankingSequence::AminoAcid(aa) => Some(aa.to_string()),
            FlankingSequence::Sequence(pep) => pep
                .sequence()
                .first()
                .and_then(|s| s.aminoacid.one_letter_code())
                .map(|c| c.to_string()),
        } {
            attributes.push(Attribute {
                name: term!(MS:1001113|c-terminal flanking residue),
                value: Value::String(string).into(),
            });
        }
        if let Some(acc) = self.species_accession {
            attributes.push(Attribute {
                name: term!(MS:1001467|taxonomy: NCBI TaxID),
                value: Value::Int(i64::from(acc)).into(),
            });
        }
        if let Some(desc) = &self.species_common_name {
            attributes.push(Attribute {
                name: term!(MS:1001468|taxonomy: common name),
                value: Value::String(desc.to_string()).into(),
            });
        }
        if let Some(desc) = &self.species_scientific_name {
            attributes.push(Attribute {
                name: term!(MS:1001469|taxonomy: scientific name),
                value: Value::String(desc.to_string()).into(),
            });
        }
        if let Some(missed) = self.missed_cleavages {
            attributes.push(Attribute {
                name: term!(MS:1003044|number of missed cleavages),
                value: Value::Int(i64::from(missed)).into(),
            });
        }
        if let Some(termini) = self.enzymatic_termini {
            attributes.push(Attribute {
                name: term!(MS:1003048|number of enzymatic termini),
                value: Value::Int(i64::from(termini)).into(),
            });
        }
        for name in &self.set_names {
            attributes.push(Attribute {
                name: term!(MS:1003212|library attribute set name),
                value: Value::String(name.to_string()).into(),
            });
        }
        attributes.extend_from_slice(&self.attributes);
        attributes
    }

    /// Parse an attribute that is part of a protein description. It saves a properly parsed attribute
    /// in the protein description, and then returns `true`. If the attribute could not be recognised
    /// as a protein description attribute it returns `false` and ignores the attribute.
    /// # Errors
    /// If the attributes contains the wrong type of data, or otherwise contains an invalid value.
    pub(crate) fn populate_from_attribute<'a>(
        &mut self,
        attribute: &Attribute,
        context: &Context<'a>,
    ) -> Result<bool, BoxedError<'a, MzSpecLibErrorKind>> {
        match attribute.name.accession {
            curie!(MS:1003048) => {
                let termini = attribute.value.scalar().to_u64().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid number of enzymatic termini",
                        e.to_string(),
                        context.clone(),
                    )
                })?;
                if termini > 2 {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid number of enzymatic termini",
                        "The number of enzymatic termini can only be 0, 1, or 2.",
                        context.clone(),
                    ));
                }
                self.enzymatic_termini = Some(termini as u8);
            }
            curie!(MS:1003044) => {
                let missed = attribute.value.scalar().to_u64().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid number of missed cleavages",
                        e.to_string(),
                        context.clone(),
                    )
                })?;
                if let Ok(missed) = missed.try_into() {
                    self.missed_cleavages = Some(missed);
                } else {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid number of missed cleavages",
                        "The number of missed cleavages has to be less than 2^16",
                        context.clone(),
                    ));
                }
            }
            curie!(MS:1001112) => {
                let string = attribute.value.scalar().to_string();
                let value = string.trim();
                if value.eq_ignore_ascii_case("n-term") || value == "-" {
                    self.flanking_sequences.0 = FlankingSequence::Terminal;
                } else if let Ok(aa) = value.parse::<AminoAcid>() {
                    self.flanking_sequences.0 = FlankingSequence::AminoAcid(aa);
                } else {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid flanking sequence",
                        "The N terminal flanking residue has to be 'N-term' or an amino acid.",
                        context.clone(),
                    ));
                }
            }
            curie!(MS:1001113) => {
                let string = attribute.value.scalar().to_string();
                let value = string.trim();
                if value.eq_ignore_ascii_case("c-term") || value == "-" {
                    self.flanking_sequences.1 = FlankingSequence::Terminal;
                } else if let Ok(aa) = value.parse::<AminoAcid>() {
                    self.flanking_sequences.1 = FlankingSequence::AminoAcid(aa);
                } else {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid flanking sequence",
                        "The C terminal flanking residue has to be 'C-term' or an amino acid.",
                        context.clone(),
                    ));
                }
            }
            curie!(MS:1003212) => self
                .set_names
                .push(attribute.value.scalar().to_string().into_boxed_str()),
            curie!(MS:1000885) => {
                let string = attribute.value.scalar().to_string();
                if let Some((acc, description)) = string.split_once(' ') {
                    // Make sure only the accession is stored here if the description is also provided
                    self.accession = Some(acc.to_string().into_boxed_str());
                    self.description = Some(description.to_string().into_boxed_str());
                } else {
                    self.accession = Some(string.into_boxed_str());
                }
            }
            curie!(MS:1000886) => {
                self.name = Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001013) => {
                self.database_name = Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001016) => {
                self.database_version = Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001088) => {
                self.description = Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001468) => {
                self.species_common_name =
                    Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001469) => {
                self.species_scientific_name =
                    Some(attribute.value.scalar().to_string().into_boxed_str());
            }
            curie!(MS:1001467) => {
                let accession = attribute.value.scalar().to_u64().map_err(|e| {
                    BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid NCBI TaxID",
                        e.to_string(),
                        context.clone(),
                    )
                })?;
                if let Ok(accession) = accession.try_into() {
                    self.species_accession = Some(accession);
                } else {
                    return Err(BoxedError::new(
                        MzSpecLibErrorKind::Attribute,
                        "Invalid NCBI TaxID",
                        "The NCBI TaxID has to be less than 2^32",
                        context.clone(),
                    ));
                }
            }
            curie!(MS:1001045) => {
                if let AttributeValue::Term(term) = &attribute.value {
                    self.cleavage_agent = CleaveAgent::Term(term.clone());
                } else {
                    self.cleavage_agent =
                        CleaveAgent::Name(attribute.value.scalar().to_string().into_boxed_str());
                }
            }
            _ => return Ok(false),
        }
        Ok(true)
    }
}

/// A cleavage agent or protease
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub enum CleaveAgent {
    /// Unknown
    #[default]
    Unknown,
    /// Named
    Name(Box<str>),
    /// With a CV term
    Term(crate::mzspeclib::Term),
}
