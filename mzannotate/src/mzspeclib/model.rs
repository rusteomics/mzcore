use std::{collections::HashMap, fmt::Display};

use mzcore::{
    prelude::{IsAminoAcid, MolecularCharge, MolecularFormula, MultiChemical, PeptidoformIon},
    sequence::FlankingSequence,
    system::isize::Charge,
};
use mzdata::{Param, params::Value};

use crate::{
    mzspeclib::{Attribute, AttributeSet, AttributeValue, EntryType},
    term,
};

/// The header for a spectral library
#[derive(Debug, Clone)]
pub struct LibraryHeader {
    /// The version of the format
    pub format_version: String,
    /// The attributes for this library
    pub attributes: Vec<Attribute>,
    /// The attribute classes for this library
    pub attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
}

impl Default for LibraryHeader {
    fn default() -> Self {
        Self {
            format_version: "1.0".into(),
            attributes: Vec::new(),
            attribute_classes: HashMap::new(),
        }
    }
}

impl Display for LibraryHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let version = Attribute::new(
            term!(MS:1003186|library format version),
            AttributeValue::Scalar(Value::String(self.format_version.clone())),
            None,
        );
        writeln!(f, "<mzSpecLib>\n{version}")?;
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl LibraryHeader {
    /// Create a new library header
    pub const fn new(
        format_version: String,
        attributes: Vec<Attribute>,
        attribute_classes: HashMap<EntryType, Vec<AttributeSet>>,
    ) -> Self {
        Self {
            format_version,
            attributes,
            attribute_classes,
        }
    }
}

/// An ID
pub type Id = u32;

/// An analyte that resulted in the spectrum.
#[derive(Default, Debug, Clone)]
pub struct Analyte {
    /// The numeric id of this analyte
    pub id: Id,
    /// The target itself
    pub target: AnalyteTarget,
    /// The matched protein(s)
    pub proteins: Vec<ProteinDescription>,
    /// Other parameters
    pub params: Vec<Param>,
}

impl Analyte {
    /// Create a new analyte
    pub const fn new(id: Id, target: AnalyteTarget) -> Self {
        Self {
            id,
            target,
            proteins: Vec::new(),
            params: Vec::new(),
        }
    }

    /// Get all attributes of this analyte
    pub fn attributes(&self) -> Vec<Attribute> {
        let mut result = self.target.terms();
        result.extend(
            self.proteins
                .iter()
                .enumerate()
                .flat_map(|(i, p)| p.attributes(i as u32)),
        );
        result
    }
}

/// The analyte itself
#[derive(Debug, Clone)]
pub enum AnalyteTarget {
    /// For an unknown class of analytes, where the charge (MS:1000041) could be known
    Unknown(Option<Charge>),
    /// A peptidoform ion (MS:1003270), with the charge (MS:1000041) stored inside if known
    PeptidoformIon(PeptidoformIon),
    /// A molecular formula (MS:1000866), with the charge (MS:1000041) stored inside if known
    MolecularFormula(MolecularFormula),
}

impl Default for AnalyteTarget {
    fn default() -> Self {
        Self::Unknown(None)
    }
}

impl MultiChemical for AnalyteTarget {
    fn formulas_inner(
        &self,
        _sequence_index: mzcore::prelude::SequencePosition,
        _peptidoform_index: usize,
        _peptidoform_ion_index: usize,
    ) -> mzcore::quantities::Multi<MolecularFormula> {
        match self {
            Self::Unknown(_) => mzcore::quantities::Multi::default(),
            Self::PeptidoformIon(pep) => pep.formulas(),
            Self::MolecularFormula(f) => f.into(),
        }
    }
}

impl AnalyteTarget {
    /// Update the charge to a new value. This does not change the definition of the peptidoform
    /// ion if that already contains a charge with the same value. This means that any complex
    /// charge carriers defined in ProForma will not be touched.
    pub fn set_charge(&mut self, charge: Charge) {
        match self {
            Self::Unknown(c) => *c = Some(charge),
            Self::PeptidoformIon(pep) => {
                if pep
                    .get_charge_carriers()
                    .is_none_or(|c| c.charge() != charge)
                {
                    pep.set_charge_carriers(Some(MolecularCharge::proton(charge)));
                }
            }
            Self::MolecularFormula(f) => f.set_charge(charge),
        }
    }

    /// Used to write out the analyte target for mzSpecLib files
    // TODO: look into other theoretical things I can easily generate
    pub fn terms(&self) -> Vec<Attribute> {
        match self {
            Self::Unknown(None) => Vec::new(),
            Self::Unknown(Some(c)) => vec![Attribute {
                name: term!(MS:1000041|charge state),
                value: AttributeValue::Scalar(Value::Int(c.value as i64)),
                group_id: None,
            }],
            Self::PeptidoformIon(pep) => {
                let mut attributes = vec![Attribute {
                    name: term!(MS:1003270|proforma peptidoform ion notation),
                    value: AttributeValue::Scalar(Value::String(pep.to_string())),
                    group_id: None,
                }];
                if let Some(c) = pep.get_charge_carriers().map(MolecularCharge::charge)
                    && c.value != 0
                {
                    attributes.push(Attribute {
                        name: term!(MS:1000041|charge state),
                        value: AttributeValue::Scalar(Value::Int(c.value as i64)),
                        group_id: None,
                    });
                }
                let f = pep.formulas();
                if f.len() == 1 && f[0].additional_mass() == 0.0 {
                    attributes.push(Attribute {
                        name: term!(MS:1000866|molecular formula),
                        value: AttributeValue::Scalar(Value::String(f[0].hill_notation_core())),
                        group_id: None,
                    });
                }
                if pep.peptidoforms().len() == 1 {
                    attributes.push(Attribute {
                        name: term!(MS:1000888|stripped peptide sequence),
                        value: AttributeValue::Scalar(Value::String(
                            pep.peptidoforms()[0]
                                .sequence()
                                .iter()
                                .map(|s| s.aminoacid.aminoacid().one_letter_code().unwrap_or('X'))
                                .collect(),
                        )),
                        group_id: None,
                    });
                }
                attributes
            }
            Self::MolecularFormula(f) => {
                let mut attributes = vec![Attribute {
                    name: term!(MS:1000866|molecular formula),
                    value: AttributeValue::Scalar(Value::String(f.hill_notation_core())),
                    group_id: None,
                }];
                if f.charge().value != 0 {
                    attributes.push(Attribute {
                        name: term!(MS:1000041|charge state),
                        value: AttributeValue::Scalar(Value::Int(f.charge().value as i64)),
                        group_id: None,
                    });
                }
                attributes
            }
        }
    }
}

/// All protein level metadata
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ProteinDescription {
    /// The id of the group / protein
    pub id: Id,
    /// Protein accession (MS:1000885) `sp|P12955|PEPD_HUMAN`
    pub accession: Option<Box<str>>,
    /// Protein name (MS:1000886) `Alpha-enolase`
    pub name: Option<Box<str>>,
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
    pub fn attributes(&self, id: u32) -> Vec<Attribute> {
        let mut attributes = Vec::new();
        if let Some(acc) = &self.accession {
            attributes.push(Attribute {
                name: term!(MS:1000885|protein accession),
                value: Value::String(acc.to_string()).into(),
                group_id: Some(id),
            });
        }
        if let Some(name) = &self.name {
            attributes.push(Attribute {
                name: term!(MS:1000886|protein name),
                value: Value::String(name.to_string()).into(),
                group_id: Some(id),
            });
        }
        match &self.cleavage_agent {
            CleaveAgent::Unknown => (),
            CleaveAgent::Name(name) => attributes.push(Attribute {
                name: term!(MS:1001045|cleavage agent name),
                value: Value::String(name.to_string()).into(),
                group_id: Some(id),
            }),
            CleaveAgent::Term(term) => attributes.push(Attribute {
                name: term!(MS:1001045|cleavage agent name),
                value: term.clone().into(),
                group_id: Some(id),
            }),
        }
        if let Some(desc) = &self.description {
            attributes.push(Attribute {
                name: term!(MS:1001088|protein description),
                value: Value::String(desc.to_string()).into(),
                group_id: Some(id),
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
                group_id: Some(id),
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
                group_id: Some(id),
            });
        }
        if let Some(acc) = self.species_accession {
            attributes.push(Attribute {
                name: term!(MS:1001467|taxonomy: NCBI TaxID),
                value: Value::Int(i64::from(acc)).into(),
                group_id: Some(id),
            });
        }
        if let Some(desc) = &self.species_common_name {
            attributes.push(Attribute {
                name: term!(MS:1001468|taxonomy: common name),
                value: Value::String(desc.to_string()).into(),
                group_id: Some(id),
            });
        }
        if let Some(desc) = &self.species_scientific_name {
            attributes.push(Attribute {
                name: term!(MS:1001469|taxonomy: scientific name),
                value: Value::String(desc.to_string()).into(),
                group_id: Some(id),
            });
        }
        if let Some(missed) = self.missed_cleavages {
            attributes.push(Attribute {
                name: term!(MS:1003044|number of missed cleavages),
                value: Value::Int(i64::from(missed)).into(),
                group_id: Some(id),
            });
        }
        if let Some(termini) = self.enzymatic_termini {
            attributes.push(Attribute {
                name: term!(MS:1003048|number of enzymatic termini),
                value: Value::Int(i64::from(termini)).into(),
                group_id: Some(id),
            });
        }
        for name in &self.set_names {
            attributes.push(Attribute {
                name: term!(MS:1003212|library attribute set name),
                value: Value::String(name.to_string()).into(),
                group_id: Some(id),
            });
        }
        attributes.extend(self.attributes.iter().map(|a| {
            let mut a = a.clone();
            a.group_id = Some(id);
            a
        }));
        attributes
    }
}

/// A cleavage agent or protease
#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum CleaveAgent {
    /// Unknown
    #[default]
    Unknown,
    /// Named
    Name(Box<str>),
    /// With a CV term
    Term(crate::mzspeclib::Term),
}

#[derive(Default, Debug, Clone)]
pub struct Interpretation {
    pub id: Id,
    pub attributes: Vec<Attribute>,
    pub analyte_refs: Vec<Id>,
    pub members: Vec<InterpretationMember>,
}

impl Interpretation {
    pub const fn new(
        id: Id,
        attributes: Vec<Attribute>,
        analyte_refs: Vec<Id>,
        members: Vec<InterpretationMember>,
    ) -> Self {
        Self {
            id,
            attributes,
            analyte_refs,
            members,
        }
    }
}

impl Display for Interpretation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Interpretation={}>", self.id)?;
        if !self.analyte_refs.is_empty() {
            let val = AttributeValue::List(
                self.analyte_refs
                    .iter()
                    .map(|v| Value::Int(i64::from(*v)))
                    .collect(),
            );
            let mixture_ids = Attribute::new(term!(MS:1003163|analyte mixture members), val, None);
            writeln!(f, "{mixture_ids}")?;
        }
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        for member in &self.members {
            write!(f, "{member}")?;
        }
        Ok(())
    }
}

#[derive(Default, Debug, Clone)]
pub struct InterpretationMember {
    pub id: Id,
    pub attributes: Vec<Attribute>,
    pub analyte_ref: Id,
}

impl Display for InterpretationMember {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<InterpretationMember={}>", self.id)?;
        for attr in &self.attributes {
            f.write_str(attr.to_string().as_str())?;
        }
        Ok(())
    }
}
