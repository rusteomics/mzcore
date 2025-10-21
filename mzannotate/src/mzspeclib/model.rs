use std::{collections::HashMap, fmt::Display};

use itertools::Itertools;
use mzcore::{
    prelude::{
        Chemical, IsAminoAcid, MolecularCharge, MolecularFormula, MultiChemical, PeptidoformIon,
    },
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
        for t in [
            &EntryType::Spectrum,
            &EntryType::Analyte,
            &EntryType::Interpretation,
            &EntryType::Cluster,
        ] {
            for group in self.attribute_classes.get(t).iter().flat_map(|s| s.iter()) {
                writeln!(f, "<AttributeSet {t}={}>", group.id)?;
                for attr in group
                    .attributes
                    .values()
                    .flat_map(|a| a.iter().map(|(a, _)| a))
                {
                    writeln!(f, "{attr}")?;
                }
            }
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

    /// Check if this attribute is already set.
    pub(crate) fn is_already_defined(
        &self,
        attribute: &Attribute,
        entry: EntryType,
        names: &[&str],
    ) -> bool {
        self.attribute_classes
            .get(&entry)
            .iter()
            .flat_map(|s| s.iter())
            .filter(|s| names.contains(&s.id.as_str()))
            .any(|s| {
                s.attributes
                    .values()
                    .flat_map(|a| a.iter())
                    .any(|a| a.0.name == attribute.name && a.0.value.equivalent(&attribute.value))
            })
    }
}

/// An ID
pub type Id = u32;

/// An analyte that resulted in the spectrum.
#[derive(Default, Debug, Clone, PartialEq, Eq)]
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
#[derive(Debug, Clone, PartialEq, Eq)]
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
    pub fn terms(&self) -> Vec<Attribute> {
        let mut attributes = Vec::new();

        if let Some(charge) = match self {
            Self::Unknown(c) => *c,
            Self::MolecularFormula(f) => (f.charge().value != 0).then_some(f.charge()),
            Self::PeptidoformIon(pep) => pep.get_charge_carriers().map(MolecularCharge::charge),
        } {
            attributes.push(Attribute {
                name: term!(MS:1000041|charge state),
                value: AttributeValue::Scalar(Value::Int(charge.value as i64)),
                group_id: None,
            });
        }

        if let Some(formula) = match self {
            Self::Unknown(_) => None,
            Self::MolecularFormula(f) => Some(f.clone()),
            Self::PeptidoformIon(pep) => (pep.formulas()
                + pep
                    .get_charge_carriers()
                    .map(Chemical::formula)
                    .unwrap_or_default())
            .to_vec()
            .into_iter()
            .exactly_one()
            .ok(),
        } {
            if formula.additional_mass() == 0.0 {
                attributes.push(Attribute {
                    name: term!(MS:1000866|molecular formula),
                    value: AttributeValue::Scalar(Value::String(formula.hill_notation_core())),
                    group_id: None,
                });
            }
            let charge = formula.charge().value as f64;
            if charge != 0.0 {
                attributes.push(Attribute {
                    name: term!(MS:1003053|theoretical monoisotopic m/z),
                    value: AttributeValue::Scalar(Value::Float(
                        formula.monoisotopic_mass().value / charge,
                    )),
                    group_id: None,
                });
                attributes.push(Attribute {
                    name: term!(MS:1003054|theoretical average m/z),
                    value: AttributeValue::Scalar(Value::Float(
                        formula.average_weight().value / charge,
                    )),
                    group_id: None,
                });
            }
            attributes.push(Attribute {
                name: term!(MS:1003243|adduct ion mass),
                value: AttributeValue::Scalar(Value::Float(formula.monoisotopic_mass().value)),
                group_id: None,
            });
            attributes.push(Attribute {
                name: term!(MS:1000224|molecular mass),
                value: AttributeValue::Scalar(Value::Float(formula.average_weight().value)),
                group_id: None,
            });
        }

        if let Some(pep) = match self {
            Self::Unknown(_) | Self::MolecularFormula(_) => None,
            Self::PeptidoformIon(pep) => Some(pep),
        } {
            if let Ok(formula) = pep.formulas().to_vec().into_iter().exactly_one() {
                attributes.push(Attribute {
                    name: term!(MS:1001117|theoretical neutral mass),
                    value: AttributeValue::Scalar(Value::Float(formula.monoisotopic_mass().value)),
                    group_id: None,
                });
            }
            attributes.push(Attribute {
                name: term!(MS:1003270|proforma peptidoform ion notation),
                value: AttributeValue::Scalar(Value::String(pep.to_string())),
                group_id: None,
            });
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
                attributes.push(Attribute {
                    name: term!(MS:1003043|number of residues),
                    value: AttributeValue::Scalar(Value::Int(
                        pep.peptidoforms()[0].sequence().len() as i64,
                    )),
                    group_id: None,
                });
            }
        }

        attributes
    }
}

/// All protein level metadata
#[derive(Default, Debug, Clone, PartialEq, Eq)]
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
        if let Some(name) = &self.database_name {
            attributes.push(Attribute {
                name: term!(MS:1001013|database name),
                value: Value::String(name.to_string()).into(),
                group_id: Some(id),
            });
        }
        if let Some(version) = &self.database_version {
            attributes.push(Attribute {
                name: term!(MS:1001016|database version),
                value: Value::String(version.to_string()).into(),
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

/// An interpreattion, this gives details on how well the analytes match to the spectrum.
#[derive(Default, Debug, Clone, PartialEq)]
pub struct Interpretation {
    /// The ID
    pub id: Id,
    /// PSM level probability (MS:1002357) in range 0..=1
    pub probability: Option<f64>,
    /// Othe rattributes
    pub attributes: Vec<Attribute>,
    /// Analyte mixture members (MS:1003163)
    pub analyte_refs: Vec<Id>,
    /// The interpretation members, this links to a specific analyte and gives additional attributes
    /// about the interpretation of that analyte.
    pub members: HashMap<Id, Vec<Attribute>>,
}

impl Interpretation {
    pub fn attributes(&self) -> Vec<Attribute> {
        // TODO: add easily calculated attributes
        let mut attributes = Vec::new();
        if !self.analyte_refs.is_empty() {
            attributes.push(Attribute::new(
                term!(MS:1003163|analyte mixture members),
                AttributeValue::List(
                    self.analyte_refs
                        .iter()
                        .map(|v| Value::Int(i64::from(*v)))
                        .collect(),
                ),
                None,
            ));
        }
        if let Some(probability) = self.probability {
            attributes.push(Attribute::new(
                term!(MS:1002357|PSM-level probability),
                Value::Float(probability),
                None,
            ));
        }
        attributes.extend_from_slice(&self.attributes);

        attributes
    }
}
