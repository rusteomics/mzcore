use std::{collections::HashMap, fmt::Display};

use mzcore::{
    prelude::{IsAminoAcid, MolecularCharge, MolecularFormula, MultiChemical, PeptidoformIon},
    sequence::FlankingSequence,
    system::isize::Charge,
};
use mzdata::{curie, params::Value};

use crate::{
    mzspeclib::{
        Attribute, AttributeSet, AttributeValue, Attributed, AttributedMut, EntryType,
        impl_attributed,
    },
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

impl_attributed!(mut LibraryHeader);

impl Display for LibraryHeader {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let version = Attribute::new(
            term!(MS:1003186|"library format version"),
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
    /// Other attributes for this analyte
    pub attributes: Vec<Attribute>,
}

impl_attributed!(mut Analyte);

impl Display for Analyte {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<Analyte={}>", self.id)?;
        for attr in self.target.terms() {
            writeln!(f, "{attr}")?;
        }
        for attr in &self.attributes {
            writeln!(f, "{attr}")?;
        }
        Ok(())
    }
}

impl Analyte {
    /// Create a new analyte
    pub const fn new(id: Id, target: AnalyteTarget, attributes: Vec<Attribute>) -> Self {
        Self {
            id,
            target,
            proteins: Vec::new(),
            attributes,
        }
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
                    pep.set_charge_carriers(Some(MolecularCharge::proton(charge)))
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
                name: term!(MS:1000041|"charge state"),
                value: AttributeValue::Scalar(Value::Int(c.value as i64)),
                group_id: None,
            }],
            Self::PeptidoformIon(pep) => {
                let mut attributes = vec![Attribute {
                    name: term!(MS:1003270|"proforma peptidoform ion notation"),
                    value: AttributeValue::Scalar(Value::String(pep.to_string())),
                    group_id: None,
                }];
                if let Some(c) = pep.get_charge_carriers().map(MolecularCharge::charge)
                    && c.value != 0
                {
                    attributes.push(Attribute {
                        name: term!(MS:1000041|"charge state"),
                        value: AttributeValue::Scalar(Value::Int(c.value as i64)),
                        group_id: None,
                    });
                }
                let f = pep.formulas();
                if f.len() == 1 && f[0].additional_mass() == 0.0 {
                    attributes.push(Attribute {
                        name: term!(MS:1000866|"molecular formula"),
                        value: AttributeValue::Scalar(Value::String(f[0].hill_notation_core())),
                        group_id: None,
                    });
                }
                if pep.peptidoforms().len() == 1 {
                    attributes.push(Attribute {
                        name: term!(MS:1000888|"stripped peptide sequence"),
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
                    name: term!(MS:1000866|"molecular formula"),
                    value: AttributeValue::Scalar(Value::String(f.hill_notation_core())),
                    group_id: None,
                }];
                if f.charge().value != 0 {
                    attributes.push(Attribute {
                        name: term!(MS:1000041|"charge state"),
                        value: AttributeValue::Scalar(Value::Int(f.charge().value as i64)),
                        group_id: None,
                    });
                }
                attributes
            }
        }
    }
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub struct ProteinDescription {
    pub id: Id,
    pub flanking_sequences: (FlankingSequence, FlankingSequence),
    /// The number of enzymatic termini (0, 1, or 2)
    pub enzymatic_termini: Option<u8>,
    /// The number of missed cleavages in the peptide
    pub missed_cleavages: Option<u16>,
    pub set_names: Vec<Box<str>>,
    pub accession: Option<Box<str>>,
    pub name: Option<Box<str>>,
    pub description: Option<Box<str>>,
    pub cleavage_agent: CleaveAgent,
    pub species_scientific_name: Option<Box<str>>,
    pub species_common_name: Option<Box<str>>,
    /// NCBI species accession
    pub species_accession: Option<u32>,
}

impl ProteinDescription {
    /// Check if this description is any if any information is present
    pub fn is_empty(&self) -> bool {
        *self == Self::default()
    }
}

#[derive(Default, Debug, Clone, PartialEq, Eq)]
pub enum CleaveAgent {
    #[default]
    Unknown,
    Name(Box<str>),
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

impl_attributed!(mut Interpretation);

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
            let mixture_ids =
                Attribute::new(term!(MS:1_003_163|"analyte mixture members"), val, None);
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

impl_attributed!(mut InterpretationMember);

impl Display for InterpretationMember {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "<InterpretationMember={}>", self.id)?;
        for attr in &self.attributes {
            f.write_str(attr.to_string().as_str())?;
        }
        Ok(())
    }
}
