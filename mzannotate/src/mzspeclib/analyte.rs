use itertools::Itertools;
use mzcore::{
    prelude::{
        Chemical, IsAminoAcid, MolecularCharge, MolecularFormula, MultiChemical, PeptidoformIon,
    },
    system::isize::Charge,
};
use mzdata::{Param, params::Value};

use crate::{
    mzspeclib::{Attribute, AttributeValue, Attributes, Id, ProteinDescription},
    term,
};

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
    pub fn attributes(&self) -> Attributes {
        let mut result = self.target.terms();
        for protein in &self.proteins {
            let res = protein.attributes();
            result.push(res);
        }
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
    pub fn terms(&self) -> Attributes {
        let mut attributes = vec![Vec::new(); 1];

        if let Some(charge) = match self {
            Self::Unknown(c) => *c,
            Self::MolecularFormula(f) => (f.charge().value != 0).then_some(f.charge()),
            Self::PeptidoformIon(pep) => pep.get_charge_carriers().map(MolecularCharge::charge),
        } {
            attributes[0].push(Attribute {
                name: term!(MS:1000041|charge state),
                value: AttributeValue::Scalar(Value::Int(charge.value as i64)),
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
                attributes[0].push(Attribute {
                    name: term!(MS:1000866|molecular formula),
                    value: AttributeValue::Scalar(Value::String(formula.hill_notation_core())),
                });
            }
            let charge = formula.charge().value as f64;
            if charge != 0.0 {
                attributes[0].push(Attribute {
                    name: term!(MS:1003053|theoretical monoisotopic m/z),
                    value: AttributeValue::Scalar(Value::Float(
                        formula.monoisotopic_mass().value / charge,
                    )),
                });
                attributes[0].push(Attribute {
                    name: term!(MS:1003054|theoretical average m/z),
                    value: AttributeValue::Scalar(Value::Float(
                        formula.average_weight().value / charge,
                    )),
                });
            }
            // TODO: when PSI-MS gets new mass terms use those
            attributes[0].push(Attribute {
                name: term!(MS:1003243|adduct ion mass),
                value: AttributeValue::Scalar(Value::Float(formula.monoisotopic_mass().value)),
            });
            attributes[0].push(Attribute {
                name: term!(MS:1000224|molecular mass),
                value: AttributeValue::Scalar(Value::Float(formula.average_weight().value)),
            });
        }

        if let Some(pep) = match self {
            Self::Unknown(_) | Self::MolecularFormula(_) => None,
            Self::PeptidoformIon(pep) => Some(pep),
        } {
            if let Ok(formula) = pep.formulas().to_vec().into_iter().exactly_one() {
                attributes[0].push(Attribute {
                    name: term!(MS:1001117|theoretical neutral mass),
                    value: AttributeValue::Scalar(Value::Float(formula.monoisotopic_mass().value)),
                });
            }
            attributes[0].push(Attribute {
                name: term!(MS:1003270|proforma peptidoform ion notation),
                value: AttributeValue::Scalar(Value::String(pep.to_string())),
            });
            if pep.peptidoforms().len() == 1 {
                attributes[0].push(Attribute {
                    name: term!(MS:1000888|stripped peptide sequence),
                    value: AttributeValue::Scalar(Value::String(
                        pep.peptidoforms()[0]
                            .sequence()
                            .iter()
                            .map(|s| s.aminoacid.aminoacid().one_letter_code().unwrap_or('X'))
                            .collect(),
                    )),
                });
                // No reason to crash on a length that is a 'bit' big
                if let Ok(v) = i64::try_from(pep.peptidoforms()[0].sequence().len()) {
                    attributes[0].push(Attribute {
                        name: term!(MS:1003043|number of residues),
                        value: AttributeValue::Scalar(Value::Int(v)),
                    });
                }
            }
        }

        attributes
    }
}
