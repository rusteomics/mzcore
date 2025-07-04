use std::{
    fmt::Display,
    ops::{Add, AddAssign},
    str::FromStr,
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    chemistry::MolecularFormula,
    error::{Context, CustomError},
    helper_functions::{explain_number_error, next_number},
    parse_json::{ParseJson, use_serde},
    quantities::Multi,
    sequence::AminoAcid,
};

/// All possible neutral losses
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum NeutralLoss {
    /// Gain of a specific formula
    Gain(u16, MolecularFormula),
    /// Loss of a specific formula
    Loss(u16, MolecularFormula),
    /// Loss of a side chain of an amino acid
    SideChainLoss(MolecularFormula, AminoAcid),
}

impl ParseJson for NeutralLoss {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        let parse_inner = |value: Value, context: &str| {
            if let Value::Array(mut arr) = value {
                if arr.len() == 2 {
                    if let Some(n) = arr[0].as_u64() {
                        Ok((
                            u16::try_from(n).map_err(|_| CustomError::error("Invalid NeutralLoss", format!("The {context} amount is too big, the number has to be below {}", u16::MAX), Context::show(n)))?,
                            MolecularFormula::from_json_value(arr.pop().unwrap())?,
                        ))
                    } else {
                        Err(CustomError::error(
                            "Invalid NeutralLoss",
                            format!("The {context} amount is not a number"),
                            Context::show(arr[0].to_string()),
                        ))
                    }
                } else {
                    Err(CustomError::error(
                        "Invalid NeutralLoss",
                        format!("The {context} is a sequence but does not have 2 children"),
                        Context::show(arr.iter().join(",")),
                    ))
                }
            } else if value.is_object() {
                Ok((1, MolecularFormula::from_json_value(value)?))
            } else {
                Err(CustomError::error(
                    "Invalid NeutralLoss",
                    format!("The {context} has to be either a map or a sequence"),
                    Context::show(value),
                ))
            }
        };
        if let Value::Object(map) = value {
            let (key, value) = map.into_iter().next().unwrap();
            match key.as_str() {
                "Loss" => parse_inner(value, "Loss").map(|(n, f)| Self::Loss(n, f)),
                "Gain" => parse_inner(value, "Gain").map(|(n, f)| Self::Gain(n, f)),
                "SideChainLoss" => {
                    if let Value::Array(mut arr) = value {
                        if arr.len() == 2 {
                            let aa = AminoAcid::from_json_value(arr.pop().unwrap())?;
                            let formula = MolecularFormula::from_json_value(arr.pop().unwrap())?;
                            Ok(Self::SideChainLoss(formula, aa))
                        } else {
                            Err(CustomError::error(
                                "Invalid NeutralLoss",
                                "The SideChainLoss is a sequence but does not have 2 children",
                                Context::show(arr.iter().join(",")),
                            ))
                        }
                    } else {
                        Err(CustomError::error(
                            "Invalid NeutralLoss",
                            "The SideChainLoss has to be a sequence",
                            Context::show(value),
                        ))
                    }
                }
                _ => Err(CustomError::error(
                    "Invalid NeutralLoss",
                    "The tag has to be Loss/Gain/SideChainLoss",
                    Context::show(key),
                )),
            }
        } else {
            Err(CustomError::error(
                "Invalid NeutralLoss",
                "The JSON value has to be a map",
                Context::show(value),
            ))
        }
    }
}

/// A diagnostic ion, defined in M (not MH+) chemical formula
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct DiagnosticIon(pub MolecularFormula);

impl ParseJson for DiagnosticIon {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl NeutralLoss {
    /// Check if this neutral loss if empty (has an empty molecular formula)
    pub fn is_empty(&self) -> bool {
        match self {
            Self::Gain(0, _) | Self::Loss(0, _) => true,
            Self::Loss(_, f) | Self::Gain(_, f) | Self::SideChainLoss(f, _) => f.is_empty(),
        }
    }

    /// Generate a nice HTML notation for this `NeutralLoss`
    pub fn hill_notation_html(&self) -> String {
        match self {
            Self::Loss(n, c) => format!("-{n}{}", c.hill_notation_html().trim_start_matches('+')),
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, c) => format!("+{n}{}", c.hill_notation_html().trim_start_matches('+')),
        }
    }

    /// Generate a nice fancy notation for this `NeutralLoss`
    pub fn hill_notation_fancy(&self) -> String {
        match self {
            Self::Loss(n, c) => format!("-{n}{}", c.hill_notation_fancy().trim_start_matches('+')),
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, c) => format!("+{n}{}", c.hill_notation_fancy().trim_start_matches('+')),
        }
    }

    /// Generate a notation for this `NeutralLoss` with pure ASCII characters
    pub fn hill_notation(&self) -> String {
        match self {
            Self::Loss(n, c) => format!("-{n}{}", c.hill_notation().trim_start_matches('+')),
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, c) => format!("+{n}{}", c.hill_notation().trim_start_matches('+')),
        }
    }
}

impl FromStr for NeutralLoss {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Allow a numeric neutral loss
        if let Ok(number) = s.parse::<f64>() {
            if number > 0.0 {
                Ok(Self::Gain(
                    1,
                    MolecularFormula::with_additional_mass(number),
                ))
            } else {
                Ok(Self::Loss(
                    1,
                    MolecularFormula::with_additional_mass(number.abs()),
                ))
            }
        } else if let Some(c) = s.chars().next() {
            // Or match a molecular formula
            let loss = match c {
                '+' => Ok(false),
                '-' => Ok(true),
                _ => Err(CustomError::error(
                    "Invalid neutral loss",
                    "A neutral loss can only start with '+' or '-'",
                    Context::line(None, s, 0, 1),
                )),
            }?;
            let (amount, start) = if let Some(amount) = next_number::<false, false, u16>(s, 1..) {
                (
                    amount.2.map_err(|err| {
                        CustomError::error(
                            "Invalid neutral loss",
                            format!("The amount specifier {}", explain_number_error(&err)),
                            Context::line(None, s, 1, amount.0),
                        )
                    })?,
                    amount.0,
                )
            } else {
                (1, 1)
            };
            let formula = MolecularFormula::from_pro_forma(s, start.., false, false, true, true)?;
            Ok(if loss {
                Self::Loss(amount, formula)
            } else {
                Self::Gain(amount, formula)
            })
        } else {
            Err(CustomError::error(
                "Invalid neutral loss",
                "A neutral loss cannot be an empty string",
                Context::None,
            ))
        }
    }
}

impl Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.hill_notation())
    }
}

impl Add<&NeutralLoss> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &NeutralLoss) -> Self::Output {
        match rhs {
            NeutralLoss::Gain(n, mol) => self + mol * n,
            NeutralLoss::Loss(n, mol) => self - mol * n,
            NeutralLoss::SideChainLoss(mol, _) => self - mol,
        }
    }
}

impl AddAssign<&NeutralLoss> for MolecularFormula {
    fn add_assign(&mut self, rhs: &NeutralLoss) {
        match rhs {
            NeutralLoss::Gain(n, mol) => *self += mol * n,
            NeutralLoss::Loss(n, mol) => *self -= mol * n,
            NeutralLoss::SideChainLoss(mol, _) => *self -= mol,
        }
    }
}

impl AddAssign<NeutralLoss> for MolecularFormula {
    fn add_assign(&mut self, rhs: NeutralLoss) {
        *self += &rhs;
    }
}

impl Add<&NeutralLoss> for &Multi<MolecularFormula> {
    type Output = Multi<MolecularFormula>;
    fn add(self, rhs: &NeutralLoss) -> Self::Output {
        match rhs {
            NeutralLoss::Gain(n, mol) => self + mol * n,
            NeutralLoss::Loss(n, mol) => self - mol * n,
            NeutralLoss::SideChainLoss(mol, _) => self - mol,
        }
    }
}

impl_binop_ref_cases!(impl Add, add for MolecularFormula, NeutralLoss, MolecularFormula);
impl_binop_ref_cases!(impl Add, add for Multi<MolecularFormula>, NeutralLoss, Multi<MolecularFormula>);

impl<'a> std::iter::Sum<&'a NeutralLoss> for MolecularFormula {
    fn sum<I: Iterator<Item = &'a NeutralLoss>>(iter: I) -> Self {
        let mut output = Self::default();
        for value in iter {
            output += value;
        }
        output
    }
}

impl std::iter::Sum<NeutralLoss> for MolecularFormula {
    fn sum<I: Iterator<Item = NeutralLoss>>(iter: I) -> Self {
        let mut output = Self::default();
        for value in iter {
            output += value;
        }
        output
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        fragment::NeutralLoss, molecular_formula, parse_json::ParseJson, prelude::AminoAcid,
    };

    #[test]
    fn deserialise_json() {
        let original = r#"{"Loss":{"elements":[["H",null,2],["O",null,1]],"additional_mass":0.0,"labels":[]}}"#;
        let new = r#"{"Loss":[1,{"elements":[["H",null,2],["O",null,1]],"additional_mass":0.0,"labels":[]}]}"#;
        let current = serde_json::to_string(&NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)))
            .expect("Could not serialise");

        let original_v = NeutralLoss::from_json(original).expect("Could not deserialise original");
        let new_v = NeutralLoss::from_json(new).expect("Could not deserialise new");
        let current_v = NeutralLoss::from_json(&current).expect("Could not deserialise current");

        assert_eq!(original_v, new_v);
        assert_eq!(current_v, new_v);

        let side_chain_loss =
            NeutralLoss::SideChainLoss(molecular_formula!(H 2 O 1), AminoAcid::Alanine);
        let side_chain_loss_json =
            serde_json::to_string(&side_chain_loss).expect("Could not serialise side chain loss");
        println!("{side_chain_loss_json}");
        let side_chain_loss_back = NeutralLoss::from_json(&side_chain_loss_json)
            .expect("Could not deserialise side chain loss");

        assert_eq!(side_chain_loss, side_chain_loss_back);
    }
}
