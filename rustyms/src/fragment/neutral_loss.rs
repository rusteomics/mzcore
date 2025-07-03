use std::{
    fmt::Display,
    ops::{Add, AddAssign},
    str::FromStr,
};

use serde::{
    Deserialize, Serialize,
    de::{Error, Visitor},
};

use crate::{
    chemistry::MolecularFormula,
    error::{Context, CustomError},
    helper_functions::{explain_number_error, next_number},
    quantities::Multi,
    sequence::AminoAcid,
};

/// All possible neutral losses
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum NeutralLoss {
    /// Gain of a specific formula
    Gain(u16, MolecularFormula),
    /// Loss of a specific formula
    Loss(u16, MolecularFormula),
    /// Loss of a side chain of an amino acid
    SideChainLoss(MolecularFormula, AminoAcid),
}

/// Possible internal values, needed to allow serde to still parse the older format which did not contain an amount value
#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum Field {
    Formula(MolecularFormula),
    Amount(u16, MolecularFormula),
    SideChain(MolecularFormula, AminoAcid),
}

impl<'de> Deserialize<'de> for NeutralLoss {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        deserializer.deserialize_map(NeutralLossVisitor)
    }
}

struct NeutralLossVisitor;

impl<'de> Visitor<'de> for NeutralLossVisitor {
    type Value = NeutralLoss;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("a neutral loss")
    }

    fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
    where
        A: serde::de::MapAccess<'de>,
    {
        if let Some((key, value)) = map.next_entry::<&str, Field>()? {
            match key {
                "Loss" => match value {
                    Field::Formula(f) => Ok(NeutralLoss::Loss(1, f)),
                    Field::Amount(n, f) => Ok(NeutralLoss::Loss(n, f)),
                    Field::SideChain(_, _) => Err(A::Error::custom(
                        "a loss neutral loss cannot contain side chain loss information",
                    )),
                },
                "Gain" => match value {
                    Field::Formula(f) => Ok(NeutralLoss::Gain(1, f)),
                    Field::Amount(n, f) => Ok(NeutralLoss::Gain(n, f)),
                    Field::SideChain(_, _) => Err(A::Error::custom(
                        "a gain neutral loss cannot contain side chain loss information",
                    )),
                },
                "SideChainLoss" => match value {
                    Field::SideChain(f, a) => Ok(NeutralLoss::SideChainLoss(f, a)),
                    _ => Err(A::Error::custom(
                        "a side chain neutral loss cannot contain basic loss/gain information",
                    )),
                },
                v => Err(A::Error::custom(format!(
                    "expected Loss/Gain/SideChainLoss not '{v}'"
                ))),
            }
        } else {
            Err(A::Error::custom("Empty neutral loss definition"))
        }
    }
}

/// A diagnostic ion, defined in M (not MH+) chemical formula
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct DiagnosticIon(pub MolecularFormula);

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
    use crate::{fragment::NeutralLoss, molecular_formula, prelude::AminoAcid};

    #[test]
    fn deserialise_json() {
        let original = r#"{"Loss":{"elements":[["H",null,2],["O",null,1]],"additional_mass":0.0,"labels":[]}}"#;
        let new = r#"{"Loss":[1,{"elements":[["H",null,2],["O",null,1]],"additional_mass":0.0,"labels":[]}]}"#;
        let current = serde_json::to_string(&NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)))
            .expect("Could not serialise");

        let original_v =
            serde_json::from_str::<NeutralLoss>(original).expect("Could not deserialise original");
        let new_v = serde_json::from_str::<NeutralLoss>(new).expect("Could not deserialise new");
        let current_v =
            serde_json::from_str::<NeutralLoss>(&current).expect("Could not deserialise current");

        assert_eq!(original_v, new_v);
        assert_eq!(current_v, new_v);

        let side_chain_loss =
            NeutralLoss::SideChainLoss(molecular_formula!(H 2 O 1), AminoAcid::Alanine);
        let side_chain_loss_json =
            serde_json::to_string(&side_chain_loss).expect("Could not serialise side chain loss");
        println!("{side_chain_loss_json}");
        let side_chain_loss_back = serde_json::from_str::<NeutralLoss>(&side_chain_loss_json)
            .expect("Could not deserialise side chain loss");

        assert_eq!(side_chain_loss, side_chain_loss_back);
    }
}
