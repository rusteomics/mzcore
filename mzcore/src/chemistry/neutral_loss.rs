use std::{
    borrow::Cow,
    fmt::Display,
    ops::{Add, AddAssign},
    str::FromStr,
};

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    chemistry::MolecularFormula,
    helper_functions::{explain_number_error, next_number},
    parse_json::{ParseJson, use_serde},
    quantities::Multi,
    sequence::AminoAcid,
    space::{Space, UsedSpace},
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

impl Space for NeutralLoss {
    fn space(&self) -> UsedSpace {
        (UsedSpace::stack(1)
            + match self {
                Self::Gain(a, f) | Self::Loss(a, f) => a.space() + f.space(),
                Self::SideChainLoss(f, a) => f.space() + a.space(),
            })
        .set_total::<Self>()
    }
}

impl ParseJson for NeutralLoss {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        let parse_inner = |value: Value, context: &str| {
            if let Value::Array(mut arr) = value {
                if arr.len() == 2 {
                    if let Some(n) = arr[0].as_u64() {
                        Ok((
                            u16::try_from(n).map_err(|_| BoxedError::new(BasicKind::Error,
                                "Invalid NeutralLoss",
                                format!("The {context} amount is too big, the number has to be below {}", u16::MAX),
                                Context::show(n.to_string())))?,
                            MolecularFormula::from_json_value(arr.pop().unwrap())?,
                        ))
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid NeutralLoss",
                            format!("The {context} amount is not a number"),
                            Context::show(arr[0].to_string()),
                        ))
                    }
                } else {
                    Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid NeutralLoss",
                        format!("The {context} is a sequence but does not have 2 children"),
                        Context::show(arr.iter().join(",")),
                    ))
                }
            } else if value.is_object() {
                Ok((1, MolecularFormula::from_json_value(value)?))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid NeutralLoss",
                    format!("The {context} has to be either a map or a sequence"),
                    Context::show(value.to_string()),
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
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid NeutralLoss",
                                "The SideChainLoss is a sequence but does not have 2 children",
                                Context::show(arr.iter().join(",")),
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid NeutralLoss",
                            "The SideChainLoss has to be a sequence",
                            Context::show(value.to_string()),
                        ))
                    }
                }
                _ => Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid NeutralLoss",
                    "The tag has to be Loss/Gain/SideChainLoss",
                    Context::show(key),
                )),
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid NeutralLoss",
                "The JSON value has to be a map",
                Context::show(value.to_string()),
            ))
        }
    }
}

/// A diagnostic ion, defined in M (not MH+) chemical formula
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct DiagnosticIon(pub MolecularFormula);

impl Space for DiagnosticIon {
    fn space(&self) -> UsedSpace {
        self.0.space()
    }
}

impl ParseJson for DiagnosticIon {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
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
            Self::Loss(n, formula) => {
                notation_helper('-', *n, formula, MolecularFormula::hill_notation_html, true)
            }
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, formula) => {
                notation_helper('+', *n, formula, MolecularFormula::hill_notation_html, true)
            }
        }
    }

    /// Generate a nice fancy notation for this `NeutralLoss`
    pub fn hill_notation_fancy(&self) -> String {
        match self {
            Self::Loss(n, formula) => notation_helper(
                '-',
                *n,
                formula,
                MolecularFormula::hill_notation_fancy,
                true,
            ),
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, formula) => notation_helper(
                '+',
                *n,
                formula,
                MolecularFormula::hill_notation_fancy,
                true,
            ),
        }
    }

    /// Generate a notation for this `NeutralLoss` with pure ASCII characters
    pub fn hill_notation(&self) -> String {
        match self {
            Self::Loss(n, formula) => {
                notation_helper('-', *n, formula, MolecularFormula::hill_notation, false)
            }
            Self::SideChainLoss(_, aa) => format!("-sidechain_{aa}"),
            Self::Gain(n, formula) => {
                notation_helper('+', *n, formula, MolecularFormula::hill_notation, false)
            }
        }
    }
}

fn notation_helper(
    sign: char,
    n: u16,
    formula: &MolecularFormula,
    f: impl Fn(&MolecularFormula) -> String,
    fancy: bool,
) -> String {
    if formula.elements().is_empty() && n != 1 {
        format!(
            "{sign}{n}{}{}",
            if fancy { '×' } else { 'x' },
            f(formula).trim_start_matches('+')
        )
    } else if n == 1 {
        format!("{sign}{}", f(formula).trim_start_matches('+'))
    } else {
        format!("{sign}{n}{}", f(formula).trim_start_matches('+'))
    }
}

impl FromStr for NeutralLoss {
    type Err = BoxedError<'static, BasicKind>;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Ok(number) = s.parse::<f64>() {
            // Allow a simple numeric neutral loss
            if number >= 0.0 {
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
        } else if (s.contains('x') && !s.contains("xe")) || s.contains('×') {
            // Allow a factor times a numeric neutral loss
            match s.split_once('x').or_else(|| s.split_once('×')) {
                None => unreachable!(),
                Some((start, end)) => {
                    if start.starts_with('-') || start.starts_with('+') {
                        let amount = start.parse::<i32>().map_err(|error|BoxedError::new(BasicKind::Error,
                            "Invalid neutral loss",
                            "The text before the times symbol should be a valid number, like: `-1x12`",
                            Context::line_with_comment(None, s, 0, start.len(), Some(Cow::Owned(explain_number_error(&error).to_string()))).to_owned(),
                        ))?;
                        let mass = end.parse::<f64>().map_err(|error|BoxedError::new(BasicKind::Error,
                            "Invalid neutral loss",
                            "The text after the times symbol should be a valid number, like: `-1x12`",
                            Context::line_with_comment(None, s, 0, start.len(), Some(Cow::Owned(error.to_string()))).to_owned(),
                        ))?;
                        if amount >= 0 {
                            Ok(Self::Gain(
                                amount as u16,
                                MolecularFormula::with_additional_mass(mass),
                            ))
                        } else {
                            Ok(Self::Loss(
                                amount.unsigned_abs() as u16,
                                MolecularFormula::with_additional_mass(mass),
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid neutral loss",
                            "A neutral loss can only start with '+' or '-'",
                            Context::line(None, s, 0, 1).to_owned(),
                        ))
                    }
                }
            }
        } else if let Some(c) = s.chars().next() {
            // Or match a molecular formula
            let loss = match c {
                '+' => Ok(false),
                '-' => Ok(true),
                _ => Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid neutral loss",
                    "A neutral loss can only start with '+' or '-'",
                    Context::line(None, s, 0, 1).to_owned(),
                )),
            }?;
            let (amount, start) = if let Some(amount) = next_number::<false, false, u16>(s, 1..) {
                (
                    amount.2.map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid neutral loss",
                            format!("The amount specifier {}", explain_number_error(&err)),
                            Context::line(None, s, 1, amount.0).to_owned(),
                        )
                    })?,
                    amount.0 + 1,
                )
            } else {
                (1, 1)
            };
            let formula = MolecularFormula::pro_forma_inner::<false, false>(
                &Context::none().lines(0, s),
                s,
                start..,
            )
            .map_err(BoxedError::to_owned)?;
            Ok(if loss {
                Self::Loss(amount, formula)
            } else {
                Self::Gain(amount, formula)
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid neutral loss",
                "A neutral loss cannot be an empty string",
                Context::none(),
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
    use std::str::FromStr;

    use crate::{
        chemistry::MolecularFormula, chemistry::NeutralLoss, parse_json::ParseJson,
        sequence::AminoAcid,
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

    #[test]
    fn parse() {
        assert_eq!(
            NeutralLoss::from_str("-12"),
            Ok(NeutralLoss::Loss(
                1,
                MolecularFormula::with_additional_mass(12.0)
            ))
        );
        assert_eq!(
            NeutralLoss::from_str("-1x12"),
            Ok(NeutralLoss::Loss(
                1,
                MolecularFormula::with_additional_mass(12.0)
            ))
        );
        assert_eq!(
            NeutralLoss::from_str("-2x12"),
            Ok(NeutralLoss::Loss(
                2,
                MolecularFormula::with_additional_mass(12.0)
            ))
        );
        assert_eq!(
            NeutralLoss::from_str("-2×12"),
            Ok(NeutralLoss::Loss(
                2,
                MolecularFormula::with_additional_mass(12.0)
            ))
        );
        assert_eq!(
            NeutralLoss::from_str("-H2O"),
            Ok(NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)))
        );
        assert_eq!(
            NeutralLoss::from_str("-2H2O"),
            Ok(NeutralLoss::Loss(2, molecular_formula!(H 2 O 1)))
        );
    }
}
