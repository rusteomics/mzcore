use std::{fmt::Display, ops::Add, str::FromStr};

use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    formula::MolecularFormula,
    Multi,
};

include!("shared/neutral_loss.rs");

impl NeutralLoss {
    /// Check if this neutral loss if empty (has an empty molecular formula)
    pub fn is_empty(&self) -> bool {
        match self {
            Self::Loss(f) | Self::Gain(f) => f.is_empty(),
        }
    }

    /// Generate a nice HTML notation for this `NeutralLoss`
    pub fn hill_notation_html(&self) -> String {
        match self {
            Self::Loss(c) => format!("-{}", c.hill_notation_html().trim_start_matches('+')),
            Self::Gain(c) => format!("+{}", c.hill_notation_html().trim_start_matches('+')),
        }
    }

    /// Generate a nice fancy notation for this `NeutralLoss`
    pub fn hill_notation_fancy(&self) -> String {
        match self {
            Self::Loss(c) => format!("-{}", c.hill_notation_fancy().trim_start_matches('+')),
            Self::Gain(c) => format!("+{}", c.hill_notation_fancy().trim_start_matches('+')),
        }
    }

    /// Generate a notation for this `NeutralLoss` with pure ASCII characters
    pub fn hill_notation(&self) -> String {
        match self {
            Self::Loss(c) => format!("-{}", c.hill_notation().trim_start_matches('+')),
            Self::Gain(c) => format!("+{}", c.hill_notation().trim_start_matches('+')),
        }
    }
}

impl FromStr for NeutralLoss {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        // Allow a numeric neutral loss
        if let Ok(number) = s.parse::<f64>() {
            if number > 0.0 {
                Ok(Self::Gain(MolecularFormula::with_additional_mass(number)))
            } else {
                Ok(Self::Loss(MolecularFormula::with_additional_mass(
                    number.abs(),
                )))
            }
        } else if let Some(c) = s.chars().next() {
            // Or match a molecular formula
            match c {
                '+' => Ok(Self::Gain(MolecularFormula::from_pro_forma(
                    s,
                    1..,
                    false,
                    false,
                    true,
                )?)),
                '-' => Ok(Self::Loss(MolecularFormula::from_pro_forma(
                    s,
                    1..,
                    false,
                    false,
                    true,
                )?)),
                _ => Err(CustomError::error(
                    "Invalid neutral loss",
                    "A neutral loss can only start with '+' or '-'",
                    Context::line(None, s, 0, 1),
                )),
            }
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
        write!(
            f,
            "{}",
            match self {
                Self::Loss(c) => format!("-{c}"),
                Self::Gain(c) => format!("+{c}"),
            }
        )
    }
}

impl std::ops::Add<&NeutralLoss> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &NeutralLoss) -> Self::Output {
        match rhs {
            NeutralLoss::Gain(mol) => self + mol,
            NeutralLoss::Loss(mol) => self - mol,
        }
    }
}

impl std::ops::Add<&NeutralLoss> for &Multi<MolecularFormula> {
    type Output = Multi<MolecularFormula>;
    fn add(self, rhs: &NeutralLoss) -> Self::Output {
        match rhs {
            NeutralLoss::Gain(mol) => self + mol,
            NeutralLoss::Loss(mol) => self - mol,
        }
    }
}

impl_binop_ref_cases!(impl Add, add for MolecularFormula, NeutralLoss, MolecularFormula);
impl_binop_ref_cases!(impl Add, add for Multi<MolecularFormula>, NeutralLoss, Multi<MolecularFormula>);
