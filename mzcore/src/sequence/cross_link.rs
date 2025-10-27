use std::{cmp::Ordering, collections::BTreeSet};

use context_error::*;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    parse_json::{ParseJson, use_serde},
    sequence::PlacementRule,
};

/// Indicate the cross-link side, it contains a set of all placement rules that apply for the placed
/// location to find all possible ways of breaking and/or neutral losses. These numbers are the
/// index into the [`LinkerSpecificity`] rules.
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
pub enum CrossLinkSide {
    /// The cross-link is symmetric, or if asymmetric it can be placed in both orientations
    Symmetric(BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'left' side
    Left(BTreeSet<usize>),
    /// The cross-link is asymmetric and this is the 'right' side
    Right(BTreeSet<usize>),
}

impl PartialOrd for CrossLinkSide {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for CrossLinkSide {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Symmetric(_), Self::Symmetric(_)) | (Self::Left(_), Self::Left(_)) => {
                Ordering::Equal
            }
            (Self::Symmetric(_), _) => Ordering::Greater,
            (_, Self::Symmetric(_)) => Ordering::Less,
            (Self::Left(_), _) => Ordering::Greater,
            (_, Self::Left(_)) => Ordering::Less,
            (Self::Right(_), Self::Right(_)) => Ordering::Equal,
        }
    }
}

impl std::hash::Hash for CrossLinkSide {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let (i, r) = match self {
            Self::Symmetric(r) => (0, r),
            Self::Left(r) => (1, r),
            Self::Right(r) => (2, r),
        };
        state.write_u8(i);
        state.write(
            &r.iter()
                .sorted()
                .flat_map(|r| r.to_ne_bytes())
                .collect_vec(),
        );
    }
}

impl ParseJson for CrossLinkSide {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

/// The name of a cross-link
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum CrossLinkName {
    /// A branch
    Branch,
    /// A cross-link
    Name(Box<str>),
}

impl ParseJson for CrossLinkName {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

/// The linker position specificities for a linker
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum LinkerSpecificity {
    /// A symmetric specificity where both ends have the same specificity.
    Symmetric {
        /// The placement rules for both ends.
        rules: Vec<PlacementRule>,
        /// All stubs that can be left after cleaving or breaking of the cross-link.
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        /// All possible neutral losses from the intact cross-linker.
        neutral_losses: Vec<NeutralLoss>,
        /// All diagnostic ions from the cross-linker
        diagnostic: Vec<DiagnosticIon>,
    },
    /// An asymmetric specificity where both ends have a different specificity.
    Asymmetric {
        /// The placement rules for both ends, these can be asymmetric thus are provided for 'right' and 'left' separately.
        rules: (Vec<PlacementRule>, Vec<PlacementRule>),
        /// All stubs that can be left after cleaving or breaking of the cross-link. The stubs are specific for right and left in the same orientation as the rules.
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        /// All possible neutral losses from the intact cross-linker.
        neutral_losses: Vec<NeutralLoss>,
        /// All diagnostic ions from the cross-linker
        diagnostic: Vec<DiagnosticIon>,
    },
}

impl ParseJson for LinkerSpecificity {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Object(map) = value {
            let (key, value) = map.into_iter().next().unwrap();
            match key.as_str() {
                "Symmetric" => {
                    if let Value::Object(mut map) = value {
                        Ok(Self::Symmetric {
                            rules: Vec::<PlacementRule>::from_json_value(
                                map.remove("rules").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'rules' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            stubs: Vec::<(MolecularFormula, MolecularFormula)>::from_json_value(
                                map.remove("stubs").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'stubs' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            neutral_losses: Vec::<NeutralLoss>::from_json_value(
                                map.remove("neutral_losses").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'neutral_losses' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            diagnostic: Vec::<DiagnosticIon>::from_json_value(
                                map.remove("diagnostic").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'diagnostic' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                        })
                    } else if let Value::Array(mut arr) = value {
                        if arr.len() == 3 {
                            let diagnostic =
                                Vec::<DiagnosticIon>::from_json_value(arr.pop().unwrap())?;
                            let stubs =
                                Vec::<(MolecularFormula, MolecularFormula)>::from_json_value(
                                    arr.pop().unwrap(),
                                )?;
                            let rules = Vec::<PlacementRule>::from_json_value(arr.pop().unwrap())?;

                            Ok(Self::Symmetric {
                                rules,
                                stubs,
                                neutral_losses: Vec::new(),
                                diagnostic,
                            })
                        } else {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid NeutralLoss",
                                "The Symmetric is a sequence but does not have 3 children",
                                Context::show(arr.iter().join(",")),
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid LinkerSpecificity",
                            "The Symmetric value has to be a map or a sequence",
                            Context::show(value.to_string()),
                        ))
                    }
                }
                "Asymmetric" => {
                    if let Value::Object(mut map) = value {
                        Ok(Self::Asymmetric {
                            rules: <(Vec<PlacementRule>, Vec<PlacementRule>)>::from_json_value(
                                map.remove("rules").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'rules' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            stubs: Vec::<(MolecularFormula, MolecularFormula)>::from_json_value(
                                map.remove("stubs").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'stubs' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            neutral_losses: Vec::<NeutralLoss>::from_json_value(
                                map.remove("neutral_losses").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'neutral_losses' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                            diagnostic: Vec::<DiagnosticIon>::from_json_value(
                                map.remove("diagnostic").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid LinkerSpecificity",
                                        "The required property 'diagnostic' is missing",
                                        Context::show(
                                            map.iter()
                                                .map(|(k, v)| format!("\"{k}\": {v}"))
                                                .join(","),
                                        ),
                                    )
                                })?,
                            )?,
                        })
                    } else if let Value::Array(mut arr) = value {
                        if arr.len() == 3 {
                            let diagnostic =
                                Vec::<DiagnosticIon>::from_json_value(arr.pop().unwrap())?;
                            let stubs =
                                Vec::<(MolecularFormula, MolecularFormula)>::from_json_value(
                                    arr.pop().unwrap(),
                                )?;
                            let rules =
                                <(Vec<PlacementRule>, Vec<PlacementRule>)>::from_json_value(
                                    arr.pop().unwrap(),
                                )?;

                            Ok(Self::Asymmetric {
                                rules,
                                stubs,
                                neutral_losses: Vec::new(),
                                diagnostic,
                            })
                        } else {
                            Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid NeutralLoss",
                                "The Asymmetric is a sequence but does not have 3 children",
                                Context::show(arr.iter().join(",")),
                            ))
                        }
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid LinkerSpecificity",
                            "The Asymmetric value has to be a map or a sequence",
                            Context::show(value.to_string()),
                        ))
                    }
                }
                _ => Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid LinkerSpecificity",
                    "The tag has to be Symmetric/Asymmetric",
                    Context::show(key),
                )),
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid LinkerSpecificity",
                "The JSON value has to be a map",
                Context::show(value.to_string()),
            ))
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        chemistry::{DiagnosticIon, NeutralLoss},
        molecular_formula,
        parse_json::ParseJson,
        sequence::{AminoAcid, LinkerSpecificity, PlacementRule, Position},
    };

    #[test]
    fn deserialise_json() {
        let old = r#"{"Asymmetric":[[[{"AminoAcid":[["Selenocysteine"],"AnyCTerm"]},{"AminoAcid":[["GlutamicAcid"],"Anywhere"]}],[{"AminoAcid":[["Selenocysteine"],"AnyNTerm"]}]],[[{"elements":[["U",null,1]],"additional_mass":0.0},{"elements":[["U",null,1]],"additional_mass":0.0}]],[{"elements":[["Te",null,1]],"additional_mass":0.0},{"elements":[["Ne",null,1]],"additional_mass":0.0},{"elements":[["H",null,2],["He",null,3]],"additional_mass":0.0},{"elements":[["H",null,1],["He",null,2]],"additional_mass":0.0},{"elements":[["I",null,1],["Er",null,1]],"additional_mass":0.0},{"elements":[["H",null,12],["C",null,12],["O",null,1]],"additional_mass":0.0}]]}"#;
        let current = LinkerSpecificity::Symmetric {
            rules: vec![
                PlacementRule::Terminal(Position::AnyNTerm),
                PlacementRule::AminoAcid(
                    vec![AminoAcid::Alanine, AminoAcid::Leucine],
                    Position::Anywhere,
                ),
            ],
            stubs: vec![(molecular_formula!(H 2 O 1), molecular_formula!(H 2 O 1))],
            neutral_losses: vec![NeutralLoss::Gain(2, molecular_formula!(H 2 O 1))],
            diagnostic: vec![DiagnosticIon(molecular_formula!(H 2 O 1))],
        };
        let current_text =
            serde_json::to_string(&current).expect("Could not serialise linker specificity");

        let _old = LinkerSpecificity::from_json(old).expect("Could not deserialise old json");
        let current_back = LinkerSpecificity::from_json(&current_text)
            .expect("Could not deserialise current json");

        assert_eq!(current, current_back);
    }
}
