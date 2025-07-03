use itertools::Itertools;
use serde::{
    Deserialize, Serialize,
    de::{Error, Visitor},
};

use std::{cmp::Ordering, collections::BTreeSet};

use crate::{chemistry::MolecularFormula, fragment::*, sequence::PlacementRule};

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

/// The name of a cross-link
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum CrossLinkName {
    /// A branch
    Branch,
    /// A cross-link
    Name(String),
}

/// The linker position specificities for a linker
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
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

/// Possible internal values, needed to allow serde to still parse the older format which did not contain a neutral loss value
#[derive(Debug, Deserialize)]
#[serde(untagged)]
enum Field {
    SymmetricOld(
        Vec<PlacementRule>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ),
    SymmetricNew {
        rules: Vec<PlacementRule>,
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        neutral_losses: Vec<NeutralLoss>,
        diagnostic: Vec<DiagnosticIon>,
    },
    AsymmetricOld(
        (Vec<PlacementRule>, Vec<PlacementRule>),
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ),
    AsymmetricNew {
        rules: (Vec<PlacementRule>, Vec<PlacementRule>),
        stubs: Vec<(MolecularFormula, MolecularFormula)>,
        neutral_losses: Vec<NeutralLoss>,
        diagnostic: Vec<DiagnosticIon>,
    },
}

impl<'de> Deserialize<'de> for LinkerSpecificity {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        deserializer.deserialize_map(LinkerSpecificityVisitor)
    }
}

struct LinkerSpecificityVisitor;

impl<'de> Visitor<'de> for LinkerSpecificityVisitor {
    type Value = LinkerSpecificity;

    fn expecting(&self, formatter: &mut std::fmt::Formatter) -> std::fmt::Result {
        formatter.write_str("a neutral loss")
    }

    fn visit_map<A>(self, mut map: A) -> Result<Self::Value, A::Error>
    where
        A: serde::de::MapAccess<'de>,
    {
        if let Some((key, value)) = map.next_entry::<&str, Field>()? {
            match key {
                "Symmetric" => match value {
                    Field::SymmetricOld(rules, stubs, diagnostic) => {
                        Ok(LinkerSpecificity::Symmetric {
                            rules,
                            stubs,
                            neutral_losses: Vec::new(),
                            diagnostic,
                        })
                    }
                    Field::SymmetricNew {
                        rules,
                        stubs,
                        neutral_losses,
                        diagnostic,
                    } => Ok(LinkerSpecificity::Symmetric {
                        rules,
                        stubs,
                        neutral_losses,
                        diagnostic,
                    }),
                    _ => Err(A::Error::custom(
                        "a symmetric linker specifity cannot contain asymmetric information",
                    )),
                },
                "Asymmetric" => match value {
                    Field::AsymmetricOld(rules, stubs, diagnostic) => {
                        Ok(LinkerSpecificity::Asymmetric {
                            rules,
                            stubs,
                            neutral_losses: Vec::new(),
                            diagnostic,
                        })
                    }
                    Field::AsymmetricNew {
                        rules,
                        stubs,
                        neutral_losses,
                        diagnostic,
                    } => Ok(LinkerSpecificity::Asymmetric {
                        rules,
                        stubs,
                        neutral_losses,
                        diagnostic,
                    }),
                    _ => Err(A::Error::custom(
                        "an asymmetric linker specifity cannot contain symmetric information",
                    )),
                },
                v => Err(A::Error::custom(format!(
                    "expected Symmetric/Asymmetric not '{v}'"
                ))),
            }
        } else {
            Err(A::Error::custom("Empty linker specificity definition"))
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        fragment::{DiagnosticIon, NeutralLoss},
        molecular_formula,
        prelude::AminoAcid,
        sequence::{LinkerSpecificity, PlacementRule, Position},
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

        let _old =
            serde_json::from_str::<LinkerSpecificity>(old).expect("Could not deserialise old json");
        let current_back = serde_json::from_str::<LinkerSpecificity>(&current_text)
            .expect("Could not deserialise current json");

        assert_eq!(current, current_back);
    }
}
