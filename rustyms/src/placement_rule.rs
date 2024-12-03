//! Rules regarding the placement of modifications

use std::str::FromStr;

use serde::{Deserialize, Serialize};

use crate::{
    error::{Context, CustomError},
    modification::{Modification, ModificationId, Ontology, SimpleModificationInner},
    AminoAcid, SequenceElement, SequencePosition,
};

include!("shared/placement_rule.rs");

impl PlacementRule {
    /// Check if this rule fits with the given location
    pub fn is_possible<T>(&self, seq: &SequenceElement<T>, position: SequencePosition) -> bool {
        match self {
            Self::AminoAcid(aa, r_pos) => {
                aa.iter().any(|a| *a == seq.aminoacid.aminoacid()) && r_pos.is_possible(position)
            }
            Self::PsiModification(mod_index, r_pos) => {
                seq.modifications.iter().any(|m| {
                    if let Modification::Simple(sim) = m {
                        if let SimpleModificationInner::Database {
                            id:
                                ModificationId {
                                    ontology: Ontology::Psimod,
                                    id,
                                    ..
                                },
                            ..
                        } = **sim
                        {
                            id.is_some_and(|i| i == *mod_index)
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                }) && r_pos.is_possible(position)
            }
            Self::Terminal(r_pos) => {
                r_pos.is_possible(position)
                    && (position == SequencePosition::NTerm || position == SequencePosition::CTerm)
            }
            Self::Anywhere => true,
        }
    }

    /// Check if this rule fits with the given location
    pub fn is_possible_aa(&self, aa: AminoAcid, position: Position) -> bool {
        match self {
            Self::AminoAcid(allowed_aa, r_pos) => {
                allowed_aa.iter().any(|a| *a == aa) && r_pos.is_possible_position(position)
            }
            Self::PsiModification(_, _) => false,
            Self::Terminal(r_pos) => {
                r_pos.is_possible_position(position) && (position != Position::Anywhere)
            }
            Self::Anywhere => true,
        }
    }

    /// Check if any of the given rules are possible
    pub fn any_possible<T>(
        rules: &[Self],
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> bool {
        rules.iter().any(|r| r.is_possible(seq, position))
    }

    /// Check if any of the given rules are possible
    pub fn any_possible_aa(rules: &[Self], aa: AminoAcid, position: Position) -> bool {
        rules.iter().any(|r| r.is_possible_aa(aa, position))
    }
}

impl FromStr for PlacementRule {
    type Err = CustomError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((head, tail)) = s.split_once('@') {
            let aa: Vec<AminoAcid> = head
                .chars()
                .enumerate()
                .map(|(i, c)| {
                    AminoAcid::try_from(c).map_err(|()| {
                        CustomError::error(
                            "Invalid amino acid",
                            "Invalid amino acid in specified amino acids in placement rule",
                            Context::line(None, s, i, 1),
                        )
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            tail.parse().map_or_else(
                |()| {
                    Err(CustomError::error(
                        "Invalid position",
                        "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                        Context::line(None, s, head.len() + 1, tail.len()),
                    ))
                },
                |position| Ok(Self::AminoAcid(aa, position)),
            )
        } else if let Ok(position) = s.parse() {
            Ok(match position {
                Position::Anywhere => Self::Anywhere,
                pos => Self::Terminal(pos),
            })
        } else {
            Err(CustomError::error(
                "Invalid position",
                "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                Context::full_line(0, s),
            ))
        }
    }
}

impl Position {
    /// See if the given peptide position is a valid position given this [`Position`] as placement rule.
    pub fn is_possible(self, position: SequencePosition) -> bool {
        match self {
            Self::Anywhere => true,
            Self::AnyNTerm | Self::ProteinNTerm => position == SequencePosition::NTerm,
            Self::AnyCTerm | Self::ProteinCTerm => position == SequencePosition::CTerm,
        }
    }

    /// See if the given position is a valid position given this [`Position`] as placement rule.
    pub fn is_possible_position(self, position: Self) -> bool {
        match self {
            Self::Anywhere => true,
            Self::AnyNTerm => position == Self::AnyNTerm || position == Self::ProteinNTerm,
            Self::ProteinNTerm => position == Self::ProteinNTerm,
            Self::AnyCTerm => position == Self::AnyCTerm || position == Self::ProteinCTerm,
            Self::ProteinCTerm => position == Self::ProteinCTerm,
        }
    }
}

impl FromStr for Position {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "anywhere" => Ok(Self::Anywhere),
            "anynterm" => Ok(Self::AnyNTerm),
            "proteinnterm" => Ok(Self::ProteinNTerm),
            "anycterm" => Ok(Self::AnyCTerm),
            "proteincterm" => Ok(Self::ProteinCTerm),
            _ => Err(()),
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {

    use crate::{checked_aminoacid::CheckedAminoAcid, modification::RulePossible};

    use super::*;
    #[test]
    fn multi_level_rule() {
        assert!(
            !PlacementRule::PsiModification(30, Position::Anywhere).is_possible(
                &SequenceElement::new(CheckedAminoAcid::Alanine, None),
                SequencePosition::Index(0)
            ),
            "Multi level mod cannot be placed if the dependent mod is not present"
        );
        let mut seq = SequenceElement::new(CheckedAminoAcid::Alanine, None);
        seq.modifications
            .push(Ontology::Psimod.find_id(30, None).unwrap().into());
        assert!(
            PlacementRule::PsiModification(30, Position::Anywhere)
                .is_possible(&seq, SequencePosition::Index(0)),
            "Multi level mod can be placed if the dependent mod is present"
        );
    }

    #[test]
    fn place_anywhere() {
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Glutamine], Position::Anywhere).is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                crate::SequencePosition::NTerm
            ),
            "start"
        );
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Glutamine], Position::Anywhere).is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                crate::SequencePosition::Index(2)
            ),
            "middle"
        );
        assert!(
            PlacementRule::AminoAcid(vec![AminoAcid::Glutamine], Position::Anywhere).is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                crate::SequencePosition::CTerm
            ),
            "end"
        );
        assert_eq!(
            Ontology::Unimod.find_id(7, None).unwrap().is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                crate::SequencePosition::CTerm
            ),
            RulePossible::Symmetric(std::collections::BTreeSet::from([0])),
            "unimod deamidated at end"
        );
    }
}
