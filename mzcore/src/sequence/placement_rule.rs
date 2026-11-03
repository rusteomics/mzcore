//! Rules regarding the placement of modifications

use std::str::FromStr;

use context_error::*;
use mzcv::AccessionCode;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    ontology::Ontology,
    parse_json::{ParseJson, use_serde},
    sequence::{
        AminoAcid, Modification, SequenceElement, SequencePosition, SimpleModificationInner,
    },
    space::{Space, UsedSpace},
};

/// A rule determining the placement of a modification
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum PlacementRule {
    /// Placed on an aminoacid on the given position, if the list is empty it matches to nothing
    AminoAcid(ThinVec<AminoAcid>, ResiduePosition, ModificationPosition),
    /// Placed on another modification on the given position
    PsiModification(u32, ResiduePosition, ModificationPosition),
    /// Placed on a specific position regardless of the amino acid there
    Position(ResiduePosition, ModificationPosition),
}

impl Space for PlacementRule {
    fn space(&self) -> UsedSpace {
        (UsedSpace::stack(1)
            + match self {
                Self::AminoAcid(a, p, m) => a.space() + p.space() + m.space(),
                Self::PsiModification(i, p, m) => i.space() + p.space() + m.space(),
                Self::Position(p, m) => p.space() + m.space(),
            })
        .set_total::<Self>()
    }
}

impl ParseJson for PlacementRule {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let serde_json::Value::Object(map) = dbg!(&value) {
            if let Some(v) = map.get("Terminal").or_else(|| map.get("Position")) {
                if let serde_json::Value::String(v) = v {
                    let p = ResiduePosition::from_str(v).map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid Position",
                            "It should be 'Anywhere', 'AnyNTerm', 'ProteinNTerm', 'AnyCTerm', or 'ProteinCTerm'",
                            Context::default().lines(0, value.to_string()),
                        )
                    })?;
                    Ok(Self::Position(
                        p,
                        if p == ResiduePosition::Anywhere {
                            ModificationPosition::Sidechain
                        } else {
                            ModificationPosition::Terminal
                        },
                    ))
                } else {
                    <(ResiduePosition, ModificationPosition)>::from_json_value(v.clone())
                        .map(|v| Self::Position(v.0, v.1))
                }
            } else if let Some(v) = map.get("PsiModification") {
                <(u32, ResiduePosition, ModificationPosition)>::from_json_value(v.clone())
                    .map(|v| Self::PsiModification(v.0, v.1, v.2))
                    .or_else(|err| {
                        <(u32, ResiduePosition)>::from_json_value(v.clone())
                            .map(|v| {
                                Self::PsiModification(v.0, v.1, ModificationPosition::Sidechain)
                            })
                            .map_err(|_| err)
                    })
            } else if let Some(v) = map.get("AminoAcid") {
                <(ThinVec<AminoAcid>, ResiduePosition, ModificationPosition)>::from_json_value(
                    v.clone(),
                )
                .map(|v| Self::AminoAcid(v.0, v.1, v.2))
                .or_else(|err| {
                    <(ThinVec<AminoAcid>, ResiduePosition)>::from_json_value(v.clone())
                        .map(|v| Self::AminoAcid(v.0, v.1, ModificationPosition::Sidechain))
                        .map_err(|_| err)
                })
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid PlacementRule",
                    "It should be a 'Position', 'AminoAcid', or 'PsiModification' object",
                    Context::default().lines(0, value.to_string()),
                ))
            }
        } else if let serde_json::Value::String(str) = &value
            && str.eq_ignore_ascii_case("Anywhere")
        {
            Ok(Self::Position(
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain,
            ))
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid PlacementRule",
                "It should be an object",
                Context::default().lines(0, value.to_string()),
            ))
        }
    }
}

/// A position where a modification can be placed
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum ResiduePosition {
    /// At any location
    #[default]
    Anywhere,
    /// At the N term of a peptide or protein
    AnyNTerm,
    /// At the C term of a peptide or protein
    AnyCTerm,
    /// At the N term of a protein
    ProteinNTerm,
    /// At the C term of a protein
    ProteinCTerm,
}

impl ParseJson for ResiduePosition {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

impl std::fmt::Display for ResiduePosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match self {
            Self::Anywhere => "Anywhere",
            Self::AnyNTerm => "AnyNTerm",
            Self::AnyCTerm => "AnyCTerm",
            Self::ProteinNTerm => "ProteinNTerm",
            Self::ProteinCTerm => "ProteinCTerm",
        })
    }
}

impl ResiduePosition {
    /// Check if this position is a subset of another position
    pub const fn is_subset(self, other: Self) -> bool {
        matches!(
            (self, other),
            (_, Self::Anywhere)
                | (Self::ProteinCTerm, Self::ProteinCTerm | Self::AnyCTerm)
                | (Self::AnyCTerm, Self::AnyCTerm)
                | (Self::ProteinNTerm, Self::ProteinNTerm | Self::AnyNTerm)
                | (Self::AnyNTerm, Self::AnyNTerm)
        )
    }
}

impl FromStr for ResiduePosition {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "" | "anywhere" => Ok(Self::Anywhere),
            "anynterm" | "anyn-term" | "any n-term" => Ok(Self::AnyNTerm),
            "proteinnterm" | "proteinn-term" | "protein n-term" => Ok(Self::ProteinNTerm),
            "anycterm" | "anyc-term" | "any c-term" => Ok(Self::AnyCTerm),
            "proteincterm" | "proteinc-term" | "protein c-term" => Ok(Self::ProteinCTerm),
            _ => Err(()),
        }
    }
}

/// A position where a modification can be placed in the residue itself. This determines if a
/// modification is considered wrong when placed with a terminal modification `[MOD]-AA` or a
/// sidechain modification `A[MOD]A` in ProForma.
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum ModificationPosition {
    /// At or in the sidechain
    #[default]
    Sidechain,
    /// At the terminal, naturally only possible at the terminal of a peptide or protein
    Terminal,
}

impl ParseJson for ModificationPosition {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

impl std::fmt::Display for ModificationPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", match self {
            Self::Sidechain => "Sidechain",
            Self::Terminal => "Terminal",
        })
    }
}

impl FromStr for ModificationPosition {
    type Err = ();

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "sidechain" => Ok(Self::Sidechain),
            "terminal" | "term" | "terminus" => Ok(Self::Terminal),
            _ => Err(()),
        }
    }
}

impl PlacementRule {
    /// Check if this rule fits with the given location
    pub fn is_possible<T>(&self, seq: &SequenceElement<T>, position: SequencePosition) -> bool {
        match self {
            Self::AminoAcid(aa, r_pos, m_pos) => {
                aa.iter().any(|a| *a == seq.aminoacid.aminoacid())
                    && RulePosition::residue_is_possible(*r_pos, position)
                    && RulePosition::modification_is_possible(*m_pos, position)
            }
            Self::PsiModification(mod_index, r_pos, m_pos) => {
                seq.modifications.iter().any(|m| {
                    if let Modification::Simple(sim) = m {
                        if let SimpleModificationInner::Database { id, .. } = &**sim
                            && id.ontology == Ontology::Psimod
                            && let AccessionCode::Numeric(i) = id.id()
                        {
                            i == *mod_index
                        } else {
                            false
                        }
                    } else {
                        false
                    }
                }) && RulePosition::residue_is_possible(*r_pos, position)
                    && RulePosition::modification_is_possible(*m_pos, position)
            }
            Self::Position(r_pos, m_pos) => {
                RulePosition::residue_is_possible(*r_pos, position)
                    && RulePosition::modification_is_possible(*m_pos, position)
            }
        }
    }

    /// Check if this rule fits with the given location
    pub fn is_possible_aa<P: RulePosition>(&self, aa: AminoAcid, position: P) -> bool {
        match self {
            Self::AminoAcid(allowed_aa, r_pos, m_pos) => {
                allowed_aa.contains(&aa)
                    && RulePosition::residue_is_possible(*r_pos, position)
                    && RulePosition::modification_is_possible(*m_pos, position)
            }
            Self::PsiModification(..) => false,
            Self::Position(r_pos, m_pos) => {
                RulePosition::residue_is_possible(*r_pos, position)
                    && RulePosition::modification_is_possible(*m_pos, position)
            }
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
    pub fn any_possible_aa<P: RulePosition + Clone>(
        rules: &[Self],
        aa: AminoAcid,
        position: P,
    ) -> bool {
        rules.iter().any(|r| r.is_possible_aa(aa, position.clone()))
    }

    /// Check if this rule is a subset of another rule
    pub fn is_subset(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Position(ps, ms), Self::Position(po, mo)) => ms == mo && ps.is_subset(*po),
            (Self::AminoAcid(aas, ps, ms), Self::AminoAcid(aao, po, mo)) => {
                ms == mo && ps.is_subset(*po) && aas.iter().all(|a| aao.contains(a))
            }
            (Self::PsiModification(is, ps, ms), Self::PsiModification(io, po, mo)) => {
                ms == mo && is == io && ps.is_subset(*po)
            }
            _ => false,
        }
    }

    /// Combine these two rules into the first. If they cannot be combined this returns false.
    pub fn combine_rules(&mut self, other: &Self) -> bool {
        match (self, other) {
            (Self::AminoAcid(new_aa, new_pos, new_m), Self::AminoAcid(aa, pos, m))
                if pos == new_pos && new_m == m =>
            {
                for a in aa {
                    if !new_aa.contains(a) {
                        new_aa.push(*a);
                    }
                }
                new_aa.sort_unstable();
                true
            }
            (Self::Position(new_pos, new_m), Self::Position(pos, m)) if new_m == m => {
                if new_pos == pos {
                    true
                } else if matches!(
                    new_pos,
                    ResiduePosition::ProteinNTerm | ResiduePosition::AnyNTerm
                ) && matches!(
                    pos,
                    ResiduePosition::ProteinNTerm | ResiduePosition::AnyNTerm
                ) {
                    *new_pos = ResiduePosition::AnyNTerm;
                    true
                } else if matches!(
                    new_pos,
                    ResiduePosition::ProteinCTerm | ResiduePosition::AnyCTerm
                ) && matches!(
                    pos,
                    ResiduePosition::ProteinCTerm | ResiduePosition::AnyCTerm
                ) {
                    *new_pos = ResiduePosition::AnyCTerm;
                    true
                } else {
                    false
                }
            }
            (a, b) if a == b => true,
            _ => false,
        }
    }
}

impl FromStr for PlacementRule {
    type Err = BoxedError<'static, BasicKind>;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let parts = s.split('@').collect::<Vec<_>>();
        let parse_aas = |range: std::ops::Range<usize>| {
            s[range.clone()]
                .char_indices()
                .map(|(i, c)| {
                    AminoAcid::try_from(c).map_err(|()| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid amino acid",
                            "Invalid amino acid in specified amino acids in placement rule",
                            Context::default()
                                .lines(0, s)
                                .add_highlight((0, range.start + i, 1))
                                .to_owned(),
                        )
                    })
                })
                .collect::<Result<Vec<AminoAcid>, _>>()
                .map(ThinVec::from)
        };

        match parts.as_slice() {
            &[aa, mpos, pos] => {
                let aa_v = parse_aas(0..aa.len())?;
                let mpos_v = mpos.parse::<ModificationPosition>().map_err(|_| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid modification position",
                        "Use any of the following for the position: Sidechain, Terminal",
                        Context::default()
                            .line_index(0)
                            .lines(0, s)
                            .add_highlight((0, aa.len() + 1..aa.len() + 1 + mpos.len()))
                            .to_owned(),
                    )
                })?;
                let pos_v = pos.parse::<ResiduePosition>().map_err(|_|
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid position",
                        "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                        Context::default()
                            .line_index(0)
                            .lines(0, s)
                            .add_highlight((0, aa.len() + mpos.len() + 2..))
                            .to_owned(),
                    ))?;
                Ok(Self::AminoAcid(aa_v, pos_v, mpos_v))
            }
            &[one, two] => {
                // Either AA@pos or mpos@pos
                let p = two.parse::<ResiduePosition>().map_err(|_|
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid position",
                        "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                        Context::default()
                            .line_index(0)
                            .lines(0, s)
                            .add_highlight((0, one.len() + 1..))
                            .to_owned(),
                    ))?;
                if let Ok(mpos) = one.parse::<ModificationPosition>() {
                    Ok(Self::Position(p, mpos))
                } else {
                    Ok(Self::AminoAcid(
                        parse_aas(0..one.len())?,
                        p,
                        ModificationPosition::Sidechain,
                    ))
                }
            }
            &[_] => {
                if let Ok(position) = s.parse() {
                    Ok(Self::Position(
                        position,
                        if position == ResiduePosition::Anywhere {
                            ModificationPosition::Sidechain
                        } else {
                            ModificationPosition::Terminal
                        },
                    ))
                } else {
                    Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid position",
                        "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                        Context::default().line_index(0).lines(0, s).to_owned(),
                    ))
                }
            }
            _ => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid position",
                "Use any of the following for the position: Anywhere, AnyNTerm, ProteinNTerm, AnyCTerm, ProteinCTerm",
                Context::default().line_index(0).lines(0, s).to_owned(),
            )),
        }
    }
}

/// Trait to be able to use any of the two position definitions when testing rules
pub trait RulePosition: Copy {
    /// Check if this position is allowed with this residue position rule
    fn residue_is_possible(rule: ResiduePosition, position: Self) -> bool;
    /// Check if this position is allowed with this modification position rule
    fn modification_is_possible(rule: ModificationPosition, position: Self) -> bool;
}

impl RulePosition for ResiduePosition {
    fn residue_is_possible(rule: ResiduePosition, position: Self) -> bool {
        match rule {
            Self::Anywhere => true,
            Self::AnyNTerm => position == Self::AnyNTerm || position == Self::ProteinNTerm,
            Self::ProteinNTerm => position == Self::ProteinNTerm,
            Self::AnyCTerm => position == Self::AnyCTerm || position == Self::ProteinCTerm,
            Self::ProteinCTerm => position == Self::ProteinCTerm,
        }
    }

    fn modification_is_possible(_rule: ModificationPosition, _position: Self) -> bool {
        true
    }
}

impl RulePosition for SequencePosition {
    fn residue_is_possible(rule: ResiduePosition, position: Self) -> bool {
        match rule {
            ResiduePosition::Anywhere => true,
            ResiduePosition::AnyNTerm | ResiduePosition::ProteinNTerm => {
                position == Self::NTerm || matches!(position, Self::Index(0, _))
            }
            ResiduePosition::AnyCTerm | ResiduePosition::ProteinCTerm => {
                position == Self::CTerm
                    || matches!(position, Self::Index(i, l) if i == l.saturating_sub(1))
            }
        }
    }

    fn modification_is_possible(rule: ModificationPosition, position: Self) -> bool {
        match rule {
            ModificationPosition::Sidechain => matches!(position, Self::Index(_, _)),
            ModificationPosition::Terminal => position == Self::NTerm || position == Self::CTerm,
        }
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {

    use super::*;
    use crate::{
        ontology::STATIC_ONTOLOGIES,
        sequence::{CheckedAminoAcid, RulePossible},
    };
    #[test]
    fn multi_level_rule() {
        let ontologies = &STATIC_ONTOLOGIES;
        assert!(
            !PlacementRule::PsiModification(
                30,
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            )
            .is_possible(
                &SequenceElement::new(CheckedAminoAcid::Alanine, None),
                SequencePosition::Index(0, 5),
            ),
            "Multi level mod cannot be placed if the dependent mod is not present"
        );
        let mut seq = SequenceElement::new(CheckedAminoAcid::Alanine, None);
        seq.modifications.push(
            ontologies
                .psimod()
                .get_by_index(&AccessionCode::Numeric(30))
                .unwrap()
                .into(),
        );
        assert!(
            PlacementRule::PsiModification(
                30,
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            )
            .is_possible(&seq, SequencePosition::Index(0, 5),),
            "Multi level mod can be placed if the dependent mod is present"
        );
    }

    #[test]
    fn place_anywhere() {
        assert!(
            // TODO: needs to be able the tell apart side chain and terminus for N and C terminal
            // AAs
            PlacementRule::AminoAcid(
                vec![AminoAcid::Glutamine].into(),
                ResiduePosition::ProteinNTerm,
                ModificationPosition::Sidechain
            )
            .is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                SequencePosition::Index(0, 5),
            ),
            "start"
        );
        assert!(
            PlacementRule::AminoAcid(
                vec![AminoAcid::Glutamine].into(),
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            )
            .is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                SequencePosition::Index(2, 5),
            ),
            "middle"
        );
        assert!(
            PlacementRule::AminoAcid(
                vec![AminoAcid::Glutamine].into(),
                ResiduePosition::ProteinCTerm,
                ModificationPosition::Sidechain
            )
            .is_possible(
                &SequenceElement::new(CheckedAminoAcid::Q, None),
                SequencePosition::Index(4, 5),
            ),
            "end"
        );
        assert_eq!(
            STATIC_ONTOLOGIES
                .unimod()
                .get_by_index(&AccessionCode::Numeric(7))
                .unwrap()
                .is_possible(
                    &SequenceElement::new(CheckedAminoAcid::Q, None),
                    SequencePosition::Index(4, 5),
                ),
            RulePossible::Symmetric(std::collections::BTreeSet::from([0])),
            "unimod deamidated at end"
        );
    }

    #[test]
    fn parse_json() {
        assert_eq!(
            PlacementRule::from_json_value(serde_json::Value::String("Anywhere".to_string())),
            Ok(PlacementRule::Position(
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            ))
        );
        assert_eq!(
            PlacementRule::from_json_value(serde_json::Value::Object(serde_json::Map::from_iter(
                [(
                    "Position".to_string(),
                    serde_json::Value::String("Anywhere".to_string())
                )]
            ))),
            Ok(PlacementRule::Position(
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            ))
        );
        assert_eq!(
            PlacementRule::from_json_value(serde_json::Value::Object(serde_json::Map::from_iter(
                [(
                    "Terminal".to_string(),
                    serde_json::Value::String("AnyNTerm".to_string())
                )]
            ))),
            Ok(PlacementRule::Position(
                ResiduePosition::AnyNTerm,
                ModificationPosition::Terminal
            ))
        );
        assert_eq!(
            PlacementRule::from_json_value(serde_json::Value::Object(serde_json::Map::from_iter(
                [(
                    "AminoAcid".to_string(),
                    serde_json::Value::Array(vec![
                        serde_json::Value::Array(vec![serde_json::Value::String(
                            "Cysteine".to_string()
                        ),]),
                        serde_json::Value::String("Anywhere".to_string())
                    ])
                )]
            ))),
            Ok(PlacementRule::AminoAcid(
                vec![AminoAcid::Cysteine].into(),
                ResiduePosition::Anywhere,
                ModificationPosition::Sidechain
            ))
        );
    }
}
