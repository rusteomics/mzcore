use std::sync::Arc;

use itertools::Itertools;
use ordered_float::OrderedFloat;
use rustyms::{
    chemistry::MolecularFormula,
    fragment::{DiagnosticIon, NeutralLoss},
    ontology::Ontology,
    sequence::{
        LinkerSpecificity, ModificationId, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};

use thin_vec::ThinVec;

#[derive(Debug, Default)]
pub(crate) struct OntologyModification {
    pub formula: MolecularFormula,
    pub name: String,
    pub ontology: Ontology,
    pub id: usize,
    pub description: String,
    pub synonyms: ThinVec<String>,
    pub cross_ids: ThinVec<(String, String)>,
    pub data: ModData,
}

#[derive(Debug)]
pub(crate) enum ModData {
    Mod {
        specificities: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
    },
    Linker {
        length: Option<OrderedFloat<f64>>,
        specificities: Vec<LinkerSpecificity>,
    },
}

impl Default for ModData {
    fn default() -> Self {
        Self::Mod {
            specificities: Vec::new(),
        }
    }
}

impl OntologyModification {
    /// Simplify the placement rules
    pub(crate) fn simplify_rules(&mut self) {
        match &mut self.data {
            ModData::Mod {
                specificities: ref mut rules,
            } => {
                let mut new = Vec::new();
                for rule in rules.iter() {
                    let rule = (
                        rule.0.clone(),
                        rule.1.iter().unique().sorted().cloned().collect(),
                        rule.2.iter().unique().sorted().cloned().collect(),
                    ); // Remove duplicate neutral losses and diagnostic ions, and sort for a better guarantee of equality
                    if new.is_empty() {
                        new.push(rule.clone());
                    } else {
                        let mut found = false;
                        for new_rule in &mut new {
                            // Check if there is a rule with the same neutral loss and diagnostic ions (these can be location specific)
                            if new_rule.1 == rule.1 && new_rule.2 == rule.2 {
                                found = true;
                                // Check if there are other rules in this set of neutral&diagnostic that also use AA placements
                                // If there are, and they are on the same position, merge the AA set
                                for position in &rule.0 {
                                    let mut pos_found = false;
                                    for new_position in &mut new_rule.0 {
                                        if let (
                                            PlacementRule::AminoAcid(new_aa, new_pos),
                                            PlacementRule::AminoAcid(aa, pos),
                                        ) = (new_position, position)
                                        {
                                            if *new_pos == *pos {
                                                for a in aa {
                                                    if !new_aa.contains(a) {
                                                        new_aa.push(*a);
                                                    }
                                                }
                                                new_aa.sort_unstable();
                                                pos_found = true;
                                                break;
                                            }
                                        }
                                    }
                                    if !pos_found {
                                        new_rule.0.push(position.clone());
                                    }
                                }
                            }
                        }
                        if !found {
                            new.push(rule.clone());
                        }
                    }
                }
                rules.clear();
                rules.extend(new);
            }
            ModData::Linker { specificities, .. } => {
                *specificities = specificities.iter().unique().sorted().cloned().collect();
            }
        }
    }

    pub(crate) fn into_mod(mut self) -> (Option<usize>, String, SimpleModification) {
        self.simplify_rules();
        let id = ModificationId {
            ontology: self.ontology,
            name: self.name.clone(),
            id: Some(self.id),
            description: self.description,
            synonyms: self.synonyms,
            cross_ids: self.cross_ids,
        };
        match self.data {
            ModData::Mod { specificities } => (
                Some(self.id),
                self.name.to_ascii_lowercase(),
                Arc::new(SimpleModificationInner::Database {
                    id,
                    formula: self.formula,
                    specificities,
                }),
            ),
            ModData::Linker {
                specificities,
                length,
            } => (
                Some(self.id),
                self.name.to_ascii_lowercase(),
                Arc::new(SimpleModificationInner::Linker {
                    specificities,
                    formula: self.formula,
                    id,
                    length,
                }),
            ),
        }
    }
}

// pub fn position_from_u8(value: u8) -> Result<Position, String> {
//     match value {
//         b'1' => Ok(Position::Anywhere),
//         b'2' => Ok(Position::Anywhere),
//         b'3' => Ok(Position::AnyNTerm),
//         b'4' => Ok(Position::AnyCTerm),
//         b'5' => Ok(Position::ProteinNTerm),
//         b'6' => Ok(Position::ProteinCTerm),
//         n => Err(format!("Outside range: {n}")),
//     }
// }

pub(crate) fn position_from_str(value: &str) -> Result<Position, String> {
    match value {
        "" | "Anywhere" => Ok(Position::Anywhere),
        "Any N-term" => Ok(Position::AnyNTerm),
        "Any C-term" => Ok(Position::AnyCTerm),
        "Protein N-term" => Ok(Position::ProteinNTerm),
        "Protein C-term" => Ok(Position::ProteinCTerm),
        n => Err(format!("Not valid position: {n}")),
    }
}
