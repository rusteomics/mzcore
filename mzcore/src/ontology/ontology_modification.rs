use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::Ontology,
    sequence::{
        LinkerLength, LinkerSpecificity, ModificationId, PlacementRule, SimpleModification,
        SimpleModificationInner,
    },
};
use itertools::Itertools;
use mzcv::SynonymScope;
use thin_vec::ThinVec;

#[derive(Debug, Default)]
pub(crate) struct OntologyModification {
    pub formula: MolecularFormula,
    pub name: Box<str>,
    pub ontology: Ontology,
    pub id: u32,
    pub description: Box<str>,
    pub synonyms: ThinVec<(SynonymScope, Box<str>)>,
    pub cross_ids: ThinVec<(Option<Box<str>>, Box<str>)>,
    pub data: ModData,
    pub obsolete: bool,
}

#[derive(Debug)]
pub(crate) enum ModData {
    Mod {
        specificities: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
    },
    Linker {
        length: LinkerLength,
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
                specificities: rules,
            } => {
                let mut new = Vec::new();
                for rule in rules.iter() {
                    let rule: (Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>) = (
                        rule.0.iter().unique().sorted().cloned().collect(),
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
                                            && *new_pos == *pos
                                        {
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
                *specificities = specificities
                    .iter()
                    .map(|rule| match rule {
                        LinkerSpecificity::Asymmetric {
                            rules,
                            stubs,
                            neutral_losses,
                            diagnostic,
                        } => {
                            let left = rules.0.iter().unique().sorted().cloned().collect();
                            let right = rules.1.iter().unique().sorted().cloned().collect();

                            if left == right {
                                LinkerSpecificity::Symmetric {
                                    rules: left,
                                    stubs: stubs.iter().unique().sorted().cloned().collect(),
                                    neutral_losses: neutral_losses
                                        .iter()
                                        .unique()
                                        .sorted()
                                        .cloned()
                                        .collect(),
                                    diagnostic: diagnostic
                                        .iter()
                                        .unique()
                                        .sorted()
                                        .cloned()
                                        .collect(),
                                }
                            } else {
                                LinkerSpecificity::Asymmetric {
                                    rules: (left, right),
                                    stubs: stubs.iter().unique().sorted().cloned().collect(),
                                    neutral_losses: neutral_losses
                                        .iter()
                                        .unique()
                                        .sorted()
                                        .cloned()
                                        .collect(),
                                    diagnostic: diagnostic
                                        .iter()
                                        .unique()
                                        .sorted()
                                        .cloned()
                                        .collect(),
                                }
                            }
                        }
                        LinkerSpecificity::Symmetric {
                            rules,
                            stubs,
                            neutral_losses,
                            diagnostic,
                        } => LinkerSpecificity::Symmetric {
                            rules: rules.iter().unique().sorted().cloned().collect(),
                            stubs: stubs.iter().unique().sorted().cloned().collect(),
                            neutral_losses: neutral_losses
                                .iter()
                                .unique()
                                .sorted()
                                .cloned()
                                .collect(),
                            diagnostic: diagnostic.iter().unique().sorted().cloned().collect(),
                        },
                    })
                    .unique()
                    .sorted()
                    .collect();
            }
        }
    }
}

impl From<OntologyModification> for SimpleModificationInner {
    fn from(mut value: OntologyModification) -> Self {
        value.simplify_rules();
        let id = ModificationId::new(
            value.ontology,
            value.name,
            Some(value.id),
            value.description,
            value.synonyms,
            value.cross_ids,
            value.obsolete,
        );
        match value.data {
            ModData::Mod { specificities } => Self::Database {
                id,
                formula: value.formula,
                specificities,
            },
            ModData::Linker {
                specificities,
                length,
            } => Self::Linker {
                specificities,
                formula: value.formula,
                id,
                length,
            },
        }
    }
}

impl From<OntologyModification> for SimpleModification {
    fn from(value: OntologyModification) -> Self {
        Self::new(value.into())
    }
}
