use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::Ontology,
    sequence::{
        CrossId, LinkerLength, LinkerSpecificity, ModificationId, PlacementRule,
        SimpleModification, SimpleModificationInner,
    },
};
use context_error::{BoxedError, Context, CreateError, combine_error};
use itertools::Itertools;
use mzcv::{
    AccessionCode, AccessionCodeParseError, CVError, Comment, Modifier, OboIdentifier,
    RelationType, SynonymScope,
};
use thin_vec::ThinVec;

#[derive(Debug, Default)]
pub(crate) struct OntologyModification {
    pub formula: MolecularFormula,
    pub name: Box<str>,
    pub ontology: Ontology,
    pub id: AccessionCode,
    pub description: Box<str>,
    pub synonyms: ThinVec<(SynonymScope, Box<str>)>,
    pub cross_ids: ThinVec<CrossId>,
    pub data: ModData,
    pub obsolete: bool,
    pub parents: ThinVec<AccessionCode>,
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
                specificities: old_rules,
            } => {
                simplify_rules(old_rules);
            }
            ModData::Linker { specificities, .. } => {
                simplify_xl_rules(specificities);
            }
        }
    }

    pub(crate) fn add_relationships(
        &mut self,
        relationships: &[(RelationType, OboIdentifier, Vec<Modifier>, Comment)],
    ) -> Vec<BoxedError<'static, CVError>> {
        let mut errors = Vec::new();

        for rel in relationships {
            if rel.0 == RelationType::IsA {
                match rel.1.1.parse() {
                    Ok(v) => self.parents.push(v),
                    Err(err) => combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid ID",
                            match err {
                                AccessionCodeParseError::Empty => {
                                    "A relationship ID cannot be empty"
                                }
                                AccessionCodeParseError::InvalidCharacters(_) => {
                                    "A relationship ID can only contain alphanumeric characters"
                                }
                                AccessionCodeParseError::TooLong(_) => {
                                    "A relationship ID can at max be 8 characters"
                                }
                            },
                            Context::default().lines(0, rel.1.1.to_string()),
                        ),
                    ),
                }
            }
        }

        errors
    }

    pub(crate) fn finish(mods: Vec<Self>) -> Vec<SimpleModification> {
        let mut links = std::collections::HashMap::new();

        for m in &mods {
            for rel in &m.parents {
                for other in &mods {
                    if other.id == *rel {
                        links
                            .entry(m.id)
                            .or_insert_with(|| (ThinVec::new(), ThinVec::new()))
                            .0
                            .push(other.id);
                        links
                            .entry(other.id)
                            .or_insert_with(|| (ThinVec::new(), ThinVec::new()))
                            .1
                            .push(m.id);
                    }
                }
            }
        }

        mods.into_iter()
            .map(|mut m| {
                m.simplify_rules();
                let mut id = ModificationId::new(
                    m.ontology,
                    m.name,
                    m.id,
                    m.description,
                    m.synonyms,
                    m.cross_ids,
                    m.obsolete,
                );
                if let Some((parents, children)) = links.remove(&m.id) {
                    id.parents = parents;
                    id.children = children;
                }
                std::sync::Arc::new(match m.data {
                    ModData::Mod { specificities } => SimpleModificationInner::Database {
                        id,
                        formula: m.formula,
                        specificities,
                    },
                    ModData::Linker {
                        specificities,
                        length,
                    } => SimpleModificationInner::Linker {
                        specificities,
                        formula: m.formula,
                        id,
                        length,
                    },
                })
            })
            .collect()
    }
}

fn simplify_rules(old_rules: &mut Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>) {
    let mut simplified_rules: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)> =
        Vec::new();
    for rule in old_rules.iter() {
        let rule: (Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>) = (
            compress_rules(&rule.0),
            rule.1.iter().unique().sorted().cloned().collect(),
            rule.2.iter().unique().sorted().cloned().collect(),
        ); // Remove duplicate neutral losses and diagnostic ions, and sort for a better guarantee of equality
        // if simplified_rules.is_empty() {
        //     simplified_rules.push(rule.clone());
        // } else {
        let mut found = false;
        for simplified_rule in &mut simplified_rules {
            // Check if there is a rule with the same neutral loss and diagnostic ions (these can be location specific)
            if simplified_rule.1 == rule.1 && simplified_rule.2 == rule.2 {
                found = true;
                combine_rules(&mut simplified_rule.0, &rule.0);
            }
        }
        if !found {
            simplified_rules.push(rule);
        }
        // }
    }
    old_rules.clear();
    old_rules.extend(simplified_rules);
}

fn simplify_xl_rules(old_rules: &mut Vec<LinkerSpecificity>) {
    let mut new_rules = Vec::new();

    for old_rule in old_rules.iter() {
        let rule = match old_rule {
            LinkerSpecificity::Asymmetric {
                rules,
                stubs,
                neutral_losses,
                diagnostic,
            } => {
                let left = compress_rules(&rules.0);
                let right = compress_rules(&rules.1);

                if left == right {
                    LinkerSpecificity::Symmetric {
                        rules: left,
                        stubs: stubs.iter().unique().sorted().cloned().collect(),
                        neutral_losses: neutral_losses.iter().unique().sorted().cloned().collect(),
                        diagnostic: diagnostic.iter().unique().sorted().cloned().collect(),
                    }
                } else {
                    LinkerSpecificity::Asymmetric {
                        rules: (left, right),
                        stubs: stubs.iter().unique().sorted().cloned().collect(),
                        neutral_losses: neutral_losses.iter().unique().sorted().cloned().collect(),
                        diagnostic: diagnostic.iter().unique().sorted().cloned().collect(),
                    }
                }
            }
            LinkerSpecificity::Symmetric {
                rules,
                stubs,
                neutral_losses,
                diagnostic,
            } => LinkerSpecificity::Symmetric {
                rules: compress_rules(rules),
                stubs: stubs.iter().unique().sorted().cloned().collect(),
                neutral_losses: neutral_losses.iter().unique().sorted().cloned().collect(),
                diagnostic: diagnostic.iter().unique().sorted().cloned().collect(),
            },
        };

        let mut found = false;
        for new_rule in &mut new_rules {
            match (&rule, new_rule) {
                (
                    LinkerSpecificity::Symmetric {
                        rules,
                        stubs,
                        neutral_losses,
                        diagnostic,
                    },
                    LinkerSpecificity::Symmetric {
                        rules: new_rules,
                        stubs: stubs2,
                        neutral_losses: losses2,
                        diagnostic: diagnostic2,
                    },
                ) if *stubs == *stubs2
                    && *neutral_losses == *losses2
                    && *diagnostic == *diagnostic2 =>
                {
                    found = true;
                    combine_rules(new_rules, rules);
                }
                _ => (),
            }
        }
        if !found {
            new_rules.push(rule);
        }
    }

    old_rules.clear();
    old_rules.extend_from_slice(&new_rules);
}

fn combine_rules(rules: &mut Vec<PlacementRule>, additional_rules: &[PlacementRule]) {
    for rule in additional_rules {
        let mut pos_found = false;
        for new_position in rules.iter_mut() {
            if new_position.combine_rules(rule) {
                pos_found = true;
                break;
            }
        }
        if !pos_found {
            rules.push(rule.clone());
        }
    }
    rules.sort_unstable();
}

fn compress_rules(rules: &[PlacementRule]) -> Vec<PlacementRule> {
    let mut new: Vec<PlacementRule> = Vec::new();
    for rule in rules.iter().unique() {
        let mut pos_found = false;
        for new_position in &mut new {
            if new_position.combine_rules(rule) {
                pos_found = true;
                break;
            }
        }
        if !pos_found {
            new.push(rule.clone());
        }
    }
    new.sort_unstable();
    new
}

#[test]
fn compress() {
    let rules = compress_rules(&[
        PlacementRule::AminoAcid(
            vec![crate::sequence::AminoAcid::Serine].into(),
            crate::sequence::Position::Anywhere,
        ),
        PlacementRule::AminoAcid(
            vec![crate::sequence::AminoAcid::Threonine].into(),
            crate::sequence::Position::Anywhere,
        ),
        PlacementRule::AminoAcid(
            vec![crate::sequence::AminoAcid::Tyrosine].into(),
            crate::sequence::Position::Anywhere,
        ),
        PlacementRule::AminoAcid(
            vec![crate::sequence::AminoAcid::Lysine].into(),
            crate::sequence::Position::Anywhere,
        ),
        PlacementRule::Position(crate::sequence::Position::ProteinNTerm),
    ]);
    assert_eq!(rules.len(), 2);
}
