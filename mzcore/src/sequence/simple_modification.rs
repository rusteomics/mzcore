use context_error::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use thin_vec::ThinVec;

use std::{
    collections::BTreeSet,
    fmt::{Display, Write},
    sync::Arc,
};

use crate::{
    chemistry::{AmbiguousLabel, Chemical, DiagnosticIon, MolecularFormula, NeutralLoss},
    glycan::{GlycanPeptideFragment, GlycanStructure, MonoSaccharide},
    ontology::Ontology,
    parse_json::ParseJson,
    quantities::Multi,
    sequence::{
        AminoAcid, GnoComposition, GnoSubsumption, LinkerSpecificity, ModificationId,
        PlacementRule, Position, RulePossible, SequenceElement, SequencePosition,
    },
    system::{Mass, OrderedMass, dalton},
};

/// A modification on an amino acid, wrapped in an [`std::sync::Arc`] to not have to clone modifications from databases.
pub type SimpleModification = Arc<SimpleModificationInner>;

/// A tag for a mass modification
#[derive(
    Clone, Copy, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize,
)]
pub enum MassTag {
    /// No tag
    #[default]
    None,
    /// An ontology
    Ontology(Ontology),
    /// Obs:
    Observed,
}

impl Display for MassTag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::None => Ok(()),
            Self::Ontology(o) => write!(f, "{}:", o.char()),
            Self::Observed => write!(f, "Obs:"),
        }
    }
}

/// A modification on an amino acid
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum SimpleModificationInner {
    /// A modification defined with a monoisotopic mass shift, the mass tag (if defined), the mass itself, and the number of digits it was defined with
    Mass(MassTag, OrderedMass, Option<u8>),
    /// A modification defined with a molecular formula
    Formula(MolecularFormula),
    /// A glycan without a defined structure
    Glycan(Vec<(MonoSaccharide, isize)>),
    /// A glycan with a defined structure
    GlycanStructure(GlycanStructure),
    /// A modification from the GNOme ontology
    Gno {
        /// The composition, weight/composition/topology
        composition: GnoComposition,
        /// The id/name
        id: ModificationId,
        /// The structure score, u16::MAX is used to indicate that no structure score is available
        structure_score: u16,
        /// The subsumption level
        subsumption_level: GnoSubsumption,
        /// The underlying glycan motif, first is the human description, the second id the GNOme ID of the motif
        motif: Option<(Box<str>, Box<str>)>,
        /// Taxonomy of the animals in which this glycan is found, defined as a list of species name with taxonomy ID
        taxonomy: ThinVec<(Box<str>, usize)>,
        /// Locations of where the glycan exists
        glycomeatlas: ThinVec<(Box<str>, ThinVec<(Box<str>, Box<str>)>)>,
    },
    /// A modification from one of the modification ontologies
    Database {
        /// The placement rules, neutral losses, and diagnostic ions
        specificities: Vec<(Vec<PlacementRule>, Vec<NeutralLoss>, Vec<DiagnosticIon>)>,
        /// The chemical formula for this modification (diff formula)
        formula: MolecularFormula,
        /// The id/name
        id: ModificationId,
    },
    /// A cross-linker
    Linker {
        /// All possible specificities for this linker
        specificities: Vec<LinkerSpecificity>,
        /// The chemical formula for this linker (diff formula)
        formula: MolecularFormula,
        /// The id/name
        id: ModificationId,
        /// The length, if known
        length: Option<OrderedFloat<f64>>,
    },
}

impl Chemical for SimpleModificationInner {
    /// Get the molecular formula for this modification.
    fn formula_inner(
        &self,
        position: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.formula_inner(
            position,
            peptidoform_index,
            GlycanPeptideFragment::FULL,
            None,
        )
        .to_vec()
        .pop()
        .unwrap()
    }
}

impl SimpleModificationInner {
    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    pub fn ontology_url(&self) -> Option<String> {
        self.description().and_then(ModificationId::url)
    }

    /// Get the description, only defined for modifications from ontologies.
    pub const fn description(&self) -> Option<&ModificationId> {
        match self {
            Self::Mass(_, _, _) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_) => {
                None
            }
            Self::Database { id, .. } | Self::Linker { id, .. } | Self::Gno { id, .. } => Some(id),
        }
    }

    /// Internal formula code with the logic to make all labels right
    pub(crate) fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        glycan_fragmentation: GlycanPeptideFragment,
        attachment: Option<AminoAcid>,
    ) -> Multi<MolecularFormula> {
        match self {
            Self::Mass(_, m, _)
            | Self::Gno {
                composition: GnoComposition::Weight(m),
                ..
            } => MolecularFormula::with_additional_mass(m.value).into(),
            Self::Gno {
                composition: GnoComposition::Composition(monosaccharides),
                ..
            }
            | Self::Glycan(monosaccharides) => {
                let mut options = Vec::new();

                if let Some(range) = glycan_fragmentation.core() {
                    for option in MonoSaccharide::composition_options(monosaccharides, range) {
                        options.push(
                            option
                                .iter()
                                .fold(MolecularFormula::default(), |acc, i| {
                                    acc + i.0.formula_inner(sequence_index, peptidoform_index)
                                        * i.1 as i32
                                })
                                .with_label(AmbiguousLabel::GlycanFragmentComposition(option)),
                        );
                    }
                }
                if glycan_fragmentation.full() {
                    options.push(monosaccharides.iter().fold(
                        MolecularFormula::default(),
                        |acc, i| {
                            acc + i.0.formula_inner(sequence_index, peptidoform_index) * i.1 as i32
                        },
                    ));
                }
                if options.is_empty() {
                    options.push(MolecularFormula::default());
                }
                options.into()
            }
            Self::GlycanStructure(glycan)
            | Self::Gno {
                composition: GnoComposition::Topology(glycan),
                ..
            } => {
                let mut options = Vec::new();

                if let Some(range) = glycan_fragmentation.core() {
                    for option in glycan.clone().determine_positions().core_options(
                        range,
                        peptidoform_index,
                        attachment.map(|a| (a, sequence_index)),
                    ) {
                        options.push(
                            option
                                .1
                                .with_label(AmbiguousLabel::GlycanFragment(option.0)),
                        );
                    }
                }
                if glycan_fragmentation.full() {
                    options.push(glycan.formula_inner(sequence_index, peptidoform_index));
                }
                if options.is_empty() {
                    options.push(MolecularFormula::default());
                }
                options.into()
            }
            Self::Formula(formula)
            | Self::Database { formula, .. }
            | Self::Linker { formula, .. } => formula.into(),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible<T>(
        &self,
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> RulePossible {
        match self {
            Self::Database { specificities, .. } if specificities.is_empty() => {
                RulePossible::Symmetric(BTreeSet::default())
            }
            Self::Database { specificities, .. } => {
                // If any of the rules match the current situation then it can be placed
                let matching: BTreeSet<usize> = specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(index, (rules, _, _))| {
                        PlacementRule::any_possible(rules, seq, position).then_some(index)
                    })
                    .collect();
                if matching.is_empty() {
                    RulePossible::No
                } else {
                    RulePossible::Symmetric(matching)
                }
            }
            Self::Linker { specificities, .. } if specificities.is_empty() => {
                RulePossible::Symmetric(BTreeSet::default())
            }
            Self::Linker { specificities, .. } => specificities
                .iter()
                .enumerate()
                .map(|(index, spec)| match spec {
                    LinkerSpecificity::Symmetric { rules, .. } => {
                        if PlacementRule::any_possible(rules, seq, position) {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric {
                        rules: (rules_left, rules_right),
                        ..
                    } => {
                        let left = PlacementRule::any_possible(rules_left, seq, position);
                        let right = PlacementRule::any_possible(rules_right, seq, position);
                        if left && right {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else if left {
                            RulePossible::AsymmetricLeft(BTreeSet::from([index]))
                        } else if right {
                            RulePossible::AsymmetricRight(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                })
                .sum::<RulePossible>(),
            _ => RulePossible::Symmetric(BTreeSet::default()),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible_aa(&self, aa: AminoAcid, position: Position) -> RulePossible {
        match self {
            Self::Database { specificities, .. } => {
                // If any of the rules match the current situation then it can be placed
                let matching: BTreeSet<usize> = specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(index, (rules, _, _))| {
                        PlacementRule::any_possible_aa(rules, aa, position).then_some(index)
                    })
                    .collect();
                if matching.is_empty() {
                    RulePossible::No
                } else {
                    RulePossible::Symmetric(matching)
                }
            }
            Self::Linker { specificities, .. } => specificities
                .iter()
                .enumerate()
                .map(|(index, spec)| match spec {
                    LinkerSpecificity::Symmetric { rules, .. } => {
                        if PlacementRule::any_possible_aa(rules, aa, position) {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric {
                        rules: (rules_left, rules_right),
                        ..
                    } => {
                        let left = PlacementRule::any_possible_aa(rules_left, aa, position);
                        let right = PlacementRule::any_possible_aa(rules_right, aa, position);
                        if left && right {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else if left {
                            RulePossible::AsymmetricLeft(BTreeSet::from([index]))
                        } else if right {
                            RulePossible::AsymmetricRight(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                })
                .sum::<RulePossible>(),
            _ => RulePossible::Symmetric(BTreeSet::default()),
        }
    }

    /// Display a modification either normalised to the internal representation or as fully valid ProForma
    /// (no glycan structure or custom modifications).
    /// # Errors
    /// When the given writer errors.
    pub fn display(&self, f: &mut impl Write, specification_compliant: bool) -> std::fmt::Result {
        match self {
            Self::Mass(tag, m, Some(digits)) => {
                write!(f, "{tag}{:+.*}", *digits as usize, m.value)?;
            }
            Self::Mass(tag, m, None) => {
                write!(f, "{tag}{:+}", m.value)?;
            }
            Self::Formula(elements) => {
                write!(f, "Formula:{}", elements.hill_notation())?;
            }
            Self::Glycan(monosaccharides) => write!(
                f,
                "Glycan:{}",
                MonoSaccharide::display_composition(monosaccharides)
            )?,
            Self::GlycanStructure(glycan) if specification_compliant => write!(
                f,
                "Glycan:{}|INFO:Structure:{glycan}",
                MonoSaccharide::display_composition(&glycan.composition())
            )?,
            Self::GlycanStructure(glycan) => write!(f, "GlycanStructure:{glycan}")?,
            Self::Database {
                formula,
                id:
                    ModificationId {
                        name,
                        ontology: Ontology::Custom,
                        ..
                    },
                ..
            } if specification_compliant => {
                write!(f, "Formula:{formula}|INFO:Custom:{name}")?;
            }
            Self::Database { id, .. } | Self::Gno { id, .. } | Self::Linker { id, .. } => {
                write!(f, "{id}")?;
            }
        }
        Ok(())
    }

    /// Get all placement rules as text
    /// # Panics
    /// When a PSI-MOD modification rule uses an non existing modification
    pub(crate) fn placement_rules(&self) -> Vec<String> {
        match self {
            Self::Database { specificities, .. } => specificities
                .iter()
                .flat_map(|set| &set.0)
                .map(|rule| match rule {
                    PlacementRule::AminoAcid(aa, pos) => {
                        format!("{}@{pos}", aa.iter().join(""))
                    }
                    PlacementRule::Terminal(pos) => pos.to_string(),
                    PlacementRule::Anywhere => "Anywhere".to_string(),
                    PlacementRule::PsiModification(index, pos) => {
                        format!("MOD:{index:07}@{pos}")
                    }
                })
                .collect_vec(),
            Self::Linker { specificities, .. } => specificities
                .iter()
                .flat_map(|set| match set {
                    LinkerSpecificity::Symmetric { rules, .. } => rules.clone(),
                    LinkerSpecificity::Asymmetric {
                        rules: (rules_a, rules_b),
                        ..
                    } => rules_a
                        .iter()
                        .cloned()
                        .chain(rules_b.iter().cloned())
                        .collect_vec(),
                })
                .map(|rule| match rule {
                    PlacementRule::AminoAcid(aa, pos) => {
                        format!("{}@{pos}", aa.iter().join(""))
                    }
                    PlacementRule::Terminal(pos) => pos.to_string(),
                    PlacementRule::Anywhere => "Anywhere".to_string(),
                    PlacementRule::PsiModification(index, pos) => {
                        format!("MOD:{index:07}@{pos}")
                    }
                })
                .collect_vec(),
            _ => Vec::new(),
        }
    }
}

impl Display for SimpleModificationInner {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true)
    }
}

impl ParseJson for SimpleModificationInner {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Object(map) = value {
            let (key, value) = map.into_iter().next().unwrap();
            match key.as_str() {
                "Mass" => f64::from_json_value(value)
                    .map(|v| Self::Mass(MassTag::None, Mass::new::<dalton>(v).into(), None)), // TODO: add parser for new format
                "Formula" => MolecularFormula::from_json_value(value).map(Self::Formula),
                "Glycan" => Vec::from_json_value(value).map(Self::Glycan),
                "GlycanStructure" => {
                    GlycanStructure::from_json_value(value).map(Self::GlycanStructure)
                }
                "Gno" => {
                    if let Value::Object(mut map) = value {
                        let context = |map: &serde_json::Map<String, Value>| {
                            Context::show(
                                map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                            )
                        };
                        Ok(Self::Gno {
                            composition: GnoComposition::from_json_value(
                                map.remove("composition").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'composition' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            id: ModificationId::from_json_value(map.remove("id").ok_or_else(
                                || {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'id' is missing",
                                        context(&map),
                                    )
                                },
                            )?)?,
                            structure_score: Option::from_json_value(
                                map.remove("structure_score").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'structure_score' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?
                            .unwrap_or(u16::MAX),
                            subsumption_level: GnoSubsumption::from_json_value(
                                map.remove("subsumption_level").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'subsumption_level' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            motif: Option::from_json_value(map.remove("motif").ok_or_else(
                                || {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'motif' is missing",
                                        context(&map),
                                    )
                                },
                            )?)?,
                            taxonomy: ThinVec::from_json_value(
                                map.remove("taxonomy").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'taxonomy' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            glycomeatlas: ThinVec::from_json_value(
                                map.remove("glycomeatlas").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'glycomeatlas' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                        })
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid Gno SimpleModification",
                            "The value has to be a map",
                            Context::show(key),
                        ))
                    }
                }
                "Database" => {
                    if let Value::Object(mut map) = value {
                        let context = |map: &serde_json::Map<String, Value>| {
                            Context::show(
                                map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                            )
                        };
                        Ok(Self::Database {
                            specificities: Vec::from_json_value(
                                map.remove("specificities").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'specificities' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            formula: MolecularFormula::from_json_value(
                                map.remove("formula").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'formula' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            id: ModificationId::from_json_value(map.remove("id").ok_or_else(
                                || {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'id' is missing",
                                        context(&map),
                                    )
                                },
                            )?)?,
                        })
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid Database SimpleModification",
                            "The value has to be a map",
                            Context::show(key),
                        ))
                    }
                }
                "Linker" => {
                    if let Value::Object(mut map) = value {
                        let context = |map: &serde_json::Map<String, Value>| {
                            Context::show(
                                map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                            )
                        };
                        Ok(Self::Linker {
                            specificities: Vec::from_json_value(
                                map.remove("specificities").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'specificities' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            formula: MolecularFormula::from_json_value(
                                map.remove("formula").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'formula' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?,
                            id: ModificationId::from_json_value(map.remove("id").ok_or_else(
                                || {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'id' is missing",
                                        context(&map),
                                    )
                                },
                            )?)?,
                            length: Option::<f64>::from_json_value(
                                map.remove("length").ok_or_else(|| {
                                    BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SimpleModification",
                                        "The required property 'length' is missing",
                                        context(&map),
                                    )
                                })?,
                            )?
                            .map(Into::into),
                        })
                    } else {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid Database SimpleModification",
                            "The value has to be a map",
                            Context::show(key),
                        ))
                    }
                }

                _ => Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SimpleModification",
                    "The tag has to be Mass/Formula/Glycan/GlycanStructure/Gno/Database/Linker",
                    Context::show(key),
                )),
            }
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid SimpleModification",
                "The JSON value has to be a map",
                Context::show(value.to_string()),
            ))
        }
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use super::*;
    use crate::molecular_formula;

    #[test]
    fn modification_masses() {
        let ((modification, _), _) = SimpleModificationInner::pro_forma(
            "G:G09675LS",
            &mut Vec::new(),
            &mut Vec::new(),
            &crate::ontology::STATIC_ONTOLOGIES,
        )
        .unwrap();
        let modification = modification.defined().unwrap();
        assert_eq!(
            modification.formula_inner(
                SequencePosition::Index(0),
                0,
                GlycanPeptideFragment::FULL,
                None,
            ),
            molecular_formula!(C 6 H 10 O 5).into()
        );
        assert_eq!(
            modification.formula_inner(
                SequencePosition::Index(0),
                0,
                GlycanPeptideFragment::CORE,
                None,
            ),
            molecular_formula!(C 6 H 10 O 5).into()
        );
        assert_eq!(
            modification.formula_inner(
                SequencePosition::Index(0),
                0,
                GlycanPeptideFragment::FREE,
                None,
            ),
            Multi::default()
        );
    }
}
