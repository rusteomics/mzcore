//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use std::{
    cmp::Ordering,
    collections::{BTreeSet, HashMap, HashSet},
    fmt::{Display, Write},
    sync::Arc,
};

use crate::{
    annotation::model::{FragmentationModel, GlycanModel, GlycanPeptideFragment},
    chemistry::{AmbiguousLabel, CachedCharge, Chemical, MolecularFormula},
    fragment::*,
    glycan::{GlycanStructure, MonoSaccharide},
    helper_functions::merge_hashmap,
    molecular_formula,
    ontology::Ontology,
    quantities::Multi,
    sequence::Linked,
    sequence::{
        AminoAcid, MUPSettings, Peptidoform, PlacementRule, Position, SequenceElement,
        SequencePosition,
    },
    system::OrderedMass,
};

/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum Modification {
    /// Any of the simple modifications
    Simple(SimpleModification),
    /// A cross link to another (or the same) peptide, a branch is also seen as a cross-link but then the name is None.
    CrossLink {
        /// The index of the peptide this cross-link is bound to (can be the index for this peptide if it is an intra link)
        peptide: usize,
        /// The sequence index where this cross-link is bound to
        sequence_index: SequencePosition,
        /// The linker that defines the chemical structure that is the actual linker
        linker: SimpleModification,
        /// The name of the cross-linker, if [`CrossLinkName::Branch`] it is a branch instead of cross-link
        name: CrossLinkName,
        /// To determine if the cross-link is placed symmetrically or if asymmetrically if this is the left or right side
        side: CrossLinkSide,
    },
    /// An ambiguous modification, that can be placed at multiple locations
    Ambiguous {
        /// The name of the group
        group: String,
        /// The id to compare be able to find the other locations where this modifications can be placed
        id: usize,
        /// The modification itself
        modification: SimpleModification,
        /// If present the localisation score, meaning the chance/ratio for this modification to show up on this exact spot
        localisation_score: Option<OrderedFloat<f64>>,
        /// If this is the preferred location or not
        preferred: bool,
    },
}

/// Indicate the cross-link side, it contains a set of all placement rules that apply for the placed
/// location to find all possible ways of breaking and/or neutral losses. These numbers are the
/// index into the [`LinkerSpecificity`] rules.
#[derive(Clone, PartialEq, Eq, Debug, Serialize, Deserialize)]
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

/// A modification on an amino acid, wrapped in an [`std::sync::Arc`] to not have to clone modifications from databases.
pub type SimpleModification = Arc<SimpleModificationInner>;

/// A modification on an amino acid
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum SimpleModificationInner {
    /// A modification defined with a monoisotopic mass shift
    Mass(OrderedMass),
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
        /// The structure score
        structure_score: Option<usize>,
        /// The subsumption level
        subsumption_level: GnoSubsumption,
        /// The underlying glycan motif, first is the human description, the second id the GNOme ID of the motif
        motif: Option<(String, String)>,
        /// Taxonomy of the animals in which this glycan is found, defined as a list of species name with taxonomy ID
        taxonomy: thin_vec::ThinVec<(String, usize)>,
        /// Locations of where the glycan exists
        glycomeatlas: thin_vec::ThinVec<(String, Vec<(String, String)>)>,
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

/// A modification id/name
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash, Default)]
pub struct ModificationId {
    /// The ontology where this linker is defined
    pub ontology: Ontology,
    /// The name
    pub name: String,
    /// The id
    pub id: Option<usize>,
    /// The description, mostly for search results
    pub description: String,
    /// Any synonyms
    pub synonyms: thin_vec::ThinVec<String>,
    /// Cross reference IDs
    pub cross_ids: thin_vec::ThinVec<(String, String)>,
}

/// The name of a cross-link
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum CrossLinkName {
    /// A branch
    Branch,
    /// A cross-link
    Name(String),
}

/// The linker position specificities for a linker
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum LinkerSpecificity {
    /// A symmetric specificity where both ends have the same specificity.
    /// The first list is all possible positions. The second list is all
    /// stubs that can be left after cleaving or breaking of the cross-link.
    Symmetric(
        Vec<PlacementRule>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ),
    /// An asymmetric specificity where both ends have a different specificity.
    /// The first list is all possible positions. The second list is all
    /// stubs that can be left after cleaving or breaking of the cross-link.
    Asymmetric(
        (Vec<PlacementRule>, Vec<PlacementRule>),
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ),
}

/// All possible compositions in the GNO ontology
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub enum GnoComposition {
    /// Only the mass is known
    Weight(OrderedMass),
    /// The composition,
    Composition(Vec<(MonoSaccharide, isize)>),
    /// The (full) structure is known
    Topology(GlycanStructure),
}

/// All possible subsumption levels in the GNOme database indicating different levels of description for a glycan species
#[derive(
    Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Default, Hash, Serialize, Deserialize,
)]
pub enum GnoSubsumption {
    /// Indicates only the average weight is defined
    #[default]
    AverageWeight,
    /// Indicates the basic composition, without isomeric information
    BaseComposition,
    /// Indicates the composition, with isomeric information
    Composition,
    /// Indicates the topology, without linkage and anomeric information
    Topology,
    /// Indicates the topology, without reducing end ring and anomeric information
    Saccharide,
}

impl Display for GnoSubsumption {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AverageWeight => write!(f, "Average weight"),
            Self::BaseComposition => write!(f, "Base composition (no isomeric information)"),
            Self::Composition => write!(f, "Composition"),
            Self::Topology => write!(f, "Topology (no linkage)"),
            Self::Saccharide => write!(f, "Saccharide"),
        }
    }
}

impl ModificationId {
    /// Get the accession number name for the ontology
    pub fn url(&self) -> Option<String> {
        match self.ontology {
            Ontology::Unimod => Some(format!(
                "https://www.unimod.org/modifications_view.php?editid1={}",
                self.id.unwrap_or_default()
            )),
            Ontology::Psimod => Some(format!(
                "https://ontobee.org/ontology/MOD?iri=http://purl.obolibrary.org/obo/MOD_{:05}",
                self.id.unwrap_or_default()
            )),
            Ontology::Gnome => Some(format!(
                "http://glytoucan.org/Structures/Glycans/{}",
                self.name.to_ascii_uppercase()
            )),
            Ontology::Resid => Some(format!(
                "https://proteininformationresource.org/cgi-bin/resid?id=AA{:04}",
                self.id.unwrap_or_default()
            )),
            Ontology::Xlmod | Ontology::Custom => None,
        }
    }
}

impl Display for ModificationId {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.ontology == Ontology::Gnome {
            write!(
                f,
                "{}:{}",
                self.ontology.char(),
                self.name.to_ascii_uppercase()
            )
        } else {
            write!(f, "{}:{}", self.ontology.char(), self.name)
        }
    }
}

/// The result of checking if a modification can be placed somewhere.
#[derive(Debug, PartialEq, Eq, Serialize, Deserialize, Clone)]
pub enum RulePossible {
    /// This modification cannot be placed
    No,
    /// This modification can be placed and if it is a cross-link it can be placed on both ends
    Symmetric(BTreeSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'left' side of the cross-link
    AsymmetricLeft(BTreeSet<usize>),
    /// This modification can be placed and if it is a cross-link it can only be placed on the 'right' side of the cross-link
    AsymmetricRight(BTreeSet<usize>),
}

impl RulePossible {
    /// Flatten this into a bool, check if the rule is not [`Self::No`]
    pub fn any_possible(self) -> bool {
        self != Self::No
    }
}

impl std::ops::Add for RulePossible {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Self::Symmetric(a), _) | (_, Self::Symmetric(a)) => Self::Symmetric(a),
            (Self::AsymmetricLeft(l), Self::AsymmetricRight(r))
            | (Self::AsymmetricRight(l), Self::AsymmetricLeft(r)) => {
                let overlap: BTreeSet<usize> = l.intersection(&r).copied().collect();
                if overlap.is_empty() {
                    Self::No
                } else {
                    Self::Symmetric(overlap)
                }
            }
            (Self::AsymmetricLeft(l), _) | (_, Self::AsymmetricLeft(l)) => Self::AsymmetricLeft(l),
            (Self::AsymmetricRight(r), _) | (_, Self::AsymmetricRight(r)) => {
                Self::AsymmetricRight(r)
            }
            _ => Self::No,
        }
    }
}

impl std::iter::Sum for RulePossible {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Self::No, |acc, i| acc + i)
    }
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
        match self {
            Self::Mass(_) | Self::Formula(_) | Self::Glycan(_) | Self::GlycanStructure(_) => None,
            Self::Database { id, .. } | Self::Linker { id, .. } | Self::Gno { id, .. } => id.url(),
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
            Self::Mass(m)
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
                    for option in MonoSaccharide::composition_options(
                        monosaccharides,
                        *range.start() as usize..=*range.end() as usize,
                    ) {
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
                    LinkerSpecificity::Symmetric(rules, _, _) => {
                        if PlacementRule::any_possible(rules, seq, position) {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric((rules_left, rules_right), _, _) => {
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
                    LinkerSpecificity::Symmetric(rules, _, _) => {
                        if PlacementRule::any_possible_aa(rules, aa, position) {
                            RulePossible::Symmetric(BTreeSet::from([index]))
                        } else {
                            RulePossible::No
                        }
                    }
                    LinkerSpecificity::Asymmetric((rules_left, rules_right), _, _) => {
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
            Self::Mass(m) => {
                write!(f, "{:+}", m.value)?;
            }
            Self::Formula(elements) => {
                write!(f, "Formula:{}", elements.hill_notation())?;
            }
            Self::Glycan(monosaccharides) => write!(
                f,
                "Glycan:{}",
                monosaccharides
                    .iter()
                    .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
            )?,
            Self::GlycanStructure(glycan) if specification_compliant => write!(
                f,
                "Glycan:{}|INFO:Structure:{glycan}",
                glycan
                    .composition()
                    .iter()
                    .fold(String::new(), |mut acc, (g, a)| {
                        write!(&mut acc, "{g}{a}").unwrap();
                        acc
                    })
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
                        format!(
                            "{}@{pos}",
                            Ontology::Psimod.find_id(*index, None).unwrap_or_else(|| {
                                panic!(
                    "Invalid PsiMod placement rule, non existing modification {index}"
                )
                            })
                        )
                    }
                })
                .collect_vec(),
            Self::Linker { specificities, .. } => specificities
                .iter()
                .flat_map(|set| match set {
                    LinkerSpecificity::Symmetric(rules, _, _) => rules.clone(),
                    LinkerSpecificity::Asymmetric((rules_a, rules_b), _, _) => rules_a
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
                        format!(
                            "{}@{pos}",
                            Ontology::Psimod.find_id(index, None).unwrap_or_else(|| {
                                panic!(
                "Invalid PsiMod placement rule, non existing modification {index}"
            )
                            })
                        )
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

impl From<SimpleModification> for Modification {
    fn from(value: SimpleModification) -> Self {
        Self::Simple(value)
    }
}

impl From<SimpleModificationInner> for Modification {
    fn from(value: SimpleModificationInner) -> Self {
        Self::Simple(Arc::new(value))
    }
}

impl CrossLinkSide {
    /// Get all allowed placement rules with all applicable neutral losses, stubs, and diagnostic ions.
    pub(crate) fn allowed_rules(
        &self,
        linker: &SimpleModification,
    ) -> (
        Vec<NeutralLoss>,
        Vec<(MolecularFormula, MolecularFormula)>,
        Vec<DiagnosticIon>,
    ) {
        let selected_rules = match self {
            Self::Left(r) | Self::Right(r) | Self::Symmetric(r) => r,
        };
        let mut neutral = Vec::new();
        let mut stubs = Vec::new();
        let mut diagnostic = Vec::new();

        match &**linker {
            SimpleModificationInner::Linker { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    match rule {
                        LinkerSpecificity::Asymmetric(_, n, d) => {
                            diagnostic.extend_from_slice(d);
                            match self {
                                Self::Left(_) => stubs.extend(n.iter().cloned()),
                                Self::Right(_) => {
                                    stubs.extend(n.iter().map(|(l, r)| (r.clone(), l.clone())));
                                }
                                Self::Symmetric(_) => stubs.extend(n.iter().flat_map(|(l, r)| {
                                    vec![(l.clone(), r.clone()), (r.clone(), l.clone())]
                                })),
                            }
                        }
                        LinkerSpecificity::Symmetric(_, n, d) => {
                            stubs.extend_from_slice(n);
                            diagnostic.extend_from_slice(d);
                        }
                    }
                }
            }
            SimpleModificationInner::Database { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    neutral.extend_from_slice(&rule.1);
                    diagnostic.extend_from_slice(&rule.2);
                }
            }
            _ => (),
        }
        (neutral, stubs, diagnostic)
    }
}

impl Modification {
    /// Check if this modification is a simple modification.
    pub const fn is_simple(&self) -> bool {
        matches!(self, Self::Simple(_))
    }
    /// Check if this modification is a cross-link.
    pub const fn is_cross_link(&self) -> bool {
        matches!(self, Self::CrossLink { .. })
    }
    /// Check if this modification is an ambiguous modification.
    pub const fn is_ambiguous(&self) -> bool {
        matches!(self, Self::Ambiguous { .. })
    }

    /// Get the formula for the whole addition (or subtraction) for this modification
    pub(crate) fn formula_inner(
        &self,
        all_peptides: &[Peptidoform<Linked>],
        visited_peptides: &[usize],
        applied_cross_links: &mut Vec<CrossLinkName>,
        allow_ms_cleavable: bool,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        glycan_model: &GlycanModel,
        attachment: Option<AminoAcid>,
    ) -> (
        Multi<MolecularFormula>,
        HashMap<FragmentKind, Multi<MolecularFormula>>,
        HashSet<CrossLinkName>,
    ) {
        match self {
            Self::Simple(modification) | Self::Ambiguous { modification, .. } => {
                match &**modification {
                    // A linker that is not cross-linked is hydrolysed
                    SimpleModificationInner::Linker { formula, .. } => (
                        (formula.clone() + molecular_formula!(H 2 O 1)).into(),
                        HashMap::new(),
                        HashSet::new(),
                    ),

                    s => {
                        let (default_rules, specific_rules) =
                            glycan_model.get_peptide_fragments(attachment);

                        let f = s.formula_inner(
                            sequence_index,
                            peptidoform_index,
                            default_rules,
                            attachment,
                        );
                        let specific = specific_rules
                            .into_iter()
                            .map(|(k, settings)| {
                                (
                                    k,
                                    s.formula_inner(
                                        sequence_index,
                                        peptidoform_index,
                                        settings,
                                        attachment,
                                    ),
                                )
                            })
                            .collect();
                        (f, specific, HashSet::new())
                    }
                }
            }
            Self::CrossLink {
                peptide: other_peptide,
                linker,
                name,
                side,
                ..
            } => {
                if applied_cross_links.contains(name) {
                    (Multi::default(), HashMap::new(), HashSet::default())
                } else if visited_peptides.contains(other_peptide) {
                    applied_cross_links.push(name.clone());
                    (
                        linker
                            .formula_inner(
                                sequence_index,
                                peptidoform_index,
                                glycan_model.default_peptide_fragment,
                                attachment,
                            )
                            .with_label(&AmbiguousLabel::CrossLinkBound(name.clone())),
                        HashMap::new(),
                        HashSet::from([name.clone()]),
                    )
                } else {
                    applied_cross_links.push(name.clone());
                    let link = linker.formula_inner(
                        sequence_index,
                        peptidoform_index,
                        glycan_model.default_peptide_fragment,
                        attachment,
                    );
                    let (_, stubs, _) = side.allowed_rules(linker);

                    if allow_ms_cleavable && !stubs.is_empty() {
                        let mut options: Vec<MolecularFormula> = stubs
                            .iter()
                            .map(|s| {
                                s.0.clone().with_label(AmbiguousLabel::CrossLinkBroken(
                                    name.clone(),
                                    s.0.clone(),
                                ))
                            })
                            .unique()
                            .collect();
                        let mut seen_peptides = HashSet::from([name.clone()]);
                        let mut specific = HashMap::new();

                        options.extend_from_slice(&{
                            let (f, f_specific, seen) = all_peptides[*other_peptide]
                                .formulas_inner(
                                    *other_peptide,
                                    all_peptides,
                                    visited_peptides,
                                    applied_cross_links,
                                    false,
                                    glycan_model,
                                );
                            seen_peptides.extend(seen);
                            specific = merge_hashmap(specific, f_specific);
                            (f * link)
                                .with_label(&AmbiguousLabel::CrossLinkBound(name.clone()))
                                .to_vec()
                        });

                        (options.into(), specific, seen_peptides)
                    } else {
                        let (f, specific, mut seen) = all_peptides[*other_peptide].formulas_inner(
                            *other_peptide,
                            all_peptides,
                            visited_peptides,
                            applied_cross_links,
                            false,
                            glycan_model,
                        );
                        seen.insert(name.clone());
                        (
                            (f * link).with_label(&AmbiguousLabel::CrossLinkBound(name.clone())),
                            specific,
                            seen,
                        )
                    }
                }
            }
        }
    }

    /// Get the formula for a modification, if it is a cross linked modification only get the cross link
    pub fn formula(&self) -> MolecularFormula {
        match self {
            Self::Simple(s) => s.formula(),
            Self::CrossLink { linker, .. } => linker.formula(),
            Self::Ambiguous { modification, .. } => modification.formula(),
        }
    }

    /// Check if this is a simple modification
    pub const fn simple(&self) -> Option<&SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } | Self::Ambiguous { .. } => None,
        }
    }

    /// Check if this is a simple modification
    pub fn into_simple(self) -> Option<SimpleModification> {
        match self {
            Self::Simple(sim) => Some(sim),
            Self::CrossLink { .. } | Self::Ambiguous { .. } => None,
        }
    }

    /// Get a url for more information on this modification. Only defined for modifications from ontologies.
    pub fn ontology_url(&self) -> Option<String> {
        match self {
            Self::Simple(modification)
            | Self::Ambiguous { modification, .. }
            | Self::CrossLink {
                linker: modification,
                ..
            } => modification.ontology_url(),
        }
    }

    /// Check to see if this modification can be placed on the specified element
    pub fn is_possible<T>(
        &self,
        seq: &SequenceElement<T>,
        position: SequencePosition,
    ) -> RulePossible {
        self.simple()
            .map_or(RulePossible::Symmetric(BTreeSet::new()), |s| {
                s.is_possible(seq, position)
            })
    }

    /// Generate theoretical fragments for side chains (glycans)
    pub(crate) fn generate_theoretical_fragments(
        &self,
        model: &FragmentationModel,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment> {
        match self {
            Self::Simple(modification) | Self::Ambiguous { modification, .. } => modification
                .generate_theoretical_fragments(
                    model,
                    peptidoform_ion_index,
                    peptidoform_index,
                    charge_carriers,
                    full_formula,
                    attachment,
                ),
            Self::CrossLink { .. } => Vec::new(),
        }
    }
}

impl SimpleModificationInner {
    /// Generate theoretical fragments for side chains (glycans)
    pub(crate) fn generate_theoretical_fragments(
        &self,
        model: &FragmentationModel,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment> {
        match self {
            Self::GlycanStructure(glycan)
            | Self::Gno {
                composition: GnoComposition::Topology(glycan),
                ..
            } => glycan
                .clone()
                .determine_positions()
                .generate_theoretical_fragments(
                    model,
                    peptidoform_ion_index,
                    peptidoform_index,
                    charge_carriers,
                    full_formula,
                    attachment,
                ),
            Self::Glycan(composition)
            | Self::Gno {
                composition: GnoComposition::Composition(composition),
                ..
            } => MonoSaccharide::theoretical_fragments(
                composition,
                model,
                peptidoform_ion_index,
                peptidoform_index,
                charge_carriers,
                full_formula,
                attachment,
            ),
            _ => Vec::new(),
        }
    }
}

/// The structure to lookup ambiguous modifications, with a list of all modifications (the order is fixed) with for each modification their name and the actual modification itself (if already defined)
pub type AmbiguousLookup = Vec<AmbiguousLookupEntry>;

/// An entry in the ambiguous lookup
#[derive(Debug, Clone)]
pub struct AmbiguousLookupEntry {
    /// The name of the modification
    pub name: String,
    /// The group of the modification
    pub group: Option<usize>,
    /// The modification itself
    pub modification: Option<SimpleModification>,
    /// The allowed locations, the actual allowed locations is the intersection of this set with the ruleset from the modification
    position: Option<Vec<PlacementRule>>,
    /// The maximal number of this modification on one place
    limit: Option<usize>,
    /// Determines if this modification can colocalise with placed modifications eg if the modification of unknown position is allowed at the second M '[Oxidation]?MM[Dioxidation]M'
    colocalise_placed_modifications: bool,
    /// Determines if this modification can colocalise with other modifications of unknown position
    colocalise_modifications_of_unknown_position: bool,
}

impl AmbiguousLookupEntry {
    /// Create a new ambiguous lookup entry
    pub const fn new(name: String, modification: Option<SimpleModification>) -> Self {
        Self {
            name,
            group: None,
            modification,
            limit: None,
            position: None,
            colocalise_placed_modifications: true,
            colocalise_modifications_of_unknown_position: true,
        }
    }

    /// Copy settings into this lookup entry
    pub fn copy_settings(&mut self, settings: &MUPSettings) {
        self.position.clone_from(&settings.position);
        self.limit = settings.limit;
        self.colocalise_placed_modifications = settings.colocalise_placed_modifications;
        self.colocalise_modifications_of_unknown_position =
            settings.colocalise_modifications_of_unknown_position;
    }

    /// Get the settings for this modification of unknown position
    pub fn as_settings(&self) -> MUPSettings {
        MUPSettings {
            position: self.position.clone(),
            limit: self.limit,
            colocalise_placed_modifications: self.colocalise_placed_modifications,
            colocalise_modifications_of_unknown_position: self
                .colocalise_modifications_of_unknown_position,
        }
    }
}

/// The structure to lookup cross-links, with a list of all linkers (the order is fixed) with for each linker their name or None if it is a branch and the actual linker itself (if already defined)
pub type CrossLinkLookup = Vec<(CrossLinkName, Option<SimpleModification>)>;

impl Modification {
    /// Display a modification either normalised to the internal representation or as fully valid ProForma
    /// (no glycan structure or custom modifications). `display_ambiguous` shows or hides the modification
    /// definition of any ambiguous modifications (eg true results in '+1#1' false in '#1').
    /// # Errors
    /// When the given writer errors.
    pub fn display(
        &self,
        f: &mut impl Write,
        specification_compliant: bool,
        display_ambiguous: bool,
    ) -> std::fmt::Result {
        match self {
            Self::Simple(sim) => sim.display(f, specification_compliant),
            Self::CrossLink { name, linker, .. } => {
                linker.display(f, specification_compliant)?;
                write!(f, "{name}")?;
                Ok(())
            }
            Self::Ambiguous {
                group,
                modification,
                localisation_score,
                ..
            } => {
                if display_ambiguous {
                    modification.display(f, specification_compliant)?;
                }
                write!(
                    f,
                    "\x23{group}{}",
                    localisation_score
                        .map(|v| format!("({v})"))
                        .unwrap_or_default()
                )?;
                Ok(())
            }
        }
    }
}

impl Display for Modification {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(f, true, true)
    }
}

impl Display for CrossLinkName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Branch => write!(f, "#BRANCH"),
            Self::Name(n) => write!(f, "#XL{n}"),
        }
    }
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_reading_custom_modifications_json() {
    use serde_json;
    let data = r#"[ [ 1, "uranium linker", { "Linker": { "specificities": [ { "Asymmetric": [ [ [ { "AminoAcid": [ [ "Selenocysteine" ], "AnyCTerm" ] }, { "AminoAcid": [ [ "GlutamicAcid" ], "Anywhere" ] } ], [ { "AminoAcid": [ [ "Selenocysteine" ], "AnyNTerm" ] } ] ], [ [ { "elements": [ [ "U", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "U", null, 1 ] ], "additional_mass": 0.0 } ] ], [ { "elements": [ [ "Te", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "Ne", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 2 ], [ "He", null, 3 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 1 ], [ "He", null, 2 ] ], "additional_mass": 0.0 }, { "elements": [ [ "I", null, 1 ], [ "Er", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 12 ], [ "C", null, 12 ], [ "O", null, 1 ] ], "additional_mass": 0.0 } ] ] } ], "formula": { "elements": [ [ "U", null, 2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Uranium linker", "id": 1, "description": "Have some uranium, its delicious!", "synonyms": [], "cross_ids": [ [ "Pubmed", "21714143" ] ] }, "length": 23.9 } } ], [ 2, "helium", { "Database": { "specificities": [ [ [ { "AminoAcid": [ [ "Alanine" ], "Anywhere" ] } ], [], [] ], [ [ { "AminoAcid": [ [ "Methionine" ], "Anywhere" ] } ], [ { "Loss": { "elements": [], "additional_mass": 12.0 } } ], [] ] ], "formula": { "elements": [ [ "He", null, 2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Helium", "id": 2, "description": "heeeelium", "synonyms": [ "heeeelium", "funny gas" ], "cross_ids": [ [ "Pubmed", "42" ], [ "Unimod", "12" ], [ "Resid", "A12" ] ] } } } ], [ 3, "db18", { "Database": { "specificities": [ [ [ { "AminoAcid": [ [ "Cysteine" ], "Anywhere" ] } ], [ { "Gain": { "elements": [], "additional_mass": 372.25 } }, { "Gain": { "elements": [], "additional_mass": 373.258 } }, { "Gain": { "elements": [], "additional_mass": 371.242 } }, { "Gain": { "elements": [], "additional_mass": 240.171 } }, { "Gain": { "elements": [], "additional_mass": 239.163 } }, { "Gain": { "elements": [], "additional_mass": 241.179 } }, { "Gain": { "elements": [], "additional_mass": 197.129 } }, { "Gain": { "elements": [], "additional_mass": 198.137 } }, { "Gain": { "elements": [], "additional_mass": 196.121 } }, { "Gain": { "elements": [], "additional_mass": 619.418 } }, { "Gain": { "elements": [], "additional_mass": 621.4343 } }, { "Gain": { "elements": [], "additional_mass": 649.465 } }, { "Gain": { "elements": [], "additional_mass": 677.497 } }, { "Gain": { "elements": [], "additional_mass": 618.41 } }, { "Gain": { "elements": [], "additional_mass": 620.426 } }, { "Gain": { "elements": [], "additional_mass": 648.457 } }, { "Gain": { "elements": [], "additional_mass": 676.489 } }, { "Gain": { "elements": [], "additional_mass": 620.426 } }, { "Gain": { "elements": [], "additional_mass": 622.442 } }, { "Gain": { "elements": [], "additional_mass": 650.473 } }, { "Gain": { "elements": [], "additional_mass": 678.504 } }, { "Gain": { "elements": [], "additional_mass": 28.006 } } ], [ { "elements": [], "additional_mass": 372.25 }, { "elements": [], "additional_mass": 240.171 }, { "elements": [], "additional_mass": 197.129 }, { "elements": [], "additional_mass": 619.418 }, { "elements": [], "additional_mass": 621.434 }, { "elements": [], "additional_mass": 649.465 }, { "elements": [], "additional_mass": 677.497 } ] ] ], "formula": { "elements": [], "additional_mass": 676.489 }, "id": { "ontology": "Custom", "name": "DB18", "id": 3, "description": "", "synonyms": [], "cross_ids": [] } } } ], [ 4, "disulfide", { "Linker": { "specificities": [ { "Symmetric": [ [ { "AminoAcid": [ [ "Cysteine" ], "Anywhere" ] } ], [ [ { "elements": [], "additional_mass": 0.0 }, { "elements": [], "additional_mass": 0.0 } ], [ { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 }, { "elements": [], "additional_mass": 0.0 } ], [ { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, -1 ] ], "additional_mass": 0.0 } ] ], [] ] } ], "formula": { "elements": [ [ "H", null, -2 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "Disulfide", "id": 4, "description": "A disulfide with all potential neutral losses", "synonyms": [], "cross_ids": [] }, "length": 0.0 } } ], [ 5, "dsso", { "Linker": { "specificities": [ { "Symmetric": [ [ { "AminoAcid": [ [ "Lysine" ], "Anywhere" ] } ], [ [ { "elements": [ [ "H", null, 1 ], [ "C", null, 3 ], [ "N", null, -1 ], [ "O", null, 3 ], [ "S", null, 1 ] ], "additional_mass": 0.0 }, { "elements": [ [ "H", null, 1 ], [ "C", null, 3 ], [ "N", null, -1 ], [ "O", null, 2 ] ], "additional_mass": 0.0 } ] ], [] ] } ], "formula": { "elements": [ [ "H", null, 2 ], [ "C", null, 6 ], [ "N", null, -2 ], [ "O", null, 5 ], [ "S", null, 1 ] ], "additional_mass": 0.0 }, "id": { "ontology": "Custom", "name": "DSSO", "id": 5, "description": "", "synonyms": [], "cross_ids": [] }, "length": 0.0 } } ]]"#;
    let mods: Vec<(usize, String, SimpleModification)> = serde_json::from_str(data).unwrap();
    assert!(mods.len() > 1);
}
