//! Handle modification related issues, access provided if you want to dive deeply into modifications in your own code.

use std::{
    collections::{BTreeSet, HashMap, HashSet},
    fmt::{Display, Write},
    path::Path,
    sync::Arc,
};

use context_error::*;
use itertools::Itertools;
use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::{FragmentationModel, GlycanModel},
    chemistry::{AmbiguousLabel, CachedCharge, Chemical, MolecularFormula},
    fragment::*,
    glycan::{GlycanStructure, MonoSaccharide},
    helper_functions::merge_hashmap,
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    parse_json::{ParseJson, use_serde},
    quantities::Multi,
    sequence::{
        AminoAcid, CrossLinkName, CrossLinkSide, Linked, LinkerSpecificity, MUPSettings,
        Peptidoform, PlacementRule, SequenceElement, SequencePosition, SimpleModification,
        SimpleModificationInner,
    },
    system::OrderedMass,
};

/// A modification on an amino acid
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
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

/// A modification id/name
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
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

impl ParseJson for ModificationId {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}

/// All possible compositions in the GNO ontology
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum GnoComposition {
    /// Only the mass is known
    Weight(OrderedMass),
    /// The composition,
    Composition(Vec<(MonoSaccharide, isize)>),
    /// The (full) structure is known
    Topology(GlycanStructure),
}

impl ParseJson for GnoComposition {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
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

impl ParseJson for GnoSubsumption {
    fn from_json_value(value: serde_json::Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
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
#[derive(Clone, Debug, Deserialize, Eq, PartialEq, Serialize)]
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
        let mut all_neutral = Vec::new();
        let mut all_stubs = Vec::new();
        let mut all_diagnostic = Vec::new();

        match &**linker {
            SimpleModificationInner::Linker { specificities, .. } => {
                for rule in specificities
                    .iter()
                    .enumerate()
                    .filter_map(|(i, r)| selected_rules.contains(&i).then_some(r))
                {
                    match rule {
                        LinkerSpecificity::Asymmetric {
                            diagnostic,
                            stubs,
                            neutral_losses,
                            ..
                        } => {
                            all_neutral.extend_from_slice(neutral_losses);
                            all_diagnostic.extend_from_slice(diagnostic);
                            match self {
                                Self::Left(_) => all_stubs.extend(stubs.iter().cloned()),
                                Self::Right(_) => {
                                    all_stubs
                                        .extend(stubs.iter().map(|(l, r)| (r.clone(), l.clone())));
                                }
                                Self::Symmetric(_) => {
                                    all_stubs.extend(stubs.iter().flat_map(|(l, r)| {
                                        vec![(l.clone(), r.clone()), (r.clone(), l.clone())]
                                    }));
                                }
                            }
                        }
                        LinkerSpecificity::Symmetric {
                            stubs,
                            neutral_losses,
                            diagnostic,
                            ..
                        } => {
                            all_stubs.extend_from_slice(stubs);
                            all_neutral.extend_from_slice(neutral_losses);
                            all_diagnostic.extend_from_slice(diagnostic);
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
                    all_neutral.extend_from_slice(&rule.1);
                    all_diagnostic.extend_from_slice(&rule.2);
                }
            }
            _ => (),
        }
        (all_neutral, all_stubs, all_diagnostic)
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
        peptidoform_ion_index: usize,
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
                                    peptidoform_ion_index,
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
                            peptidoform_ion_index,
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

    /// Get the name if this is a Unimod modification (for use in mzPAF)
    pub(crate) fn unimod_name(&self) -> Option<&str> {
        match self {
            Self::Simple(s)
            | Self::CrossLink { linker: s, .. }
            | Self::Ambiguous {
                modification: s, ..
            } => match &**s {
                SimpleModificationInner::Database {
                    id:
                        ModificationId {
                            ontology: Ontology::Unimod,
                            name,
                            ..
                        },
                    ..
                } => Some(name),
                _ => None,
            },
        }
    }
}

/// The structure to lookup ambiguous modifications, with a list of all modifications (the order is fixed) with for each modification their name and the actual modification itself (if already defined)
pub type AmbiguousLookup = Vec<AmbiguousLookupEntry>;

/// An entry in the ambiguous lookup
#[derive(Clone, Debug)]
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

impl Display for CrossLinkName {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Branch => write!(f, "#BRANCH"),
            Self::Name(n) => write!(f, "#XL{n}"),
        }
    }
}

/// Parse a custom modifications JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom database in previous version of
/// the library.
/// # Errors
/// If the string could not be parsed.
pub fn parse_custom_modifications(
    path: &Path,
) -> Result<CustomDatabase, BoxedError<'static, BasicKind>> {
    let string = std::fs::read_to_string(path).map_err(|err| {
        BoxedError::new(
            BasicKind::Error,
            "Could not parse custom modifications file",
            err.to_string(),
            Context::default().source(path.to_string_lossy()).to_owned(),
        )
    })?;
    CustomDatabase::from_json(&string)
}

/// Parse a custom modifications JSON string. The parser is guaranteed to be backwards compatible
/// with any JSON made by the serde serialisation of the custom database in previous version of
/// the library.
/// # Errors
/// If the string could not be parsed.
pub fn parse_custom_modifications_str(
    value: &str,
) -> Result<CustomDatabase, BoxedError<'static, BasicKind>> {
    CustomDatabase::from_json(value)
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_reading_custom_modifications_json_2024() {
    let data = include_str!("custom_modifications_2024.json");
    let mods = parse_custom_modifications_str(data).unwrap();
    assert!(mods.len() > 1);
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_reading_custom_modifications_json_2025() {
    let data = include_str!("custom_modifications_20250207.json");
    let mods = parse_custom_modifications_str(data).unwrap();
    assert!(mods.len() > 1);
}
