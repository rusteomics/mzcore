//! Handle positioned glycan structures
use std::{hash::Hash, ops::RangeInclusive};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    annotation::model::{FragmentationModel, GlycanModel},
    chemistry::{CachedCharge, Chemical, MolecularFormula},
    fragment::{Fragment, FragmentType, GlycanBreakPos},
    quantities::Multi,
    sequence::{AminoAcid, SequencePosition},
    system::usize::Charge,
};

use crate::uom::num_traits::Zero;

use super::{glycan::MonoSaccharide, position::GlycanPosition};

/// The index in the branches as stored in the structure
pub type GlycanBranchIndex = usize;
/// The index in the branches when the branches are sorted on mass, this is used to properly render the names of the branches for human consumption
pub type GlycanBranchMassIndex = usize;

/// Rose tree representation of glycan structure
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct PositionedGlycanStructure {
    pub(super) sugar: MonoSaccharide,
    pub(super) branches: Vec<PositionedGlycanStructure>,
    pub(super) inner_depth: usize,
    pub(super) outer_depth: usize,
    /// The branches taken to get to this location (from the root) as the index in the branches and the index in the branches when sorted by mass.
    /// For a general glycan with a fucose on the first hexnac and a bisection after the core double
    /// hexnac + hex, this variable will contain an empty list for the root hexnac. For the fucose
    /// this variable will contain `[(0, 1)]` indicating it is the first branch in the structure but
    /// the second branch if the branches are sorted by mass. For the monosaccharides in the left
    /// bisection this variable will contain `[(1, 0), (0, 0)]`, indicating that it took the main
    /// branch (and not the fucose) and that it took the left branch for the second bisection which
    /// is heavier than the right branch.
    pub(super) branch: Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>,
}

impl Chemical for PositionedGlycanStructure {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.sugar.formula_inner(sequence_index, peptidoform_index)
            + self
                .branches
                .iter()
                .map(|f| f.formula_inner(sequence_index, peptidoform_index))
                .sum::<MolecularFormula>()
    }
}

impl PositionedGlycanStructure {
    /// All core options, with the Y breakage positions leading to this fragment
    pub fn core_options(
        &self,
        range: RangeInclusive<u8>,
        peptidoform_index: usize,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<(Vec<GlycanPosition>, MolecularFormula)> {
        self.internal_break_points(0, peptidoform_index, attachment)
            .iter()
            .filter(|(_, _, depth)| range.contains(depth))
            .map(|(f, pos, _)| {
                (
                    pos.iter()
                        .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                        .map(GlycanBreakPos::position)
                        .cloned()
                        .collect(),
                    f.clone(),
                )
            })
            .collect()
    }

    /// All possible bonds that can be broken and the molecular formula that would be held over if these bonds all broke and the broken off parts are lost.
    fn internal_break_points(
        &self,
        depth: u8,
        peptidoform_index: usize,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<(MolecularFormula, Vec<GlycanBreakPos>, u8)> {
        // Find every internal fragment ending at this bond (in a B breakage) (all bonds found are Y breakages and endings)
        // Walk through all branches and determine all possible breakages
        if self.branches.is_empty() {
            vec![
                (
                    self.formula_inner(SequencePosition::default(), peptidoform_index),
                    vec![GlycanBreakPos::End(self.position(attachment, false))],
                    depth + u8::from(!self.sugar.is_fucose()),
                ),
                (
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(self.position(attachment, false))],
                    depth,
                ),
            ]
        } else {
            self.branches
                .iter()
                .map(|b| {
                    b.internal_break_points(
                        depth + u8::from(!b.sugar.is_fucose()),
                        peptidoform_index,
                        attachment,
                    )
                }) // get all previous options
                .fold(Vec::new(), |accumulator, branch_options| {
                    if accumulator.is_empty() {
                        branch_options
                    } else {
                        let mut new_accumulator = Vec::new();
                        for base in &accumulator {
                            for option in &branch_options {
                                new_accumulator.push((
                                    &option.0 + &base.0,
                                    [option.1.clone(), base.1.clone()].concat(),
                                    option.2.max(base.2),
                                ));
                            }
                        }
                        new_accumulator
                    }
                })
                .into_iter()
                .map(|(m, b, d)| {
                    (
                        m + self
                            .sugar
                            .formula_inner(SequencePosition::default(), peptidoform_index),
                        b,
                        d,
                    )
                })
                .chain(std::iter::once((
                    // add the option of it breaking here
                    MolecularFormula::default(),
                    vec![GlycanBreakPos::Y(self.position(attachment, false))],
                    depth,
                )))
                .collect()
        }
    }

    /// Get uncharged diagnostic ions from all positions
    fn diagnostic_ions(
        &self,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        attachment: Option<(AminoAcid, SequencePosition)>,
        model: &GlycanModel,
    ) -> Vec<Fragment> {
        let mut output = self.sugar.diagnostic_ions(
            peptidoform_ion_index,
            peptidoform_index,
            crate::fragment::DiagnosticPosition::Glycan(
                self.position(attachment, false),
                self.sugar.clone(),
            ),
            true,
            model,
        );
        output.extend(self.branches.iter().flat_map(|b| {
            b.diagnostic_ions(peptidoform_ion_index, peptidoform_index, attachment, model)
        }));

        output
    }

    /// Generate all theoretical fragments for this glycan
    /// * `full_formula` the total formula of the whole peptide + glycan
    pub fn generate_theoretical_fragments(
        &self,
        model: &FragmentationModel,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment> {
        let charges_other = charge_carriers.range(model.glycan.other_charge_range);
        let charges_oxonium = charge_carriers.range(model.glycan.oxonium_charge_range);
        model
            .glycan
            .allow_structural
            .then(|| {
                // Get all base fragments from this node and all its children
                let mut base_fragments = self
                    .oxonium_fragments(peptidoform_ion_index, peptidoform_index, attachment)
                    .into_iter()
                    .flat_map(|f| f.with_charge_range_slice(&charges_oxonium))
                    .flat_map(|f| f.with_neutral_losses(&model.glycan.neutral_losses))
                    .collect_vec();
                // Generate all Y fragments
                base_fragments.extend(
                    self.internal_break_points(0, peptidoform_index, attachment)
                        .iter()
                        .filter(|(_, bonds, _)| {
                            bonds.iter().all(|b| !matches!(b, GlycanBreakPos::B(_)))
                                && !bonds.iter().all(|b| matches!(b, GlycanBreakPos::End(_)))
                        })
                        .flat_map(move |(f, bonds, _)| {
                            full_formula.iter().map(move |full| {
                                Fragment::new(
                                    full - self.formula_inner(
                                        SequencePosition::default(),
                                        peptidoform_index,
                                    ) + f,
                                    Charge::zero(),
                                    peptidoform_ion_index,
                                    peptidoform_index,
                                    FragmentType::Y(
                                        bonds
                                            .iter()
                                            .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                                            .map(GlycanBreakPos::position)
                                            .cloned()
                                            .collect(),
                                    ),
                                )
                            })
                        })
                        .flat_map(|f| f.with_charge_range_slice(&charges_other))
                        .flat_map(|f| f.with_neutral_losses(&model.glycan.neutral_losses)),
                );
                // Generate all diagnostic ions
                base_fragments.extend(
                    self.diagnostic_ions(
                        peptidoform_ion_index,
                        peptidoform_index,
                        attachment,
                        &model.glycan,
                    )
                    .into_iter()
                    .flat_map(|f| f.with_charge_range_slice(&charges_oxonium)),
                );
                base_fragments
            })
            .unwrap_or_default()
    }

    /// Generate all fragments without charge and neutral loss options
    fn oxonium_fragments(
        &self,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment> {
        // Find all B type fragments (with and without Y breakage)
        let mut base_fragments = self
            .internal_break_points(0, peptidoform_index, attachment)
            .into_iter()
            .filter(|(m, _, _)| *m != MolecularFormula::default())
            .map(|(formula, breakages, _)| {
                Fragment::new(
                    formula,
                    Charge::zero(),
                    peptidoform_ion_index,
                    peptidoform_index,
                    FragmentType::B {
                        b: self.position(attachment, true),
                        y: breakages
                            .iter()
                            .filter(|b| matches!(b, GlycanBreakPos::Y(_)))
                            .map(GlycanBreakPos::position)
                            .cloned()
                            .collect(),
                        end: breakages
                            .iter()
                            .filter(|b| matches!(b, GlycanBreakPos::End(_)))
                            .map(GlycanBreakPos::position)
                            .cloned()
                            .collect(),
                    },
                )
            })
            .collect_vec();

        // Extend with the theoretical fragments for all branches of this position
        base_fragments.extend(self.branches.iter().flat_map(|b| {
            b.oxonium_fragments(peptidoform_ion_index, peptidoform_index, attachment)
        }));
        base_fragments
    }

    fn position(
        &self,
        attachment: Option<(AminoAcid, SequencePosition)>,
        outer: bool,
    ) -> GlycanPosition {
        GlycanPosition {
            inner_depth: self.inner_depth,
            series_number: if outer {
                self.outer_depth + 1
            } else {
                self.inner_depth
            },
            branch: self.branch.clone(),
            attachment,
        }
    }
}
