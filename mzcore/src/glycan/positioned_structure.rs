//! Handle positioned glycan structures
use std::hash::Hash;

use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{Chemical, MolecularFormula},
    glycan::position::GlycanBreakPos,
    sequence::{AminoAcid, SequencePosition},
};

use super::{glycan::MonoSaccharide, position::GlycanPosition};

/// The index in the branches as stored in the structure
pub type GlycanBranchIndex = usize;
/// The index in the branches when the branches are sorted on mass, this is used to properly render the names of the branches for human consumption
pub type GlycanBranchMassIndex = usize;

/// Rose tree representation of glycan structure
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct PositionedGlycanStructure {
    /// The sugar at this level
    pub sugar: MonoSaccharide,
    /// All branches
    pub branches: Vec<Self>,
    /// The inner depth, the number of steps needed to get to the root / attached amino acid
    pub inner_depth: u16,
    /// The outer depth, the number of steps needed to get to the closest leaf node
    pub outer_depth: u16,
    /// The branches taken to get to this location (from the root) as the index in the branches and the index in the branches when sorted by mass.
    /// For a general glycan with a fucose on the first hexnac and a bisection after the core double
    /// hexnac + hex, this variable will contain an empty list for the root hexnac. For the fucose
    /// this variable will contain `[(0, 1)]` indicating it is the first branch in the structure but
    /// the second branch if the branches are sorted by mass. For the monosaccharides in the left
    /// bisection this variable will contain `[(1, 0), (0, 0)]`, indicating that it took the main
    /// branch (and not the fucose) and that it took the left branch for the second bisection which
    /// is heavier than the right branch.
    pub branch: Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>,
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
        range: (Option<usize>, Option<usize>),
        peptidoform_index: usize,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<(Vec<GlycanPosition>, MolecularFormula)> {
        self.internal_break_points(0, peptidoform_index, attachment)
            .iter()
            .filter(|(_, _, depth)| {
                range.0.is_none_or(|s| *depth as usize >= s)
                    && range.1.is_none_or(|e| *depth as usize <= e)
            })
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
    pub fn internal_break_points(
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

    /// Get the glycan position for this level
    pub fn position(
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
            branch: self.branch.clone().into(),
            attachment,
        }
    }
}
