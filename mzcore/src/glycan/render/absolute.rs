use itertools::Itertools;

use crate::{
    chemistry::Chemical,
    glycan::{
        GlycanBranchIndex, GlycanBranchMassIndex, GlycanPosition, GlycanSelection, GlycanStructure,
        RenderedGlycan,
        render::{
            element::GlycanRoot,
            shape::{Colour, Shape},
        },
    },
};

impl GlycanStructure {
    /// Render this glycan to the internal representation. This can then be rendered to SVG or a bitmap.
    ///  * `basis`: the text or symbol to draw at the root of the tree.
    ///  * `column_size`: the size (in pixels) of one block in the glycan, the full size with the padding and sugar size included.
    ///  * `sugar_size`: the size (in pixels) of a monosaccharide.
    ///  * `stroke_size`: the size (in pixels) of the strokes in the graphic.
    ///  * `direction`: the direction the draw the image in.
    ///  * `selection`: the selection of the glycan to draw, used to render fragments.
    ///  * `foreground`: the colour to be used for the foreground, in RGB order.
    ///  * `background`: the colour to be used for the background, in RGB order, this is used to fill 'empty' sugars if the isomeric state is unknown.
    ///  * `footnotes`: used to gather modification texts that are too big to place in line. The caller will have to find their own way of displaying this to the user.
    ///
    /// # Errors
    /// If the underlying buffer errors the error is returned. Otherwise `Ok(false)` is returned if the given `root_break` is not valid, and `Ok(true)` is returned if the rendering was fully successful.
    pub fn render<'a>(
        &'a self,
        basis: GlycanRoot,
        column_size: f32,
        sugar_size: f32,
        stroke_size: f32,
        direction: GlycanDirection,
        selection: GlycanSelection<'a>,
        foreground: [u8; 3],
        background: [u8; 3],
        footnotes: &'a mut Vec<String>,
    ) -> Option<RenderedGlycan> {
        self.position_absolute(0, &[], footnotes).render(
            basis,
            column_size,
            sugar_size,
            stroke_size,
            direction,
            selection,
            foreground,
            background,
            footnotes,
        )
    }

    /// Build the rendered glycan.
    /// # Panics
    /// If any monosaccharide does not have a defined mass.
    fn position_absolute(
        &self,
        depth: u16,
        path: &[(GlycanBranchIndex, GlycanBranchMassIndex)],
        footnotes: &mut Vec<String>,
    ) -> AbsolutePositionedGlycan {
        let (shape, colour, inner_modifications, outer_modifications) = self.sugar.get_shape();
        // Automatically make footnotes out of long outer modification texts
        let outer_modifications = if outer_modifications.len() > 6 {
            let index = footnotes.iter().position(|e| *e == outer_modifications);
            index.map_or_else(
                || {
                    let index = footnotes.len();
                    footnotes.push(outer_modifications);
                    OuterModifications::Footnote(index)
                },
                OuterModifications::Footnote,
            )
        } else if !outer_modifications.is_empty() {
            OuterModifications::Text(outer_modifications)
        } else {
            OuterModifications::Empty
        };

        if self.branches.is_empty() {
            AbsolutePositionedGlycan {
                y: 0,
                x: 0.0,
                mid_point: 0.5,
                width: 1.0,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                position: GlycanPosition {
                    inner_depth: depth,
                    series_number: depth,
                    branch: path.to_vec().into(),
                    attachment: None,
                },
                title: self.sugar.to_string(),
                branch_index: 0,
                branches: Vec::new(),
                sides: Vec::new(),
            }
        } else {
            let mut y_depth = 0;
            let mut branches = Vec::new();
            let mut sides = Vec::new();
            for (mass_index, (branch_index, branch)) in self
                .branches
                .iter()
                .enumerate()
                .sorted_unstable_by(|(_, a), (_, b)| {
                    b.formula()
                        .monoisotopic_mass()
                        .partial_cmp(&a.formula().monoisotopic_mass())
                        .unwrap()
                })
                .enumerate()
                .sorted_unstable_by(|a, b| a.1.0.cmp(&b.1.0))
            {
                let mut new_path = path.to_vec();
                new_path.push((branch_index, mass_index));
                let mut rendered = branch.position_absolute(
                    depth + 1,
                    if self.branches.len() > 1 {
                        &new_path
                    } else {
                        path
                    },
                    footnotes,
                );
                rendered.branch_index = branch_index;
                if rendered.is_sideways() && sides.len() < 2 {
                    if sides.is_empty() && rendered.shape == Shape::Triangle {
                        rendered.shape = Shape::LeftPointingTriangle;
                    } else if sides.len() == 1 && rendered.shape == Shape::Triangle {
                        rendered.shape = Shape::RightPointingTriangle;
                    }
                    sides.push(rendered);
                } else {
                    y_depth = y_depth.max(rendered.y);
                    branches.push(rendered);
                }
            }
            // Update all branch placements
            let mut displacement = 0.0;
            for branch in &mut branches {
                branch.transpose(y_depth - branch.y, displacement);
                displacement += branch.width;
            }
            if !branches.is_empty() {
                y_depth += 1;
            }
            // Determine the center point for this sugar
            let mut center = match branches.len() {
                0 => 0.5,
                1 => branches[0].mid_point,
                n => {
                    // Find the median midpoint of the branches
                    (branches[n / 2 - (n + 1) % 2].x
                        + branches[n / 2 - (n + 1) % 2].mid_point
                        + branches[n / 2].x
                        + branches[n / 2].mid_point)
                        / 2.0
                }
            };
            let mut width = branches.last().map_or(1.0, |b| b.x + b.width);
            if !sides.is_empty() {
                sides[0].transpose(y_depth, center + 0.5);
                width = width.max(center + 0.5 + sides[0].width);
            }
            if sides.len() == 2 {
                let mut x = center - 0.5 - sides[1].width;
                if x < 0.0 {
                    let shift = -x;
                    center += shift;
                    for branch in &mut branches {
                        branch.transpose(0, shift);
                    }
                    sides[0].transpose(0, shift);
                    width += shift;
                    x = 0.0;
                }
                sides[1].transpose(y_depth, x);
            }
            AbsolutePositionedGlycan {
                y: y_depth,
                x: 0.0,
                mid_point: center,
                width,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                position: GlycanPosition {
                    inner_depth: depth,
                    series_number: depth,
                    branch: path.to_vec().into(),
                    attachment: None,
                },
                title: self.sugar.to_string(),
                branch_index: 0,
                branches,
                sides,
            }
        }
    }
}

/// An absolute positioned glycan.
#[derive(Clone, Debug)]
pub(super) struct AbsolutePositionedGlycan {
    /// The depth of this sugar along the main axis of the glycan, starting at 0 at the top (in the leaves)
    pub(super) y: usize,
    /// The sideways placement of this whole tree starting at 0 at the leftmost monosaccharide, 1.0 is the width of one monosaccharide
    pub(super) x: f32,
    /// The sideways placement of this sugar within this tree, for the absolute sideways placement of this sugar add this to `x`
    pub(super) mid_point: f32,
    /// The total width of the (sub)tree with all of its branches and sides
    pub(super) width: f32,
    /// The shape of the monosaccharide
    pub(super) shape: Shape,
    /// The colour of the monosaccharide
    pub(super) colour: Colour,
    /// Text to be shown inside the monosaccharide
    pub(super) inner_modifications: String,
    /// Text to be shown outside the monosaccharide
    pub(super) outer_modifications: OuterModifications,
    /// The position of this sugar
    pub(super) position: GlycanPosition,
    /// Full name of the glycan
    pub(super) title: String,
    /// The index into the branches of the parent monosaccharide
    pub(super) branch_index: usize,
    /// All branches that go up the tree
    pub(super) branches: Vec<AbsolutePositionedGlycan>,
    /// All branches that go to the side (Fucoses)
    pub(super) sides: Vec<AbsolutePositionedGlycan>,
}

#[derive(Clone, Debug)]
/// Modifications that are to be shown outside of the monosaccharide
pub(super) enum OuterModifications {
    /// Too long of a text, or it did not fit, so show as a footnote
    Footnote(usize),
    /// Text
    Text(String),
    /// No modification
    Empty,
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
/// The direction of the rendered glycan
pub enum GlycanDirection {
    /// A top down tree, with the root at the bottom
    TopDown,
    /// A left to right tree, with the root at the right
    LeftToRight,
}

/// A subtree of a rendered glycan, used to restrict the canvas for glycan fragments
#[derive(Clone, Debug)]
pub(super) struct SubTree<'a> {
    /// The root for this sub tree
    pub(super) tree: &'a AbsolutePositionedGlycan,
    /// Total depth of the glycans with the breaks applied
    pub(super) depth: usize,
    /// The horizontal offset from the left
    pub(super) left_offset: f32,
    /// The horizontal offset from the right
    pub(super) right_offset: f32,
    /// If this fragment is topped by a breaking symbol, needed to calculate the correct height for the canvas
    pub(super) break_top: bool,
    /// If this fragment is bottomed by a breaking symbol, needed to calculate the correct height for the canvas
    pub(super) break_bottom: bool,
    /// All breaking branches, standardised to the linked root
    pub(super) branch_breaks: Vec<(u16, Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>)>,
}

impl AbsolutePositionedGlycan {
    /// Transpose this glycan and all of its branches
    fn transpose(&mut self, y: usize, x: f32) {
        self.y += y;
        self.x += x;
        for branch in &mut self.branches {
            branch.transpose(y, x);
        }
        for side in &mut self.sides {
            side.transpose(y, x);
        }
    }

    /// Check if this sugar should be rendered to the side of the parent sugar
    fn is_sideways(&self) -> bool {
        self.colour == Colour::Red
            && self.shape == Shape::Triangle
            && self.branches.is_empty()
            && self.sides.is_empty()
    }

    /// Get the subtree starting on the given position, return None if the starting position is not valid, it also indicates the depth of this subtree for the given branch breakages and if a break tops the structure
    pub(super) fn get_subtree<'a>(&'a self, selection: GlycanSelection<'a>) -> Option<SubTree<'a>> {
        /// Calculate the maximal depth, break top, left and right offset
        fn canvas_size(
            tree: &AbsolutePositionedGlycan,
            breakages: &[(u16, Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>)],
        ) -> (usize, bool, f32, f32) {
            let lx = (tree.x + tree.mid_point - 0.5).max(0.0);
            let rx = (tree.width - tree.mid_point - 0.5).max(0.0);
            // The tree is cut here
            if breakages.iter().any(|b| b.0 == 0) {
                return (0, true, lx, rx);
            }

            let total_branches = tree.branches.len() + tree.sides.len();
            let (depth, break_top, left_offset, right_offset) = match total_branches {
                0 => (0, false, lx, rx),
                1 => tree.branches.first().map_or((0, false, lx, rx), |branch| {
                    canvas_size(
                        branch,
                        &breakages
                            .iter()
                            .map(|b| (b.0 - 1, b.1.clone()))
                            .collect_vec(),
                    )
                }),
                _ => tree
                    .branches
                    .iter()
                    .enumerate()
                    .map(|(i, branch)| {
                        (
                            i,
                            canvas_size(
                                branch,
                                &breakages
                                    .iter()
                                    .filter(|b| {
                                        b.1.first().map(|b| b.0) == Some(branch.branch_index)
                                    })
                                    .map(|b| (b.0 - 1, b.1[1..].to_vec()))
                                    .collect_vec(),
                            ),
                        )
                    })
                    .fold((0, false, lx, rx), |acc, (i, v)| {
                        (
                            acc.0.max(v.0),
                            if v.0 >= acc.0 { v.1 } else { acc.1 },
                            if i == 0 { v.2 } else { acc.2 },
                            if i == tree.branches.len() - 1 {
                                v.3
                            } else {
                                acc.3
                            },
                        )
                    }),
            };
            (
                depth + 1,
                break_top,
                if tree.sides.len() == 2 {
                    left_offset.min(tree.x + tree.mid_point - 1.5).max(0.0)
                } else {
                    left_offset
                },
                if tree.sides.is_empty() {
                    right_offset
                } else {
                    right_offset.min(tree.width - tree.mid_point - 1.5).max(0.0)
                },
            )
        }

        let (tree, rules, break_bottom) = match selection {
            GlycanSelection::Subtree(root, branch_breaks) => {
                let start = root.unwrap_or(&self.position);
                let mut tree = self;
                let mut depth = 0;
                let mut branch_choices = start.branch.clone();
                branch_choices.reverse();
                while depth < start.inner_depth {
                    depth += 1;

                    let total_branches = tree.branches.len() + tree.sides.len();
                    match total_branches {
                        0 => return None,
                        1 => tree = tree.branches.first().or_else(|| tree.sides.first())?,
                        _ => {
                            let index = branch_choices.pop()?;
                            tree = tree
                                .branches
                                .iter()
                                .chain(tree.sides.iter())
                                .find(|b| b.branch_index == index.0)?;
                        }
                    }
                }

                let rules = branch_breaks
                    .iter()
                    .filter(|b| {
                        b.inner_depth >= start.inner_depth && b.branch.starts_with(&start.branch)
                    })
                    .map(|b| {
                        (
                            b.inner_depth - start.inner_depth,
                            b.branch[start.branch.len()..].to_vec(),
                        )
                    })
                    .collect_vec();
                (tree, rules, root.is_some())
            }
            GlycanSelection::SingleSugar(position) => {
                let mut tree = self;
                let mut depth = 0;
                let mut branch_choices = position.branch.clone();
                branch_choices.reverse();
                while depth < position.inner_depth {
                    depth += 1;

                    let total_branches = tree.branches.len() + tree.sides.len();
                    match total_branches {
                        0 => return None,
                        1 => tree = tree.branches.first().or_else(|| tree.sides.first())?,
                        _ => {
                            let index = branch_choices.pop()?;
                            tree = tree
                                .branches
                                .iter()
                                .find(|b| b.branch_index == index.0)
                                .or_else(|| {
                                    tree.sides.iter().find(|b| b.branch_index == index.0)
                                })?;
                        }
                    }
                }

                let rules = tree
                    .branches
                    .iter()
                    .chain(tree.sides.iter())
                    .map(|b| {
                        (1, vec![(b.branch_index, b.branch_index)])
                        // TODO: the mass_index should be stored here, but currently that is unused so for now this does not introduce incorrect behaviour
                    })
                    .collect_vec();
                (tree, rules, true)
            }
        };
        let (depth, break_top, left_offset, right_offset) = canvas_size(tree, &rules);
        Some(SubTree {
            tree,
            depth,
            left_offset,
            right_offset,
            break_top,
            break_bottom,
            branch_breaks: rules,
        })
    }
}
