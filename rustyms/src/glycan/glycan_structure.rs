//! Handle glycan structures
use std::{
    ops::Range,
    str::FromStr,
    {fmt::Display, hash::Hash},
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{Chemical, MolecularFormula},
    error::{Context, CustomError},
    glycan::{
        glycan::{BaseSugar, MonoSaccharide},
        lists::GLYCAN_PARSE_LIST,
        {GlycanBranchIndex, GlycanBranchMassIndex, GlycanPosition, PositionedGlycanStructure},
    },
    helper_functions::{end_of_enclosure, next_char},
    sequence::SequencePosition,
};

/// Rose tree representation of glycan structure
#[allow(dead_code)]
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct GlycanStructure {
    pub(super) sugar: MonoSaccharide,
    pub(super) branches: Vec<GlycanStructure>,
}

impl GlycanStructure {
    /// Create a new glycan structure
    #[allow(dead_code)]
    pub const fn new(sugar: MonoSaccharide, branches: Vec<Self>) -> Self {
        Self { sugar, branches }
    }

    /// Parse a short IUPAC glycan structure.
    /// # Panics
    /// If there is no single sugar found
    /// # Errors
    /// When the format is not correct, could be unknown monosaccharide, or an open brace
    pub fn from_short_iupac(
        line: &str,
        range: Range<usize>,
        line_index: usize,
    ) -> Result<Self, CustomError> {
        let mut offset = range.start;
        let mut branch = Self {
            sugar: MonoSaccharide::new(BaseSugar::Decose, &[]),
            branches: Vec::new(),
        }; // Starting sugar, will be removed
        let mut last_branch: &mut Self = &mut branch;
        let bytes = line.as_bytes();

        while offset < range.end {
            while bytes[offset] == b'[' {
                let end = end_of_enclosure(line, offset + 1, b'[', b']').ok_or_else(|| {
                    CustomError::error(
                        "Invalid iupac short glycan",
                        "No closing brace for branch",
                        Context::line(Some(line_index), line, offset, range.end - offset),
                    )
                })?;
                last_branch.branches.push(Self::from_short_iupac(
                    line,
                    offset + 1..end,
                    line_index,
                )?);

                offset = end + 1;
            }
            let (sugar, new_offset) = MonoSaccharide::from_short_iupac(line, offset, line_index)?;
            offset = new_offset;

            last_branch.branches.push(Self {
                sugar: sugar.clone(),
                branches: Vec::new(),
            });
            last_branch = last_branch.branches.last_mut().unwrap();

            offset = Self::ignore_linking_information(bytes, offset, &range);
        }
        branch
            .branches
            .pop()
            .map_or_else(
                || {
                    Err(CustomError::error(
                        "Invalid iupac short glycan",
                        "No glycan found",
                        Context::line(Some(line_index), line.to_string(), range.start, range.len()),
                    ))
                },
                Ok,
            )
            .map(Self::reroot)
    }

    /// # Panics
    /// It panics if a brace was not closed that was not close to the end of the input (more then 10 bytes from the end).
    fn ignore_linking_information(bytes: &[u8], mut offset: usize, range: &Range<usize>) -> usize {
        if offset < bytes.len() && bytes[offset] == b'(' {
            if let Some(end) = next_char(bytes, offset + 1, b')') {
                offset = end + 1; // just ignore all linking stuff I do not care
            } else {
                // This only happens for incomplete branches where the last parts of the branch are unknown.
                assert!(range.end - offset < 10); // make sure it is the last part
                offset = range.end; // assume it is the last not closed brace
            }
        }
        offset
    }

    /// Inverts the tree, gets a tree where the an outer branch is chosen as root.
    /// It inverts it by choosing the last (rightmost) branch as new root.
    /// # Panics
    /// If there is no sugar in the starting structure.
    fn reroot(self) -> Self {
        let mut new_structure: Option<Vec<Self>> = None;
        let mut old_structure = Some(self);

        while let Some(mut old) = old_structure.take() {
            // Define new sugar
            let mut new_sugar = Self {
                sugar: old.sugar,
                branches: Vec::new(),
            };
            // If there is already some info in the new structure add that as a branch
            if let Some(new_structure) = new_structure {
                for branch in new_structure {
                    new_sugar.branches.push(branch);
                }
            }
            let mut new_branches = vec![new_sugar];
            // Take the last branch from the old sugar
            if let Some(last) = old.branches.pop() {
                old_structure = Some(last);
            }
            // Put all the other old branches on the new sugar
            for branch in old.branches {
                new_branches.push(branch);
            }
            new_structure = Some(new_branches);
        }

        let mut new = new_structure.unwrap();
        assert_eq!(new.len(), 1);
        new.pop().unwrap()
    }

    /// Recursively show the structure of this glycan
    fn display_tree(&self) -> String {
        if self.branches.is_empty() {
            self.sugar.to_string()
        } else {
            format!(
                "{}({})",
                self.sugar,
                self.branches.iter().map(Self::display_tree).join(",")
            )
        }
    }

    /// Check if this structure contains the given monosaccharide, see [`MonoSaccharide::equivalent`]
    /// for details on how the matching happens. It is possible to only select a small part of this
    /// glycan using the `root_break` and `branch_break` parameters. See
    /// [`Subtree`](crate::glycan::render::GlycanSelection::Subtree) for more info on how this
    /// selection works.
    pub fn contains(
        &self,
        monosaccharide: &MonoSaccharide,
        precise: bool,
        root_break: Option<&GlycanPosition>,
        branch_breaks: &[GlycanPosition],
    ) -> Option<bool> {
        fn check_subtree(
            target: &MonoSaccharide,
            precise: bool,
            element: &GlycanStructure,
            breaks: &[(usize, Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>)],
        ) -> bool {
            let total_branches = element.branches.len();
            if element.sugar.equivalent(target, precise) {
                return true;
            }
            for (index, branch) in element.branches.iter().enumerate() {
                if !((total_branches == 1 && breaks.iter().any(|b| b.0 == 1))
                    || breaks
                        .iter()
                        .any(|b| b.0 == 1 && b.1.first().map(|b| b.0) == Some(index)))
                    && check_subtree(
                        target,
                        precise,
                        branch,
                        &breaks
                            .iter()
                            .filter(|b| {
                                (total_branches > 1 && b.1.first().map(|b| b.0) == Some(index)
                                    || total_branches == 1)
                                    && b.0 > 0
                            })
                            .map(|b| (b.0 - 1, b.1[usize::from(total_branches > 1)..].to_vec()))
                            .collect_vec(),
                    )
                {
                    return true;
                }
            }
            false
        }

        let base = GlycanPosition {
            inner_depth: 0,
            series_number: 0,
            branch: Vec::new(),
            attachment: None,
        };
        let start = root_break.unwrap_or(&base);
        let mut tree = self;
        let mut depth = 0;
        let mut branch_choices = start.branch.clone();
        branch_choices.reverse();
        while depth < start.inner_depth {
            depth += 1;

            let total_branches = tree.branches.len();
            match total_branches {
                0 => return None,
                1 => tree = tree.branches.first()?,
                _ => {
                    let index = branch_choices.pop()?;
                    tree = tree
                        .branches
                        .iter()
                        .enumerate()
                        .find(|(branch_index, _)| *branch_index == index.0)?
                        .1;
                }
            }
        }

        let rules = branch_breaks
            .iter()
            .filter(|b| b.inner_depth >= start.inner_depth && b.branch.starts_with(&start.branch))
            .map(|b| {
                (
                    b.inner_depth - start.inner_depth,
                    b.branch[start.branch.len()..].to_vec(),
                )
            })
            .collect_vec();
        Some(check_subtree(monosaccharide, precise, tree, &rules))
    }
}

impl Display for GlycanStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display_tree())
    }
}

impl Chemical for GlycanStructure {
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

impl FromStr for GlycanStructure {
    type Err = CustomError;
    /// Parse a textual structure representation of a glycan (outside ProForma format)
    /// Example: Hex(Hex(HexNAc)) => Hex-Hex-HexNAc (linear)
    /// Example: Hex(Fuc,Hex(HexNAc,Hex(HexNAc)))
    ///          =>  Hex-Hex-HexNAc
    ///              └Fuc  └Hex-HexNAc
    /// # Errors
    /// Return an Err if the format is not correct
    fn from_str(line: &str) -> Result<Self, CustomError> {
        Self::parse(line, 0..line.len())
    }
}

impl GlycanStructure {
    /// Parse a textual structure representation of a glycan (outside ProForma format)
    /// Example: Hex(Hex(HexNAc)) => Hex-Hex-HexNAc (linear)
    /// Example: Hex(Fuc,Hex(HexNAc,Hex(HexNAc)))
    ///          =>  Hex-Hex-HexNAc
    ///              └Fuc  └Hex-HexNAc
    /// # Errors
    /// Return an Err if the format is not correct
    pub fn parse(line: &str, range: Range<usize>) -> Result<Self, CustomError> {
        Self::parse_internal(line, range).map(|(g, _)| g)
    }

    /// # Errors
    /// Return an Err if the format is not correct
    fn parse_internal(line: &str, range: Range<usize>) -> Result<(Self, usize), CustomError> {
        // Parse at the start the first recognised glycan name
        if let Some(name) = GLYCAN_PARSE_LIST
            .iter()
            .find(|name| line[range.clone()].starts_with(&name.0))
        {
            // If the name is followed by a bracket parse a list of branches
            let index = range.start + name.0.len();
            if line.as_bytes()[index] == b'(' {
                // Find the end of this list
                let end = end_of_enclosure(line, index + 1, b'(', b')').ok_or_else(|| {
                    CustomError::error(
                        "Invalid glycan branch",
                        "No valid closing delimiter",
                        Context::line(None, line, index, 1),
                    )
                })?;
                // Parse the first branch
                let mut index = index + 1;
                let mut branches = Vec::new();
                let (glycan, pos) = Self::parse_internal(line, index..end)?;
                index = pos;
                branches.push(glycan);
                // Keep parsing until the end of this branch level (until the ')' is reached)
                while index < end {
                    if line.as_bytes()[index] != b',' {
                        return Err(CustomError::error(
                            "Invalid glycan structure",
                            "Branches should be separated by commas ','",
                            Context::line(None, line, index, 1),
                        ));
                    }
                    index += 1;
                    let (glycan, pos) = Self::parse_internal(line, index..end)?;
                    branches.push(glycan);
                    index = pos;
                }
                Ok((
                    Self {
                        sugar: name.1.clone(),
                        branches,
                    },
                    end + 1,
                ))
            } else {
                Ok((
                    Self {
                        sugar: name.1.clone(),
                        branches: Vec::new(),
                    },
                    range.start + name.0.len(),
                ))
            }
        } else {
            Err(CustomError::error(
                "Could not parse glycan structure",
                "Could not parse the following part",
                Context::line(None, line, range.start, range.len()),
            ))
        }
    }

    /// Annotate all positions in this tree with all positions
    pub fn determine_positions(self) -> PositionedGlycanStructure {
        self.internal_pos(0, &[]).0
    }

    /// Given the inner depth determine the correct positions and branch ordering
    /// Return the positioned tree and the outer depth.
    /// # Panics
    /// When any of the masses in this glycan cannot be compared see [`f64::partial_cmp`].
    fn internal_pos(
        self,
        inner_depth: usize,
        branch: &[(GlycanBranchIndex, GlycanBranchMassIndex)],
    ) -> (PositionedGlycanStructure, usize) {
        // Sort the branches on decreasing molecular weight
        let branches = self
            .branches
            .into_iter()
            .enumerate()
            .sorted_unstable_by(|(_, a), (_, b)| {
                b.formula()
                    .monoisotopic_mass()
                    .partial_cmp(&a.formula().monoisotopic_mass())
                    .unwrap()
            })
            .collect_vec();

        // Get the correct branch indices adding a new layer of indices when needed
        let branches: Vec<(PositionedGlycanStructure, usize)> = if branches.len() == 1 {
            branches
                .into_iter()
                .map(|(_, b)| b.internal_pos(inner_depth + 1, branch))
                .collect()
        } else {
            branches
                .into_iter()
                .enumerate()
                .map(|(mass_index, (index, b))| {
                    let mut new_branch = branch.to_vec();
                    new_branch.push((index, mass_index));
                    b.internal_pos(inner_depth + 1, &new_branch)
                })
                .collect()
        };

        let outer_depth = branches.iter().map(|b| b.1).max().unwrap_or(0);
        (
            PositionedGlycanStructure {
                sugar: self.sugar,
                branches: branches.into_iter().map(|b| b.0).collect(),
                branch: branch.to_vec(),
                inner_depth,
                outer_depth,
            },
            outer_depth + 1,
        )
    }

    /// Get the composition of a `GlycanStructure`. The result is normalised (sorted and deduplicated).
    /// # Panics
    /// If one monosaccharide species has occurrence outside the range of [`isize::MIN`] to [`isize::MAX`].
    pub fn composition(&self) -> Vec<(MonoSaccharide, isize)> {
        let composition = self.composition_inner();
        MonoSaccharide::simplify_composition(composition)
            .expect("One monosaccharide species has a number outside of the range of isize")
    }

    /// Get the composition in monosaccharides of this glycan
    fn composition_inner(&self) -> Vec<(MonoSaccharide, isize)> {
        let mut output = vec![(self.sugar.clone(), 1)];
        output.extend(self.branches.iter().flat_map(Self::composition_inner));
        output
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    use super::*;
    use crate::molecular_formula;

    #[test]
    fn parse_glycan_structure_01() {
        assert_eq!(
            GlycanStructure::from_str("hep(hex)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                    branches: Vec::new()
                }],
            }
        );
    }

    #[test]
    fn parse_glycan_structure_02() {
        assert_eq!(
            GlycanStructure::from_str("hex(hex,hep)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: Vec::new()
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                        branches: Vec::new()
                    }
                ],
            }
        );
    }

    #[test]
    fn parse_glycan_structure_03() {
        assert_eq!(
            GlycanStructure::from_str("hex(hex(hex),hep)").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                branches: vec![
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: vec![GlycanStructure {
                            sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                .with_name("Hex"),
                            branches: Vec::new()
                        }]
                    },
                    GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                        branches: Vec::new()
                    }
                ],
            }
        );
    }

    #[test]
    fn parse_glycan_structure_04() {
        assert_eq!(
            GlycanStructure::from_str("hep(hex(hex(hex(hep),hex)))").unwrap(),
            GlycanStructure {
                sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[]).with_name("Hep"),
                branches: vec![GlycanStructure {
                    sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                    branches: vec![GlycanStructure {
                        sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[]).with_name("Hex"),
                        branches: vec![
                            GlycanStructure {
                                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                    .with_name("Hex"),
                                branches: vec![GlycanStructure {
                                    sugar: MonoSaccharide::new(BaseSugar::Heptose(None), &[])
                                        .with_name("Hep"),
                                    branches: Vec::new(),
                                }],
                            },
                            GlycanStructure {
                                sugar: MonoSaccharide::new(BaseSugar::Hexose(None), &[])
                                    .with_name("Hex"),
                                branches: Vec::new(),
                            },
                        ],
                    }],
                }],
            }
        );
    }

    #[test]
    fn correct_masses() {
        let (sugar, _) = MonoSaccharide::from_short_iupac("Neu5Ac", 0, 0).unwrap();
        dbg!(&sugar);

        assert_eq!(sugar.formula(), molecular_formula!(C 11 H 17 N 1 O 8));
    }

    #[test]
    fn correct_structure_g43728nl() {
        // Furanoses added for error detection
        let structure = GlycanStructure::from_short_iupac(
            "Neu5Ac(?2-?)Galf(?1-?)GlcNAc(?1-?)Man(?1-?)[Galf(?1-?)GlcNAc(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc", 
            0..101,
            0
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(HexNAc(Hex(Hex(HexNAc(Hexf(NonAcAN))),Hex(HexNAc(Hexf)))))"
        );
    }

    #[test]
    fn correct_structure_g36564am() {
        let structure = GlycanStructure::from_short_iupac(
            "Gal(?1-?)GlcNAc(?1-?)Man(?1-?)[GlcNAc(?1-?)Man(?1-?)][GlcNAc(?1-?)]Man(?1-?)GlcNAc",
            0..82,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(Hex(Hex(HexNAc(Hex)),Hex(HexNAc),HexNAc))"
        );
    }
    #[test]
    fn correct_structure_g04605kt() {
        let structure = GlycanStructure::from_short_iupac(
            "L-GlcNAc(b1-2)L-Man(a1-3)[GlcNAc(b1-4)][L-Gal(b1-4)GlcNAc(b1-2)L-Man(a1-6)]Man(b1-4)L-GlcNAc(b1-4)GlcNAc(b1-",
            0..108,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(HexNAc(Hex(Hex(HexNAc),HexNAc,Hex(HexNAc(Hex)))))"
        );
    }

    #[test]
    fn correct_structure_g67881ee() {
        // Fully specified version of g36564am
        // Furanoses added for error detection
        let structure = GlycanStructure::from_short_iupac(
            "GlcNAc(b1-2)Man(a1-3)[GlcNAc(b1-4)][Galf(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-",
            0..87,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(Hex(Hex(HexNAc),HexNAc,Hex(HexNAc(Hexf))))"
        );
    }

    #[test]
    fn correct_structure_g11771hd() {
        let structure = GlycanStructure::from_short_iupac(
            "GlcNAc(?1-?)[GlcNAc(?1-?)]Man(?1-?)[Man(?1-?)Man(?1-?)]Man(?1-?)GlcNAc(?1-?)GlcNAc(?1-",
            0..86,
            0,
        )
        .unwrap();

        assert_eq!(
            structure.to_string(),
            "HexNAc(HexNAc(Hex(Hex(HexNAc,HexNAc),Hex(Hex))))"
        );
    }
}
