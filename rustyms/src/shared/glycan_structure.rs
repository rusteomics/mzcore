use crate::helper_functions::{end_of_enclosure, next_char};
use std::ops::Range;

/// Rose tree representation of glycan structure
#[allow(dead_code)]
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
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

    /// Parse a short IUPAC glycan structure
    /// # Panics
    /// Panics if there is no single sugar found
    /// # Errors
    /// Errors when the format is not correct, could be unknown monosaccharide, or an open brace
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
}

impl Display for GlycanStructure {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.display_tree())
    }
}

impl Chemical for GlycanStructure {
    fn formula_inner(
        &self,
        sequence_index: crate::SequencePosition,
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
