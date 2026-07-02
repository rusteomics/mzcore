use std::{fmt::Write, num::NonZeroU16};

use crate::chemistry::{Element, MolecularFormula};

#[derive(Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(super) struct StructuralFormula {
    pub(super) elements: Vec<(Element, Option<NonZeroU16>)>, // TODO: handle *
    pub(super) connections: Vec<(usize, usize, Connection)>,
    // TODO: handle ambiguous connections
    // TODO: chimeric things
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub(super) enum Connection {
    SingleCovalent,
    DoubleCovalent,
    TripleCovalent,
}

impl Connection {
    fn covalent_bonds(self) -> usize {
        match self {
            Self::SingleCovalent => 1,
            Self::DoubleCovalent => 2,
            Self::TripleCovalent => 3,
        }
    }
}

impl StructuralFormula {
    // TODO: should be the default to_string as well
    fn to_smiles(&self) -> Option<String> {
        todo!()
    }

    pub(super) fn to_dot(&self) -> Option<String> {
        let mut res = String::new();
        writeln!(&mut res, "graph MOL {{").unwrap();

        for (index, (element, isotope)) in self.elements.iter().enumerate() {
            writeln!(
                &mut res,
                "n{index} [label=\"{}{element}\", shape=none, margin=0, fontcolor={}]",
                isotope.map(|i| i.to_string()).unwrap_or_default(),
                match element {
                    Element::O => "red",
                    Element::N => "blue",
                    _ => "black",
                }
            )
            .unwrap();
        }

        for (i1, i2, ty) in &self.connections {
            writeln!(&mut res, "n{i1} -- n{i2} [color=\"{}\"]", match ty {
                Connection::SingleCovalent => "black",
                Connection::DoubleCovalent => "black:invis:black",
                Connection::TripleCovalent => "black:invis:black:invis:black",
            })
            .unwrap();
        }

        writeln!(&mut res, "}}").unwrap();
        Some(res)
    }

    /// Get the composition of this structure. It returns `None` if any isotope is invalid.
    pub(super) fn composition(&self) -> Option<MolecularFormula> {
        self.elements
            .iter()
            .fold(Some(MolecularFormula::default()), |acc, (e, i)| {
                acc.and_then(|acc| MolecularFormula::new(&[(*e, *i, 1)], &[]).map(|s| acc + s))
            })
    }

    /// Infer missing elements in this graph, given a lookup list of groups for each amount of
    /// covalent bonds missing (1, 2, & 3). The given groups are each structural formulas with the
    /// element at the 0 index needing the correct covalent bond.
    // TODO: needs to be able to switch the group from OH on C to H on anything else
    pub(super) fn infer(&mut self, group: [(Vec<(usize, Connection)>, Self); 3]) {
        let mut added = Vec::new();
        for (i, (e, _)) in self.elements.iter().enumerate() {
            let sum: usize = self
                .connections
                .iter()
                .filter(|c| c.0 == i || c.1 == i)
                .map(|c| c.2.covalent_bonds())
                .sum();
            let missing = missing_bonds(*e, sum);
            if let Some(missing) = missing.checked_sub(1) {
                let (connections, g) = &group[missing];
                let offset = self.elements.len() + added.len();
                added.extend_from_slice(&g.elements);
                for c in &g.connections {
                    self.connections.push((c.0 + offset, c.1 + offset, c.2));
                }
                for c in connections {
                    self.connections.push((c.0 + offset, i, c.1));
                }

                // TODO: handle if there are more empty spots than defined in the group, maybe also
                // allow only having H as fill group
            }
        }
        self.elements.extend_from_slice(&added);
    }
}

fn missing_bonds(element: Element, in_place: usize) -> usize {
    match element {
        Element::H => 1_usize.saturating_sub(in_place),
        Element::C => 4_usize.saturating_sub(in_place),
        Element::O => 2_usize.saturating_sub(in_place),
        Element::N => 3_usize.saturating_sub(in_place),
        _ => todo!(),
    }
}
