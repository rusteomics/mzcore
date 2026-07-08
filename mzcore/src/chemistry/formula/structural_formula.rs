use std::{fmt::Write, num::NonZeroU16};

use crate::{
    chemistry::{Element, MolecularFormula},
    system::i8::Charge,
};

/// A structural chemical formula. This takes a graph based approach with separate nodes and edges.
/// Because of this approach a single structural formula can be used to describe multiple structures
/// that are not covalently bonded, as is expressed in SMILES using the dot `.` bond.
///
/// Chimeric information is currently not stored so will be deleted if read in from SMILES and
/// exported again.
#[derive(Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct StructuralFormula {
    /// The atoms (or nodes in the graph) with the element, isotope, and charge
    pub atoms: Vec<(Option<Element>, Option<NonZeroU16>, Charge)>,
    /// The bonds (or edges in the graph)
    pub connections: Vec<(usize, usize, Connection)>,
    // TODO: handle ambiguous connections
}

/// A bond between atoms
#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum Connection {
    /// A single covalent bond (σ)
    #[default]
    SingleCovalent,
    /// A double covalent bond (σ + π)
    DoubleCovalent,
    /// A triple covalent bond (σ + 2π)
    TripleCovalent,
    /// A quadruple bond (σ + 2π + δ)
    QuadrupleCovalent,
    /// An aromatic bond, every atom needs to have an even number of these (0 is even)
    Aromatic,
}

impl Connection {
    fn covalent_bonds(self) -> usize {
        match self {
            Self::SingleCovalent => 1,
            Self::DoubleCovalent => 2,
            Self::TripleCovalent => 3,
            Self::QuadrupleCovalent => 4,
            Self::Aromatic => 1, // Yes the correct answer is 1.5 but that is handled somewhere else
        }
    }
}

impl StructuralFormula {
    // TODO: should be the default to_string as well
    fn to_smiles(&self) -> Option<String> {
        todo!()
    }

    /// Get a graph in dot language to display this structure for debug purposes.
    pub fn to_dot(&self) -> String {
        let mut res = String::new();
        writeln!(&mut res, "graph MOL {{").unwrap();

        for (index, (element, isotope, c)) in self.atoms.iter().enumerate() {
            writeln!(
                &mut res,
                "n{index} [label=\"{}{}{}\", shape=none, margin=0, fontcolor={}]",
                isotope.map(|i| i.to_string()).unwrap_or_default(),
                element.map_or("*", |e| e.symbol()),
                if c.value == 0 {
                    String::new()
                } else {
                    format!("{:+}", c.value)
                },
                match element {
                    Some(Element::O) => "red",
                    Some(Element::N) => "blue",
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
                Connection::QuadrupleCovalent => "black:invis:black:invis:black:invis:black",
                Connection::Aromatic => "black:invis:grey",
            })
            .unwrap();
        }

        writeln!(&mut res, "}}").unwrap();
        res
    }

    /// Get the composition of this structure. It returns `None` if any isotope is invalid. Any
    /// unknown elements are ignored.
    pub fn composition(&self) -> Option<MolecularFormula> {
        self.atoms
            .iter()
            .fold(Some(MolecularFormula::default()), |acc, (e, i, c)| {
                e.map_or(acc.clone(), |e| {
                    acc.and_then(|acc| {
                        MolecularFormula::new(
                            &[(e, *i, 1), (Element::Electron, None, -c.value as i32)],
                            &[],
                        )
                        .map(|s| acc + s)
                    })
                })
            })
    }

    /// Infer missing elements in this graph, given a lookup list of groups for each amount of
    /// covalent bonds missing (1, 2, & 3). The given groups are each structural formulas with the
    /// element at the 0 index needing the correct covalent bond.
    // TODO: needs to be able to switch the group from OH on C to H on anything else
    pub fn infer(&mut self, group: [(Vec<(usize, Connection)>, Self); 3]) {
        let mut added = Vec::new();
        for (i, (e, ..)) in self.atoms.iter().enumerate() {
            // TODO: does the charge need to be handled?
            let sum: usize = self
                .connections
                .iter()
                .filter(|c| c.0 == i || c.1 == i)
                .map(|c| c.2.covalent_bonds())
                .sum();
            let missing = missing_bonds(*e, sum);
            if let Some(missing) = missing.checked_sub(1) {
                let (connections, g) = &group[missing];
                let offset = self.atoms.len() + added.len();
                added.extend_from_slice(&g.atoms);
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
        self.atoms.extend_from_slice(&added);
    }

    /// Infer missing hydrogens in this graph. But only on the given selection of atoms (to allow
    /// for radicals in SMILES).
    pub fn infer_hydrogens(&mut self, selection: &[usize]) {
        for index in selection {
            let element = self.atoms[*index].0;
            // Just assume that the number of aromatic bonds works out (i.e. is even), if this ever
            // is an odd number this is a violation of the assumptions and just plain
            // does not make sense.
            let sum: usize = self
                .connections
                .iter()
                .filter(|c| c.0 == *index || c.1 == *index)
                .fold((0, false), |(acc, prev_arom), (_, _, c)| {
                    if prev_arom && *c == Connection::Aromatic {
                        (
                            acc + c.covalent_bonds()
                                + usize::from(
                                    element != Some(Element::S),
                                    // Specifically S needs to be detected here to not have it add
                                    // uneccessary hydrogens, the bonds are not really 1.5 on bpth
                                    // sides but just 1
                                ),
                            false,
                        )
                    } else {
                        (
                            acc + c.covalent_bonds(),
                            prev_arom || *c == Connection::Aromatic,
                        )
                    }
                })
                .0;
            let missing = missing_bonds(element, sum);
            for _ in 0..missing {
                let new_index = self.atoms.len();
                self.atoms.push((Some(Element::H), None, Charge::default()));
                self.connections.push((*index, new_index, Connection::SingleCovalent));
            }
        }
    }

    /// Sort the connections to have the lowest index as the first index and all connections sorted
    /// on the first index then the second index.
    pub(super) fn normalise_connections(&mut self) {
        for (i1, i2, _) in self.connections.iter_mut() {
            *i1 = *i1.min(i2);
            *i2 = *i1.max(i2);
        }
        self.connections.sort_unstable();
    }
}

/// Follows the OpenSMILES specification for the valence
fn missing_bonds(element: Option<Element>, in_place: usize) -> usize {
    match element {
        Some(
            Element::H
            | Element::F
            | Element::Cl
            | Element::Br
            | Element::I
            | Element::At
            | Element::Ts,
        ) => 1_usize.saturating_sub(in_place),
        Some(Element::O) => 2_usize.saturating_sub(in_place),
        Some(Element::B) => 3_usize.saturating_sub(in_place),
        Some(Element::C) => 4_usize.saturating_sub(in_place),
        Some(Element::N | Element::P) => {
            if in_place > 3 { 5_usize } else { 3 }.saturating_sub(in_place)
        }
        Some(Element::S) => if in_place > 4 {
            6_usize
        } else if in_place > 2 {
            4
        } else {
            2
        }
        .saturating_sub(in_place),
        _ => 0,
    }
}
