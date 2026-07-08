use std::{fmt::Write, num::NonZeroU16};

use crate::{
    chemistry::{Element, MolecularFormula},
    system::i8::Charge,
};

#[derive(Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct StructuralFormula {
    pub elements: Vec<(Option<Element>, Option<NonZeroU16>, Charge)>,
    pub connections: Vec<(usize, usize, Connection)>,
    // TODO: handle ambiguous connections
    // TODO: chimeric things
}

#[derive(Clone, Copy, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum Connection {
    #[default]
    SingleCovalent,
    DoubleCovalent,
    TripleCovalent,
    QuadrupleCovalent,
}

impl Connection {
    fn covalent_bonds(self) -> usize {
        match self {
            Self::SingleCovalent => 1,
            Self::DoubleCovalent => 2,
            Self::TripleCovalent => 3,
            Self::QuadrupleCovalent => 4,
        }
    }
}

impl StructuralFormula {
    // TODO: should be the default to_string as well
    fn to_smiles(&self) -> Option<String> {
        todo!()
    }

    pub fn to_dot(&self) -> Option<String> {
        let mut res = String::new();
        writeln!(&mut res, "graph MOL {{").unwrap();

        for (index, (element, isotope, c)) in self.elements.iter().enumerate() {
            writeln!(
                &mut res,
                "n{index} [label=\"{}{}{:+}\", shape=none, margin=0, fontcolor={}]",
                isotope.map(|i| i.to_string()).unwrap_or_default(),
                element.map_or("*", |e| e.symbol()),
                c.value,
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
            })
            .unwrap();
        }

        writeln!(&mut res, "}}").unwrap();
        Some(res)
    }

    /// Get the composition of this structure. It returns `None` if any isotope is invalid. Any
    /// unknown elements are ignored.
    pub fn composition(&self) -> Option<MolecularFormula> {
        self.elements
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
        for (i, (e, ..)) in self.elements.iter().enumerate() {
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

    /// Infer missing hydrogens in this graph.
    pub fn infer_hydrogens(&mut self) {
        let mut added = Vec::new();
        for (i, (e, isotope, charge)) in self.elements.iter().enumerate() {
            if charge.value != 0 {
                continue; // Poor way of detecting square bracket atoms in smiles
            }
            // TODO: does the charge need to be handled?
            let sum: usize = self
                .connections
                .iter()
                .filter(|c| c.0 == i || c.1 == i)
                .map(|c| c.2.covalent_bonds())
                .sum();
            let missing = missing_bonds(*e, sum);
            for _ in 0..missing {
                let offset = self.elements.len() + added.len();
                added.push((Some(Element::H), None, Charge::default()));
                self.connections.push((i, offset, Connection::SingleCovalent));
            }
        }
        self.elements.extend_from_slice(&added);
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
