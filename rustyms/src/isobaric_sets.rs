use std::cmp::Ordering;

use itertools::Itertools;

use crate::{
    annotation::model::GlycanModel,
    chemistry::Chemical,
    quantities::Tolerance,
    sequence::{
        AminoAcid, CheckedAminoAcid, Modification, Peptidoform, PlacementRule, Position,
        SemiAmbiguous, SequenceElement, SequencePosition, SimpleLinear, SimpleModification,
        SimpleModificationInner,
    },
    system::{Mass, Ratio, fraction},
};

/// A list of building blocks for a sequence defined by its sequence elements and its mass.
pub type BuildingBlocks = Vec<(SequenceElement<SemiAmbiguous>, Mass)>;
/// A list of all combinations of terminal modifications and their accompanying amino acid
pub type TerminalBuildingBlocks = Vec<(SequenceElement<SemiAmbiguous>, SimpleModification, Mass)>;
/// Get the possible building blocks for sequences based on the given modifications.
/// Useful for any automated sequence generation, like isobaric set generation or de novo sequencing.
/// The result is for each location (N term, center, C term) the list of all possible building blocks with its mass, sorted on mass.
/// # Panics
/// Panics if any of the modifications does not have a defined mass.
pub fn building_blocks(
    amino_acids: &[AminoAcid],
    fixed: &[(SimpleModification, Option<PlacementRule>)],
    variable: &[(SimpleModification, Option<PlacementRule>)],
) -> (
    TerminalBuildingBlocks,
    BuildingBlocks,
    TerminalBuildingBlocks,
) {
    /// Enforce the placement rules of predefined modifications.
    fn can_be_placed(
        modification: &SimpleModification,
        seq: &SequenceElement<SemiAmbiguous>,
        position: SequencePosition,
    ) -> bool {
        if let SimpleModificationInner::Database { specificities, .. } = &**modification {
            specificities.is_empty()
                || specificities
                    .iter()
                    .any(|(rules, _, _)| PlacementRule::any_possible(rules, seq, position))
        } else {
            true
        }
    }
    fn n_term_options(amino_acids: &[AminoAcid], rule: &PlacementRule) -> Vec<AminoAcid> {
        match rule {
            PlacementRule::AminoAcid(aa, Position::AnyNTerm | Position::ProteinNTerm) => {
                amino_acids
                    .iter()
                    .filter(|a| aa.contains(a))
                    .copied()
                    .collect_vec()
            }
            PlacementRule::Terminal(Position::AnyNTerm | Position::ProteinNTerm) => {
                amino_acids.iter().copied().collect_vec()
            }
            _ => Vec::new(),
        }
    }
    fn c_term_options(amino_acids: &[AminoAcid], rule: &PlacementRule) -> Vec<AminoAcid> {
        match rule {
            PlacementRule::AminoAcid(aa, Position::AnyCTerm | Position::ProteinCTerm) => {
                amino_acids
                    .iter()
                    .filter(|a| aa.contains(a))
                    .copied()
                    .collect_vec()
            }
            PlacementRule::Terminal(Position::AnyCTerm | Position::ProteinCTerm) => {
                amino_acids.iter().copied().collect_vec()
            }
            _ => Vec::new(),
        }
    }
    fn generate_terminal(
        position: &impl Fn(&PlacementRule) -> Vec<AminoAcid>,
        fixed: &[(SimpleModification, Option<PlacementRule>)],
        variable: &[(SimpleModification, Option<PlacementRule>)],
    ) -> TerminalBuildingBlocks {
        let mut options = fixed
            .iter()
            .chain(variable.iter())
            .flat_map(|(modification, rule)| {
                rule.as_ref().map_or_else(
                    || {
                        if let SimpleModificationInner::Database { specificities, .. } =
                            &**modification
                        {
                            specificities
                                .iter()
                                .flat_map(|(rules, _, _)| {
                                    rules.iter().flat_map(position).collect_vec()
                                })
                                .unique()
                                .map(|a| {
                                    (
                                        SequenceElement::new(CheckedAminoAcid::new(a), None),
                                        modification.clone(),
                                    )
                                })
                                .collect_vec()
                        } else {
                            Vec::new()
                        }
                    },
                    |rule| {
                        position(rule)
                            .into_iter()
                            .map(|a| {
                                (
                                    SequenceElement::new(CheckedAminoAcid::new(a), None),
                                    modification.clone(),
                                )
                            })
                            .collect_vec()
                    },
                )
            })
            .flat_map(|(a, m)| {
                #[expect(clippy::redundant_clone)] // not redundant
                let mc = m.clone();
                a.formulas_all(
                    &[],
                    &[],
                    &mut Vec::new(),
                    false,
                    SequencePosition::default(),
                    0,
                    &GlycanModel::DISALLOW,
                )
                .0
                .iter()
                .map(|f| f.monoisotopic_mass() + m.formula().monoisotopic_mass())
                .map(|mass| (a.clone(), mc.clone(), mass))
                .collect_vec()
            })
            .collect_vec();
        options.sort_unstable_by(|a, b| a.2.value.total_cmp(&b.2.value));
        options
    }

    let generate = |position| {
        let mut options: Vec<(SequenceElement<SemiAmbiguous>, Mass)> = amino_acids
            .iter()
            .flat_map(|aa| {
                let mut options = Vec::new();
                options.extend(
                    fixed
                        .iter()
                        .filter(|&m| {
                            m.1.as_ref().map_or_else(
                                || {
                                    can_be_placed(
                                        &m.0,
                                        &SequenceElement::new(aa.into(), None),
                                        position,
                                    )
                                },
                                |rule| {
                                    rule.is_possible(
                                        &SequenceElement::new(aa.into(), None),
                                        position,
                                    )
                                },
                            )
                        })
                        .map(|m| {
                            SequenceElement::new(aa.into(), None)
                                .with_simple_modification(m.0.clone())
                        }),
                );
                if options.is_empty() {
                    vec![SequenceElement::new(aa.into(), None)]
                } else {
                    options
                }
            })
            .flat_map(|seq| {
                let mut options = vec![seq.clone()];
                options.extend(
                    variable
                        .iter()
                        .filter(|&m| {
                            m.1.as_ref().map_or_else(
                                || can_be_placed(&m.0, &seq, position),
                                |rule| rule.is_possible(&seq, position),
                            )
                        })
                        .map(|m| {
                            let mut modifications = seq.modifications.clone();
                            modifications.push(Modification::Simple(m.0.clone()));
                            let mut result = SequenceElement::new(seq.aminoacid, None);
                            result.modifications = modifications;
                            result
                        }),
                );
                options
            })
            .flat_map(|s| {
                s.formulas_all(
                    &[],
                    &[],
                    &mut Vec::new(),
                    false,
                    position,
                    0,
                    &GlycanModel::DISALLOW,
                )
                .0
                .iter()
                .map(|f| (s.clone(), f.monoisotopic_mass()))
                .collect_vec()
            })
            .collect();
        options.sort_unstable_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        options
    };

    // Create the building blocks
    (
        generate_terminal(&|rule| n_term_options(amino_acids, rule), fixed, variable),
        generate(SequencePosition::Index(0)),
        generate_terminal(&|rule| c_term_options(amino_acids, rule), fixed, variable),
    )
}

/// Find the isobaric sets for the given mass with the given modifications and ppm error.
/// The modifications are placed on any location they are allowed based on the given placement
/// rules, so using any modifications which provide those is advised. If the provided [`Peptidoform`]
/// has multiple formulas, it uses the formula with the lowest monoisotopic mass.
/// # Panics
/// Panics if any of the modifications does not have a defined mass. Or if the weight of the
/// base selection is already in the tolerance of the given mass.
pub fn find_isobaric_sets(
    mass: Mass,
    tolerance: Tolerance<Mass>,
    amino_acids: &[AminoAcid],
    fixed: &[(SimpleModification, Option<PlacementRule>)],
    variable: &[(SimpleModification, Option<PlacementRule>)],
    base: Option<&Peptidoform<SimpleLinear>>,
) -> IsobaricSetIterator {
    let bounds = tolerance.bounds(mass);
    let base_mass = base
        .and_then(|b| {
            b.formulas()
                .mass_bounds()
                .into_option()
                .map(|(f, _)| f.monoisotopic_mass())
        })
        .unwrap_or_default();
    let bounds = (bounds.0 - base_mass, bounds.1 - base_mass);
    assert!(
        bounds.0.value > 0.0,
        "Cannot have a base selection that has a weight within the tolerance of the intended final mass for isobaric search."
    );
    let (n_term, center, c_term) = building_blocks(amino_acids, fixed, variable);

    IsobaricSetIterator::new(n_term, c_term, center, bounds, base)
}

/// Iteratively generate isobaric sets based on the given settings.
#[derive(Debug)]
pub struct IsobaricSetIterator {
    n_term: Vec<(SequenceElement<SemiAmbiguous>, SimpleModification, Mass)>,
    c_term: Vec<(SequenceElement<SemiAmbiguous>, SimpleModification, Mass)>,
    center: Vec<(SequenceElement<SemiAmbiguous>, Mass)>,
    sizes: (Mass, Mass),
    bounds: (Mass, Mass),
    state: (Option<usize>, Option<usize>, Vec<usize>),
    base: Option<Peptidoform<SimpleLinear>>,
}

impl IsobaricSetIterator {
    /// `n_term` & `c_term` are the possible combinations of terminal modifications with their valid placements and the full mass of this combo
    /// # Panics
    /// If there is not at least one element in the `center` list.
    fn new(
        n_term: Vec<(SequenceElement<SemiAmbiguous>, SimpleModification, Mass)>,
        c_term: Vec<(SequenceElement<SemiAmbiguous>, SimpleModification, Mass)>,
        center: Vec<(SequenceElement<SemiAmbiguous>, Mass)>,
        bounds: (Mass, Mass),
        base: Option<&Peptidoform<SimpleLinear>>,
    ) -> Self {
        let sizes = (center.first().unwrap().1, center.last().unwrap().1);
        let mut iter = Self {
            n_term,
            c_term,
            center,
            sizes,
            bounds,
            state: (None, None, Vec::new()),
            base: base.cloned(),
        };
        while iter.current_mass() < iter.bounds.0 - iter.sizes.0 {
            iter.state.2.push(0);
        }
        iter
    }

    fn current_mass(&self) -> Mass {
        self.state.0.map(|i| self.n_term[i].2).unwrap_or_default()
            + self.state.1.map(|i| self.c_term[i].2).unwrap_or_default()
            + self
                .state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].1)
                .sum::<Mass>()
    }

    fn mass_fits(&self) -> Ordering {
        let mass = self.current_mass();
        if mass < self.bounds.0 {
            Ordering::Less
        } else if mass > self.bounds.1 {
            Ordering::Greater
        } else {
            Ordering::Equal
        }
    }

    /// # Panics
    /// If the base sequence is empty.
    fn peptide(&self) -> Peptidoform<SimpleLinear> {
        let mut sequence = Vec::with_capacity(
            self.base.as_ref().map(Peptidoform::len).unwrap_or_default()
                + self.state.2.len()
                + usize::from(self.state.0.is_some())
                + usize::from(self.state.1.is_some()),
        );
        if self
            .base
            .as_ref()
            .is_some_and(|b| !b.get_n_term().is_empty())
        {
            sequence.push(self.base.as_ref().unwrap().sequence()[0].clone());
        } else if let Some(n) = self.state.0.map(|i| self.n_term[i].clone()) {
            sequence.push(n.0.into());
        }
        if let Some(base) = &self.base {
            let n = usize::from(base.get_n_term().is_empty());
            let c = usize::from(base.get_c_term().is_empty());
            sequence.extend(base.sequence()[n..base.len() - n - c].iter().cloned());
        }
        sequence.extend(
            self.state
                .2
                .iter()
                .copied()
                .map(|i| self.center[i].0.clone().into()),
        );
        if self
            .base
            .as_ref()
            .is_some_and(|b| !b.get_c_term().is_empty())
        {
            sequence.push(
                self.base
                    .as_ref()
                    .unwrap()
                    .sequence()
                    .last()
                    .unwrap()
                    .clone(),
            );
        } else if let Some(c) = self.state.1.map(|i| self.c_term[i].clone()) {
            sequence.push(c.0.into());
        }
        Peptidoform::new(sequence)
            .n_term(self.base.as_ref().map_or_else(
                || {
                    self.state.0.map_or(Vec::new(), |i| {
                        vec![Modification::Simple(self.n_term[i].1.clone())]
                    })
                },
                |b| b.get_n_term().to_vec(),
            ))
            .c_term(self.base.as_ref().map_or_else(
                || {
                    self.state.0.map_or(Vec::new(), |i| {
                        vec![Modification::Simple(self.c_term[i].1.clone())]
                    })
                },
                |b| b.get_c_term().to_vec(),
            ))
    }

    /// Reset the state for the center selection
    fn reset_center_state(&mut self) {
        self.state.2.clear();
        while self.current_mass() < self.bounds.0 - self.sizes.0 {
            self.state.2.push(0);
        }
    }
}

impl Iterator for IsobaricSetIterator {
    type Item = Peptidoform<SimpleLinear>;
    fn next(&mut self) -> Option<Self::Item> {
        loop {
            // N terminal loop
            loop {
                // C terminal loop

                // Main loop
                while !self.state.2.is_empty() {
                    // Do state + 1 at the highest level where this is still possible and check if that one fits the bounds
                    // Until every level is full then pop and try with one fewer number of amino acids
                    while !self.state.2.iter().all(|s| *s == self.center.len() - 1) {
                        let mut level = self.state.2.len() - 1;
                        loop {
                            if self.state.2[level] == self.center.len() - 1 {
                                if level == 0 {
                                    break;
                                }
                                level -= 1;
                            } else {
                                // Update this level
                                self.state.2[level] += 1;
                                // Reset the levels above, has to start at minimal at the index of this level to prevent 'rotations' of the set to show up
                                for l in level + 1..self.state.2.len() {
                                    self.state.2[l] = self.state.2[level];
                                }
                                match self.mass_fits() {
                                    Ordering::Greater => {
                                        // If the mass is too great the level below will have the be changed, otherwise it could only be getting heavier with the next iteration(s)
                                        if level == 0 {
                                            break;
                                        }
                                        level -= 1;
                                    }
                                    Ordering::Equal => {
                                        return Some(self.peptide());
                                    }
                                    Ordering::Less => {
                                        // If there a way to reach at least the lower limit by having all the heaviest options selected try and reach them.
                                        // Otherwise this will increase this level again next iteration.
                                        if self.state.2[0..level]
                                            .iter()
                                            .map(|i| self.center[*i].1)
                                            .sum::<Mass>()
                                            + Ratio::new::<fraction>(
                                                (self.state.2.len() - level) as f64,
                                            ) * self.sizes.1
                                            > self.bounds.0
                                        {
                                            level = self.state.2.len() - 1;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    self.state.2.pop();
                    // Stop the search when there is no possibility for a fitting answer
                    if self.sizes.1 * Ratio::new::<fraction>(self.state.2.len() as f64)
                        < self.bounds.0
                    {
                        break;
                    }
                    // Reset the levels to be all 0s again
                    for level in 0..self.state.2.len() {
                        self.state.2[level] = 0;
                    }
                }

                // Try the next C terminal option
                if let Some(c) = self.state.1 {
                    if c + 1 >= self.c_term.len() {
                        break;
                    }
                    self.state.1 = Some(c + 1);
                } else if self.c_term.is_empty()
                    || self
                        .base
                        .as_ref()
                        .is_some_and(|b| !b.get_c_term().is_empty())
                {
                    // If the base piece has a defined C term mod do not try possible C term mods in the isobaric generation
                    break;
                } else {
                    self.state.1 = Some(0);
                }
                self.reset_center_state();
            }
            // Try the next N terminal option
            if let Some(n) = self.state.0 {
                if n + 1 >= self.n_term.len() {
                    break;
                }
                self.state.0 = Some(n + 1);
            } else if self.n_term.is_empty()
                || self
                    .base
                    .as_ref()
                    .is_some_and(|b| !b.get_n_term().is_empty())
            {
                // If the base piece has a defined N term mod do not try possible N term mods in the isobaric generation
                break;
            } else {
                self.state.1 = Some(0);
            }
            self.reset_center_state();
        }
        None
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {

    use super::*;
    #[test]
    fn simple_isobaric_sets() {
        let pep = Peptidoform::pro_forma("AG", None)
            .unwrap()
            .into_unambiguous()
            .unwrap();
        let sets: Vec<Peptidoform<SimpleLinear>> = find_isobaric_sets(
            pep.bare_formula().monoisotopic_mass(),
            Tolerance::new_ppm(10.0),
            AminoAcid::UNIQUE_MASS_AMINO_ACIDS,
            &[],
            &[],
            None,
        )
        .collect();
        assert_eq!(
            &sets,
            &[
                Peptidoform::pro_forma("GA", None)
                    .unwrap()
                    .into_simple_linear()
                    .unwrap(),
                Peptidoform::pro_forma("Q", None)
                    .unwrap()
                    .into_simple_linear()
                    .unwrap(),
            ]
        );
    }
}
