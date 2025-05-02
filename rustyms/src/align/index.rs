use std::{
    collections::{HashMap, HashSet},
    ops::RangeInclusive,
    sync::Arc,
};

use itertools::Itertools;

use crate::{
    align::diagonal_array::DiagonalArray,
    prelude::{MassMode, Peptidoform},
    quantities::{Multi, Tolerance},
    sequence::{AtMax, SimpleLinear, UnAmbiguous},
    system::{Mass, OrderedMass},
};

use super::{
    AlignScoring, AlignType, Alignment,
    mass_alignment::{align_cached, calculate_masses},
};

#[derive(Debug)]
pub struct OneToManyIndex<'lifetime, const STEPS: u16, Complexity> {
    sequences: Vec<(Arc<Peptidoform<Complexity>>, DiagonalArray<Multi<Mass>>)>,
    mass_index: HashMap<OrderedMass, HashSet<usize>>,
    mode: MassMode,
    scoring: AlignScoring<'lifetime>,
    align_type: AlignType,
    matches: Option<(u8, RangeInclusive<Mass>)>,
}

impl<'lifetime, const STEPS: u16, Complexity: AtMax<SimpleLinear>>
    OneToManyIndex<'lifetime, STEPS, Complexity>
{
    pub fn new<Source: HasPeptidoform<Complexity>>(
        sequences: impl IntoIterator<Item = Source>,
        mode: MassMode,
        scoring: AlignScoring<'lifetime>,
        align_type: AlignType,
        matches: Option<(u8, RangeInclusive<Mass>)>,
    ) -> Self {
        let sequences = sequences
            .into_iter()
            .filter_map(|p| p.peptidoform())
            .map(|p| {
                let masses = calculate_masses::<STEPS>(&p, mode);
                (p, masses)
            })
            .collect::<Vec<_>>();
        let mass_index = if let Some((_matches, range)) = matches.clone() {
            let mut map: HashMap<OrderedMass, HashSet<usize>> = HashMap::new();
            let tolerance = scoring.tolerance;
            let start = *range.start();
            for (index, mass) in sequences
                .iter()
                .enumerate()
                .flat_map(|(index, (_, masses))| {
                    masses
                        .iter()
                        .flat_map(|masses| masses.iter())
                        .filter(|mass| range.contains(mass))
                        .map(move |mass| (index, standardise_mass(start, *mass, tolerance)))
                })
            {
                map.entry(mass).or_default().insert(index);
            }
            map
        } else {
            HashMap::new()
        };

        Self {
            sequences,
            mass_index,
            mode,
            scoring,
            align_type,
            matches,
        }
    }

    pub fn align_one<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        other: &'alignment Peptidoform<OtherComplexity>,
    ) -> Vec<Alignment<'alignment, Complexity, OtherComplexity>> {
        let other_masses = calculate_masses::<STEPS>(other, self.mode);
        if let Some((matches, range)) = self.matches.clone() {
            let mut interesting: HashMap<usize, u8> = HashMap::new();

            for mass in other_masses
                .iter()
                .flat_map(|masses| masses.iter())
                .filter(|mass| range.contains(mass))
                .map(|mass| standardise_mass(*range.start(), *mass, self.scoring.tolerance))
            {
                // Do I need to check multiple bins?
                if let Some(indices) = self.mass_index.get(&mass) {
                    for index in indices {
                        *interesting.entry(*index).or_default() += 1;
                    }
                }
            }

            interesting
                .iter()
                .filter_map(|(index, times)| (*times >= matches).then_some(*index))
                .map(|index| &self.sequences[index])
                .map(|(template, template_masses)| {
                    align_cached::<STEPS, Complexity, OtherComplexity>(
                        template,
                        template_masses,
                        other,
                        &other_masses,
                        self.scoring,
                        self.align_type,
                    )
                })
                .collect()
        } else {
            self.sequences
                .iter()
                .map(|(template, template_masses)| {
                    align_cached::<STEPS, Complexity, OtherComplexity>(
                        template,
                        template_masses,
                        other,
                        &other_masses,
                        self.scoring,
                        self.align_type,
                    )
                })
                .collect()
        }
    }

    // TODO: parallel, Arc<pep> input, Arc<pep> alignment, any kind of input type that could give a pep
    pub fn align<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        others: &'alignment [Peptidoform<OtherComplexity>],
    ) -> Vec<Vec<Alignment<'alignment, Complexity, OtherComplexity>>> {
        others.iter().map(|o| self.align_one(o)).collect()
    }
}

fn standardise_mass(start: Mass, mass: Mass, tolerance: Tolerance<OrderedMass>) -> OrderedMass {
    match tolerance {
        Tolerance::Absolute(da) => Mass::new::<crate::system::dalton>(
            ((mass - start) / *da)
                .value
                .floor()
                .mul_add(da.value, start.value),
        )
        .into(),
        Tolerance::Relative(ppm) => {
            todo!()
        }
    }
}

pub trait HasPeptidoform<Complexity> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<Complexity>>>;
}

impl HasPeptidoform<UnAmbiguous> for crate::imgt::Allele<'_> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<UnAmbiguous>>> {
        Some(Arc::new(self.sequence.clone()))
    }
}

impl HasPeptidoform<UnAmbiguous> for &crate::imgt::Allele<'_> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<UnAmbiguous>>> {
        Some(Arc::new(self.sequence.clone()))
    }
}
