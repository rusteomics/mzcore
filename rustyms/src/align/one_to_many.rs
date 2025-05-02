use std::sync::Arc;

use crate::{
    align::diagonal_array::DiagonalArray,
    prelude::{MassMode, Peptidoform},
    quantities::Multi,
    sequence::{AtMax, SimpleLinear, UnAmbiguous},
    system::Mass,
};

use super::{
    AlignScoring, AlignType, Alignment,
    mass_alignment::{align_cached, calculate_masses},
};

#[derive(Debug)]
pub struct OneToManyIndex<'lifetime, const STEPS: u16, Complexity> {
    sequences: Vec<(Arc<Peptidoform<Complexity>>, DiagonalArray<Multi<Mass>>)>,
    mode: MassMode,
    scoring: AlignScoring<'lifetime>,
    align_type: AlignType,
}

impl<'lifetime, const STEPS: u16, Complexity: AtMax<SimpleLinear>>
    OneToManyIndex<'lifetime, STEPS, Complexity>
{
    pub fn new<Source: HasPeptidoform<Complexity>>(
        sequences: impl IntoIterator<Item = Source>,
        mode: MassMode,
        scoring: AlignScoring<'lifetime>,
        align_type: AlignType,
    ) -> Self {
        Self {
            sequences: sequences
                .into_iter()
                .filter_map(|p| p.peptidoform())
                .map(|p| {
                    let masses = calculate_masses::<STEPS>(&p, mode);
                    (p, masses)
                })
                .collect(),
            mode,
            scoring,
            align_type,
        }
    }

    pub fn align_one<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        other: &'alignment Peptidoform<OtherComplexity>,
    ) -> Vec<Alignment<'alignment, Complexity, OtherComplexity>> {
        let other_masses = calculate_masses::<STEPS>(other, self.mode);

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

    // TODO: parallel, Arc<pep> input, Arc<pep> alignment, any kind of input type that could give a pep
    pub fn align<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        others: &'alignment [Peptidoform<OtherComplexity>],
    ) -> Vec<Vec<Alignment<'alignment, Complexity, OtherComplexity>>> {
        others.iter().map(|o| self.align_one(o)).collect()
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
