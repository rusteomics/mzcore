#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::collections::HashSet;

use crate::peptidoform::{AnnotatedPeptide, Annotation, Region, UnAmbiguous};
use crate::Peptidoform;

pub use super::fancy::FancyDisplay;
pub use super::shared::*;

/// Get a specific germline
pub fn get_germline(
    species: Species,
    gene: Gene,
    allele: Option<usize>,
) -> Option<Allele<'static>> {
    super::germlines(species).and_then(|g| g.find(species, gene, allele))
}

/// The selection rules for iterating over a selection of germlines.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Selection<S1: std::hash::BuildHasher, S2: std::hash::BuildHasher> {
    /// The species you want, None allows all, otherwise only the species specified will be returned
    pub species: Option<HashSet<Species, S1>>,
    /// The chain of genes you want, None allows all, otherwise only the chains specified will be returned
    pub chains: Option<HashSet<ChainType, S2>>,
    /// The kind of genes you want, None allows all, otherwise only the genes specified will be returned
    pub genes: Option<HashSet<GeneType>>,
    /// The way of handling alleles you want
    pub allele: AlleleSelection,
}

impl<S1: std::hash::BuildHasher, S2: std::hash::BuildHasher> Selection<S1, S2> {
    /// Builder pattern method to add a species selection, will replace any previously set species selection
    #[must_use]
    pub fn species(self, species: impl Into<HashSet<Species, S1>>) -> Self {
        Self {
            species: Some(species.into()),
            ..self
        }
    }

    /// Builder pattern method to add a chain selection, will replace any previously set chain selection
    #[must_use]
    pub fn chain(self, chains: impl Into<HashSet<ChainType, S2>>) -> Self {
        Self {
            chains: Some(chains.into()),
            ..self
        }
    }

    /// Builder pattern method to add a gene selection, will replace any previously set gene selection
    #[must_use]
    pub fn gene(self, genes: impl Into<HashSet<GeneType>>) -> Self {
        Self {
            genes: Some(genes.into()),
            ..self
        }
    }

    /// Builder pattern method to add an allele selection, will replace any previously set allele selection
    #[must_use]
    pub fn allele(self, allele: AlleleSelection) -> Self {
        Self { allele, ..self }
    }
}

impl<
        S1: std::hash::BuildHasher + Clone + Send + Sync,
        S2: std::hash::BuildHasher + Clone + Send + Sync,
    > Selection<S1, S2>
{
    /// Get the selected alleles
    pub fn germlines(self) -> impl Iterator<Item = Allele<'static>> {
        super::all_germlines()
            .filter(move |g| self.species.as_ref().is_none_or(|s| s.contains(&g.species)))
            .flat_map(|g| g.into_iter().map(|c| (g.species, c.0, c.1)))
            .filter(move |(_, kind, _)| self.chains.as_ref().is_none_or(|k| k.contains(kind)))
            .flat_map(|(species, _, c)| c.into_iter().map(move |g| (species, g.0, g.1)))
            .filter(move |(_, gene, _)| self.genes.as_ref().is_none_or(|s| contains_gene(s, *gene)))
            .flat_map(|(species, _, germlines)| germlines.iter().map(move |a| (species, a)))
            .flat_map(move |(species, germline)| {
                germline
                    .into_iter()
                    .take(self.allele.take_num())
                    .map(move |(a, seq)| (species, &germline.name, *a, seq))
            })
            .map(Into::into)
    }

    #[cfg(feature = "rayon")]
    /// Get the selected alleles in parallel fashion, only available if you enable the feature "rayon" (on by default)
    pub fn par_germlines(self) -> impl ParallelIterator<Item = Allele<'static>> {
        super::par_germlines()
            .filter(move |g| self.species.as_ref().is_none_or(|s| s.contains(&g.species)))
            .flat_map(|g| g.into_par_iter().map(|c| (g.species, c.0, c.1)))
            .filter(move |(_, kind, _)| self.chains.as_ref().is_none_or(|k| k.contains(kind)))
            .flat_map(|(species, _, c)| c.into_par_iter().map(move |g| (species, g.0, g.1)))
            .filter(move |(_, gene, _)| self.genes.as_ref().is_none_or(|s| contains_gene(s, *gene)))
            .flat_map(|(species, _, germlines)| {
                germlines.into_par_iter().map(move |a| (species, a))
            })
            .flat_map(move |(species, germline)| {
                germline
                    .into_par_iter()
                    .take(self.allele.take_num())
                    .map(move |(a, seq)| (species, &germline.name, *a, seq))
            })
            .map(Into::into)
    }
}

fn contains_gene(s: &HashSet<GeneType>, gene: GeneType) -> bool {
    s.contains(&gene) || matches!(gene, GeneType::C(_)) && s.contains(&GeneType::C(None))
}

impl<S1: std::hash::BuildHasher, S2: std::hash::BuildHasher> Default for Selection<S1, S2> {
    /// Get a default selection, which gives all kinds and genes but only returns the first allele
    fn default() -> Self {
        Self {
            species: None,
            chains: None,
            genes: None,
            allele: AlleleSelection::First,
        }
    }
}

/// The allele handling strategy
#[derive(Debug, Clone, Copy, Hash, PartialEq, Eq)]
pub enum AlleleSelection {
    /// Return all alleles
    All,
    /// Only return the first allele. It can have a number higher than 1 if the previous alleles are not functional.
    First,
}

impl AlleleSelection {
    const fn take_num(self) -> usize {
        match self {
            Self::First => 1,
            Self::All => usize::MAX,
        }
    }
}

/// A returned allele
#[non_exhaustive] // Do not let anyone build it themselves
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Allele<'a> {
    /// The species where this gene originates from
    pub species: Species,
    /// The gene where this is the sequence for, eg `IGHV3-23`
    pub gene: std::borrow::Cow<'a, Gene>,
    /// The allele number, in IMGT this follows the name, eg `*01` is the allele in `IGHV3-23*01`
    pub number: usize,
    /// The actual sequence, the sequences present in the database are pure amino acids, no modifications are to be expected
    pub sequence: &'a Peptidoform<UnAmbiguous>,
    /// The regions in the sequence, every region has an annotation and a length, all lengths together are the same length as the full sequence
    pub regions: &'a [(Region, usize)],
    /// Any additional annotations, every annotation has beside the kind it is also its location, as index in the sequence
    pub annotations: &'a [(Annotation, usize)],
}

impl Allele<'_> {
    /// Get the IMGT name for this allele
    pub fn name(&self) -> String {
        format!("{}*{:02}", self.gene, self.number)
    }

    /// Get the biologists name for this allele with fancy non ASCII characters
    pub fn fancy_name(&self) -> String {
        format!("{}*{:02}", self.gene.to_fancy_string(), self.number)
    }
}

impl AnnotatedPeptide for Allele<'_> {
    type Complexity = UnAmbiguous;
    fn peptide(&self) -> &Peptidoform<Self::Complexity> {
        self.sequence
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        self.annotations
    }
    fn regions(&self) -> &[(Region, usize)] {
        self.regions
    }
}

impl<'a> From<(Species, &'a Gene, usize, &'a AnnotatedSequence)> for Allele<'a> {
    fn from(value: (Species, &'a Gene, usize, &'a AnnotatedSequence)) -> Self {
        Self {
            species: value.0,
            gene: std::borrow::Cow::Borrowed(value.1),
            number: value.2,
            sequence: &value.3.sequence,
            regions: &value.3.regions,
            annotations: &value.3.annotations,
        }
    }
}

impl Germlines {
    /// Find a specific germline.
    pub fn find(&self, species: Species, gene: Gene, allele: Option<usize>) -> Option<Allele<'_>> {
        let chain = match gene.chain {
            ChainType::Heavy => &self.h,
            ChainType::LightKappa => &self.k,
            ChainType::LightLambda => &self.l,
            ChainType::Iota => &self.i,
        };
        let genes = match gene.kind {
            GeneType::V => &chain.variable,
            GeneType::J => &chain.joining,
            GeneType::C(None) => &chain.c,
            GeneType::C(Some(Constant::A)) => &chain.a,
            GeneType::C(Some(Constant::D)) => &chain.d,
            GeneType::C(Some(Constant::E)) => &chain.e,
            GeneType::C(Some(Constant::G)) => &chain.g,
            GeneType::C(Some(Constant::M)) => &chain.m,
            GeneType::C(Some(Constant::O)) => &chain.o,
            GeneType::C(Some(Constant::T)) => &chain.t,
        };
        genes
            .binary_search_by(|g| g.name.cmp(&gene))
            .ok()
            .and_then(|g| {
                let g = &genes[g];
                allele.map_or(g.alleles.first(), |a| {
                    g.alleles.iter().find(|(ga, _)| a == *ga)
                })
            })
            .map(move |(a, seq)| Allele {
                species,
                gene: std::borrow::Cow::Owned(gene),
                number: *a,
                sequence: &seq.sequence,
                regions: &seq.regions,
                annotations: &seq.annotations,
            })
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use std::collections::HashSet;

    use crate::imgt::select::contains_gene;

    use super::Selection;
    use super::{ChainType, GeneType, Species};

    #[test]
    fn try_first_human() {
        let selection = Selection::default()
            .species([Species::HomoSapiens])
            .chain([ChainType::Heavy])
            .gene([GeneType::V]);
        let first = selection.germlines().next().unwrap();
        assert_eq!(first.name(), "IGHV1-2*01");
    }

    #[test]
    fn try_first_g_human() {
        let selection = Selection::default()
            .species([Species::HomoSapiens])
            .chain([ChainType::Heavy])
            .gene([GeneType::C(Some(crate::imgt::Constant::G))]);
        let first = selection.germlines().next().unwrap();
        assert_eq!(first.name(), "IGHGP*01");
    }

    #[test]
    fn gene_selections() {
        let constant = HashSet::from([GeneType::C(None)]);
        assert!(contains_gene(&constant, GeneType::C(None)));
        assert!(contains_gene(
            &constant,
            GeneType::C(Some(crate::imgt::Constant::G))
        ));
        assert!(contains_gene(
            &constant,
            GeneType::C(Some(crate::imgt::Constant::A))
        ));
        let constant_g = HashSet::from([GeneType::C(Some(crate::imgt::Constant::G))]);
        assert!(!contains_gene(&constant_g, GeneType::C(None)));
        assert!(contains_gene(
            &constant_g,
            GeneType::C(Some(crate::imgt::Constant::G))
        ));
        assert!(!contains_gene(
            &constant_g,
            GeneType::C(Some(crate::imgt::Constant::A))
        ));
    }
}
