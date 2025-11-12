use mzcv::CVIndex;
#[cfg(feature = "rayon")]
use rayon::prelude::*;
use std::collections::HashSet;

use mzcore::sequence::{
    AnnotatedPeptide, Annotation, HasPeptidoformImpl, Peptidoform, Region, UnAmbiguous,
};

pub(super) use super::*;

/// The selection rules for iterating over a selection of germlines.
#[derive(Clone, Debug, Eq, PartialEq)]
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
    'a,
    S1: std::hash::BuildHasher + Clone + Send + Sync + 'a,
    S2: std::hash::BuildHasher + Clone + Send + Sync + 'a,
> Selection<S1, S2>
{
    /// Get the selected alleles
    pub fn germlines(self, cv: &'a CVIndex<IMGT>) -> impl Iterator<Item = Allele<'a>> {
        cv.data()
            .iter()
            .filter(move |(s, _)| {
                self.species
                    .as_ref()
                    .is_none_or(|species| species.contains(s))
            })
            .flat_map(|(s, g)| g.into_iter().map(|c| (*s, c.0, c.1)))
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

    /// Get the selected alleles in parallel fashion, only available if you enable the feature "rayon" (on by default)
    #[cfg(feature = "rayon")]
    pub fn par_germlines(self, cv: &'a CVIndex<IMGT>) -> impl ParallelIterator<Item = Allele<'a>> {
        cv.data()
            .par_iter()
            .filter(move |(s, _)| {
                self.species
                    .as_ref()
                    .is_none_or(|species| species.contains(s))
            })
            .flat_map(|(s, g)| g.into_par_iter().map(|c| (*s, c.0, c.1)))
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
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq)]
pub enum AlleleSelection {
    /// Return all alleles
    All,
    /// Only return the first allele. It can have a number higher than 1 if the previous alleles are not functional.
    First,
}

impl AlleleSelection {
    #[allow(dead_code)] // If internal no data is turned on
    const fn take_num(self) -> usize {
        match self {
            Self::First => 1,
            Self::All => usize::MAX,
        }
    }
}

/// A returned allele
#[non_exhaustive] // Do not let anyone build it themselves
#[derive(Clone, Debug, Eq, PartialEq)]
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

impl HasPeptidoformImpl for Allele<'_> {
    type Complexity = UnAmbiguous;
    fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
        self.sequence
    }
}

impl AnnotatedPeptide for Allele<'_> {
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
    /// Find a specific allele.
    pub fn find_allele(&self, gene: Gene, allele: Option<usize>) -> Option<Allele<'_>> {
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
                allele
                    .map_or_else(
                        || g.alleles.first(),
                        |a| g.alleles.iter().find(|(ga, _)| a == *ga),
                    )
                    .map(|res| (g, res))
            })
            .map(move |(g, (a, seq))| Allele {
                species: g.species,
                gene: std::borrow::Cow::Owned(gene),
                number: *a,
                sequence: &seq.sequence,
                regions: &seq.regions,
                annotations: &seq.annotations,
            })
    }

    /// Find a specific germline.
    pub fn find_germline(&self, gene: Gene) -> Option<std::sync::Arc<Germline>> {
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
            .map(|g| genes[g].clone())
    }

    /// Remove a specific germline.
    pub fn remove_germline(&mut self, gene: Gene) {
        let chain = match gene.chain {
            ChainType::Heavy => &mut self.h,
            ChainType::LightKappa => &mut self.k,
            ChainType::LightLambda => &mut self.l,
            ChainType::Iota => &mut self.i,
        };
        let genes = match gene.kind {
            GeneType::V => &mut chain.variable,
            GeneType::J => &mut chain.joining,
            GeneType::C(None) => &mut chain.c,
            GeneType::C(Some(Constant::A)) => &mut chain.a,
            GeneType::C(Some(Constant::D)) => &mut chain.d,
            GeneType::C(Some(Constant::E)) => &mut chain.e,
            GeneType::C(Some(Constant::G)) => &mut chain.g,
            GeneType::C(Some(Constant::M)) => &mut chain.m,
            GeneType::C(Some(Constant::O)) => &mut chain.o,
            GeneType::C(Some(Constant::T)) => &mut chain.t,
        };
        genes
            .binary_search_by(|g| g.name.cmp(&gene))
            .ok()
            .map(|g| genes.remove(g));
    }
}

#[cfg(all(test, not(feature = "internal-no-data")))]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use std::collections::HashSet;

    use crate::{ChainType, GeneType, STATIC_IMGT, Selection, Species, select::contains_gene};

    #[test]
    fn try_first_human() {
        let selection = Selection::default()
            .species([Species::HomoSapiens])
            .chain([ChainType::Heavy])
            .gene([GeneType::V]);
        let first = selection.germlines(&STATIC_IMGT).next().unwrap();
        assert_eq!(first.name(), "IGHV1-2*01");
    }

    #[test]
    fn try_first_g_human() {
        let selection = Selection::default()
            .species([Species::HomoSapiens])
            .chain([ChainType::Heavy])
            .gene([GeneType::C(Some(crate::Constant::G))]);
        let first = selection.germlines(&STATIC_IMGT).next().unwrap();
        assert_eq!(first.name(), "IGHGP*01");
    }

    #[test]
    fn gene_selections() {
        let constant = HashSet::from([GeneType::C(None)]);
        assert!(contains_gene(&constant, GeneType::C(None)));
        assert!(contains_gene(
            &constant,
            GeneType::C(Some(crate::Constant::G))
        ));
        assert!(contains_gene(
            &constant,
            GeneType::C(Some(crate::Constant::A))
        ));
        let constant_g = HashSet::from([GeneType::C(Some(crate::Constant::G))]);
        assert!(!contains_gene(&constant_g, GeneType::C(None)));
        assert!(contains_gene(
            &constant_g,
            GeneType::C(Some(crate::Constant::G))
        ));
        assert!(!contains_gene(
            &constant_g,
            GeneType::C(Some(crate::Constant::A))
        ));
    }
}
