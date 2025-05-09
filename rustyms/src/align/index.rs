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
    sequence::{AtMax, SemiAmbiguous, SimpleLinear, UnAmbiguous},
    system::{Mass, OrderedMass},
};

use super::{
    AlignScoring, AlignType, Alignment,
    mass_alignment::{align_cached, calculate_masses},
};

#[derive(Debug)]
pub struct OneToManyIndex<'lifetime, const STEPS: u16, Complexity> {
    sequences: Vec<(Arc<Peptidoform<Complexity>>, DiagonalArray<Multi<Mass>>)>,
    mass_index: Vec<(Mass, Vec<usize>)>,
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
            let mut map: Vec<(Mass, Vec<usize>)> = Vec::new();
            for (index, mass) in sequences
                .iter()
                .enumerate()
                .flat_map(|(index, (_, masses))| {
                    masses
                        .iter()
                        .flat_map(|masses| masses.iter())
                        .filter(|mass| range.contains(mass))
                        .map(move |mass| (index, *mass))
                })
            {
                match map.binary_search_by(|value| value.0.value.total_cmp(&mass.value)) {
                    Ok(i) => match map[i].1.binary_search(&index) {
                        Ok(_) => (),
                        Err(i2) => map[i].1.insert(i2, index),
                    },
                    Err(i) => map.insert(i, (mass, vec![index])),
                }
            }
            map
        } else {
            Vec::new()
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
    ) -> Vec<(usize, Alignment<'alignment, Complexity, OtherComplexity>)> {
        let other_masses = calculate_masses::<STEPS>(other, self.mode);
        if let Some((matches, range)) = self.matches.clone() {
            let mut interesting: HashMap<usize, u8> = HashMap::new();

            for mass in other_masses
                .iter()
                .flat_map(|masses| masses.iter())
                .filter(|mass| range.contains(mass))
            {
                let (low, high) = self.scoring.tolerance.bounds(*mass);
                let first_index = self
                    .mass_index
                    .binary_search_by(|value| value.0.value.total_cmp(&low.value))
                    .unwrap_or_else(|i| (i + 1))
                    .min(self.mass_index.len() - 1);
                let last_index = self
                    .mass_index
                    .binary_search_by(|value| value.0.value.total_cmp(&high.value))
                    .map_or_else(|i| i, |i| i + 1)
                    .min(self.mass_index.len() - 1)
                    .max(first_index);

                for index in self.mass_index[first_index..last_index]
                    .iter()
                    .flat_map(|(_, indices)| indices)
                {
                    let value = interesting.entry(*index).or_default();
                    *value = value.saturating_add(1);
                }
            }

            interesting
                .iter()
                .filter_map(|(index, times)| (*times >= matches).then_some(*index))
                .map(|index| (index, &self.sequences[index]))
                .map(|(index, (template, template_masses))| {
                    (
                        index,
                        align_cached::<STEPS, Complexity, OtherComplexity>(
                            template,
                            template_masses,
                            other,
                            &other_masses,
                            self.scoring,
                            self.align_type,
                        ),
                    )
                })
                .collect()
        } else {
            self.sequences
                .iter()
                .enumerate()
                .map(|(index, (template, template_masses))| {
                    (
                        index,
                        align_cached::<STEPS, Complexity, OtherComplexity>(
                            template,
                            template_masses,
                            other,
                            &other_masses,
                            self.scoring,
                            self.align_type,
                        ),
                    )
                })
                .collect()
        }
    }

    pub fn debug_align_one<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        other: &'alignment Peptidoform<OtherComplexity>,
    ) -> Vec<(
        usize,
        Alignment<'alignment, Complexity, OtherComplexity>,
        u8,
    )> {
        let other_masses = calculate_masses::<STEPS>(other, self.mode);
        if let Some((_matches, range)) = self.matches.clone() {
            let mut interesting: HashMap<usize, u8> = HashMap::new();

            for mass in other_masses
                .iter()
                .flat_map(|masses| masses.iter())
                .filter(|mass| range.contains(mass))
            {
                let (low, high) = self.scoring.tolerance.bounds(*mass);
                let first_index = self
                    .mass_index
                    .binary_search_by(|value| value.0.value.total_cmp(&low.value))
                    .unwrap_or_else(|i| (i + 1))
                    .min(self.mass_index.len() - 1);
                let last_index = self
                    .mass_index
                    .binary_search_by(|value| value.0.value.total_cmp(&high.value))
                    .unwrap_or_else(|i| i.saturating_sub(1))
                    .min(self.mass_index.len() - 1)
                    .max(first_index);

                for index in self.mass_index[first_index..=last_index]
                    .iter()
                    .flat_map(|(_, indices)| indices)
                {
                    let value = interesting.entry(*index).or_default();
                    *value = value.saturating_add(1);
                }
            }

            self.sequences
                .iter()
                .enumerate()
                .map(|(index, (template, template_masses))| {
                    (
                        index,
                        align_cached::<STEPS, Complexity, OtherComplexity>(
                            template,
                            template_masses,
                            other,
                            &other_masses,
                            self.scoring,
                            self.align_type,
                        ),
                    )
                })
                .map(|(i, a)| (i, a, interesting.get(&i).map_or(0, |v| *v)))
                .collect()
        } else {
            self.sequences
                .iter()
                .enumerate()
                .map(|(index, (template, template_masses))| {
                    (
                        index,
                        align_cached::<STEPS, Complexity, OtherComplexity>(
                            template,
                            template_masses,
                            other,
                            &other_masses,
                            self.scoring,
                            self.align_type,
                        ),
                        0,
                    )
                })
                .collect()
        }
    }

    pub fn debug_align<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        others: &'alignment [Peptidoform<OtherComplexity>],
    ) -> Vec<
        Vec<(
            usize,
            Alignment<'alignment, Complexity, OtherComplexity>,
            u8,
        )>,
    > {
        others.iter().map(|o| self.debug_align_one(o)).collect()
    }

    // TODO: parallel, Arc<pep> input, Arc<pep> alignment, any kind of input type that could give a pep
    pub fn align<'alignment, OtherComplexity: AtMax<SimpleLinear>>(
        &'alignment self,
        others: &'alignment [Peptidoform<OtherComplexity>],
    ) -> Vec<Vec<(usize, Alignment<'alignment, Complexity, OtherComplexity>)>> {
        others.iter().map(|o| self.align_one(o)).collect()
    }
}

pub trait HasPeptidoform<Complexity> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<Complexity>>>;
}

#[cfg(feature = "imgt")]
impl HasPeptidoform<UnAmbiguous> for crate::imgt::Allele<'_> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<UnAmbiguous>>> {
        Some(Arc::new(self.sequence.clone()))
    }
}

#[cfg(feature = "imgt")]
impl HasPeptidoform<UnAmbiguous> for &crate::imgt::Allele<'_> {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<UnAmbiguous>>> {
        Some(Arc::new(self.sequence.clone()))
    }
}

#[cfg(feature = "identification")]
impl HasPeptidoform<SemiAmbiguous> for crate::identification::FastaData {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<SemiAmbiguous>>> {
        Some(Arc::new(self.peptide().clone()))
    }
}

#[cfg(feature = "identification")]
impl HasPeptidoform<SemiAmbiguous> for &crate::identification::FastaData {
    fn peptidoform(&self) -> Option<Arc<Peptidoform<SemiAmbiguous>>> {
        Some(Arc::new(self.peptide().clone()))
    }
}

#[cfg(all(
    feature = "imgt",
    feature = "identification",
    not(feature = "internal-no-data")
))]
#[test]
fn test_index() {
    use crate::{
        identification::{IdentifiedPeptidoformSource, PeaksData, csv},
        imgt::{AlleleSelection, ChainType, Selection, Species},
    };
    use std::{fs::File, io::BufReader};

    std::fs::create_dir_all("dump").unwrap();
    let file = File::create("dump/test_align_index_results.csv").unwrap();

    // Setup the data needed
    let peptides = PeaksData::parse_reader(
        BufReader::new(
            include_bytes!("../../data/200305_HER_test_04_DENOVO_excerpt.csv").as_slice(),
        ),
        None,
        false,
        None,
    )
    .unwrap()
    .filter_map(|p| p.ok().and_then(|mut p| p.peptide.1.pop()))
    .take(100)
    .collect::<Vec<_>>();
    let germlines = Selection::default()
        .species([Species::HomoSapiens])
        .chain([ChainType::Heavy])
        .allele(AlleleSelection::First)
        .germlines()
        .collect::<Vec<_>>();

    let mut results = Vec::new();

    let range_very_very_low =
        Mass::new::<crate::system::dalton>(0.0)..=Mass::new::<crate::system::dalton>(200.0);
    let range_very_low =
        Mass::new::<crate::system::dalton>(100.0)..=Mass::new::<crate::system::dalton>(300.0);
    let range_low =
        Mass::new::<crate::system::dalton>(200.0)..=Mass::new::<crate::system::dalton>(400.0);
    let range_med =
        Mass::new::<crate::system::dalton>(300.0)..=Mass::new::<crate::system::dalton>(500.0);
    let range_high =
        Mass::new::<crate::system::dalton>(400.0)..=Mass::new::<crate::system::dalton>(600.0);

    for (name, range) in [
        ("very_very_low", range_very_very_low.clone()),
        ("very_low", range_very_low.clone()),
        ("low", range_low.clone()),
        ("med", range_med.clone()),
        ("high", range_high.clone()),
    ] {
        let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
            &germlines,
            MassMode::Monoisotopic,
            AlignScoring::default(),
            AlignType::EITHER_GLOBAL,
            Some((1, range.clone())),
        );

        results.extend(
            index
                .debug_align(&peptides)
                .iter()
                .enumerate()
                .flat_map(|(index, row)| row.iter().map(move |cell| (index, cell)))
                .map(|(index, cell)| {
                    [
                        ("Peptide index".to_string(), index.to_string()),
                        ("Germline index".to_string(), cell.0.to_string()),
                        ("Peptide".to_string(), cell.1.seq_b().to_string()),
                        ("Germline".to_string(), cell.1.seq_a().to_string()),
                        (
                            "Alignment score".to_string(),
                            cell.1.normalised_score().to_string(),
                        ),
                        ("Mass matches".to_string(), cell.2.to_string()),
                        ("Mass range".to_string(), name.to_string()),
                    ]
                }),
        );
    }

    csv::write_csv(file, results, ',').unwrap()
}
