use std::collections::HashMap;
use std::fmt::Write;

use itertools::Itertools;
use rustyms::{
    align::AlignScoring,
    imgt::{AnnotatedSequence, Gene, Germline, Germlines, Species},
    sequence::{Annotation, Peptidoform, Region, UnAmbiguous},
};

use crate::imgt_gene::IMGTGene;
use crate::structs::DataItem;

use crate::structs::SingleSeq;

pub(crate) fn combine(
    data: impl Iterator<Item = Result<DataItem, String>>,
) -> (
    HashMap<Species, Germlines>,
    Vec<(Species, IMGTGene, String)>,
) {
    let mut grouped = HashMap::new();
    let mut errors = Vec::new();
    let mut temp: Vec<(Species, SingleSeq)> = Vec::new();

    for element in data.flatten() {
        let species = element.species;
        for gene in element.genes {
            match gene.clone().finish() {
                Ok(gene) => temp.push((species, gene)),
                Err(err) => {
                    errors.push((species, gene, err));
                }
            }
        }
    }

    // Combine temp seqs
    let mut deduped_temp: Vec<(Species, TemporaryGermline)> = Vec::new();
    'sequences: for (species, seq) in temp {
        for (dspecies, dseq) in &mut deduped_temp {
            if *dspecies == species && dseq.name == seq.name {
                dseq.add(seq);
                continue 'sequences;
            }
        }
        // If not found
        deduped_temp.push((
            species,
            TemporaryGermline {
                name: seq.name.clone(),
                alleles: vec![(seq.allele, vec![TemporarySequence::from_single(seq)])],
            },
        ));
    }

    // Save temp seqs in final data structure
    for (species, entry) in deduped_temp {
        grouped
            .entry(species)
            .or_insert_with(|| Germlines::new(species))
            .insert(entry.finalise());
    }
    (grouped, errors)
}

struct TemporaryGermline {
    name: Gene,
    alleles: Vec<(usize, Vec<TemporarySequence>)>,
}

impl TemporaryGermline {
    fn add(&mut self, single: SingleSeq) {
        if let Some(al) = self.alleles.iter_mut().find(|al| al.0 == single.allele) {
            if let Some(s) =
                al.1.iter_mut()
                    .find(|s| s.sequence == single.sequence.sequence)
            {
                s.add_single(single);
            } else {
                al.1.push(TemporarySequence::from_single(single));
            }
            al.1.sort();
        } else {
            // If not found
            self.alleles
                .push((single.allele, vec![TemporarySequence::from_single(single)]));
            self.alleles.sort_unstable_by_key(|a| a.0); // Maybe do the fancy insert at the right place trick
        }
    }

    fn finalise(self) -> Germline {
        Germline {
            name: self.name,
            alleles: self
                .alleles
                .into_iter()
                .filter_map(|(a, seqs)| {
                    Some((
                        a,
                        seqs.iter()
                            .find(|s| !s.sequence.is_empty())
                            .map(TemporarySequence::annotated_sequence)?,
                    ))
                })
                .collect(),
        }
    }
}

impl std::fmt::Display for TemporaryGermline {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        const MAX_WIDTH: usize = 100; // Pure sequence disregards any spaces in front
        const SHOW_DNA: bool = false;
        writeln!(f, "█ GENE: {}", self.name)?;
        let mut first_allele = None;
        for allele in &self.alleles {
            writeln!(f, "╻ *{:02}", allele.0)?;
            let mut reference = None;
            for (index, seq) in allele.1.iter().enumerate() {
                let main_branch = if index == allele.1.len() - 1 {
                    ' '
                } else {
                    '│'
                };
                writeln!(
                    f,
                    "{}┬╸ACC: {}",
                    if index == allele.1.len() - 1 {
                        '└'
                    } else {
                        '├'
                    },
                    seq.acc.iter().join(" ")
                )?;
                write!(f, "{main_branch}├─SEQ:")?;
                let seq_str = seq.sequence.to_string();
                if seq_str.chars().count() < 90 {
                    writeln!(f, " {seq_str}")?;
                } else {
                    writeln!(f)?;
                    let lines = seq_str
                        .chars()
                        .collect_vec()
                        .chunks(MAX_WIDTH)
                        .map(|c| c.iter().collect::<String>())
                        .collect_vec();
                    for line in lines {
                        writeln!(f, "{main_branch}│ {line}")?;
                    }
                }

                write!(f, "{main_branch}├─REG: ")?;
                let regions = seq.regions();
                if regions.len() > 1 {
                    writeln!(f)?;
                }
                for region in &regions {
                    if regions.len() > 1 {
                        write!(f, "{main_branch}│    ⊕ ")?;
                    }
                    writeln!(
                        f,
                        "{} supported by {}",
                        region.0.iter().fold(String::new(), |mut acc, (r, l)| {
                            write!(&mut acc, "[{r},{l}]").unwrap();
                            acc
                        }),
                        region.1.iter().map(|i| seq.acc[*i].clone()).join(" "),
                    )?;
                }
                write!(f, "{main_branch}├─ANN: ")?;
                let conserved = seq.conserved();
                if conserved.len() > 1 {
                    writeln!(f)?;
                }
                for cons in &conserved {
                    if conserved.len() > 1 {
                        write!(f, "{main_branch}│    ⊕ ")?;
                    }
                    writeln!(
                        f,
                        "{} supported by {}",
                        cons.0.iter().fold(String::new(), |mut acc, (r, l)| {
                            write!(&mut acc, "[{r},{l}]").unwrap();
                            acc
                        }),
                        cons.1.iter().map(|i| seq.acc[*i].clone()).join(" "),
                    )?;
                }
                let scoring = AlignScoring::<'_> {
                    matrix: rustyms::align::matrix::BLOSUM90,
                    ..Default::default()
                };
                if let Some(first_allele) = first_allele {
                    let alignment = rustyms::align::align::<1, UnAmbiguous, UnAmbiguous>(
                        first_allele,
                        &seq.sequence,
                        scoring,
                        rustyms::align::AlignType::GLOBAL,
                    )
                    .stats();
                    writeln!(
                        f,
                        "{main_branch}├─ALLELE DIF: {} AAs",
                        alignment.length - alignment.identical,
                    )?;
                }
                if let Some(reference) = reference {
                    let alignment = rustyms::align::align::<1, UnAmbiguous, UnAmbiguous>(
                        reference,
                        &seq.sequence,
                        scoring,
                        rustyms::align::AlignType::GLOBAL,
                    )
                    .stats();
                    writeln!(
                        f,
                        "{main_branch}├─OPTION DIF: {} AAs",
                        alignment.length - alignment.identical,
                    )?;
                }
                writeln!(
                    f,
                    "{main_branch}└{}DNA: {} sequence(s)",
                    if SHOW_DNA { "┬" } else { "─" },
                    seq.dna.len()
                )?;
                if SHOW_DNA {
                    for (di, dna) in seq.dna.iter().enumerate() {
                        let lines = dna
                            .0
                            .chars()
                            .collect_vec()
                            .chunks(MAX_WIDTH)
                            .map(|c| c.iter().collect::<String>())
                            .collect_vec();
                        for (line_i, line) in lines.iter().enumerate() {
                            writeln!(
                                f,
                                "{main_branch} {}{}",
                                if line_i == 0 {
                                    if di == seq.dna.len() - 1 {
                                        '└'
                                    } else {
                                        '├'
                                    }
                                } else if di == seq.dna.len() - 1 {
                                    ' '
                                } else {
                                    '│'
                                },
                                line
                            )?;
                        }
                        writeln!(
                            f,
                            "{main_branch} {} supported by {}",
                            if di == seq.dna.len() - 1 { ' ' } else { '│' },
                            dna.1.iter().map(|i| seq.acc[*i].clone()).join(" "),
                        )?;
                    }
                }
                if reference.is_none() {
                    reference = Some(&seq.sequence);
                }
                if first_allele.is_none() {
                    first_allele = Some(&seq.sequence);
                }
            }
        }
        writeln!(f)
    }
}

#[derive(Debug, Eq, PartialEq)]
struct TemporarySequence {
    acc: Vec<String>,
    sequence: Peptidoform<UnAmbiguous>,
    regions: HashMap<Vec<(Region, usize)>, Vec<usize>>,
    annotations: HashMap<Vec<(Annotation, usize)>, Vec<usize>>,
    dna: HashMap<String, Vec<usize>>,
}

impl TemporarySequence {
    fn from_single(single: SingleSeq) -> Self {
        Self {
            acc: vec![single.acc],
            sequence: single.sequence.sequence,
            dna: [(single.dna, vec![0])].into(),
            regions: [(single.sequence.regions, vec![0])].into(),
            annotations: [(single.sequence.annotations, vec![0])].into(),
        }
    }

    fn add_single(&mut self, single: SingleSeq) {
        let index = self.acc.len();
        self.acc.push(single.acc);
        self.dna.entry(single.dna).or_default().push(index);
        self.regions
            .entry(single.sequence.regions)
            .or_default()
            .push(index);
        self.annotations
            .entry(single.sequence.annotations)
            .or_default()
            .push(index);
    }

    fn annotated_sequence(&self) -> AnnotatedSequence {
        AnnotatedSequence {
            sequence: self.sequence.clone(),
            regions: self.regions()[0].0.clone(),
            annotations: self.conserved()[0].0.clone(),
        }
    }

    fn regions(&self) -> Vec<(Vec<(Region, usize)>, Vec<usize>)> {
        let mut vec = self
            .regions
            .iter()
            .map(|(r, a)| (r.to_owned(), a.to_owned()))
            .collect_vec();
        vec.sort_by_key(|s| -(s.1.len() as isize));
        vec.sort_by_key(|s| -(s.0.len() as isize));
        vec
    }

    fn conserved(&self) -> Vec<(Vec<(Annotation, usize)>, Vec<usize>)> {
        let mut vec = self
            .annotations
            .iter()
            .map(|(r, a)| (r.to_owned(), a.to_owned()))
            .collect_vec();
        vec.sort_by_key(|s| -(s.1.len() as isize));
        vec.sort_by_key(|s| -(s.0.len() as isize));
        vec
    }
}

impl PartialOrd for TemporarySequence {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for TemporarySequence {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        other.acc.len().cmp(&self.acc.len()).then(
            (other.regions()[0].0.len() + other.conserved()[0].0.len())
                .cmp(&(self.regions()[0].0.len() + self.conserved()[0].0.len())),
        )
    }
}
