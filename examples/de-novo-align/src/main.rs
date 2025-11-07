//! Align many peptides to a database using mass based alignment
#![allow(clippy::type_complexity)]
use std::{collections::BTreeMap, fs::File, io::BufWriter, sync::Arc};

use clap::Parser;
use itertools::{Itertools, MinMaxResult};
use mzalign::prelude::*;
use mzcore::{csv::write_csv, ontology::Ontologies, prelude::*, sequence::SemiAmbiguous};
use mzident::{FastaData, PeptidoformPresent, SpectrumIds, prelude::*};
use rayon::prelude::*;

#[derive(Debug, Parser)]
struct Cli {
    /// The input identified peptidoforms file
    #[arg(short, long)]
    peptides: String,
    /// The fasta database of known proteins
    #[arg(short, long)]
    database: String,
    /// Where to store the results
    #[arg(long)]
    out_path: String,
    /// Run alignment sequentially instead of in parallel (useful for profiling)
    #[arg(long)]
    no_parallel: bool,
}

fn process_alignment_group<I>(
    alignments: I,
) -> Vec<(
    Alignment<Arc<FastaData>, IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>>,
    bool,
)>
where
    I: Iterator<
        Item = Alignment<Arc<FastaData>, IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>>,
    >,
{
    let alignments = alignments.collect_vec();
    let max = alignments
        .iter()
        .max_by(|a, b| a.normalised_score().total_cmp(&b.normalised_score()))
        .map(Alignment::normalised_score)
        .unwrap_or_default();
    let mut alignments = alignments
        .into_iter()
        .filter(|a| (a.normalised_score() - max).abs() < 1E-9)
        .collect_vec();
    if alignments.len() == 1
        && let Some(a) = alignments.pop()
    {
        vec![(a, true)]
    } else {
        alignments.into_iter().map(|a| (a, false)).collect_vec()
    }
}

fn run_alignments(
    index: &AlignIndex<4, Arc<FastaData>>,
    peptides: Vec<IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>>,
    scoring: AlignScoring,
    sequential: bool,
) -> Vec<(
    Alignment<Arc<FastaData>, IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>>,
    bool,
)> {
    if sequential {
        index
            .align(peptides, scoring, AlignType::EITHER_GLOBAL)
            .flat_map(process_alignment_group)
            .collect()
    } else {
        index
            .par_align(peptides, scoring, AlignType::EITHER_GLOBAL)
            .flat_map(process_alignment_group)
            .collect()
    }
}

#[allow(clippy::missing_panics_doc)]
fn main() {
    let args = Cli::parse();
    let out_file = BufWriter::new(File::create(args.out_path).expect("Could not create out file"));
    let peptides = open_identified_peptidoforms_file(args.peptides, &Ontologies::init().0, false)
        .unwrap()
        .filter_map(Result::ok)
        .filter_map(IdentifiedPeptidoform::into_semi_ambiguous)
        .collect_vec();
    let database = FastaData::parse_file(args.database).expect("Could not open database");
    let index = AlignIndex::<4, Arc<FastaData>>::new(
        database.into_iter().map(Arc::new),
        MassMode::Monoisotopic,
    );

    let scoring = AlignScoring {
        pair: PairMode::DatabaseToPeptidoform,
        ..Default::default()
    };

    let alignments: Vec<_> = run_alignments(&index, peptides, scoring, args.no_parallel)
        .into_iter()
        .map(|(alignment, unique)| {
            BTreeMap::from([
                (
                    "Peptide".to_string(),
                    alignment.seq_b().peptidoform().to_string(),
                ),
                (
                    "Spectra ref".to_string(),
                    match alignment.seq_b().scans() {
                        SpectrumIds::None => String::new(),
                        SpectrumIds::FileNotKnown(scans) => scans.iter().join(";"),
                        SpectrumIds::FileKnown(scans) => scans
                            .iter()
                            .map(|(file, scans)| {
                                format!("{}:{}", file.to_string_lossy(), scans.iter().join(";"))
                            })
                            .join("|"),
                    },
                ),
                (
                    "De novo score".to_string(),
                    alignment
                        .seq_b()
                        .score
                        .map_or(String::new(), |s| s.to_string()),
                ),
                (
                    "Protein".to_string(),
                    alignment.seq_a().identifier().to_string(),
                ),
                (
                    "Alignment score".to_string(),
                    alignment.normalised_score().to_string(),
                ),
                ("Unique".to_string(), unique.to_string()),
                ("Start".to_string(), alignment.start_a().to_string()),
                (
                    "End".to_string(),
                    (alignment.start_a() + alignment.len_a()).to_string(),
                ),
                ("Path".to_string(), alignment.short()),
                (
                    "Mass".to_string(),
                    match alignment.seq_b().peptidoform().formulas().mass_bounds() {
                        MinMaxResult::NoElements => "-".to_string(),
                        MinMaxResult::OneElement(m) => m.monoisotopic_mass().value.to_string(),
                        MinMaxResult::MinMax(min, max) => format!(
                            "{} - {}",
                            min.monoisotopic_mass().value,
                            max.monoisotopic_mass().value
                        ),
                    },
                ),
                (
                    "Z".to_string(),
                    alignment
                        .seq_b()
                        .charge()
                        .map_or(0, |c| c.value)
                        .to_string(),
                ),
                (
                    "Peptide length".to_string(),
                    alignment.seq_b().peptidoform().len().to_string(),
                ),
                (
                    "Retention time".to_string(),
                    alignment
                        .seq_b()
                        .retention_time()
                        .map_or(f64::NAN, |t| t.value)
                        .to_string(),
                ),
            ])
        })
        .collect();

    write_csv(out_file, alignments, ',').unwrap();
}
