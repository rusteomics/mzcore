//! Align many peptides to a database using mass based alignment
use std::{collections::HashMap, fs::File, io::BufWriter};

use align::AlignScoring;
use clap::Parser;
use identification::SpectrumIds;
use itertools::{Itertools, MinMaxResult};
use rayon::prelude::*;
use rustyms::{
    align::{AlignType, align},
    identification::{
        FastaData, IdentifiedPeptidoform, PeptidoformPresent, csv::write_csv,
        open_identified_peptides_file,
    },
    prelude::Peptidoform,
    sequence::{HasPeptidoform, HasPeptidoformImpl, SemiAmbiguous},
    *,
};

#[derive(Debug, Parser)]
struct Cli {
    /// The input identified peptides file
    #[arg(short, long)]
    peptides: String,
    /// The fasta database of known proteins
    #[arg(short, long)]
    database: String,
    /// Where to store the results
    #[arg(long)]
    out_path: String,
}

fn main() {
    let args = Cli::parse();
    let out_file = BufWriter::new(File::create(args.out_path).unwrap());
    let peptides = open_identified_peptides_file(args.peptides, None, false)
        .unwrap()
        .filter_map(Result::ok)
        .filter_map(IdentifiedPeptidoform::into_semi_ambiguous)
        .collect_vec();
    let database = FastaData::parse_file(args.database).unwrap();

    let alignments: Vec<_> = peptides
        .par_iter()
        .flat_map(|peptide| {
            let alignments = database
                .iter()
                .map(|db| {
                    (
                        db,
                        peptide,
                        align::<
                            4,
                            &Peptidoform<SemiAmbiguous>,
                            &IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent>,
                        >(
                            db.peptide(),
                            peptide,
                            AlignScoring::default(),
                            AlignType::EITHER_GLOBAL,
                        ),
                    )
                })
                .collect_vec();
            let max = alignments
                .iter()
                .max_by(|a, b| a.2.normalised_score().total_cmp(&b.2.normalised_score()))
                .unwrap()
                .2
                .normalised_score();
            let mut alignments = alignments
                .into_iter()
                .filter(|a| a.2.normalised_score() == max)
                .collect_vec();
            if alignments.len() == 1 {
                let (d, p, a) = alignments.pop().unwrap();
                vec![(d, p, a, true)]
            } else {
                alignments
                    .into_iter()
                    .map(|(d, p, a)| (d, p, a, false))
                    .collect_vec()
            }
        })
        .map(|(db, peptide, alignment, unique)| {
            HashMap::from([
                (
                    "Peptide".to_string(),
                    alignment.seq_b().peptidoform().to_string(),
                ),
                (
                    "Spectra ref".to_string(),
                    match peptide.scans() {
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
                    peptide.score.map_or(String::new(), |s| s.to_string()),
                ),
                ("Protein".to_string(), db.identifier().to_string()),
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
                    peptide.charge().map_or(0, |c| c.value).to_string(),
                ),
                (
                    "Peptide length".to_string(),
                    peptide.peptidoform().len().to_string(),
                ),
                (
                    "Retention time".to_string(),
                    peptide
                        .retention_time()
                        .map_or(f64::NAN, |t| t.value)
                        .to_string(),
                ),
            ])
        })
        .collect();

    write_csv(out_file, alignments, ',').unwrap();
}
