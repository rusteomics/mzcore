//! Benchmark a simple peptides to germline alignment

use std::{
    fs::File,
    io::{BufReader, BufWriter, prelude::*},
    ops::RangeInclusive,
    time::{Duration, Instant},
};

use rustyms::{
    align::{AlignScoring, AlignType, OneToManyIndex, align},
    identification::{IdentifiedPeptidoformSource, PeaksData},
    imgt::{Allele, AlleleSelection, ChainType, Selection, Species},
    prelude::*,
    quantities::Tolerance,
    sequence::{SemiAmbiguous, UnAmbiguous},
    system::Mass,
};

fn main() {
    std::fs::create_dir_all("dump").unwrap();
    let file = File::create("dump/benchmark_results.csv").unwrap();

    // Setup the data needed
    let peptides = PeaksData::parse_reader(
        BufReader::new(include_bytes!("../data/200305_HER_test_04_DENOVO_excerpt.csv").as_slice()),
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

    let range_low =
        Mass::new::<rustyms::system::dalton>(200.0)..=Mass::new::<rustyms::system::dalton>(300.0);
    let range_med =
        Mass::new::<rustyms::system::dalton>(300.0)..=Mass::new::<rustyms::system::dalton>(400.0);
    let range_high =
        Mass::new::<rustyms::system::dalton>(400.0)..=Mass::new::<rustyms::system::dalton>(500.0);

    // Do the benchmarking
    results.push(measure(bench_align, &(&germlines, &peptides), "Align"));
    results.push(measure(
        bench_index,
        &(&germlines, &peptides, None),
        "Index - none",
    ));
    for matches in 1..10 {
        for (name, range) in [
            ("low", range_low.clone()),
            ("med", range_med.clone()),
            ("high", range_high.clone()),
        ] {
            results.push(measure(
                bench_index,
                &(&germlines, &peptides, Some((matches, range.clone()))),
                &format!("Index - {matches} - {name}"),
            ));
        }
    }

    // Save the results to a csv
    let mut sink = BufWriter::new(file);
    sink.write_all("Name,Average(ns),StandardDeviation(ns),Runs\n".as_bytes())
        .unwrap();
    for item in results {
        sink.write_fmt(format_args!(
            "{},{},{},{}\n",
            item.0, item.1, item.2, item.3
        ))
        .unwrap();
    }
    sink.flush().unwrap();
}

fn bench_align(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let results: Vec<Vec<_>> = data
        .0
        .iter()
        .map(|g| {
            data.1
                .iter()
                .map(|p| {
                    align::<4, UnAmbiguous, SemiAmbiguous>(
                        g.sequence,
                        p,
                        AlignScoring::default(),
                        AlignType::EITHER_GLOBAL,
                    )
                })
                .collect()
        })
        .collect();
}

fn bench_index(
    data: &(
        &[Allele<'_>],
        &[Peptidoform<SemiAmbiguous>],
        Option<(u8, RangeInclusive<Mass>)>,
    ),
) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        data.2.clone(),
    );

    index.align(data.1);
}

fn measure_multiple<T>(
    function: fn(&T),
    subjects: &[(&str, T)],
    description: &str,
) -> Vec<(String, u128, u128, usize)> {
    let mut output = Vec::with_capacity(subjects.len());
    for (name, item) in subjects {
        output.push(measure(function, item, &format!("{description} - {name}")));
    }
    output
}

fn measure<T>(function: fn(&T), subject: &T, description: &str) -> (String, u128, u128, usize) {
    let mut times = Vec::new();
    function(subject);
    let start = Instant::now();
    let mut now;

    for _ in 0..5 {
        now = Instant::now();
        function(subject);
        times.push(now.elapsed());
    }

    let average = start.elapsed().checked_div(5).unwrap();

    // Lets run for 3 more seconds, including cloning of the subject
    for _ in 0..3_000_000_000 / average.as_nanos() {
        now = Instant::now();
        function(subject);
        times.push(now.elapsed());
    }

    let average = times
        .iter()
        .fold(Duration::new(0, 0), |total, item| {
            total.checked_add(*item).unwrap()
        })
        .checked_div(times.len() as u32)
        .unwrap();

    let mut deviation = 0.0;
    for run in &times {
        let difference = run.as_nanos() as f64 - average.as_nanos() as f64;
        deviation += difference * difference;
    }
    deviation /= times.len() as f64;
    deviation = deviation.sqrt();

    (
        description.to_string(),
        average.as_nanos(),
        deviation as u128,
        times.len(),
    )
}
