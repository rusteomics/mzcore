//! Benchmark a simple peptides to germline alignment

use std::{
    fs::File,
    io::{BufReader, BufWriter, prelude::*},
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

    // Do the benchmarking
    results.push(measure(bench_align, &(&germlines, &peptides), "Align"));
    results.push(measure(
        bench_index,
        &(&germlines, &peptides),
        "Index - none",
    ));
    results.push(measure(
        bench_index_1,
        &(&germlines, &peptides),
        "Index - 1",
    ));
    results.push(measure(
        bench_index_2,
        &(&germlines, &peptides),
        "Index - 2",
    ));
    results.push(measure(
        bench_index_3,
        &(&germlines, &peptides),
        "Index - 3",
    ));
    results.push(measure(
        bench_index_4,
        &(&germlines, &peptides),
        "Index - 4",
    ));
    results.push(measure(
        bench_index_5,
        &(&germlines, &peptides),
        "Index - 5",
    ));
    results.push(measure(
        bench_index_6,
        &(&germlines, &peptides),
        "Index - 6",
    ));
    results.push(measure(
        bench_index_7,
        &(&germlines, &peptides),
        "Index - 7",
    ));
    results.push(measure(
        bench_index_8,
        &(&germlines, &peptides),
        "Index - 8",
    ));
    results.push(measure(
        bench_index_9,
        &(&germlines, &peptides),
        "Index - 9",
    ));

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

fn scoring() -> AlignScoring<'static> {
    AlignScoring::<'_> {
        tolerance: Tolerance::Absolute(Mass::new::<rustyms::system::dalton>(0.02).into()),
        ..Default::default()
    }
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
                        scoring(),
                        AlignType::EITHER_GLOBAL,
                    )
                })
                .collect()
        })
        .collect();
}

fn bench_index(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        None,
    );

    index.align(data.1);
}

fn bench_index_1(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            1,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_2(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            2,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_3(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            3,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_4(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            4,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_5(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            5,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_6(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            6,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_7(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            7,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_8(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            8,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
    );

    index.align(data.1);
}

fn bench_index_9(data: &(&[Allele<'_>], &[Peptidoform<SemiAmbiguous>])) {
    let index: OneToManyIndex<'_, 4, UnAmbiguous> = OneToManyIndex::new(
        data.0,
        MassMode::Monoisotopic,
        scoring(),
        AlignType::EITHER_GLOBAL,
        Some((
            9,
            Mass::new::<rustyms::system::dalton>(200.0)
                ..=Mass::new::<rustyms::system::dalton>(300.0),
        )),
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
