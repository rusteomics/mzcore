#![allow(clippy::missing_panics_doc, clippy::float_cmp)]
//! Integration tests for alignments
use mzalign::prelude::*;
use mzcore::{prelude::*, sequence::SimpleLinear};
use mzident::{PeptidoformPresent, csv::parse_csv, prelude::*};

const MAB_SEQUENCE: &str = "GPGGGGSLPEVREKHEFLNRLKQLPLLESQIATIEQSAPSQSDQEQLFSNVQYFAHYCRKYAPLYAAEAKRVFSLEKKMSNYIQFKSKCRIEPVCLLLHGSPGAGKSVATNLIGRSLAEKLNSSVYSLPPDPDHFDGYKQQAVVIMDDLCQNPDGKDVSLFCQMVSSVDFVPPMAALEEKGILFTSPFVLASTNAGSINAPTVSDSRALARRFHFDMNIEVISMYSQNGKINMPMSVKTCDDECCPVNFKKCCPLVCGKAIQFIDRRTQVRYSLDMLVTEMFREYNHRHSVGTTLEALFQ";

/// Test on many peptidoforms if alignment to itself always results in a perfect match
#[test]
fn identity() {
    let peptidoforms: Vec<_> =
        open_identified_peptidoforms_file("tests/example_peptidoforms.csv", None, false)
            .unwrap()
            .filter_map(Result::ok)
            .filter_map(IdentifiedPeptidoform::into_simple_linear)
            .collect();

    for peptidoform in &peptidoforms {
        let alignment = align::<
            4,
            &IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent>,
            &IdentifiedPeptidoform<SimpleLinear, PeptidoformPresent>,
        >(
            peptidoform,
            peptidoform,
            AlignScoring::default(),
            AlignType::GLOBAL,
        );
        assert_eq!(alignment.score().normalised, 1.0);
        assert_eq!(alignment.stats().identity(), 1.0);
    }
}

/// Test for a set of peptidoforms from a database matcher if the match to the database is indeed perfect
#[test]
fn alignment_to_database() {
    let peptidoforms: Vec<_> =
        open_identified_peptidoforms_file("tests/example_peptidoforms.csv", None, false)
            .unwrap()
            .filter_map(Result::ok)
            .filter_map(IdentifiedPeptidoform::into_simple_linear)
            .collect();
    let database = Peptidoform::pro_forma(MAB_SEQUENCE, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let index = AlignIndex::<4, Peptidoform<SimpleLinear>>::new([database], MassMode::Monoisotopic);

    for (id, peptidoform) in peptidoforms.iter().enumerate() {
        let alignment = index
            .align_one(
                peptidoform,
                AlignScoring {
                    pair: PairMode::DatabaseToPeptidoform,
                    ..Default::default()
                },
                AlignType::GLOBAL_B,
            )
            .next()
            .unwrap();
        assert!(
            !(alignment.score().normalised != 1.0 || alignment.stats().identity() != 1.0),
            "Alignment of ({id}) {} did not work out: {}",
            peptidoform.peptidoform(),
            alignment.short()
        );
    }
}

/// Test for a set of randomly paired peptidoforms if the alignment algorithm still gives the same result
#[test]
fn pairwise() {
    use std::io::Write;

    let lines: Vec<_> = parse_csv("tests/pairwise_examples.csv", b',', None)
        .unwrap()
        .filter_map(Result::ok)
        .collect();
    let mut out =
        std::io::BufWriter::new(std::fs::File::create("tests/pairwise_examples_new.csv").unwrap());
    writeln!(
        &mut out,
        "a,b,path,score,absolute score,maximal score,identical,mass similar,gaps,length"
    )
    .unwrap();
    let mut different_paths = Vec::new();

    for (id, line) in lines.iter().enumerate() {
        let a = Peptidoform::pro_forma(line.index_column("a").unwrap().0, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let b = Peptidoform::pro_forma(line.index_column("b").unwrap().0, None)
            .unwrap()
            .into_simple_linear()
            .unwrap();
        let path = line.index_column("path").unwrap().0;
        let absolute_score: isize = line
            .index_column("absolute score")
            .unwrap()
            .0
            .parse()
            .unwrap();
        let maximal_score: isize = line
            .index_column("maximal score")
            .unwrap()
            .0
            .parse()
            .unwrap();
        let alignment = align::<4, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
            &a,
            &b,
            AlignScoring::default(),
            AlignType::GLOBAL,
        );
        let line_nm = id + 2;
        let score = alignment.score();
        let stats = alignment.stats();
        assert_eq!(score.absolute, absolute_score, "Pair at line {line_nm}",);
        assert_eq!(score.max, maximal_score, "Pair at line {line_nm}",);
        if alignment.short() != path {
            different_paths.push((alignment.short(), path, line_nm));
        }
        writeln!(
            &mut out,
            "{},{},{},{},{},{},{},{},{},{}",
            line.index_column("a").unwrap().0,
            line.index_column("b").unwrap().0,
            alignment.short(),
            score.normalised,
            score.absolute,
            score.max,
            stats.identical,
            stats.mass_similar,
            stats.gaps,
            stats.length
        )
        .unwrap();
    }

    if !different_paths.is_empty() {
        for (found, expected, line) in &different_paths {
            println!("Line {line}: found '{found}' instead of '{expected}'");
        }

        println!("See the newly created 'pairwise_examples_new.csv' for the full new data");
        panic!(
            "{} pairs had a different path (but are otherwise similarly optimal)",
            different_paths.len()
        );
    }
}
