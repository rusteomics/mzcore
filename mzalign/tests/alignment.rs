#![allow(clippy::missing_panics_doc, clippy::float_cmp)]
//! Integration tests for alignments
use mzalign::prelude::*;
use mzcore::{csv::parse_csv, prelude::*, sequence::SimpleLinear};
use mzident::{PeptidoformPresent, prelude::*};

const MAB_SEQUENCE: &str = "GPGGGGSLPEVREKHEFLNRLKQLPLLESQIATIEQSAPSQSDQEQLFSNVQYFAHYCRKYAPLYAAEAKRVFSLEKKMSNYIQFKSKCRIEPVCLLLHGSPGAGKSVATNLIGRSLAEKLNSSVYSLPPDPDHFDGYKQQAVVIMDDLCQNPDGKDVSLFCQMVSSVDFVPPMAALEEKGILFTSPFVLASTNAGSINAPTVSDSRALARRFHFDMNIEVISMYSQNGKINMPMSVKTCDDECCPVNFKKCCPLVCGKAIQFIDRRTQVRYSLDMLVTEMFREYNHRHSVGTTLEALFQ";

/// Test on many peptidoforms if alignment to itself always results in a perfect match
#[test]
fn identity() {
    let peptidoforms: Vec<_> = open_identified_peptidoforms_file(
        "tests/example_peptidoforms.csv",
        &mzcore::ontology::STATIC_ONTOLOGIES,
        false,
    )
    .unwrap()
    .filter_map(Result::ok)
    .filter_map(PSM::into_simple_linear)
    .collect();

    for peptidoform in &peptidoforms {
        let alignment = align::<
            4,
            &PSM<SimpleLinear, PeptidoformPresent>,
            &PSM<SimpleLinear, PeptidoformPresent>,
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
    let peptidoforms: Vec<_> = open_identified_peptidoforms_file(
        "tests/example_peptidoforms.csv",
        &mzcore::ontology::STATIC_ONTOLOGIES,
        false,
    )
    .unwrap()
    .filter_map(Result::ok)
    .filter_map(PSM::into_simple_linear)
    .collect();
    let database = Peptidoform::pro_forma(MAB_SEQUENCE, &mzcore::ontology::STATIC_ONTOLOGIES)
        .unwrap()
        .0
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
    let mut differences = Vec::new();

    for (id, line) in lines.iter().enumerate() {
        let a = Peptidoform::pro_forma(
            line.index_column("a").unwrap().0,
            &mzcore::ontology::STATIC_ONTOLOGIES,
        )
        .unwrap()
        .0
        .into_simple_linear()
        .unwrap();
        let b = Peptidoform::pro_forma(
            line.index_column("b").unwrap().0,
            &mzcore::ontology::STATIC_ONTOLOGIES,
        )
        .unwrap()
        .0
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
        if score.absolute != absolute_score {
            differences.push((
                "absolute score",
                score.absolute.to_string(),
                absolute_score.to_string(),
                line_nm,
            ));
        }
        if score.max != maximal_score {
            differences.push((
                "max score",
                score.max.to_string(),
                maximal_score.to_string(),
                line_nm,
            ));
        }
        if alignment.short() != path {
            differences.push(("path", alignment.short(), path.to_string(), line_nm));
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

    if !differences.is_empty() {
        for (what, found, expected, line) in &differences {
            println!("Line {line}: found {what} '{found}' instead of '{expected}'");
        }

        println!("See the newly created 'pairwise_examples_new.csv' for the full new data");
        panic!("{} pairs had differences", differences.len());
    }
}
