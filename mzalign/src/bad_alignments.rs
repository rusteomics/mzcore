#![allow(clippy::missing_panics_doc)]

use crate::{AlignScoring, AlignType, align};
use mzcore::sequence::{Peptidoform, SimpleLinear};

#[test]
fn overextended_rotation() {
    test_alignment("IVQEVS", "LEVQVES", "1i1I2=2r1=");
}

fn test_alignment(peptide_one: &str, peptide_two: &str, path: &str) {
    let first_peptide = Peptidoform::pro_forma(peptide_one, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let second_peptide = Peptidoform::pro_forma(peptide_two, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let alignment = align::<4, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
        &first_peptide,
        &second_peptide,
        AlignScoring::default(),
        AlignType::GLOBAL,
    );
    assert_eq!(
        alignment.short(),
        path,
        "Alignment of {peptide_one} vs {peptide_two} did not go to plan"
    );
}
