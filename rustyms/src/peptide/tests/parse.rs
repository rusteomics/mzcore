use std::{num::NonZeroU16, sync::Arc};

use crate::{
    model::PrimaryIonSeries,
    modification::{self, ModificationId, SimpleModificationInner},
    peptide::{
        parse::{global_modifications, parse_charge_state},
        GlobalModification,
    },
    placement_rule::{self, PlacementRule},
    system::{da, usize::Charge},
    AminoAcid, CompoundPeptidoform, Element, LinearPeptide, Model, MolecularCharge, MultiChemical,
    Peptidoform,
};

#[test]
fn parse_global_modifications() {
    let parse = |str: &str| global_modifications(str, 0, None);
    assert_eq!(
        parse("<[+5]@D>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::Anywhere,
                Some(AminoAcid::AsparticAcid),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@d>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::Anywhere,
                Some(AminoAcid::AsparticAcid),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@N-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyNTerm,
                Some(AminoAcid::AsparticAcid),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@n-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyNTerm,
                Some(AminoAcid::AsparticAcid),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@C-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                crate::placement_rule::Position::AnyCTerm,
                Some(AminoAcid::AsparticAcid),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<D>"),
        Ok((
            3,
            vec![GlobalModification::Isotope(Element::H, NonZeroU16::new(2))]
        ))
    );
    assert_eq!(
        parse("<12C>"),
        Ok((
            5,
            vec![GlobalModification::Isotope(Element::C, NonZeroU16::new(12))]
        ))
    );
    assert!(parse("<D").is_err());
    assert!(parse("<[+5]>").is_err());
    assert!(parse("<[+5]@DD>").is_err());
    assert!(parse("<[5+]@D>").is_err());
    assert!(parse("<[+5@D>").is_err());
    assert!(parse("<+5]@D>").is_err());
    assert!(parse("<[+5#g1]@D>").is_err());
    assert!(parse("<[+5#g1>").is_err());
    assert!(parse("<C12>").is_err());
    assert!(parse("<>").is_err());
    assert!(parse("<@>").is_err());
    assert!(parse("<@D,E,R,T>").is_err());
    assert!(parse("<[+5]@D,E,R,Te>").is_err());
    assert!(parse("<[+5]@D,E,R,N-term:OO>").is_err());
}

#[test]
fn charge_state_positive() {
    let parse = |str: &str| {
        parse_charge_state(str, 0).map(|(len, res)| {
            assert_eq!(
                len,
                str.len(),
                "Not full parsed: '{str}', amount parsed: {len} as '{res}'"
            );
            res
        })
    };
    assert_eq!(parse("/1"), Ok(MolecularCharge::proton(1)));
    assert_eq!(parse("/5"), Ok(MolecularCharge::proton(5)));
    assert_eq!(parse("/-5"), Ok(MolecularCharge::proton(-5)));
    assert_eq!(parse("/1[+H+]"), Ok(MolecularCharge::proton(1)));
    assert_eq!(parse("/2[+H+,+H+]"), Ok(MolecularCharge::proton(2)));
    assert_eq!(
        parse("/1[+Na+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Na 1 Electron -1)
        )]))
    );
    assert_eq!(
        parse("/3[2Na+1,1H1+1]"),
        Ok(MolecularCharge::new(&[
            (2, molecular_formula!(Na 1 Electron -1)),
            (1, molecular_formula!(H 1 Electron -1))
        ]))
    );
    assert_eq!(
        parse("/1[-OH-]"),
        Ok(MolecularCharge::new(&[(
            -1,
            molecular_formula!(O 1 H 1 Electron 1)
        ),]))
    );
    assert_eq!(
        parse("/1[+N1H3+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(N 1 H 3 Electron -1)
        ),]))
    );
    assert_eq!(
        parse("/1[+[15N1]+]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!([15 N 1] Electron -1)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+ Fe +3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/-1[+e-]"),
        Ok(MolecularCharge::new(
            &[(1, molecular_formula!(Electron 1)),]
        ))
    );
    assert_eq!(parse("/1[+H1e-1+]"), Ok(MolecularCharge::proton(1)));
    assert_eq!(
        parse("/3[+Fe1e0+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-1+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-2+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
    assert_eq!(
        parse("/3[+Fe1e-3+3]"),
        Ok(MolecularCharge::new(&[(
            1,
            molecular_formula!(Fe 1 Electron -3)
        ),]))
    );
}

#[test]
fn charge_state_negative() {
    let parse = |str: &str| parse_charge_state(str, 0);
    assert!(parse("/3[+Fe+]").is_err());
    assert!(parse("/3[+Fe]").is_err());
    assert!(parse("/3[+Fe 1]").is_err());
    assert!(parse("/3[+[54Fe1+3]").is_err());
    assert!(parse("/3[+54Fe1]+3]").is_err());
    assert!(parse("/1[1H1-1]").is_err());
    assert!(parse("/1[1H1+1").is_err());
    assert!(parse("/1[1+1]").is_err());
    assert!(parse("/1[H+1]").is_err());
    assert!(parse("/1[1H]").is_err());
    assert!(parse("/1[1H1]").is_err());
    assert!(parse("/ 1 [ 1 H 1]").is_err());
}

#[test]
fn parse_glycan() {
    let glycan = LinearPeptide::pro_forma("A[Glycan:Hex]", None).unwrap();
    let spaces = LinearPeptide::pro_forma("A[Glycan:    Hex    ]", None).unwrap();
    assert_eq!(glycan.len(), 1);
    assert_eq!(spaces.len(), 1);
    assert_eq!(glycan, spaces);
    let incorrect = CompoundPeptidoform::pro_forma("A[Glycan:Hec]", None);
    assert!(incorrect.is_err());
}

#[test]
fn parse_formula() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma("A[Glycan:Hex]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    assert_eq!(peptide.len(), 1);
    assert_eq!(glycan.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_labile() {
    let with = LinearPeptide::pro_forma("{Formula:C6H10O5}A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    assert_eq!(with.len(), 1);
    assert_eq!(without.len(), 1);
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(
        with.get_labile()[0].to_string(),
        "Formula:C6H10O5".to_string()
    );
}

#[test]
fn parse_ambiguous_modification() {
    let with = LinearPeptide::pro_forma("A[Phospho#g0]A[#g0]", None).unwrap();
    let without = LinearPeptide::pro_forma("AA", None).unwrap();
    assert_eq!(with.len(), 2);
    assert_eq!(without.len(), 2);
    assert_eq!(with[0].modifications.len(), 1);
    assert_eq!(with[1].modifications.len(), 1);
    assert!(CompoundPeptidoform::pro_forma("A[#g0]A[#g0]", None).is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[Phospho#g0]", None).is_err());
    assert!(CompoundPeptidoform::pro_forma("A[Phospho#g0]A[#g0(0.o1)]", None).is_err());
    assert_eq!(
        LinearPeptide::pro_forma("A[+12#g0]A[#g0]", None)
            .unwrap()
            .to_string(),
        "A[+12#g0]A[#g0]".to_string()
    );
    assert_eq!(
        LinearPeptide::pro_forma("A[#g0]A[+12#g0]", None)
            .unwrap()
            .to_string(),
        "A[#g0]A[+12#g0]".to_string()
    );
}

#[test]
fn parse_ambiguous_aminoacid() {
    let with = LinearPeptide::pro_forma("(?AA)C(?A)(?A)", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let without = LinearPeptide::pro_forma("AACAA", None)
        .unwrap()
        .into_linear()
        .unwrap();
    assert_eq!(with.len(), 5);
    assert_eq!(without.len(), 5);
    assert!(with[0].ambiguous.is_some());
    assert!(with[1].ambiguous.is_some());
    assert_eq!(with.formulas(), without.formulas());
    assert_eq!(with.to_string(), "(?AA)C(?A)(?A)".to_string());
}

#[test]
fn parse_hard_tags() {
    let peptide = LinearPeptide::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let glycan = LinearPeptide::pro_forma(
        "A[info:you can define a tag multiple times|Glycan:Hex|Formula:C6H10O5]",
        None,
    )
    .unwrap()
    .into_linear()
    .unwrap();
    assert_eq!(peptide.len(), 1);
    assert_eq!(glycan.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_global() {
    let deuterium = LinearPeptide::pro_forma("<D>A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let nitrogen_15 = LinearPeptide::pro_forma("<15N>A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    assert_eq!(deuterium.len(), 1);
    assert_eq!(nitrogen_15.len(), 1);
    // Formula: A + H2O
    assert_eq!(
        deuterium.formulas(),
        molecular_formula!([2 H 7] C 3 O 2 N 1).into()
    );
    assert_eq!(
        nitrogen_15.formulas(),
        molecular_formula!(H 7 C 3 O 2 [15 N 1]).into()
    );
}

#[test]
fn parse_chimeric() {
    let dimeric = CompoundPeptidoform::pro_forma("A+AA", None).unwrap();
    let trimeric = dbg!(CompoundPeptidoform::pro_forma("A+AA-[+2]+AAA", None).unwrap());
    assert_eq!(dimeric.peptidoforms().len(), 2);
    assert_eq!(dimeric.peptidoforms()[0].peptides()[0].len(), 1);
    assert_eq!(dimeric.peptidoforms()[1].peptides()[0].len(), 2);
    assert_eq!(trimeric.peptidoforms().len(), 3);
    assert_eq!(trimeric.peptidoforms()[0].peptides()[0].len(), 1);
    assert_eq!(trimeric.peptidoforms()[1].peptides()[0].len(), 2);
    assert_eq!(trimeric.peptidoforms()[2].peptides()[0].len(), 3);
    assert!(trimeric.peptidoforms()[1].peptides()[0]
        .get_c_term()
        .is_some());
}

#[test]
fn parse_unimod() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma(
        "[U:Gln->pyro-Glu]-QE[Cation:Na]AA",
        None
    ));
    assert!(peptide.is_ok());
}

#[test]
fn parse_custom() {
    let peptide = dbg!(CompoundPeptidoform::pro_forma(
        "A[C:WEEE]",
        Some(&vec![(
            Some(0),
            "weee".to_string(),
            SimpleModificationInner::Database {
                formula: molecular_formula!(U 1),
                specificities: vec![(
                    vec![PlacementRule::AminoAcid(
                        AminoAcid::CANONICAL_AMINO_ACIDS.to_vec(),
                        placement_rule::Position::Anywhere
                    )],
                    Vec::new(),
                    Vec::new()
                )],
                id: ModificationId {
                    ontology: modification::Ontology::Custom,
                    name: "WEEE".to_string(),
                    id: Some(0),
                    ..ModificationId::default()
                }
            }
            .into()
        )])
    ));
    assert!(peptide.is_ok());
    assert_eq!(
        peptide.as_ref().unwrap().to_string(),
        "A[Formula:U1|INFO:Custom:WEEE]"
    );
    assert_eq!(
        peptide.unwrap().formulas(),
        molecular_formula!(C 3 H 7 N 1 O 2 U 1).into()
    );
}

#[test]
fn parse_xl_intra() {
    let peptide = Peptidoform::pro_forma("A[XLMOD:02001#XLTEST]A[#XLTEST]", None).unwrap();
    println!("{peptide}");
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        peptide.formulas().to_vec()[0],
        (AminoAcid::Alanine.single_formula().unwrap() * 2)
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1).with_label(crate::AmbiguousLabel::CrossLinkBound(
                crate::CrossLinkName::Name("test".to_string())
            ))
    );
}

#[test]
fn parse_xl_inter() {
    let peptide = Peptidoform::pro_forma("A[XLMOD:02001#XLTEST]//A[#XLTEST]", None).unwrap();
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        peptide.formulas().to_vec()[0],
        (AminoAcid::Alanine.single_formula().unwrap() * 2
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1) * 2)
            .with_label(crate::AmbiguousLabel::CrossLinkBound(
                crate::CrossLinkName::Name("test".to_string())
            ))
    );
}

#[test]
fn dimeric_peptide() {
    // Only generate a single series, easier to reason about
    let test_model = Model::none().a(PrimaryIonSeries::default());

    // With two different sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+CC", None).unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, aC, pAA, pCC

    // With two identical sequences
    let dimeric = CompoundPeptidoform::pro_forma("AA+AA", None).unwrap();
    let fragments = dbg!(dimeric
        .generate_theoretical_fragments(Charge::new::<crate::system::charge::e>(1), &test_model));
    assert_eq!(fragments.len(), 4); // aA, pAA (both twice once for each peptide)
}

#[test]
fn parse_adduct_ions_01() {
    let peptide = CompoundPeptidoform::pro_forma("A/2[2Na+]+A", None).unwrap();
    assert_eq!(peptide.peptidoforms().len(), 2);
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0]
            .get_charge_carriers()
            .unwrap()
            .charge_carriers,
        vec![(2, molecular_formula!(Na 1 Electron -1))]
    );
    assert_eq!(
        peptide.peptidoforms()[0].peptides()[0].sequence(),
        peptide.peptidoforms()[1].peptides()[0].sequence()
    );
}

#[test]
fn hydrolysed_xl() {
    let peptide_xl = LinearPeptide::pro_forma("EMEVTK[XLMOD:02001]SESPEK", None)
        .unwrap()
        .into_unambiguous()
        .unwrap();
    let peptide_mod = LinearPeptide::pro_forma("EMEVTK[Formula:C8H12O3]SESPEK", None)
        .unwrap()
        .into_unambiguous()
        .unwrap();

    assert_eq!(peptide_xl.formula(), peptide_mod.formula());
}
