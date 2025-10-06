use std::{num::NonZeroU16, sync::Arc};

use context_error::*;

use crate::{
    chemistry::{AmbiguousLabel, Element, MolecularCharge, MultiChemical},
    molecular_formula,
    ontology::Ontology,
    sequence::{
        AminoAcid, CompoundPeptidoformIon, CrossLinkName, GlobalModification, ModificationId,
        Peptidoform, PeptidoformIon, PlacementRule, Position, SimpleModificationInner,
        peptidoform::parse::{global_modifications, parse_charge_state},
    },
    system::da,
};

#[test]
fn parse_global_modifications() {
    let parse = |str: &str| global_modifications(str, 0, None).map_err(BoxedError::to_owned);
    assert_eq!(
        parse("<[+5]@D>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                PlacementRule::AminoAcid(vec![AminoAcid::AsparticAcid], Position::Anywhere),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@d>"),
        Ok((
            8,
            vec![GlobalModification::Fixed(
                PlacementRule::AminoAcid(vec![AminoAcid::AsparticAcid], Position::Anywhere),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@N-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                PlacementRule::AminoAcid(vec![AminoAcid::AsparticAcid], Position::AnyNTerm),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@n-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                PlacementRule::AminoAcid(vec![AminoAcid::AsparticAcid], Position::AnyNTerm),
                Arc::new(SimpleModificationInner::Mass(da(5.0).into()))
            )]
        ))
    );
    assert_eq!(
        parse("<[+5]@C-term:D>"),
        Ok((
            15,
            vec![GlobalModification::Fixed(
                PlacementRule::AminoAcid(vec![AminoAcid::AsparticAcid], Position::AnyCTerm),
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
        parse_charge_state(str, 0)
            .map(|(len, res)| {
                assert_eq!(
                    len,
                    str.len(),
                    "Not full parsed: '{str}', amount parsed: {len} as '{res}'"
                );
                res
            })
            .map_err(BoxedError::to_owned)
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
}

#[test]
fn charge_state_negative() {
    let parse = |str: &str| parse_charge_state(str, 0).map_err(BoxedError::to_owned);
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
    let glycan = Peptidoform::pro_forma("A[Glycan:Hex]", None).unwrap();
    let spaces = Peptidoform::pro_forma("A[Glycan:    Hex    ]", None).unwrap();
    assert_eq!(glycan.len(), 1);
    assert_eq!(spaces.len(), 1);
    assert_eq!(glycan, spaces);
    let incorrect = CompoundPeptidoformIon::pro_forma("A[Glycan:Hec]", None);
    assert!(incorrect.is_err());
}

#[test]
fn parse_formula() {
    let peptide = Peptidoform::pro_forma("A[Formula:C6H10O5]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let glycan = Peptidoform::pro_forma("A[Glycan:Hex]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    assert_eq!(peptide.len(), 1);
    assert_eq!(glycan.len(), 1);
    assert_eq!(glycan.formulas(), peptide.formulas());
}

#[test]
fn parse_labile() {
    let with = Peptidoform::pro_forma("{Formula:C6H10O5}A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let without = Peptidoform::pro_forma("A", None)
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
    let with = Peptidoform::pro_forma("A[Phospho#g0]A[#g0]", None).unwrap();
    let without = Peptidoform::pro_forma("AA", None).unwrap();
    assert_eq!(with.len(), 2);
    assert_eq!(without.len(), 2);
    assert_eq!(with[0].modifications.len(), 1);
    assert_eq!(with[1].modifications.len(), 1);
    assert!(CompoundPeptidoformIon::pro_forma("A[#g0]A[#g0]", None).is_err());
    assert!(CompoundPeptidoformIon::pro_forma("A[Phospho#g0]A[Phospho#g0]", None).is_err());
    assert!(CompoundPeptidoformIon::pro_forma("A[Phospho#g0]A[#g0(0.o1)]", None).is_err());
    assert_eq!(
        Peptidoform::pro_forma("A[+12#g0]A[#g0]", None)
            .unwrap()
            .to_string(),
        "A[+12#g0]A[#g0]".to_string()
    );
    assert_eq!(
        Peptidoform::pro_forma("A[#g0]A[+12#g0]", None)
            .unwrap()
            .to_string(),
        "A[#g0]A[+12#g0]".to_string()
    );
}

#[test]
fn parse_terminal_ambiguous_modification() {
    // N-term
    let unplaced_n = Peptidoform::pro_forma("[deamidated]?FAAQAA", None).unwrap();
    assert!(unplaced_n.get_n_term()[0].is_ambiguous());
    assert_eq!(unplaced_n.sequence()[3].modifications.len(), 1);
    assert!(unplaced_n.sequence()[3].modifications[0].is_ambiguous());
    let placed_n = Peptidoform::pro_forma("[deamidated#u1]-FAAQ[#u1]AA", None).unwrap();
    assert!(placed_n.get_n_term()[0].is_ambiguous());
    assert_eq!(placed_n.sequence()[3].modifications.len(), 1);
    assert!(placed_n.sequence()[3].modifications[0].is_ambiguous());
    // C-term
    let unplaced_c = Peptidoform::pro_forma("[oxidation]?AHAMTEG", None).unwrap();
    assert!(unplaced_c.get_c_term()[0].is_ambiguous());
    assert_eq!(unplaced_c.sequence()[3].modifications.len(), 1);
    assert!(unplaced_c.sequence()[3].modifications[0].is_ambiguous());
    let placed_c = Peptidoform::pro_forma("AHAM[oxidation#u1]TEG-[#u1]", None).unwrap();
    assert!(placed_c.get_c_term()[0].is_ambiguous());
    assert_eq!(placed_c.sequence()[3].modifications.len(), 1);
    assert!(placed_c.sequence()[3].modifications[0].is_ambiguous());
}

#[test]
fn parse_ambiguous_aminoacid() {
    let with = Peptidoform::pro_forma("(?AA)C(?A)(?A)", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let without = Peptidoform::pro_forma("AACAA", None)
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
    let peptide = Peptidoform::pro_forma("A[Formula:C6H10O5|INFO:hello world 🦀]", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let glycan = Peptidoform::pro_forma(
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
    let deuterium = Peptidoform::pro_forma("<D>A", None)
        .unwrap()
        .into_linear()
        .unwrap();
    let nitrogen_15 = Peptidoform::pro_forma("<15N>A", None)
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
    let dimeric = CompoundPeptidoformIon::pro_forma("A+AA", None).unwrap();
    let trimeric = dbg!(CompoundPeptidoformIon::pro_forma("A+AA-[+2]+AAA", None).unwrap());
    assert_eq!(dimeric.peptidoform_ions().len(), 2);
    assert_eq!(dimeric.peptidoform_ions()[0].peptidoforms()[0].len(), 1);
    assert_eq!(dimeric.peptidoform_ions()[1].peptidoforms()[0].len(), 2);
    assert_eq!(trimeric.peptidoform_ions().len(), 3);
    assert_eq!(trimeric.peptidoform_ions()[0].peptidoforms()[0].len(), 1);
    assert_eq!(trimeric.peptidoform_ions()[1].peptidoforms()[0].len(), 2);
    assert_eq!(trimeric.peptidoform_ions()[2].peptidoforms()[0].len(), 3);
    assert_eq!(
        trimeric.peptidoform_ions()[1].peptidoforms()[0]
            .get_c_term()
            .len(),
        1
    );
}

#[test]
fn parse_unimod() {
    let peptide = dbg!(CompoundPeptidoformIon::pro_forma(
        "[U:Gln->pyro-Glu]-QE[Cation:Na]AA",
        None
    ));
    assert!(peptide.is_ok());
}

#[test]
fn parse_custom() {
    let peptide = dbg!(CompoundPeptidoformIon::pro_forma(
        "A[C:WEEE]",
        Some(&vec![(
            Some(0),
            "weee".to_string(),
            SimpleModificationInner::Database {
                formula: molecular_formula!(U 1),
                specificities: vec![(
                    vec![PlacementRule::AminoAcid(
                        AminoAcid::CANONICAL_AMINO_ACIDS.to_vec(),
                        Position::Anywhere
                    )],
                    Vec::new(),
                    Vec::new()
                )],
                id: ModificationId {
                    ontology: Ontology::Custom,
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
    let peptide = PeptidoformIon::pro_forma("A[XLMOD:02001#XLTEST]A[#XLTEST]", None).unwrap();
    println!("{peptide}");
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        peptide.formulas().to_vec()[0],
        (AminoAcid::Alanine.single_formula().unwrap() * 2)
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1).with_label(AmbiguousLabel::CrossLinkBound(
                CrossLinkName::Name("test".to_string())
            ))
    );
}

#[test]
fn parse_xl_inter() {
    let peptide = PeptidoformIon::pro_forma("A[XLMOD:02001#XLTEST]//A[#XLTEST]", None).unwrap();
    //dbg!(&singular.sequence[0].modifications);
    assert_eq!(
        peptide.formulas().to_vec()[0],
        (AminoAcid::Alanine.single_formula().unwrap() * 2_i32
            + molecular_formula!(C 8 H 10 O 2)
            + molecular_formula!(H 2 O 1) * 2_i32)
            .with_label(AmbiguousLabel::CrossLinkBound(CrossLinkName::Name(
                "test".to_string()
            )))
    );
}

#[test]
fn parse_adduct_ions_01() {
    let peptide = CompoundPeptidoformIon::pro_forma("A/2[2Na+]+A", None).unwrap();
    assert_eq!(peptide.peptidoform_ions().len(), 2);
    assert_eq!(
        peptide.peptidoform_ions()[0].peptidoforms()[0]
            .get_charge_carriers()
            .unwrap()
            .charge_carriers,
        vec![(2, molecular_formula!(Na 1 Electron -1))]
    );
    assert_eq!(
        peptide.peptidoform_ions()[0].peptidoforms()[0].sequence(),
        peptide.peptidoform_ions()[1].peptidoforms()[0].sequence()
    );
}

#[test]
fn hydrolysed_xl() {
    let peptide_xl = Peptidoform::pro_forma("EMEVTK[XLMOD:02001]SESPEK", None)
        .unwrap()
        .into_unambiguous()
        .unwrap();
    let peptide_mod = Peptidoform::pro_forma("EMEVTK[Formula:C8H12O3]SESPEK", None)
        .unwrap()
        .into_unambiguous()
        .unwrap();

    assert_eq!(peptide_xl.formula(), peptide_mod.formula());
}

#[test]
fn multiple_n_mods() {
    assert!(Peptidoform::pro_forma("[U:Carbamyl][+17.027]-EGEEEK", None).is_ok());
}

#[test]
fn multiple_c_mods() {
    assert!(Peptidoform::pro_forma("PEPTIDEG-[Methyl][Amidated]", None).is_ok());
}

#[test]
fn ambiguous_aas() {
    assert!(Peptidoform::pro_forma("[+320.85]-(?N[U:Deamidated])-[+432.85]", None).is_ok());
    assert!(Peptidoform::pro_forma("[+1550.7049484934485]-T(?HQ[U:Deamidated])", None).is_ok());
}

#[ignore = "A known bug"]
#[test]
fn ambiguous_mods() {
    let sequence = "[U:Oxidation#u0]?FTSPFVLASTNAGSINAPTVSDSRALARRFHFDM[#u0]NIEVISM[#u0]YSQNGKINM[#u0]PM[#u0]SVKTCDDE";
    let parsed = Peptidoform::pro_forma(sequence, None).unwrap();
    assert_eq!(sequence, parsed.to_string());
    assert_eq!(parsed.get_ambiguous_modifications().len(), 1);
}
