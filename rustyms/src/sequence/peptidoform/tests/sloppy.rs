use std::sync::Arc;

use crate::{
    molecular_formula,
    ontology::Ontology,
    parse_sloppy_test,
    sequence::{
        Modification, Peptidoform, SemiAmbiguous, SimpleModificationInner, SloppyParsingParameters,
    },
};

#[test]
fn sloppy_names() {
    assert_eq!(
        Modification::sloppy_modification("Deamidation (NQ)", 0..16, None, None),
        Ok(Ontology::Unimod.find_name("deamidated", None).unwrap())
    );
    assert_eq!(
        Modification::sloppy_modification("Pyro-glu from Q", 0..15, None, None),
        Ok(Ontology::Unimod.find_name("gln->pyro-glu", None).unwrap())
    );
}

#[test]
fn sloppy_names_custom() {
    let db = Some(vec![(
        Some(0),
        "test".to_string(),
        Arc::new(SimpleModificationInner::Formula(molecular_formula!(O 1))),
    )]);
    assert_eq!(
        Modification::sloppy_modification("test", 0..4, None, db.as_ref()),
        Ok(Arc::new(SimpleModificationInner::Formula(
            molecular_formula!(O 1)
        )))
    );
    assert_eq!(
        Modification::sloppy_modification("Test", 0..4, None, db.as_ref()),
        Ok(Arc::new(SimpleModificationInner::Formula(
            molecular_formula!(O 1)
        )))
    );
    assert_eq!(
        Modification::sloppy_modification("C:Test", 0..6, None, db.as_ref()),
        Ok(Arc::new(SimpleModificationInner::Formula(
            molecular_formula!(O 1)
        )))
    );
}

#[test]
fn sloppy_msfragger() {
    assert_eq!(
        Peptidoform::<SemiAmbiguous>::sloppy_pro_forma(
            "n[211]GC[779]RQSSEEK",
            0..20,
            None,
            &SloppyParsingParameters {
                ignore_prefix_lowercase_n: true,
                ..Default::default()
            }
        )
        .unwrap(),
        Peptidoform::pro_forma("[211]-GC[779]RQSSEEK", None)
            .unwrap()
            .into_semi_ambiguous()
            .unwrap()
    );
}

parse_sloppy_test!(ne "_", fuzz_01);
parse_sloppy_test!(ne "ffffffff[gln->|yro-glu]SC2N:iTRAQ4pleeeeeB]", hang_01);
parse_sloppy_test!(ne "SEQUEN[Formula:[13B2YC2][12Cu2]HKKKyro-g|||||||||||||@@||||||||||||||lmmmmmm|||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||o-glu]n[13YEQUEeedISEQU9SEmmmm]SBSE-@CSE->pyro-glm]n`n->pyrogl>pyro-gl", hang_02);
