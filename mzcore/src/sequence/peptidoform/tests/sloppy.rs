use std::sync::Arc;

use crate::{
    molecular_formula,
    ontology::Ontology,
    parse_sloppy_test,
    sequence::{
        Modification, ModificationId, Peptidoform, SemiAmbiguous, SimpleModificationInner,
        SloppyParsingParameters,
    },
};

#[test]
fn sloppy_names() {
    assert_eq!(
        Modification::sloppy_modification(
            "Deamidation (NQ)",
            0..16,
            None,
            &crate::ontology::STATIC_ONTOLOGIES
        ),
        Ok(crate::ontology::STATIC_ONTOLOGIES
            .unimod()
            .get_by_name("deamidated")
            .unwrap())
    );
    assert_eq!(
        Modification::sloppy_modification(
            "Pyro-glu from Q",
            0..15,
            None,
            &crate::ontology::STATIC_ONTOLOGIES
        ),
        Ok(crate::ontology::STATIC_ONTOLOGIES
            .unimod()
            .get_by_name("gln->pyro-glu")
            .unwrap())
    );
}

#[test]
fn sloppy_names_custom() {
    let ontologies = crate::ontology::Ontologies::init_static().with_custom([Arc::new(
        SimpleModificationInner::Database {
            formula: molecular_formula!(O 1),
            id: ModificationId {
                ontology: Ontology::Custom,
                name: "Test".to_string(),
                id: Some(0),
                ..Default::default()
            },
            specificities: Vec::new(),
        },
    )]);
    assert!(Modification::sloppy_modification("test", 0..4, None, &ontologies).is_ok());
    assert!(Modification::sloppy_modification("Test", 0..4, None, &ontologies).is_ok());
    assert!(Modification::sloppy_modification("C:Test", 0..6, None, &ontologies).is_ok());
}

#[test]
fn sloppy_msfragger() {
    assert_eq!(
        Peptidoform::<SemiAmbiguous>::sloppy_pro_forma(
            "n[211]GC[779]RQSSEEK",
            0..20,
            &crate::ontology::STATIC_ONTOLOGIES,
            &SloppyParsingParameters {
                ignore_prefix_lowercase_n: true,
                ..Default::default()
            }
        )
        .unwrap(),
        Peptidoform::pro_forma("[211]-GC[779]RQSSEEK", &crate::ontology::STATIC_ONTOLOGIES)
            .unwrap()
            .0
            .into_semi_ambiguous()
            .unwrap()
    );
}

parse_sloppy_test!(ne "_", fuzz_01);
parse_sloppy_test!(ne "ffffffff[gln->|yro-glu]SC2N:iTRAQ4pleeeeeB]", hang_01);
parse_sloppy_test!(ne "SEQUEN[Formula:[13B2YC2][12Cu2]HKKKyro-g|||||||||||||@@||||||||||||||lmmmmmm|||| |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||o-glu]n[13YEQUEeedISEQU9SEmmmm]SBSE-@CSE->pyro-glm]n`n->pyrogl>pyro-gl", hang_02);
