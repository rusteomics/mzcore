#![allow(clippy::missing_panics_doc)]
mod fuzz_crash;
mod fuzz_hang;
mod parse;
mod pro_forma_negative;
mod pro_forma_positive;
mod sloppy;

/// Create a parse test based on a given case and its name.
#[macro_export]
macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            let upper = $case.to_ascii_uppercase();
            let lower = $case.to_ascii_lowercase();
            let res_upper = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &upper,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            let res_lower = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &lower,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let (res_back, back_warnings) =
                $crate::sequence::CompoundPeptidoformIon::pro_forma_strict(
                    &back,
                    &$crate::ontology::STATIC_ONTOLOGIES,
                )
                .unwrap();
            assert_eq!(res.unwrap(), res_back, "{} != {back}", $case);
            assert!(back_warnings.len() == 0, "{back_warnings:?}");
        }
    };
    (just_parse $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
        }
    };
    (ignore $case:literal, $name:ident) => {
        #[test]
        #[allow(clippy::ignore_without_reason)]
        #[ignore]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            let upper = $case.to_ascii_uppercase();
            let lower = $case.to_ascii_lowercase();
            let res_upper = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &upper,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            let res_lower = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &lower,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &back,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (casing_specific $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            let back = res.as_ref().unwrap().to_string();
            let res_back = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                &back,
                &$crate::ontology::STATIC_ONTOLOGIES,
            )
            .map(|(a, _)| a);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
            );
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

/// Create a sloppy parse test based on a given case and its name.
#[macro_export]
macro_rules! parse_sloppy_test {
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::Peptidoform::sloppy_pro_forma(
                $case,
                0..$case.len(),
                &$crate::ontology::STATIC_ONTOLOGIES,
                &SloppyParsingParameters::default(),
            );
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}
