#![allow(clippy::missing_panics_doc)]
mod negative;
mod parse;
mod positive;
mod sloppy;

/// Create a sloppy parse test based on a given case and its name.
#[macro_export]
macro_rules! parse_sloppy_test {
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::Peptidoform::sloppy_pro_forma(
                $case,
                &$crate::ontology::STATIC_ONTOLOGIES,
                &SloppyParsingParameters::default(),
            );
            assert!(res.is_err(), "{}\n{}", $case, res.unwrap());
        }
    };
}

/// Test a batch of tests from a CSV file (the file has to be pasted inline for compiler reasons)
#[macro_export]
macro_rules! negative_tests {
    (
        Id,Example,Source,Key,Notes
        $($id:literal,$case:literal,$source:ident,$key:tt,$note:literal)*
    ) => {
            $($crate::negative_tests!($id,$case,$source,$key,$note);)*
        };
    ($id:literal,$case:literal,$source:ident,"ignore",$note:literal) => {
        paste::paste!{
            #[ignore]
            #[test]
            fn [<$source _ $id>]() {
                // Empty but still generated to show the correct number of ignored tests in the output
            }
    }};
    ($id:literal,$case:literal,$source:ident,$key:literal,$note:literal) => {
        paste::paste!{
            #[test]
            fn [<$source _ $id>]() {
                let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                    $case,
                    &$crate::ontology::STATIC_ONTOLOGIES,
                );
                println!("{}\n{:?}", $case, res);
                assert!(res.is_err() || res.is_ok_and(|(_, errs)| !errs.is_empty()));
            }
        }
    }
}

/// Test a batch of tests from a CSV file (the file has to be pasted inline for compiler reasons)
#[macro_export]
macro_rules! positive_tests {
    (
        Id,Example,Source,Key,Notes
        $($id:literal,$case:literal,$source:ident,$key:tt,$note:literal)*
    ) => {
            $($crate::positive_tests!($id,$case,$source,$key,$note);)*
        };
    ($id:literal,$case:literal,$source:ident,"ignore",$note:literal) => {
        paste::paste!{
            #[ignore]
            #[test]
            fn [<$source _ $id>]() {
                // Empty but still generated to show the correct number of ignored tests in the output
            }
    }};
    ($id:literal,$case:literal,$source:ident,"casing_specific",$note:literal) => {
        paste::paste!{
            #[test]
            fn [<$source _ $id>]() {
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
    }};
    ($id:literal,$case:literal,$source:ident,"just_parse",$note:literal) => {
        paste::paste!{
            #[test]
            fn [<$source _ $id>]() {
                use itertools::Itertools;
                let res = $crate::sequence::CompoundPeptidoformIon::pro_forma(
                    $case,
                    &$crate::ontology::STATIC_ONTOLOGIES,
                )
                .map(|(a, _)| a);
                println!("{}", $case);
                assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            }
    }};
    ($id:literal,$case:literal,$source:ident,"ignore_warnings",$note:literal) => {
        paste::paste!{
            #[test]
            fn [<$source _ $id>]() {
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
                let (res_back, _) =
                    $crate::sequence::CompoundPeptidoformIon::pro_forma_strict(
                        &back,
                        &$crate::ontology::STATIC_ONTOLOGIES,
                    )
                    .unwrap();
                assert_eq!(res.unwrap(), res_back, "{} != {back}", $case);
            }
        }
    };
    ($id:literal,$case:literal,$source:ident,$key:literal,$note:literal) => {
        paste::paste!{
            #[test]
            fn [<$source _ $id>]() {
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
        }
    }
}
