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
            let res =
                $crate::sequence::CompoundPeptidoformIon::pro_forma($case, None).map(|(a, _)| a);
            let upper = $case.to_ascii_uppercase();
            let lower = $case.to_ascii_lowercase();
            let res_upper =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&upper, None).map(|(a, _)| a);
            let res_lower =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&lower, None).map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&back, None).map(|(a, _)| a);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (just_parse $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res =
                $crate::sequence::CompoundPeptidoformIon::pro_forma($case, None).map(|(a, _)| a);
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
            let res =
                $crate::sequence::CompoundPeptidoformIon::pro_forma($case, None).map(|(a, _)| a);
            let upper = $case.to_ascii_uppercase();
            let lower = $case.to_ascii_lowercase();
            let res_upper =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&upper, None).map(|(a, _)| a);
            let res_lower =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&lower, None).map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&back, None).map(|(a, _)| a);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (casing_specific $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res =
                $crate::sequence::CompoundPeptidoformIon::pro_forma($case, None).map(|(a, _)| a);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap().into_iter().join("\n"));
            let back = res.as_ref().unwrap().to_string();
            let res_back =
                $crate::sequence::CompoundPeptidoformIon::pro_forma(&back, None).map(|(a, _)| a);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::CompoundPeptidoformIon::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

/// Create a sloppy parse test based on a given case and its name.
#[macro_export]
macro_rules! parse_sloppy_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::Peptidoform::sloppy_pro_forma(
                $case,
                0..$case.len(),
                None,
                SloppyParsingParameters::default(),
            );
            let res_leading_n = $crate::sequence::Peptidoform::sloppy_pro_forma(
                $case,
                0..$case.len(),
                None,
                SloppyParsingParameters {
                    ignore_prefix_lowercase_n: true,
                },
            );
            let res_upper = $crate::sequence::Peptidoform::sloppy_pro_forma(
                &$case.to_ascii_uppercase(),
                0..$case.len(),
                None,
                SloppyParsingParameters::default(),
            );
            let res_lower = $crate::sequence::Peptidoform::sloppy_pro_forma(
                &$case.to_ascii_lowercase(),
                0..$case.len(),
                None,
                SloppyParsingParameters::default(),
            );
            println!("{}", $case);
            dbg!(&res);
            assert!(res.is_ok());
            assert_eq!(res, res_leading_n);
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::sequence::Peptidoform::sloppy_pro_forma(
                $case,
                0..$case.len(),
                None,
                &SloppyParsingParameters::default(),
            );
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}
