#![allow(clippy::missing_panics_doc)]
mod fuzz_crash;
mod fuzz_hang;
mod parse;
mod pro_forma_negative;
mod pro_forma_positive;
mod sloppy;

#[macro_export]
macro_rules! parse_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::CompoundPeptidoform::pro_forma($case, None);
            let res_upper =
                $crate::CompoundPeptidoform::pro_forma(&$case.to_ascii_uppercase(), None);
            let res_lower =
                $crate::CompoundPeptidoform::pro_forma(&$case.to_ascii_lowercase(), None);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap());
            assert_eq!(res, res_upper);
            assert_eq!(res, res_lower);
            let back = res.as_ref().unwrap().to_string();
            let res_back = $crate::CompoundPeptidoform::pro_forma(&back, None);
            assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::CompoundPeptidoform::pro_forma($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

#[macro_export]
macro_rules! parse_sloppy_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::LinearPeptide::sloppy_pro_forma(
                $case,
                0..$case.len(),
                None,
                SloppyParsingParameters::default(),
            );
            let res_leading_n = $crate::LinearPeptide::sloppy_pro_forma(
                $case,
                0..$case.len(),
                None,
                SloppyParsingParameters {
                    ignore_prefix_lowercase_n: true,
                },
            );
            let res_upper = $crate::LinearPeptide::sloppy_pro_forma(
                &$case.to_ascii_uppercase(),
                0..$case.len(),
                None,
                SloppyParsingParameters::default(),
            );
            let res_lower = $crate::LinearPeptide::sloppy_pro_forma(
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
            let res = $crate::LinearPeptide::sloppy_pro_forma(
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
