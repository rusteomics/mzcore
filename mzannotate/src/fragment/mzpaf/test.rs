use std::sync::LazyLock;

use crate::fragment::mzpaf::write::ToMzPAF;

static BASIC_ANALYTES: LazyLock<[(std::num::NonZeroU32, crate::mzspeclib::AnalyteTarget); 2]> =
    LazyLock::new(|| {
        [
            (
                std::num::NonZeroU32::new(1).unwrap(),
                crate::mzspeclib::AnalyteTarget::PeptidoformIon(
                    mzcore::sequence::PeptidoformIon::pro_forma(
                        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                        &mzcore::ontology::STATIC_ONTOLOGIES,
                    )
                    .unwrap()
                    .0,
                ),
            ),
            (
                std::num::NonZeroU32::new(2).unwrap(),
                crate::mzspeclib::AnalyteTarget::PeptidoformIon(
                    mzcore::sequence::PeptidoformIon::pro_forma(
                        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                        &mzcore::ontology::STATIC_ONTOLOGIES,
                    )
                    .unwrap()
                    .0,
                ),
            ),
        ]
    });

/// Create a parse test based on a given case and its name.
#[macro_export]
macro_rules! mzpaf_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::fragment::Fragment::mz_paf(
                $case,
                &mzcore::ontology::STATIC_ONTOLOGIES,
                BASIC_ANALYTES.as_slice(),
            );
            match res {
                Err(err) => {
                    println!("Failed: '{}'", $case);
                    println!("{err:?}");
                    panic!("Failed test")
                }
                Ok((res, _)) => {
                    let back = res.iter().map(|a| a.to_mz_paf_string()).join(",");
                    let res_back = $crate::fragment::Fragment::mz_paf_strict(
                        &back,
                        &mzcore::ontology::STATIC_ONTOLOGIES,
                        BASIC_ANALYTES.as_slice(),
                    );
                    match res_back {
                        Ok((res_back, warnings)) => {
                            let back_back = res_back.iter().map(|a| a.to_mz_paf_string()).join(",");
                            assert_eq!(
                                back, back_back,
                                "{back} != {back_back} (from input: {})",
                                $case
                            );
                            if !warnings.is_empty() {
                                println!("{warnings:?}");
                                panic!("Output was not fully spec compliant")
                            }
                        }
                        Err(err) => {
                            println!("Failed: '{}' was exported as '{back}'", $case);
                            println!("{err:?}");
                            panic!("Failed test")
                        }
                    }
                }
            };
        }
    };
    (strict $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let res = $crate::fragment::Fragment::mz_paf_strict(
                $case,
                &mzcore::ontology::STATIC_ONTOLOGIES,
                BASIC_ANALYTES.as_slice(),
            );
            match res {
                Err(err) => {
                    println!("Failed: '{}'", $case);
                    println!("{err:?}");
                    panic!("Failed test")
                }
                Ok((res, warnings)) => {
                    let back = res.iter().map(|a| a.to_mz_paf_string()).join(",");
                    let res_back = $crate::fragment::Fragment::mz_paf_strict(
                        &back,
                        &mzcore::ontology::STATIC_ONTOLOGIES,
                        BASIC_ANALYTES.as_slice(),
                    );
                    match res_back {
                        Ok((res_back, warnings)) => {
                            let back_back = res_back.iter().map(|a| a.to_mz_paf_string()).join(",");
                            assert_eq!(
                                back, back_back,
                                "{back} != {back_back} (from input: {})",
                                $case
                            );
                            if !warnings.is_empty() {
                                println!("{warnings:?}");
                                panic!("Output was not fully spec compliant")
                            }
                        }
                        Err(err) => {
                            println!("Failed: '{}' was exported as '{back}'", $case);
                            println!("{err:?}");
                            panic!("Failed test")
                        }
                    }
                    if !warnings.is_empty() {
                        println!("{warnings:?}");
                        panic!("Failed test")
                    }
                }
            };
        }
    };
    (ne strict $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::fragment::Fragment::mz_paf_strict(
                $case,
                &mzcore::ontology::STATIC_ONTOLOGIES,
                BASIC_ANALYTES.as_slice(),
            );
            match res {
                Err(err) => {
                    println!("Failed: '{}'", $case);
                    println!("{err:?}");
                    panic!("Failed test")
                }
                Ok((_res, warnings)) => {
                    if warnings.is_empty() {
                        panic!("Example should have failed strict parsing")
                    }
                }
            };
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res =
                $crate::fragment::Fragment::mz_paf($case, &mzcore::ontology::STATIC_ONTOLOGIES, &[
                ]);
            //println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

mzpaf_test!(strict "b2-H2O/3.2ppm,b4-H2O^2/3.2ppm", spec_positive_1);
mzpaf_test!(strict "b2-H2O/3.2ppm*0.75,b4-H2O^2/3.2ppm*0.25", spec_positive_2);
mzpaf_test!(strict "1@y12/0.13,2@b9-NH3/0.23", spec_positive_3);
mzpaf_test!(strict "0@y1{K}", spec_positive_4);
mzpaf_test!(strict "0@y1{K}-NH3", spec_positive_5);
mzpaf_test!(strict "y1/-1.4ppm", spec_positive_6);
mzpaf_test!(strict "y1/-0.0002", spec_positive_7);
mzpaf_test!(strict "y4-H2O+2i[M+H+Na]^2", spec_positive_8);
mzpaf_test!(strict "?", spec_positive_9);
mzpaf_test!(strict "?^3", spec_positive_10);
mzpaf_test!(strict "?+2i^4", spec_positive_11);
mzpaf_test!(strict "?17", spec_positive_12);
mzpaf_test!(strict "?17+i/1.45ppm", spec_positive_13);
mzpaf_test!(strict "?17-H2O/-0.87ppm", spec_positive_14);
mzpaf_test!(strict "0@b2{LL}", spec_positive_15);
mzpaf_test!(strict "0@y1{K}", spec_positive_16);
mzpaf_test!(strict "0@b2{LC[Carbamidomethyl]}", spec_positive_17);
mzpaf_test!(strict "0@b1{[Acetyl]-M}", spec_positive_18);
mzpaf_test!(strict "0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19a);
mzpaf_test!(strict
    "0@y14{M[Oxidation]ACKAAAAAAAAAA}-CH4OS[M+H+Na]^2",
    spec_positive_19b
);
mzpaf_test!(strict "m3:6", spec_positive_20);
mzpaf_test!(strict "b3-C2H3NO", spec_positive_21);
mzpaf_test!(strict "m3:6-CO", spec_positive_22);
mzpaf_test!(strict "m3:6-CO-H2O^2", spec_positive_23);
mzpaf_test!(strict "m3:5/1.1ppm,m4:6/1.1ppm", spec_positive_24);
mzpaf_test!(strict "m3:5", spec_positive_25);
mzpaf_test!(strict "IY", spec_positive_26);
mzpaf_test!(strict "IH", spec_positive_27);
mzpaf_test!(strict "IL-CH2", spec_positive_28);
mzpaf_test!(strict "IC[Carbamidomethyl]", spec_positive_29);
mzpaf_test!(strict "IY[Phospho]", spec_positive_30);
mzpaf_test!(strict "IC[+58.005]", spec_positive_31);
mzpaf_test!(strict "p^2", spec_positive_32a);
mzpaf_test!(strict "p^-2", spec_positive_32b);
mzpaf_test!(strict "p-H3PO4^2", spec_positive_33);
mzpaf_test!(strict "p^4", spec_positive_34);
mzpaf_test!(strict "p+H^3", spec_positive_35);
mzpaf_test!(strict "p^3", spec_positive_36);
mzpaf_test!(strict "p+2H^2", spec_positive_37);
mzpaf_test!(strict "p^2", spec_positive_38);
mzpaf_test!(strict "p+H^2", spec_positive_39);
mzpaf_test!(strict "p+3H", spec_positive_40);
mzpaf_test!(strict "p+2H", spec_positive_41);
mzpaf_test!(strict "p+H", spec_positive_42);
mzpaf_test!(strict "p", spec_positive_43);
mzpaf_test!(strict "r[TMT127N]", spec_positive_44);
mzpaf_test!(strict "r[iTRAQ114]", spec_positive_45);
mzpaf_test!(strict "r[TMT6plex]", spec_positive_46);
mzpaf_test!(strict "r[Hex]", spec_positive_47);
mzpaf_test!(strict "r[Adenine]", spec_positive_48);
mzpaf_test!(strict "0@_{Urocanic Acid}", spec_positive_49);
mzpaf_test!(strict "f{C13H9}/-0.55ppm", spec_positive_50);
mzpaf_test!(strict "f{C12H9N}/0.06ppm", spec_positive_51);
mzpaf_test!(strict "f{C13H9N}/-2.01ppm", spec_positive_52);
mzpaf_test!(strict "f{C13H10N}/-0.11ppm", spec_positive_53);
mzpaf_test!(strict "f{C13H11N}/-0.09ppm", spec_positive_54);
mzpaf_test!(strict "f{C13H12N}/0.26ppm", spec_positive_55);
mzpaf_test!(strict "f{C14H10N}/0.19ppm", spec_positive_56);
mzpaf_test!(strict "f{C14H11N}/0.45ppm", spec_positive_57);
mzpaf_test!(strict "f{C14H10NO}/0.03ppm", spec_positive_58);
mzpaf_test!(strict "f{C16H22O}+i^3", spec_positive_59);
mzpaf_test!(strict "f{C15[13C1]H22O}^3", spec_positive_60);
mzpaf_test!(strict "s{CN=C=O}[M+H]/-0.55ppm", spec_positive_61);
mzpaf_test!(strict "s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm", spec_positive_62);
mzpaf_test!(strict "p-[Hex]", spec_positive_63);
mzpaf_test!(strict "y2+CO-H2O", spec_positive_64);
mzpaf_test!(strict "y2-H2O-NH3", spec_positive_65);
mzpaf_test!(strict "p-2H2O-HPO3-[TMT6plex]", spec_positive_66);
mzpaf_test!(strict "p-2[iTRAQ115]", spec_positive_67);
mzpaf_test!(strict "p-CO-H2O-HPO3-[iTRAQ116]", spec_positive_68);
mzpaf_test!(strict "y2-NH3-[2H1]", spec_positive_69);
mzpaf_test!(strict "y5-H2[18O1][M+Na]", spec_positive_70);
mzpaf_test!(strict "y12+i", spec_positive_71);
mzpaf_test!(strict "y12+i13C", spec_positive_72);
mzpaf_test!(strict "y12+i15N", spec_positive_73);
mzpaf_test!(strict "y12+2i13C+i15N", spec_positive_74);
mzpaf_test!(strict "y12+iA", spec_positive_75);
mzpaf_test!(strict "y12+2iA", spec_positive_76);
mzpaf_test!(strict "y4[M+Na]", spec_positive_77);
mzpaf_test!(strict "y5-H2O[M+H+Na]^2", spec_positive_78);
mzpaf_test!(strict "y6[M+[2H2]]", spec_positive_79);
mzpaf_test!(strict "y5[M+[15N1]H4]", spec_positive_80);
mzpaf_test!(strict "&1@y7/-0.002", spec_positive_81);
mzpaf_test!(strict "&y7/-0.001", spec_positive_82);
mzpaf_test!(strict "y7/0.000*0.95", spec_positive_83);
mzpaf_test!(strict "&y7/0.001", spec_positive_84);
mzpaf_test!(strict "&y7/0.002", spec_positive_85);
mzpaf_test!(strict "b6-H2O/-0.005,&y7/0.003", spec_positive_86);
mzpaf_test!(strict "y12-H2O^2/7.4ppm*0.70", spec_positive_87);
mzpaf_test!(strict "y12/3.4ppm*0.85,b9-NH3/5.2ppm*0.05", spec_positive_88);

mzpaf_test!(ne r"0@y4{Mar^2dation]ACK}-CH4OS", fuzz_0);
mzpaf_test!(ne r"y4-4", fuzz_1);
mzpaf_test!(ne r"0@y4{Mar^2dation]ACK", fuzz_2);
mzpaf_test!(ne r"0@y4{", fuzz_3);
mzpaf_test!(ne r"IM[", fuzz_4);
mzpaf_test!(ne r"IM[Carboxymethyl", fuzz_5);
mzpaf_test!(ne r"r[]", fuzz_6);
mzpaf_test!(ne r"r[Adensosine", fuzz_7);
mzpaf_test!(ne r"f{}", fuzz_8);
mzpaf_test!(ne r"f{C2H6", fuzz_9);
mzpaf_test!(ne r"m3:4/1.1ppm,m4:0/1.1ppm", fuzz_10);
mzpaf_test!(ne r"m0:4/1.1ppm,m4:5/1.1ppm", fuzz_11);
mzpaf_test!(ne r"0@w7{M[]},x5{M[]},x7{M[]}-H", fuzz_12);
mzpaf_test!(ne r"a0", fuzz_13);
mzpaf_test!(ne r"y5-H2[18O1][M+Na],y5[M+Na],y5[M+[15N1]H[15N1]H4],y12/3.4ppm*0.85,b0-NH3/5.2ppm*0.05", fuzz_14);
mzpaf_test!(ne r"IX[ethyl],p-H95O9^8,IX[ethyl],p-H95O9^8,IX[ethyl],p-H95O9^8,IX[ethyl],c0", fuzz_15);
mzpaf_test!(ne r"da0000000000000000000000000000000000", fuzz_16);
mzpaf_test!(ne r"1@y66{/4743[07kli-,657kli-,6666666666666666660[13N1]Ca-000Ca657kli-,6666666i-,6666666666666666660[13N1]Ca-000Ca657kli-,666666666666666660[13N1]Ca-030Ca-00000Ca-00000Ca-0000Ca[13N1]Ca-0000Ca4m*0.7+,b dation]ACCa-001-,65[13N1]C666666666660[-0Ca-00[13N1]Ca-001Ca-+Ki-]}", fuzz_17);
mzpaf_test!(ne r"1@f{H666666660H666666660H666666660H666666660}", fuzz_18);
mzpaf_test!(ne r"0@IG,0@IP,0@I　[N8𴴴C8nl+[<5NO]H+H2PO3^29<5NO]H+H2PO3^290@_{M0@_{M]b　[", fuzz_19);
mzpaf_test!(ne r"0@y4{A/[]}", fuzz_20);
mzpaf_test!(ne r"0@y4{A/[Na߳]}}", fuzz_21);
mzpaf_test!(ne r"y5[M+H+Na]^2", fuzz_22);

mzpaf_test!("IC[Carbamidomethyl]/-0.0008", hand_test_01);
mzpaf_test!(
    "1@p-[sidechain_Y]-[sidechain_M]^3,1@c26+2H-H2O1^3",
    hand_test_02
);
mzpaf_test!(ne strict "1@p-OH2", hand_test_03);
mzpaf_test!(ne strict "1@p-PH3O4", hand_test_04);
mzpaf_test!(ne strict "1@p-NH3-H2O", hand_test_05);
mzpaf_test!(ne strict "1@p-1H2O", hand_test_06);
mzpaf_test!(ne strict "1@p+1i", hand_test_07);
mzpaf_test!(ne strict "1@p[M+1Na]", hand_test_08);
mzpaf_test!(ne strict "0@y4{AAAAA}", hand_test_09);
mzpaf_test!(ne strict "0@m4:5", hand_test_10);
mzpaf_test!(ne strict "y4*0.2,y5*0.5", hand_test_11);
mzpaf_test!(ne strict "y4*0.2,y5", hand_test_12);
mzpaf_test!(ne strict "y5^0", hand_test_13);
mzpaf_test!(ne strict "y5[M+Na+H]", hand_test_14);
mzpaf_test!(ne "y5[M+na]", hand_test_15);
mzpaf_test!(ne "y5+h", hand_test_16);
mzpaf_test!(ne strict "y5+H2O+H2O", hand_test_17);
mzpaf_test!(ne strict "p^1", hand_test_18);
