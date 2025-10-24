use crate::fragment::mzpaf::write::ToMzPAF;

/// Create a parse test based on a given case and its name.
#[macro_export]
macro_rules! mzpaf_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            use itertools::Itertools;
            let basic_analytes = [
                (
                    1,
                    $crate::mzspeclib::AnalyteTarget::PeptidoformIon(
                        mzcore::sequence::PeptidoformIon::pro_forma(
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            None,
                        )
                        .unwrap(),
                    ),
                ),
                (
                    2,
                    $crate::mzspeclib::AnalyteTarget::PeptidoformIon(
                        mzcore::sequence::PeptidoformIon::pro_forma(
                            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
                            None,
                        )
                        .unwrap(),
                    ),
                ),
            ];
            let res = $crate::fragment::Fragment::mz_paf($case, None, &basic_analytes);
            match res {
                Err(err) => {
                    println!("Failed: '{}'", $case);
                    println!("{err}");
                    panic!("Failed test")
                }
                Ok(res) => {
                    let back = res.iter().map(|a| a.to_mz_paf_string()).join(",");
                    let res_back = $crate::fragment::Fragment::mz_paf(&back, None, &basic_analytes);
                    match res_back {
                        Ok(res_back) => {
                            let back_back = res_back.iter().map(|a| a.to_mz_paf_string()).join(",");
                            assert_eq!(
                                back, back_back,
                                "{back} != {back_back} (from input: {})",
                                $case
                            )
                        }
                        Err(err) => {
                            println!("Failed: '{}' was exported as '{back}'", $case);
                            println!("{err}");
                            panic!("Failed test")
                        }
                    }
                }
            };
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::fragment::Fragment::mz_paf($case, None, &[]);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

mzpaf_test!("b2-H2O/3.2ppm,b4-H2O^2/3.2ppm", spec_positive_1);
mzpaf_test!("b2-H2O/3.2ppm*0.75,b4-H2O^2/3.2ppm*0.25", spec_positive_2);
mzpaf_test!("1@y12/0.13,2@b9-NH3/0.23", spec_positive_3);
mzpaf_test!("0@y1{K}", spec_positive_4);
mzpaf_test!("0@y1{K}-NH3", spec_positive_5);
mzpaf_test!("y1/-1.4ppm", spec_positive_6);
mzpaf_test!("y1/-0.0002", spec_positive_7);
mzpaf_test!("y4-H2O+2i[M+H+Na]^2", spec_positive_8);
mzpaf_test!("?", spec_positive_9);
mzpaf_test!("?^3", spec_positive_10);
mzpaf_test!("?+2i^4", spec_positive_11);
mzpaf_test!("?17", spec_positive_12);
mzpaf_test!("?17+i/1.45ppm", spec_positive_13);
mzpaf_test!("?17-H2O/-0.87ppm", spec_positive_14);
mzpaf_test!("0@b2{LL}", spec_positive_15);
mzpaf_test!("0@y1{K}", spec_positive_16);
mzpaf_test!("0@b2{LC[Carbamidomethyl]}", spec_positive_17);
mzpaf_test!("0@b1{[Acetyl]-M}", spec_positive_18);
mzpaf_test!("0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19a);
mzpaf_test!("0@y44{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19b);
mzpaf_test!("0@y444{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19c);
mzpaf_test!("m3:6", spec_positive_20);
mzpaf_test!("b3-C2H3NO", spec_positive_21);
mzpaf_test!("m3:6-CO", spec_positive_22);
mzpaf_test!("m3:6-CO-H2O^2", spec_positive_23);
mzpaf_test!("m3:4/1.1ppm,m4:5/1.1ppm", spec_positive24);
mzpaf_test!("m3:5", spec_positive_25);
mzpaf_test!("IY", spec_positive_26);
mzpaf_test!("IH", spec_positive_27);
mzpaf_test!("IL-CH2", spec_positive_28);
mzpaf_test!("IC[Carbamidomethyl]", spec_positive_29);
mzpaf_test!("IY[Phospho]", spec_positive_30);
mzpaf_test!("IC[+58.005]", spec_positive_31);
mzpaf_test!("p^2", spec_positive_32a);
mzpaf_test!("p^-2", spec_positive_32b);
mzpaf_test!("p-H3PO4^2", spec_positive_33);
mzpaf_test!("p^4", spec_positive_34);
mzpaf_test!("p+H^3", spec_positive_35);
mzpaf_test!("p^3", spec_positive_36);
mzpaf_test!("p+2H^2", spec_positive_37);
mzpaf_test!("p^2", spec_positive_38);
mzpaf_test!("p+H^2", spec_positive_39);
mzpaf_test!("p+3H", spec_positive_40);
mzpaf_test!("p+2H", spec_positive_41);
mzpaf_test!("p+H", spec_positive_42);
mzpaf_test!("p", spec_positive_43);
mzpaf_test!("r[TMT127N]", spec_positive_44);
mzpaf_test!("r[iTRAQ114]", spec_positive_45);
mzpaf_test!("r[TMT6plex]", spec_positive_46);
mzpaf_test!("r[Hex]", spec_positive_47);
mzpaf_test!("r[Adenine]", spec_positive_48);
mzpaf_test!("0@_{Urocanic Acid}", spec_positive_49);
mzpaf_test!("f{C13H9}/-0.55ppm", spec_positive_50);
mzpaf_test!("f{C12H9N}/0.06ppm", spec_positive_51);
mzpaf_test!("f{C13H9N}/-2.01ppm", spec_positive_52);
mzpaf_test!("f{C13H10N}/-0.11ppm", spec_positive_53);
mzpaf_test!("f{C13H11N}/-0.09ppm", spec_positive_54);
mzpaf_test!("f{C13H12N}/0.26ppm", spec_positive_55);
mzpaf_test!("f{C14H10N}/0.19ppm", spec_positive_56);
mzpaf_test!("f{C14H11N}/0.45ppm", spec_positive_57);
mzpaf_test!("f{C14H10NO}/0.03ppm", spec_positive_58);
mzpaf_test!("f{C16H22O}+i^3", spec_positive_59);
mzpaf_test!("f{C15[13C1]H22O}^3", spec_positive_60);
// mzpaf_test!("s{CN=C=O}[M+H]/-0.55ppm", spec_positive_61); TODO: SMILES not supported
// mzpaf_test!("s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm", spec_positive_62); TODO: SMILES not supported
mzpaf_test!("p-[Hex]", spec_positive_63);
mzpaf_test!("y2+CO-H2O", spec_positive_64);
mzpaf_test!("y2-H2O-NH3", spec_positive_65);
mzpaf_test!("p-[TMT6plex]-2H2O-HPO3", spec_positive_66);
mzpaf_test!("p-2[iTRAQ115]", spec_positive_67);
mzpaf_test!("p-[iTRAQ116]-CO-H2O-HPO3", spec_positive_68);
mzpaf_test!("y2-[2H1]-NH3", spec_positive_69);
mzpaf_test!("y5-H2[18O1][M+Na]", spec_positive_70);
mzpaf_test!("y12+i", spec_positive_71);
mzpaf_test!("y12+i13C", spec_positive_72);
mzpaf_test!("y12+i15N", spec_positive_73);
mzpaf_test!("y12+2i13C+i15N", spec_positive_74);
mzpaf_test!("y12+iA", spec_positive_75);
mzpaf_test!("y12+2iA", spec_positive_76);
mzpaf_test!("y4[M+Na]", spec_positive_77);
mzpaf_test!("y5-H2O[M+H+Na]^2", spec_positive_78);
mzpaf_test!("y6[M+[2H2]]", spec_positive_79);
mzpaf_test!("y5[M+[15N1]H4]", spec_positive_80);
mzpaf_test!("&1@y7/-0.002", spec_positive_81);
mzpaf_test!("&y7/-0.001", spec_positive_82);
mzpaf_test!("y7/0.000*0.95", spec_positive_83);
mzpaf_test!("&y7/0.001", spec_positive_84);
mzpaf_test!("&y7/0.002", spec_positive_85);
mzpaf_test!("b6-H2O/-0.005,&y7/0.003", spec_positive_86);
mzpaf_test!("y12-H2O^2/7.4ppm*0.70", spec_positive_87);
mzpaf_test!("y12/3.4ppm*0.85,b9-NH3/5.2ppm*0.05", spec_positive_88);

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

mzpaf_test!("IC[Carbamidomethyl]/-0.0008", hand_test_01);
