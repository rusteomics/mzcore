#![allow(clippy::missing_panics_doc)]
use crate::{
    annotation::model::*,
    chemistry::{AmbiguousLabel, MassMode},
    ontology::{CustomDatabase, Ontology},
    sequence::{
        AminoAcid, CompoundPeptidoformIon, LinkerSpecificity, ModificationId, Peptidoform,
        PlacementRule, Position, SimpleModificationInner,
    },
    system::{MassOverCharge, Ratio, isize::Charge, ratio::ppm},
    *,
};

use itertools::Itertools;
use std::sync::Arc;

#[test]
fn triple_a() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (44.049476, "a+1"),
        (72.044390, "b+1"),
        (89.070939, "c+1"),
        (115.086589, "a+2"),
        (143.081504, "b+2"),
        (160.108053, "c+2"),
        (187.071333, "x+2"),
        (161.092069, "y+2"),
        (144.065520, "z+2"),
        (145.073345, "z·+2"),
        (116.034220, "x+1"),
        (90.054955, "y+1"),
        (73.028406, "z+1"),
        (74.036231, "z·+1"),
        (232.129183, "precursor"),
    ];
    let model = FragmentationModel::none()
        .clone()
        .a(PrimaryIonSeries::default())
        .b(PrimaryIonSeries::default())
        .c(PrimaryIonSeries::default())
        .x(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .z(PrimaryIonSeries::default().variants(vec![0, 1]));
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("AAA", None)
            .unwrap()
            .into_linear()
            .unwrap(),
        &model,
        1,
        true,
        false,
    );
}

#[test]
fn with_modifications() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi mods: -17.02655@[ 15.99491@
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (84.044389, "a1"),
        (155.081503, "a2"),
        (226.118617, "a3"),
        (112.039304, "b1"),
        (183.076418, "b2"),
        (254.113532, "b3"),
        (129.065853, "c1"),
        (200.102967, "c2"),
        (271.140081, "c3"),
        (334.106728, "x3"),
        (263.069614, "x2"),
        (192.0325, "x1"),
        (308.127463, "y3"),
        (237.09035, "y2"),
        (166.053236, "y1"),
        (291.100914, "z3"),
        (220.063801, "z2"),
        (149.026687, "z1"),
        (292.108739, "z·3"),
        (221.071626, "z·2"),
        (150.034512, "z·1"),
        (419.159491, "precursor"),
    ];
    let model = FragmentationModel::none()
        .clone()
        .a(PrimaryIonSeries::default())
        .b(PrimaryIonSeries::default())
        .c(PrimaryIonSeries::default())
        .x(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .z(PrimaryIonSeries::default().variants(vec![0, 1]));
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("[Gln->pyro-Glu]-QAAM[Oxidation]", None).unwrap(),
        &model,
        1,
        true,
        false,
    );
}

// #[test]
// fn with_possible_modifications() {
//     // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi mods: 15.99491@1 and separately 15.99491@2
//     #[expect(clippy::unreadable_literal)]
//     let theoretical_fragments = &[
//         (132.047761, "M2-b1"),
//         (148.042671, "M1-b1"),
//         (150.058326, "M1-y1"),
//         (166.053236, "M2-y1"),
//         (297.093720, "precursor"),
//     ];
//     let model = Model::none()
//         .b(PrimaryIonSeries::default())
//         .y(PrimaryIonSeries::default());
//     test(
//         theoretical_fragments,
//         LinearPeptide::pro_forma("M[Oxidation#a1]M[#a1]", None)
//             .unwrap()
//             .into_linear()
//             .expect("Not a linear peptide"),
//         &model,
//         1,
//         true,
//         false,
//     );
// }

// #[test]
// fn with_double_possible_modifications() {
//     // Manual addition over previous example
//     #[expect(clippy::unreadable_literal)]
//     let theoretical_fragments = &[
//         (132.047761, "b1"),
//         (148.042671, "b1+O"),
//         (164.038, "b1+O+O"),
//         (150.058326, "y1"),
//         (166.053236, "y1+O"),
//         (182.048, "y1+O+O"),
//         (313.089, "precursor"),
//     ];
//     let model = Model::none()
//         .b(PrimaryIonSeries::default())
//         .y(PrimaryIonSeries::default());
//     test(
//         theoretical_fragments,
//         LinearPeptide::pro_forma("M[Oxidation#a1][Oxidation#a2]M[#a1][#a2]", None)
//             .unwrap()
//             .into_linear()
//             .unwrap(),
//         &model,
//         1,
//         true,
//         false,
//     );
// }

#[test]
fn higher_charges() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (9.615716, "a1+5"),
        (30.217553, "a2+5"),
        (11.767826, "a1+4"),
        (37.520122, "a2+4"),
        (15.354676, "a1+3"),
        (49.691071, "a2+3"),
        (22.528376, "a1+2"),
        (74.032968, "a2+2"),
        (44.049476, "a1+1"),
        (147.058660, "a2+1"),
        (62.424038, "precursor"),
    ];
    let model = FragmentationModel::none()
        .clone()
        .a(PrimaryIonSeries::default());
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("ACD", None)
            .unwrap()
            .into_linear()
            .unwrap(),
        &model,
        5,
        false,
        false,
    );
}

#[test]
fn all_aminoacids() {
    // Compare rustyms with https://proteomicsresource.washington.edu/cgi-bin/fragment.cgi
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (44.049476, "a1"),
        (200.150587, "a2"),
        (314.193514, "a3"),
        (429.220457, "a4"),
        (532.229642, "a5"),
        (660.288219, "a6"),
        (789.330812, "a7"),
        (846.352276, "a8"),
        (983.411188, "a9"),
        (1096.495252, "a10"),
        (1209.579316, "a11"),
        (1337.674279, "a12"),
        (1468.714764, "a13"),
        (1615.783178, "a14"),
        (1712.835942, "a15"),
        (1799.86797, "a16"),
        (1900.915649, "a17"),
        (2086.994962, "a18"),
        (2250.05829, "a19"),
        (72.04439, "b1"),
        (228.145501, "b2"),
        (342.188429, "b3"),
        (457.215372, "b4"),
        (560.224556, "b5"),
        (688.283134, "b6"),
        (817.325727, "b7"),
        (874.347191, "b8"),
        (1011.406103, "b9"),
        (1124.490167, "b10"),
        (1237.574231, "b11"),
        (1365.669194, "b12"),
        (1496.709678, "b13"),
        (1643.778092, "b14"),
        (1740.830856, "b15"),
        (1827.862885, "b16"),
        (1928.910563, "b17"),
        (2114.989876, "b18"),
        (2278.053205, "b19"),
        (89.070939, "c1"),
        (245.17205, "c2"),
        (359.214978, "c3"),
        (474.241921, "c4"),
        (577.251105, "c5"),
        (705.309683, "c6"),
        (834.352276, "c7"),
        (891.37374, "c8"),
        (1028.432652, "c9"),
        (1141.516716, "c10"),
        (1254.60078, "c11"),
        (1382.695743, "c12"),
        (1513.736227, "c13"),
        (1660.804641, "c14"),
        (1757.857405, "c15"),
        (1844.889434, "c16"),
        (1945.937112, "c17"),
        (2132.016425, "c18"),
        (2295.079754, "c19"),
        (2350.074334, "x19"),
        (2193.973223, "x18"),
        (2079.930296, "x17"),
        (1964.903353, "x16"),
        (1861.894168, "x15"),
        (1733.83559, "x14"),
        (1604.792997, "x13"),
        (1547.771534, "x12"),
        (1410.712622, "x11"),
        (1297.628558, "x10"),
        (1184.544494, "x9"),
        (1056.449531, "x8"),
        (925.409046, "x7"),
        (778.340632, "x6"),
        (681.287868, "x5"),
        (594.25584, "x4"),
        (493.208161, "x3"),
        (307.128848, "x2"),
        (144.06552, "x1"),
        (2324.09507, "y19"),
        (2167.993959, "y18"),
        (2053.951031, "y17"),
        (1938.924088, "y16"),
        (1835.914903, "y15"),
        (1707.856326, "y14"),
        (1578.813733, "y13"),
        (1521.792269, "y12"),
        (1384.733357, "y11"),
        (1271.649293, "y10"),
        (1158.565229, "y9"),
        (1030.470266, "y8"),
        (899.429781, "y7"),
        (752.361367, "y6"),
        (655.308604, "y5"),
        (568.276575, "y4"),
        (467.228897, "y3"),
        (281.149584, "y2"),
        (118.086255, "y1"),
        (2307.06852, "z19"),
        (2150.967409, "z18"),
        (2036.924482, "z17"),
        (1921.897539, "z16"),
        (1818.888354, "z15"),
        (1690.829777, "z14"),
        (1561.787184, "z13"),
        (1504.76572, "z12"),
        (1367.706808, "z11"),
        (1254.622744, "z10"),
        (1141.53868, "z9"),
        (1013.443717, "z8"),
        (882.403232, "z7"),
        (735.334818, "z6"),
        (638.282055, "z5"),
        (551.250026, "z4"),
        (450.202348, "z3"),
        (264.123035, "z2"),
        (101.059706, "z1"),
        (2308.076345, "z·19"),
        (2151.975234, "z·18"),
        (2037.932307, "z·17"),
        (1922.905364, "z·16"),
        (1819.896179, "z·15"),
        (1691.837602, "z·14"),
        (1562.795009, "z·13"),
        (1505.773545, "z·12"),
        (1368.714633, "z·11"),
        (1255.630569, "z·10"),
        (1142.546505, "z·9"),
        (1014.451542, "z·8"),
        (883.411057, "z·7"),
        (736.342643, "z·6"),
        (639.28988, "z·5"),
        (552.257851, "z·4"),
        (451.210173, "z·3"),
        (265.13086, "z·2"),
        (102.067531, "z·1"),
        (2395.132183, "precursor"),
    ];
    let model = FragmentationModel::none()
        .clone()
        .a(PrimaryIonSeries::default())
        .b(PrimaryIonSeries::default())
        .c(PrimaryIonSeries::default())
        .x(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .z(PrimaryIonSeries::default().variants(vec![0, 1]));
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("ARNDCQEGHILKMFPSTWYV", None)
            .unwrap()
            .into_linear()
            .unwrap(),
        &model,
        1,
        false,
        false,
    );
}

#[test]
fn glycan_structure_fragmentation() {
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (4593.06932015166, "pep+N4H5S1"),
        (4301.97390015166, "pep+N4H5"),
        (4139.92108015166, "pep+N4H4"),
        (3977.86826015166, "pep+N4H3"),
        (3936.84171015166, "pep+N3H4"),
        (3774.78889015166, "pep+N3H3"),
        (3612.73607015166, "pep+N3H2"),
        (3571.70952015166, "pep+N2H3"),
        (3409.65670015166, "pep+N2H2"),
        (3247.60388015166, "pep+N2H1"),
        (3085.55106015166, "pep+N2"),
        (2882.47169015166, "pep+N"),
        (2679.39232015166, "pep"),
        (1914.68430008988, "B7:N4H5S1"),
        (1711.60492757466, "B6:N3H5S1"),
        (1508.52555505943, "B5:N2H5S1"),
        (819.28771229837, "B4α:N1H2S1"),
        (657.23488888309, "B3α:N1H1S1"),
        (454.15551636786, "B2α:H1S1"),
        (292.10269295258, "B1α:S1"),
        (528.19229579777, "B3β:N1H2"),
        (366.13947238249, "B2β:N1H1"),
        (163.06009986727, "B1β:H1"),
        (1752.63147667460, "OxY1βEnd1αB7:N4H4S1"),
        (1623.58888358929, "OxEnd1βY1αB7:N4H5"),
        (1549.55210415938, "OxY2βEnd1αB7/Y1βEnd1αB6:N3H4S1"),
        (1461.53606017401, "OxY1βY1αB7/End1βY2αB7:N3H4S1"),
        (1420.50951107406, "OxEnd1βY1αB6:N3H5"),
        (1387.49928074410, "OxY3βEnd1αB7:N3H3S1"),
        (1346.47273164415, "OxY1βEnd1αB5/Y2βEnd1αB6:N2H4S1"),
        (1299.48323675873, "OxY1βY2αB7:N4H3"),
        (
            1258.45668765878,
            "OxY1βY1αB6/Y2βY1αB7/End1βY2αB6/End1βY3αB7:N3H4",
        ),
        (1217.43013855884, "OxEnd1βY1αB5:N2H4"),
        (1184.41990822887, "OxY3βEnd1αB6:N2H2S1"),
        (1143.39335912893, "OxY2βEnd1αB5:N1H3S1"),
        (
            1096.40386424350,
            "OxY1βY2αB6/Y1βY3αB7/Y2βY2αB7/Y3βY1αB7/End1βY4αB7:N3H3",
        ),
        (
            1055.37731514356,
            "OxY1βY1αB5/Y2βY1αB6/End1βY2αB5/End1βY3αB6:N2H4",
        ),
        (981.34053571365, "OxY3βEnd1αB5:N1H3S1"),
        (934.35104082822, "OxY1βY4αB7/Y3βY2αB7:N3H2"),
        (
            893.32449172828,
            "OxY1βY2αB5/Y1βY3αB6/Y2βY2αB6/Y2βY3αB7/Y3βY1αB6/End1βY4αB6:N2H3",
        ),
        (852.29794262833, "OxY2βY1αB5/End1βY3αB5:N1H4"),
        (
            731.27166831300,
            "OxY1βY4αB6/Y2βY4αB7/Y3βY2αB6/Y3βY3αB7:N2H2",
        ),
        (
            690.24511921305,
            "OxY1βY3αB5/Y2βY2αB5/Y2βY3αB6/Y3βY1αB5/End1βY4αB5:N1H3",
        ),
    ];

    let model = FragmentationModel::none()
        .clone()
        .glycan(GlycanModel::DISALLOW.allow_structural(true));
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("MVSHHN[GNO:G43728NL]LTTGATLINEQWLLTTAK", None)
            .unwrap()
            .into_linear()
            .unwrap(),
        &model,
        1,
        true,
        true,
    );
}

#[test]
fn glycan_composition_fragmentation() {
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (4593.06932015166, "pep+N4H5S1"),
        (4301.97390015166, "pep+N4H5"),
        (4139.92108015166, "pep+N4H4"),
        (3977.86826015166, "pep+N4H3"),
        (3936.84171015166, "pep+N3H4"),
        (3774.78889015166, "pep+N3H3"),
        (3612.73607015166, "pep+N3H2"),
        (3571.70952015166, "pep+N2H3"),
        (3409.65670015166, "pep+N2H2"),
        (3247.60388015166, "pep+N2H1"),
        (3085.55106015166, "pep+N2"),
        (2882.47169015166, "pep+N"),
        (2679.39232015166, "pep"),
        (1914.68430008988, "B7:N4H5S1"),
        (1711.60492757466, "B6:N3H5S1"),
        (1508.52555505943, "B5:N2H5S1"),
        (819.28771229837, "B4α:N1H2S1"),
        (657.23488888309, "B3α:N1H1S1"),
        (454.15551636786, "B2α:H1S1"),
        (292.10269295258, "B1α:S1"),
        (528.19229579777, "B3β:N1H2"),
        (366.13947238249, "B2β:N1H1"),
        (163.06009986727, "B1β:H1"),
        (1752.63147667460, "OxY1βEnd1αB7:N4H4S1"),
        (1623.58888358929, "OxEnd1βY1αB7:N4H5"),
        (1549.55210415938, "OxY2βEnd1αB7/Y1βEnd1αB6:N3H4S1"),
        (1461.53606017401, "OxY1βY1αB7/End1βY2αB7:N3H4S1"),
        (1420.50951107406, "OxEnd1βY1αB6:N3H5"),
        (1387.49928074410, "OxY3βEnd1αB7:N3H3S1"),
        (1346.47273164415, "OxY1βEnd1αB5/Y2βEnd1αB6:N2H4S1"),
        (1299.48323675873, "OxY1βY2αB7:N4H3"),
        (
            1258.45668765878,
            "OxY1βY1αB6/Y2βY1αB7/End1βY2αB6/End1βY3αB7:N3H4",
        ),
        (1217.43013855884, "OxEnd1βY1αB5:N2H4"),
        (1184.41990822887, "OxY3βEnd1αB6:N2H2S1"),
        (1143.39335912893, "OxY2βEnd1αB5:N1H3S1"),
        (
            1096.40386424350,
            "OxY1βY2αB6/Y1βY3αB7/Y2βY2αB7/Y3βY1αB7/End1βY4αB7:N3H3",
        ),
        (
            1055.37731514356,
            "OxY1βY1αB5/Y2βY1αB6/End1βY2αB5/End1βY3αB6:N2H4",
        ),
        (981.34053571365, "OxY3βEnd1αB5:N1H3S1"),
        (934.35104082822, "OxY1βY4αB7/Y3βY2αB7:N3H2"),
        (
            893.32449172828,
            "OxY1βY2αB5/Y1βY3αB6/Y2βY2αB6/Y2βY3αB7/Y3βY1αB6/End1βY4αB6:N2H3",
        ),
        (852.29794262833, "OxY2βY1αB5/End1βY3αB5:N1H4"),
        (
            731.27166831300,
            "OxY1βY4αB6/Y2βY4αB7/Y3βY2αB6/Y3βY3αB7:N2H2",
        ),
        (
            690.24511921305,
            "OxY1βY3αB5/Y2βY2αB5/Y2βY3αB6/Y3βY1αB5/End1βY4αB5:N1H3",
        ),
    ];
    let model = FragmentationModel::none()
        .clone()
        .glycan(GlycanModel::DISALLOW.compositional_range_generic(0..=10));
    test(
        theoretical_fragments,
        Peptidoform::pro_forma("MVSHHN[Glycan:N4H5S1]LTTGATLINEQWLLTTAK", None)
            .unwrap()
            .into_linear()
            .unwrap(),
        &model,
        1,
        true,
        false,
    );
}

fn custom_database() -> CustomDatabase {
    vec![
        (
            Some(0),
            "dsso".to_string(),
            Arc::new(SimpleModificationInner::Linker {
                specificities: vec![LinkerSpecificity::Symmetric {
                    rules: vec![PlacementRule::AminoAcid(
                        vec![AminoAcid::Lysine],
                        Position::Anywhere,
                    )],
                    stubs: vec![(
                        molecular_formula!(C 3 O 2 H 1 N -1),
                        molecular_formula!(C 3 O 3 H 1 N -1 S 1),
                    )],
                    neutral_losses: Vec::new(),
                    diagnostic: Vec::new(),
                }],
                formula: molecular_formula!(C 6 O 5 H 2 N -2 S 1),
                id: ModificationId {
                    name: "DSSO".to_string(),
                    id: Some(0),
                    ontology: Ontology::Custom,
                    ..ModificationId::default()
                },
                length: None,
            }),
        ),
        (
            Some(1),
            "disulfide".to_string(),
            Arc::new(SimpleModificationInner::Linker {
                specificities: vec![LinkerSpecificity::Symmetric {
                    rules: vec![PlacementRule::AminoAcid(
                        vec![AminoAcid::Cysteine],
                        Position::Anywhere,
                    )],
                    stubs: vec![(molecular_formula!(H - 1), molecular_formula!(H - 1))],
                    neutral_losses: Vec::new(),
                    diagnostic: Vec::new(),
                }],
                formula: molecular_formula!(C 6 O 5 H 2 N -2 S 1),
                id: ModificationId {
                    name: "Disulfide".to_string(),
                    id: Some(1),
                    ontology: Ontology::Custom,
                    ..ModificationId::default()
                },
                length: None,
            }),
        ),
    ]
}

#[test]
fn intra_link() {
    #[expect(clippy::unreadable_literal)]
    let theoretical_fragments = &[
        (260.196, "y2+"),
        (407.265, "y3+"),
        (439.7200472739809, "precursor/loop#XL1++"),
        (472.17481251740094, "b3+"),
        (619.2432264279929, "b4+"),
        (732.327290402381, "b5+"),
    ];
    let peptide =
        CompoundPeptidoformIon::pro_forma("K[C:DSSO#XL1]GK[#XL1]FLK", Some(&custom_database()))
            .unwrap();
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .allow_cross_link_cleavage(true);
    test(
        theoretical_fragments,
        peptide.clone(),
        &model,
        2,
        true,
        false,
    );
    let fragments = peptide.generate_theoretical_fragments(Charge::new::<system::e>(2), &model);
    let doubly_annotated = dbg!(
        fragments
            .iter()
            .filter(|f| f.formula.as_ref().unwrap().labels().len() > 2)
            .collect_vec()
    );
    assert_eq!(doubly_annotated.len(), 0);
}

#[test]
fn ensure_no_double_xl_labels_breaking() {
    let peptide =
        CompoundPeptidoformIon::pro_forma("EVQLVESGGGLVQPGGSLRLSC[C:Disulfide#XL1]AASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC[#XL1]SRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGC[C:Disulfide#XL2]LVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYIC[#XL2]NVNHKPSNTKVDKKVEPKSC[C:Disulfide#XL3]DKT//DIQMTQSPSSLSASVGDRVTITC[C:Disulfide#XL4]RASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC[#XL4]QQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVC[C:Disulfide#XL5]LLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYAC[#XL5]EVTHQGLSSPVTKSFNRGEC[#XL3]", Some(&custom_database()))
            .unwrap();
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .allow_cross_link_cleavage(true);
    let fragments = peptide.generate_theoretical_fragments(Charge::new::<system::e>(2), &model);
    let doubly_annotated = dbg!(
        fragments
            .iter()
            .filter(|f| f.formula.as_ref().unwrap().labels().len()
                > f.formula
                    .as_ref()
                    .unwrap()
                    .labels()
                    .iter()
                    .map(|l| match l {
                        AmbiguousLabel::CrossLinkBound(n)
                        | AmbiguousLabel::CrossLinkBroken(n, _) => n.to_string(),
                        _ => String::new(),
                    })
                    .unique()
                    .count())
            .collect_vec()
    );
    assert_eq!(doubly_annotated.len(), 0);
}

#[test]
fn ensure_no_double_xl_labels_non_breaking() {
    let peptide =
        CompoundPeptidoformIon::pro_forma("EVQLVESGGGLVQPGGSLRLSC[C:Disulfide#XL1]AASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYC[#XL1]SRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGC[C:Disulfide#XL2]LVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYIC[#XL2]NVNHKPSNTKVDKKVEPKSC[C:Disulfide#XL3]DKT//DIQMTQSPSSLSASVGDRVTITC[C:Disulfide#XL4]RASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYC[#XL4]QQHYTTPPTFGQGTKVEIKRTVAAPSVFIFPPSDEQLKSGTASVVC[C:Disulfide#XL5]LLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYAC[#XL5]EVTHQGLSSPVTKSFNRGEC[#XL3]", Some(&custom_database()))
            .unwrap();
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .allow_cross_link_cleavage(false);
    let fragments = peptide.generate_theoretical_fragments(Charge::new::<system::e>(2), &model);
    let doubly_annotated = dbg!(
        fragments
            .iter()
            .filter(|f| f.formula.as_ref().unwrap().labels().len()
                > f.formula
                    .as_ref()
                    .unwrap()
                    .labels()
                    .iter()
                    .map(|l| match l {
                        AmbiguousLabel::CrossLinkBound(n)
                        | AmbiguousLabel::CrossLinkBroken(n, _) => n.to_string(),
                        _ => String::new(),
                    })
                    .unique()
                    .count())
            .collect_vec()
    );
    assert_eq!(doubly_annotated.len(), 0);
}

#[test]
fn ensure_no_double_xl_labels_small_breaking() {
    let peptide = CompoundPeptidoformIon::pro_forma(
        "EC[C:Disulfide#XL1]AC[#XL1]SC[C:Disulfide#XL3]D//DC[#XL3]",
        Some(&custom_database()),
    )
    .unwrap();
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .allow_cross_link_cleavage(true);
    let fragments = peptide.generate_theoretical_fragments(Charge::new::<system::e>(2), &model);
    let doubly_annotated = dbg!(
        fragments
            .iter()
            .filter(|f| f.formula.as_ref().unwrap().labels().len() > 2)
            .collect_vec()
    );
    assert_eq!(doubly_annotated.len(), 0);
}

#[test]
fn ensure_no_double_xl_labels_small_non_breaking() {
    let peptide = CompoundPeptidoformIon::pro_forma(
        "EC[C:Disulfide#XL1]AC[#XL1]SC[C:Disulfide#XL3]D//DC[#XL3]",
        Some(&custom_database()),
    )
    .unwrap();
    let model = FragmentationModel::none()
        .clone()
        .b(PrimaryIonSeries::default())
        .y(PrimaryIonSeries::default())
        .allow_cross_link_cleavage(false);
    let fragments = peptide.generate_theoretical_fragments(Charge::new::<system::e>(2), &model);
    let doubly_annotated = dbg!(
        fragments
            .iter()
            .filter(|f| f.formula.as_ref().unwrap().labels().len() > 2)
            .collect_vec()
    );
    assert_eq!(doubly_annotated.len(), 0);
}

fn test(
    theoretical_fragments: &[(f64, &str)],
    peptide: impl Into<CompoundPeptidoformIon>,
    model: &FragmentationModel,
    charge: isize,
    allow_left_over_generated: bool,
    allow_double_theoretical: bool,
) {
    let peptide = peptide.into();
    let mut calculated_fragments =
        peptide.generate_theoretical_fragments(Charge::new::<system::e>(charge), model);
    let mut found = Vec::new();
    let mut this_found;
    for goal in theoretical_fragments {
        this_found = false;
        let mut index = 0;
        while index < calculated_fragments.len() {
            if calculated_fragments[index]
                .mz(MassMode::Monoisotopic)
                .map(|v| v.ppm(MassOverCharge::new::<system::mz>(goal.0)))
                .is_some_and(|v| v < Ratio::new::<ppm>(20.0))
            {
                println!(
                    "Match: {}@{} with {} (labels: {})",
                    goal.1,
                    goal.0,
                    calculated_fragments[index],
                    calculated_fragments[index]
                        .formula
                        .as_ref()
                        .unwrap()
                        .labels()
                        .iter()
                        .join(",")
                );
                calculated_fragments.remove(index);
                index = index.saturating_sub(1);
                this_found = true;
                if !allow_double_theoretical {
                    break;
                }
            }
            index += 1;
        }
        found.push(this_found);
    }

    for left in calculated_fragments.iter().sorted_by(|a, b| b.cmp(a)) {
        println!(
            "Excess fragments: {left} (labels: {})",
            left.formula.as_ref().unwrap().labels().iter().join(",")
        );
    }
    let not_found: Vec<usize> = found
        .iter()
        .enumerate()
        .filter_map(|(index, is_found)| (!is_found).then_some(index))
        .collect();

    for index in &not_found {
        let (mass, name) = theoretical_fragments[*index];
        println!("Not found: {mass} {name}");
    }
    assert_eq!(not_found.len(), 0, "Not all needed fragments are found");
    if !allow_left_over_generated {
        assert_eq!(
            calculated_fragments.len(),
            0,
            "Not all generated fragments are accounted for"
        );
    }
}
