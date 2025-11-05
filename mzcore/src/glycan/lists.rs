use std::sync::LazyLock;

use thin_vec::ThinVec;

use crate::prelude::MolecularFormula;

use super::glycan::{
    BaseSugar, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, MonoSaccharide, NonoseIsomer,
    PentoseIsomer, TetroseIsomer,
};

pub(crate) const BASE_SUGARS: &[(&str, BaseSugar, &[GlycanSubstituent])] = &[
    ("sug", BaseSugar::Sugar, &[]),
    ("tri", BaseSugar::Triose, &[]),
    ("tet", BaseSugar::Tetrose(None), &[]),
    (
        "ery",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[],
    ),
    ("tho", BaseSugar::Tetrose(Some(TetroseIsomer::Threose)), &[]),
    ("pen", BaseSugar::Pentose(None), &[]),
    ("rib", BaseSugar::Pentose(Some(PentoseIsomer::Ribose)), &[]),
    (
        "ara",
        BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)),
        &[],
    ),
    ("xyl", BaseSugar::Pentose(Some(PentoseIsomer::Xylose)), &[]),
    ("lyx", BaseSugar::Pentose(Some(PentoseIsomer::Lyxose)), &[]),
    ("hex", BaseSugar::Hexose(None), &[]),
    ("glc", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[]),
    ("gal", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[]),
    ("man", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[]),
    ("all", BaseSugar::Hexose(Some(HexoseIsomer::Allose)), &[]),
    ("alt", BaseSugar::Hexose(Some(HexoseIsomer::Altrose)), &[]),
    ("gul", BaseSugar::Hexose(Some(HexoseIsomer::Gulose)), &[]),
    ("ido", BaseSugar::Hexose(Some(HexoseIsomer::Idose)), &[]),
    ("tal", BaseSugar::Hexose(Some(HexoseIsomer::Talose)), &[]),
    ("hep", BaseSugar::Heptose(None), &[]),
    (
        "gro-manhep",
        BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)),
        &[],
    ),
    (
        "neu",
        BaseSugar::Nonose(None),
        &[GlycanSubstituent::Acid, GlycanSubstituent::Amino],
    ),
    (
        "sia",
        BaseSugar::Nonose(None),
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "kdn",
        BaseSugar::Nonose(Some(NonoseIsomer::Kdn)),
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Acid,
        ],
    ),
    (
        "kdo",
        BaseSugar::Octose,
        &[GlycanSubstituent::Acid, GlycanSubstituent::Deoxy],
    ),
    (
        "fuc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "rha",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "qui",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy],
    ),
    (
        "oli",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "tyv",
        BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "asc",
        BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "abe",
        BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "par",
        BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "dig",
        BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    (
        "col",
        BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
        &[GlycanSubstituent::Deoxy, GlycanSubstituent::Deoxy],
    ),
    ("psi", BaseSugar::Hexose(Some(HexoseIsomer::Psicose)), &[]),
    ("fru", BaseSugar::Hexose(Some(HexoseIsomer::Fructose)), &[]),
    ("sor", BaseSugar::Hexose(Some(HexoseIsomer::Sorbose)), &[]),
    ("tag", BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[]),
    (
        "xul",
        BaseSugar::Pentose(Some(PentoseIsomer::Xylulose)),
        &[],
    ),
    (
        "sed",
        BaseSugar::Heptose(Some(HeptoseIsomer::Sedoheptulose)),
        &[],
    ),
    (
        "murnac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::NAcetyl, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "murngc",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::NGlycolyl,
            GlycanSubstituent::OCarboxyEthyl,
        ],
    ),
    (
        "mur",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[GlycanSubstituent::Amino, GlycanSubstituent::OCarboxyEthyl],
    ),
    (
        "api",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[GlycanSubstituent::HydroxyMethyl],
    ),
    (
        "dha",
        BaseSugar::Heptose(None),
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Acid,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "bac",
        BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
        &[
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "pse",
        BaseSugar::Nonose(Some(NonoseIsomer::Pse)),
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "leg",
        BaseSugar::Nonose(Some(NonoseIsomer::Leg)),
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
    (
        "aci",
        BaseSugar::Nonose(Some(NonoseIsomer::Aci)),
        &[
            GlycanSubstituent::Acid,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Amino,
            GlycanSubstituent::Deoxy,
            GlycanSubstituent::Deoxy,
        ],
    ),
];

// TODO: Points from mobiusklein in rusteomics/mzcore/pull/2
// * Add an additional level which defines the leaving group, to make the chemical formula difference easier
pub(crate) const POSTFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("ac", GlycanSubstituent::Acetyl),
    ("ala", GlycanSubstituent::Alanyl),
    ("am", GlycanSubstituent::Acetimidoyl),
    ("en", GlycanSubstituent::Didehydro),
    ("fo", GlycanSubstituent::Formyl),
    ("gc", GlycanSubstituent::Glycolyl),
    ("gly", GlycanSubstituent::Glycyl),
    ("gr", GlycanSubstituent::Glyceryl),
    ("4hb", GlycanSubstituent::HydroxyButyryl),
    ("3rhb", GlycanSubstituent::HydroxyButyryl),
    ("3shb", GlycanSubstituent::HydroxyButyryl),
    ("lt", GlycanSubstituent::Lactyl),
    ("lac", GlycanSubstituent::Lac),
    ("me", GlycanSubstituent::Methyl),
    ("cme", GlycanSubstituent::Methyl), // unsure about the difference with Me
    ("nac", GlycanSubstituent::NAcetyl),
    ("pyr", GlycanSubstituent::CargoxyEthylidene),
    ("tau", GlycanSubstituent::Tauryl),
    ("onic", GlycanSubstituent::Acid),
    ("uronic", GlycanSubstituent::Acid),
    ("aric", GlycanSubstituent::Aric),
    ("ol", GlycanSubstituent::Alcohol),
    ("etn", GlycanSubstituent::Ethanolamine),
    ("etoh", GlycanSubstituent::Ethanolamine),
    ("ulof", GlycanSubstituent::Ulof),
    ("ulo", GlycanSubstituent::Ulo),
    ("n2dime", GlycanSubstituent::NDiMe),
    ("ndime", GlycanSubstituent::NDiMe),
    ("pcho", GlycanSubstituent::PCholine),
    ("ce", GlycanSubstituent::Glycyl), // Same molecular formula
    ("suc", GlycanSubstituent::Suc),
    ("nfo", GlycanSubstituent::NFo),
    ("dime", GlycanSubstituent::DiMethyl),
    ("a", GlycanSubstituent::Acid),
    ("p", GlycanSubstituent::Phosphate),
    ("s", GlycanSubstituent::Sulfate),
    ("n", GlycanSubstituent::Amino),
];

pub(crate) const DOUBLE_LINKED_POSTFIX_SUBSTITUENTS: &[(&str, &[GlycanSubstituent])] = &[
    ("py", &[GlycanSubstituent::Pyruvyl]),
    ("n", &[GlycanSubstituent::Water]),
    (
        "p",
        &[GlycanSubstituent::Phosphate, GlycanSubstituent::Water],
    ),
];

pub(crate) const PREFIX_SUBSTITUENTS: &[(&str, GlycanSubstituent)] = &[
    ("deoxy", GlycanSubstituent::Deoxy),
    ("anhydro", GlycanSubstituent::Deoxy),
    ("d", GlycanSubstituent::Deoxy),
];

/// All monosaccharides ordered to be able to parse glycans by matching them from the top
pub(crate) static GLYCAN_PARSE_LIST: LazyLock<Vec<(Vec<&'static str>, MonoSaccharide)>> =
    LazyLock::new(|| {
        vec![
            (
                vec!["Phosphate", "Phospho"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Phosphate].into(),
                    proforma_name: Some("Phosphate".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Sulfate"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Sulfate].into(),
                    proforma_name: Some("Sulfate".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Sug"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Sugar,
                    substituents: ThinVec::new(),
                    proforma_name: Some("Sug".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Tri"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Triose,
                    substituents: ThinVec::new(),
                    proforma_name: Some("Tri".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Tet"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Tetrose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Tet".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Pent", "Pen", "Rib", "Ara", "Xyl", "Lyx"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Pen".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec![
                    "a-Hex", "aHex", "HexA", "GlcA", "GalA", "ManA", "AllA", "AltA", "GulA",
                    "IdoA", "TalA",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Acid].into(),
                    proforma_name: Some("aHex".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["en,a-Hex", "en,aHex"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Didehydro,
                    ]
                    .into(),
                    proforma_name: Some("en,aHex".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["d-Hex", "dHex"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("dHex".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["HexNAc(S)", "HexNAcS"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate]
                        .into(),
                    proforma_name: Some("HexNAcS".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec![
                    "HexNAc", "GlcNAc", "GalNAc", "ManNAc", "AllNAc", "AltNAc", "GulNAc", "IdoNAc",
                    "TalNAc",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl].into(),
                    proforma_name: Some("HexNAc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["HexNS"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino, GlycanSubstituent::Sulfate].into(),
                    proforma_name: Some("HexNS".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["HexN"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino].into(),
                    proforma_name: Some("HexN".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["HexS"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Sulfate].into(),
                    proforma_name: Some("HexS".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["HexP"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Phosphate].into(),
                    proforma_name: Some("HexP".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec![
                    "Hex", "Glc", "Gal", "Man", "All", "Alt", "Gul", "Ido", "Tal",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Hex".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Hep"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Heptose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Hep".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Oct"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: ThinVec::new(),
                    proforma_name: Some("Oct".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Kdo"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![GlycanSubstituent::Acid, GlycanSubstituent::Deoxy].into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Non"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Non".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Kdn"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(Some(NonoseIsomer::Kdn)),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ]
                    .into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Sia"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ]
                    .into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Dec"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Decose,
                    substituents: ThinVec::new(),
                    proforma_name: Some("Dec".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Neu5Ac", "NeuAc"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ]
                    .into(),
                    proforma_name: Some("NeuAc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Neu5Gc", "NeuGc"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                    ]
                    .into(),
                    proforma_name: Some("NeuGc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Neu"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ]
                    .into(),
                    proforma_name: Some("Neu".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Fuc"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy].into(),
                    proforma_name: Some("Fuc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Xxx"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: ThinVec::new(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Aldi"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Alcohol].into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["Me"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Methyl].into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
            (
                vec!["P"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Phosphate].into(),
                    proforma_name: Some("Phosphate".to_string().into_boxed_str()), // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements, both getting counts after them
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["H"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: ThinVec::new(),
                    proforma_name: Some("Hex".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["N"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl].into(),
                    proforma_name: Some("HexNAc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["F"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy].into(),
                    proforma_name: Some("Fuc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["S", "A"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ]
                    .into(),
                    proforma_name: Some("NeuAc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["G"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                    ]
                    .into(),
                    proforma_name: Some("NeuGc".to_string().into_boxed_str()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                vec!["X"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Custom(MolecularFormula::default()),
                    substituents: vec![GlycanSubstituent::Acetyl].into(),
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
        ]
    });
