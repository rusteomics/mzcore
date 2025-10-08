use std::sync::LazyLock;

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
pub(crate) static GLYCAN_PARSE_LIST: LazyLock<Vec<(String, MonoSaccharide)>> =
    LazyLock::new(|| {
        vec![
            (
                "phosphate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("phosphate".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "phospho".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("phosphate".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "sulfate".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("sulfate".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "sug".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Sugar,
                    substituents: vec![],
                    proforma_name: Some("Sug".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "tri".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Triose,
                    substituents: vec![],
                    proforma_name: Some("Tri".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "tet".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Tetrose(None),
                    substituents: vec![],
                    proforma_name: Some("Tet".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "pent".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "pen".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "rib".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(Some(PentoseIsomer::Ribose)),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "ara".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "xyl".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(Some(PentoseIsomer::Xylose)),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "lyx".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(Some(PentoseIsomer::Lyxose)),
                    substituents: vec![],
                    proforma_name: Some("Pen".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "a-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Acid],
                    proforma_name: Some("a-Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "en,a-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituent::Didehydro,
                    ],
                    proforma_name: Some("en,a-Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "d-hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("d-Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "dhex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("d-Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexnac(s)".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNAc(S)".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "glcnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "galnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "mannac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "allnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "altnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "idonac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Idose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "talnac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexns".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino, GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexNS".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexn".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Amino],
                    proforma_name: Some("HexN".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexs".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Sulfate],
                    proforma_name: Some("HexS".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hexp".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("HexP".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hex".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "glc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "gal".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "man".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "all".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "alt".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "gul".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "ido".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Idose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "tal".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "hep".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Heptose(None),
                    substituents: vec![],
                    proforma_name: Some("Hep".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "oct".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![],
                    proforma_name: Some("Oct".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "kdo".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Octose,
                    substituents: vec![GlycanSubstituent::Acid, GlycanSubstituent::Deoxy],
                    proforma_name: Some("Oct".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "non".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "kdn".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(Some(NonoseIsomer::Kdn)),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "sia".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ],
                    proforma_name: Some("Non".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "dec".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Decose,
                    substituents: vec![],
                    proforma_name: Some("Dec".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "neu5ac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ],
                    proforma_name: Some("NeuAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "neuac".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ],
                    proforma_name: Some("NeuAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "neu5gc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                    ],
                    proforma_name: Some("NeuGc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "neugc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                    ],
                    proforma_name: Some("NeuGc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "neu".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Deoxy,
                    ],
                    proforma_name: Some("Neu".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "fuc".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "xxx".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![],
                    proforma_name: Some("Xxx".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "aldi".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Alcohol],
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "me".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Methyl],
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            // Single letter codes, by defining them like this they will be read but exported to the standard ProForma codes
            (
                "x".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Acetyl],
                    proforma_name: None,
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "p".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::None,
                    substituents: vec![GlycanSubstituent::Phosphate],
                    proforma_name: Some("phosphate".to_string()), // TODO: technically maybe not working when multiple are in there, think it through, should be two different elements,  both getting counts after them
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "h".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![],
                    proforma_name: Some("Hex".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "n".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![GlycanSubstituent::NAcetyl],
                    proforma_name: Some("HexNAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "f".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                    substituents: vec![GlycanSubstituent::Deoxy],
                    proforma_name: Some("Fuc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "s".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ],
                    proforma_name: Some("NeuAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "a".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                    ],
                    proforma_name: Some("NeuAc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
            (
                "g".to_string(),
                MonoSaccharide {
                    base_sugar: BaseSugar::Nonose(None),
                    substituents: vec![
                        GlycanSubstituent::Acid,
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                    ],
                    proforma_name: Some("NeuGc".to_string()),
                    furanose: false,
                    configuration: None,
                },
            ),
        ]
    });
