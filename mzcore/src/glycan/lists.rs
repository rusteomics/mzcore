use std::{num::NonZeroU8, sync::LazyLock};

use thin_vec::ThinVec;

use super::{
    MonoSaccharide,
    glycan::{
        BaseSugar, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, NonoseIsomer, PentoseIsomer,
        TetroseIsomer,
    },
};
use crate::glycan::{GlycanSubstituentLeavingGroup, GlycanSubstituentLocation};

pub(crate) const BASE_SUGARS: &[(
    &str,
    BaseSugar,
    &[(GlycanSubstituent, GlycanSubstituentLocation)],
)] = &[
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
    ("neu", BaseSugar::Nonose(None), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("sia", BaseSugar::Nonose(None), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("kdn", BaseSugar::Nonose(Some(NonoseIsomer::Kdn)), &[
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("kdo", BaseSugar::Octose, &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("fuc", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("rha", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[(
        GlycanSubstituent::Deoxy,
        GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
    )]),
    ("qui", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[(
        GlycanSubstituent::Deoxy,
        GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
    )]),
    ("oli", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("tyv", BaseSugar::Hexose(Some(HexoseIsomer::Mannose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("asc", BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("abe", BaseSugar::Hexose(Some(HexoseIsomer::Gulose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("par", BaseSugar::Hexose(Some(HexoseIsomer::Altrose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("dig", BaseSugar::Hexose(Some(HexoseIsomer::Allose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("col", BaseSugar::Hexose(Some(HexoseIsomer::Talose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("psi", BaseSugar::Hexose(Some(HexoseIsomer::Psicose)), &[]),
    ("fru", BaseSugar::Hexose(Some(HexoseIsomer::Fructose)), &[]),
    ("sor", BaseSugar::Hexose(Some(HexoseIsomer::Sorbose)), &[]),
    ("tag", BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[]),
    ("vionac", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(4),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(6),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
        (
            GlycanSubstituent::NAcetyl,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(4),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
    ]),
    ("vio", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(4),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(4),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(
                NonZeroU8::new(6),
                GlycanSubstituentLeavingGroup::Hydroxy,
            ),
        ),
    ]),
    ("xul", BaseSugar::Pentose(Some(PentoseIsomer::Xylulose)), &[
    ]),
    (
        "sed",
        BaseSugar::Heptose(Some(HeptoseIsomer::Sedoheptulose)),
        &[],
    ),
    ("murnac", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::NAcetyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::OCarboxyEthyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("murngc", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::NGlycolyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::OCarboxyEthyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("mur", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::OCarboxyEthyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    (
        "api",
        BaseSugar::Tetrose(Some(TetroseIsomer::Erythrose)),
        &[(
            GlycanSubstituent::HydroxyMethyl,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        )],
    ),
    ("dha", BaseSugar::Heptose(None), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("bac", BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("pse", BaseSugar::Nonose(Some(NonoseIsomer::Pse)), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("leg", BaseSugar::Nonose(Some(NonoseIsomer::Leg)), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
    ("aci", BaseSugar::Nonose(Some(NonoseIsomer::Aci)), &[
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
        (
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
        ),
    ]),
];

pub(crate) const POSTFIX_SUBSTITUENTS: &[(
    &str,
    (GlycanSubstituent, GlycanSubstituentLeavingGroup),
)] = &[
    (
        "ac",
        (
            GlycanSubstituent::Acetyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "ala",
        (
            GlycanSubstituent::Alanyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "am",
        (
            GlycanSubstituent::Acetimidoyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "en",
        (
            GlycanSubstituent::Didehydro,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "fo",
        (
            GlycanSubstituent::Formyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "gc",
        (
            GlycanSubstituent::Glycolyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "gly",
        (
            GlycanSubstituent::Glycyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "gr",
        (
            GlycanSubstituent::Glyceryl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "4hb",
        (
            GlycanSubstituent::HydroxyButyryl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "3rhb",
        (
            GlycanSubstituent::HydroxyButyryl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "3shb",
        (
            GlycanSubstituent::HydroxyButyryl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "lt",
        (
            GlycanSubstituent::Lactyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "lac",
        (
            GlycanSubstituent::Lac,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "me",
        (
            GlycanSubstituent::Methyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "cme",
        (
            GlycanSubstituent::Methyl,
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ), // unsure about the difference with Me
    (
        "nac",
        (
            GlycanSubstituent::NAcetyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "pyr",
        (
            GlycanSubstituent::CargoxyEthylidene,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "tau",
        (
            GlycanSubstituent::Tauryl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "onic",
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "uronic",
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "aric",
        (
            GlycanSubstituent::Aric,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "ol",
        (
            GlycanSubstituent::Alcohol,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "etn",
        (
            GlycanSubstituent::Ethanolamine,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "etoh",
        (
            GlycanSubstituent::EtOH,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "ulof",
        (
            GlycanSubstituent::Ulof,
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ),
    (
        "ulo",
        (
            GlycanSubstituent::Ulo,
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ),
    (
        "n2dime",
        (
            GlycanSubstituent::NDiMe,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "ndime",
        (
            GlycanSubstituent::NDiMe,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "pcho",
        (
            GlycanSubstituent::PCholine,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "ce",
        (
            GlycanSubstituent::Glycyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ), // Same molecular formula
    (
        "suc",
        (
            GlycanSubstituent::Suc,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "nfo",
        (
            GlycanSubstituent::NFo,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "cfo",
        (
            GlycanSubstituent::Formyl,
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ),
    (
        "dime",
        (
            GlycanSubstituent::DiMethyl,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "cf",
        (
            GlycanSubstituent::Element(crate::chemistry::Element::F),
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ),
    (
        "cn",
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLeavingGroup::Hydrogen,
        ),
    ),
    (
        "sh",
        (
            GlycanSubstituent::Element(crate::chemistry::Element::S), // TODO: it has one H bound
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "a",
        (
            GlycanSubstituent::Acid,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "p",
        (
            GlycanSubstituent::Phosphate,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "s",
        (
            GlycanSubstituent::Sulfate,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
    (
        "n",
        (
            GlycanSubstituent::Amino,
            GlycanSubstituentLeavingGroup::Hydroxy,
        ),
    ),
];

pub(crate) const DOUBLE_LINKED_POSTFIX_SUBSTITUENTS: &[(&str, &[GlycanSubstituent])] = &[
    ("py", &[GlycanSubstituent::Pyruvyl]),
    ("n", &[GlycanSubstituent::Water]),
    ("p", &[
        GlycanSubstituent::Phosphate,
        GlycanSubstituent::Water,
    ]),
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
            (vec!["Phosphate", "Phospho"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Sulfate"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Sug"], MonoSaccharide {
                base_sugar: BaseSugar::Sugar,
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Tri"], MonoSaccharide {
                base_sugar: BaseSugar::Triose,
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Tet"], MonoSaccharide {
                base_sugar: BaseSugar::Tetrose(None),
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (
                vec!["Pent", "Pen", "Rib", "Ara", "Xyl", "Lyx"],
                MonoSaccharide {
                    base_sugar: BaseSugar::Pentose(None),
                    substituents: ThinVec::new(),
                    furanose: false,
                    configuration: None,
                    #[cfg(debug_assertions)]
                    source_field: String::new(),
                },
            ),
            (
                vec![
                    "a-Hex", "aHex", "HexA", "GlcA", "GalA", "ManA", "AllA", "AltA", "GulA",
                    "IdoA", "TalA",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![(
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    )]
                    .into(),
                    furanose: false,
                    configuration: None,
                    #[cfg(debug_assertions)]
                    source_field: String::new(),
                },
            ),
            (vec!["en,a-Hex", "en,aHex"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Didehydro,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["d-Hex", "dHex"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![(
                    GlycanSubstituent::Deoxy,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["HexNAc(S)", "HexNAcS"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::NAcetyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Sulfate,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (
                vec![
                    "HexNAc", "GlcNAc", "GalNAc", "ManNAc", "AllNAc", "AltNAc", "GulNAc", "IdoNAc",
                    "TalNAc",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: vec![(
                        GlycanSubstituent::NAcetyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    )]
                    .into(),
                    furanose: false,
                    configuration: None,
                    #[cfg(debug_assertions)]
                    source_field: String::new(),
                },
            ),
            (vec!["VioNAc"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                substituents: vec![
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(4),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(6),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::NAcetyl,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(4),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Vio"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                substituents: vec![
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(4),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(4),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            NonZeroU8::new(6),
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["HexNS"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Sulfate,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["HexN"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![(
                    GlycanSubstituent::Amino,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["HexS"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![(
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["HexP"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![(
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (
                vec![
                    "Hex", "Glc", "Gal", "Man", "All", "Alt", "Gul", "Ido", "Tal",
                ],
                MonoSaccharide {
                    base_sugar: BaseSugar::Hexose(None),
                    substituents: ThinVec::new(),
                    furanose: false,
                    configuration: None,
                    #[cfg(debug_assertions)]
                    source_field: String::new(),
                },
            ),
            (vec!["Hep"], MonoSaccharide {
                base_sugar: BaseSugar::Heptose(None),
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Oct"], MonoSaccharide {
                base_sugar: BaseSugar::Octose,
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Kdo"], MonoSaccharide {
                base_sugar: BaseSugar::Octose,
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Non"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Kdn"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(Some(NonoseIsomer::Kdn)),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Sia"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Dec"], MonoSaccharide {
                base_sugar: BaseSugar::Decose,
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Neu5Ac", "NeuAc"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Neu5Gc", "NeuGc"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Neu"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Deoxy,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Fuc"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                substituents: vec![(
                    GlycanSubstituent::Deoxy,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Xxx"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Aldi"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Alcohol,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["Me"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Methyl,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            // Single letter codes, by defining them like this they will be read but exported to
            // the standard ProForma codes
            (vec!["P"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["H"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: ThinVec::new(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["N"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(None),
                substituents: vec![(
                    GlycanSubstituent::NAcetyl,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["F"], MonoSaccharide {
                base_sugar: BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                substituents: vec![(
                    GlycanSubstituent::Deoxy,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["S", "A"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["G"], MonoSaccharide {
                base_sugar: BaseSugar::Nonose(None),
                substituents: vec![
                    (
                        GlycanSubstituent::Acid,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Amino,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                    (
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    ),
                ]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
            (vec!["X"], MonoSaccharide {
                base_sugar: BaseSugar::Custom(Box::default()),
                substituents: vec![(
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy),
                )]
                .into(),
                furanose: false,
                configuration: None,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            }),
        ]
    });
