use std::sync::LazyLock;

use crate::{
    annotation::model::{
        ChargePoint, ChargeRange, FragmentationModel, GlycanModel, Location, PrimaryIonSeries,
        SatelliteIonSeries, SatelliteLocation,
    },
    fragment::{FragmentKind, NeutralLoss},
    glycan::{BaseSugar, GlycanSubstituent, MonoSaccharide},
    molecular_formula,
    sequence::AminoAcid,
};

use super::GlycanPeptideFragment;

static MODEL_ALL: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (1, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, IMMONIUM_LOSSES.clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    allow_cross_link_cleavage: true,
});

static MODEL_NONE: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::none(),
    b: PrimaryIonSeries::none(),
    c: PrimaryIonSeries::none(),
    d: SatelliteIonSeries::default(),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::default(),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::none(),
    z: PrimaryIonSeries::none(),
    precursor: (Vec::new(), Vec::new(), (0, None), ChargeRange::PRECURSOR),
    immonium: None,
    modification_specific_neutral_losses: false,
    modification_specific_diagnostic_ions: None,
    glycan: GlycanModel::DISALLOW,
    allow_cross_link_cleavage: false,
});

static MODEL_UVPD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default().variants(vec![0, 1, 2]),
    b: PrimaryIonSeries::default().variants(vec![0, 2]),
    c: PrimaryIonSeries::default(),
    d: SatelliteIonSeries::base(),
    v: SatelliteIonSeries::base(),
    w: SatelliteIonSeries::base(),
    x: PrimaryIonSeries::default().variants(vec![-1, 0, 1, 2]),
    y: PrimaryIonSeries::default().variants(vec![-2, -1, 0]),
    z: PrimaryIonSeries::default(),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, IMMONIUM_LOSSES.clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::DISALLOW,
    allow_cross_link_cleavage: false,
});

static MODEL_ETHCD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default().location(Location::TakeN { skip: 0, take: 1 }),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    d: SatelliteIonSeries::base(),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .default_peptide_fragment(GlycanPeptideFragment::FULL)
        .peptide_fragment_rules(vec![
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::CORE,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::FREE,
            ),
        ]),
    allow_cross_link_cleavage: true,
});

static MODEL_EAD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, IMMONIUM_LOSSES.clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .peptide_fragment_rules(vec![
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::CORE,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::FREE,
            ),
        ]),
    allow_cross_link_cleavage: true,
});

static MODEL_EACID: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .peptide_fragment_rules(vec![
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::CORE,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::c, FragmentKind::z, FragmentKind::w],
                GlycanPeptideFragment::FULL,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                vec![FragmentKind::b, FragmentKind::y, FragmentKind::v],
                GlycanPeptideFragment::FREE,
            ),
        ]),
    allow_cross_link_cleavage: true,
});

static MODEL_CID_HCD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .location(Location::TakeN { skip: 0, take: 1 })
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::none(),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::default(),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::none(),
    precursor: (
        vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .default_peptide_fragment(GlycanPeptideFragment::CORE)
        .peptide_fragment_rules(vec![
            (
                vec![AminoAcid::Asparagine, AminoAcid::Tryptophan],
                Vec::new(),
                GlycanPeptideFragment::CORE,
            ),
            (
                vec![AminoAcid::Serine, AminoAcid::Threonine],
                Vec::new(),
                GlycanPeptideFragment::FREE,
            ),
        ]),
    allow_cross_link_cleavage: true,
});

static MODEL_ETD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::none(),
    b: PrimaryIonSeries::none(),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![-1, 0, 2]),
    d: SatelliteIonSeries::default(),
    v: SatelliteIonSeries::base(),
    w: SatelliteIonSeries::default()
        .location(SatelliteLocation {
            rules: vec![
                (vec![AminoAcid::Methionine], 5),
                (vec![AminoAcid::Leucine, AminoAcid::GlutamicAcid], 2),
                (vec![AminoAcid::Isoleucine, AminoAcid::AsparticAcid], 1),
                (
                    vec![
                        AminoAcid::Valine,
                        AminoAcid::Asparagine,
                        AminoAcid::Threonine,
                        AminoAcid::Serine,
                        AminoAcid::Tryptophan,
                        AminoAcid::Histidine,
                        AminoAcid::Phenylalanine,
                        AminoAcid::Tyrosine,
                    ],
                    0,
                ),
            ],
            base: None,
        })
        .variants(vec![-1, 0, 1]),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))])
        .variants(vec![-1, 0, 1, 2]),
    precursor: (
        vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
        ],
        vec![
            (
                vec![AminoAcid::AsparticAcid],
                vec![NeutralLoss::Loss(1, molecular_formula!(C 1 H 1 O 2))],
            ),
            (
                vec![AminoAcid::GlutamicAcid],
                vec![NeutralLoss::Loss(1, molecular_formula!(C 2 H 3 O 2))],
            ),
        ],
        (2, None),
        ChargeRange {
            start: ChargePoint::Relative(-2),
            end: ChargePoint::Relative(0),
        },
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow().default_peptide_fragment(GlycanPeptideFragment::FREE),
    allow_cross_link_cleavage: true,
});

static MODEL_TD_ETD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::none(),
    b: PrimaryIonSeries::none(),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
        ])
        .variants(vec![0, 1, 2]),
    d: SatelliteIonSeries::default(),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::default(),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::none(),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
        ])
        .variants(vec![-1, 0, 1, 2]),
    precursor: (
        vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
            NeutralLoss::Gain(1, molecular_formula!(H 1)),
            NeutralLoss::Gain(2, molecular_formula!(H 1)),
            NeutralLoss::Gain(3, molecular_formula!(H 1)),
        ],
        vec![
            (
                vec![AminoAcid::AsparticAcid],
                vec![NeutralLoss::Loss(1, molecular_formula!(C 1 H 1 O 2))],
            ),
            (
                vec![AminoAcid::GlutamicAcid],
                vec![NeutralLoss::Loss(1, molecular_formula!(C 2 H 3 O 2))],
            ),
        ],
        (0, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::DISALLOW,
    allow_cross_link_cleavage: true,
});

impl FragmentationModel {
    /// Generate all possible fragments
    pub fn all() -> &'static Self {
        LazyLock::force(&MODEL_ALL)
    }

    /// Generate no fragments (except for precursor)
    pub fn none() -> &'static Self {
        LazyLock::force(&MODEL_NONE)
    }

    /// UVPD
    /// 10.1021/acs.chemrev.9b00440 and 10.1021/jacs.6b05147
    pub fn uvpd() -> &'static Self {
        LazyLock::force(&MODEL_UVPD)
    }

    /// electron-transfer/higher-energy collisional dissociation
    pub fn ethcd() -> &'static Self {
        LazyLock::force(&MODEL_ETHCD)
    }

    /// EAD
    pub fn ead() -> &'static Self {
        LazyLock::force(&MODEL_EAD)
    }

    /// EAciD
    pub fn eacid() -> &'static Self {
        LazyLock::force(&MODEL_EACID)
    }

    /// CID Hcd
    pub fn cid_hcd() -> &'static Self {
        LazyLock::force(&MODEL_CID_HCD)
    }

    /// ETD 10.1002/jms.3919
    pub fn etd() -> &'static Self {
        LazyLock::force(&MODEL_ETD)
    }

    /// Top Down ETD
    pub fn td_etd() -> &'static Self {
        LazyLock::force(&MODEL_TD_ETD)
    }
}

/// All losses from the base immonium ions. Compiled from the sources below.
#[doc = include_str!("immonium_losses.md")]
pub(super) static IMMONIUM_LOSSES: LazyLock<Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>> =
    LazyLock::new(||
    // TODO: For B/Z there are common immonium ions, but the mass is the same (meaning the loss is different), find a way of representing that
    vec![(
        vec![AminoAcid::Arginine], vec![
            NeutralLoss::Gain(1, molecular_formula!(C 2 O 2)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 2)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 2 H 2 N 2)),
            NeutralLoss::Loss(1, molecular_formula!(C 3 H 6 N 2)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 5 N 3)),
            NeutralLoss::Loss(1, molecular_formula!(C 3 H 4 N 2 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(C 4 H 8 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 4 H 10 N 2)),
        ]),
        (vec![AminoAcid::Asparagine], vec![NeutralLoss::Loss(1, molecular_formula!(H 3 N 1))]),
        (vec![AminoAcid::AsparticAcid, AminoAcid::GlutamicAcid, AminoAcid::Serine], vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]),
        (vec![AminoAcid::Glutamine], vec![
            NeutralLoss::Gain(1, molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 3 N 1 O 1)),
        ]),
        (vec![AminoAcid::Histidine], vec![
            NeutralLoss::Gain(1, molecular_formula!(C 2 O 2)),
            NeutralLoss::Gain(1, molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 2 N 1)),
        ]),
        (vec![AminoAcid::Leucine, AminoAcid::Isoleucine, AminoAcid::AmbiguousLeucine], vec![
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 2)),
            NeutralLoss::Loss(1, molecular_formula!(C 3 H 6)),
        ]),
        (vec![AminoAcid::Lysine], vec![
            NeutralLoss::Gain(1, molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(C -2 H 1 N 1 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 5 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 2 H 7 N 1)),
        ]),
        (vec![AminoAcid::Methionine], vec![
            NeutralLoss::Loss(1, molecular_formula!(H 2 S 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 2 H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 4 S 1)),
        ]),
        (vec![AminoAcid::Phenylalanine], vec![NeutralLoss::Gain(1, molecular_formula!(C 2 O 2))]),
        (vec![AminoAcid::Threonine], vec![NeutralLoss::Loss(1, molecular_formula!(H 2 N 1))]),
        (vec![AminoAcid::Tryptophan], vec![
            NeutralLoss::Loss(1, molecular_formula!(H 4 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 1 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 2 H 4 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 4 H 6 N 2)),
        ]),
        (vec![AminoAcid::Tyrosine], vec![
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 3 N 1 O 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 5 H 7 N 1)),
        ]),
        (vec![AminoAcid::Valine], vec![
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 1 O -1)),
            NeutralLoss::Loss(1, molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 2 N 1)),
            NeutralLoss::Loss(1, molecular_formula!(C 1 H 5 N 1)),
        ]),
    ]);

/// Generate all uncharged diagnostic ions for this monosaccharide.
/// According to: <https://doi.org/10.1016/j.trac.2018.09.007>.
pub(super) static GLYCAN_LOSSES: LazyLock<Vec<(MonoSaccharide, bool, Vec<NeutralLoss>)>> =
    LazyLock::new(|| {
        vec![
            (
                MonoSaccharide::new(BaseSugar::Hexose(None), &[]),
                false,
                vec![
                    NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(1, molecular_formula!(H 4 O 2)),
                    NeutralLoss::Loss(1, molecular_formula!(C 1 H 6 O 3)),
                    NeutralLoss::Loss(1, molecular_formula!(C 2 H 6 O 3)),
                ],
            ),
            (
                MonoSaccharide::new(BaseSugar::Hexose(None), &[GlycanSubstituent::NAcetyl]),
                false,
                vec![
                    NeutralLoss::Loss(1, molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(1, molecular_formula!(H 4 O 2)),
                    NeutralLoss::Loss(1, molecular_formula!(C 2 H 4 O 2)),
                    NeutralLoss::Loss(1, molecular_formula!(C 1 H 6 O 3)),
                    NeutralLoss::Loss(1, molecular_formula!(C 2 H 6 O 3)),
                    NeutralLoss::Loss(1, molecular_formula!(C 4 H 8 O 4)),
                ],
            ),
            (
                MonoSaccharide::new(
                    BaseSugar::Nonose(None),
                    &[
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Acetyl,
                        GlycanSubstituent::Acid,
                    ],
                )
                .with_name("NeuAc"),
                false,
                vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
            ),
            (
                MonoSaccharide::new(
                    BaseSugar::Nonose(None),
                    &[
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                )
                .with_name("NeuGc"),
                false,
                vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))],
            ),
        ]
    });
