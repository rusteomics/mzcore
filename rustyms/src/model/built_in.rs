use std::sync::LazyLock;

use crate::{
    glycan::{BaseSugar, GlycanSubstituent, MonoSaccharide},
    model::{
        ChargePoint, ChargeRange, FragmentationModel, GlycanModel, Location, PrimaryIonSeries,
        SatelliteIonSeries, SatelliteLocation,
    },
    AminoAcid, NeutralLoss,
};

static MODEL_ALL: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    precursor: (
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (1, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, immonium_losses().clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
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
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, immonium_losses().clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::DISALLOW,
    allow_cross_link_cleavage: false,
});

static MODEL_ETHCD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default().location(Location::TakeN { skip: 0, take: 1 }),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    d: SatelliteIonSeries::base(),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
        .variants(vec![0, 1]),
    precursor: (
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    allow_cross_link_cleavage: true,
});

static MODEL_EAD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    precursor: (
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: Some((ChargeRange::ONE, immonium_losses().clone())),
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    allow_cross_link_cleavage: true,
});

static MODEL_EACID: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    w: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    x: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    precursor: (
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::ONE_TO_PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    allow_cross_link_cleavage: true,
});

static MODEL_CID_HCD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::default()
        .location(Location::TakeN { skip: 0, take: 1 })
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    b: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    c: PrimaryIonSeries::none(),
    d: SatelliteIonSeries::base()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::default(),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::none(),
    precursor: (
        vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
        Vec::new(),
        (0, None),
        ChargeRange::PRECURSOR,
    ),
    immonium: None,
    modification_specific_neutral_losses: true,
    modification_specific_diagnostic_ions: Some(ChargeRange::ONE),
    glycan: GlycanModel::default_allow(),
    allow_cross_link_cleavage: true,
});

static MODEL_ETD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::none(),
    b: PrimaryIonSeries::none(),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
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
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))])
        .variants(vec![-1, 0, 1, 2]),
    precursor: (
        vec![
            NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
        ],
        vec![
            (
                vec![AminoAcid::AsparticAcid],
                vec![NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2))],
            ),
            (
                vec![AminoAcid::GlutamicAcid],
                vec![NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2))],
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
    glycan: GlycanModel::DISALLOW,
    allow_cross_link_cleavage: true,
});

static MODEL_TD_ETD: LazyLock<FragmentationModel> = LazyLock::new(|| FragmentationModel {
    a: PrimaryIonSeries::none(),
    b: PrimaryIonSeries::none(),
    c: PrimaryIonSeries::default()
        .neutral_losses(vec![
            NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
        ])
        .variants(vec![0, 1, 2]),
    d: SatelliteIonSeries::default(),
    v: SatelliteIonSeries::default(),
    w: SatelliteIonSeries::default(),
    x: PrimaryIonSeries::none(),
    y: PrimaryIonSeries::none(),
    z: PrimaryIonSeries::default()
        .neutral_losses(vec![
            NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
        ])
        .variants(vec![-1, 0, 1, 2]),
    precursor: (
        vec![
            NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
            NeutralLoss::Gain(molecular_formula!(H 1)),
            NeutralLoss::Gain(molecular_formula!(H 2)),
            NeutralLoss::Gain(molecular_formula!(H 3)),
        ],
        vec![
            (
                vec![AminoAcid::AsparticAcid],
                vec![NeutralLoss::Loss(molecular_formula!(C 1 H 1 O 2))],
            ),
            (
                vec![AminoAcid::GlutamicAcid],
                vec![NeutralLoss::Loss(molecular_formula!(C 2 H 3 O 2))],
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
///
/// | AA    | [Wikipedia](https://upload.wikimedia.org/wikipedia/commons/thumb/0/01/Amino_acid_fragment_ions.png/400px-Amino_acid_fragment_ions.png) |  0.1016/1044-0305(93)87006-X  | [ionsource](https://www.ionsource.com/Card/immon/immon.htm) | [10.1002/chin.199624319](http://dx.doi.org/10.1002/chin.199624319) | [Prospector   (MS-Comp)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=mscomp) | [10.1186/1477-5956-9-2](http://dx.doi.org/10.1186/1477-5956-9-2) |  10.1016/j.ymeth.2004.08.013   | 10.1385/1597452750 (table 5)  | 10.1021/ac902712f  | [Prospector   (MS-Product)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msproduct) | [ThermoFisher](https://tools.thermofisher.com/content/sfs/brochures/cms_040030.pdf) | 10.1074/mcp.O113.035915 | 10.1074/mcp.O113.035915 | 10.1021/ac902712f | [Prospector   (MS-Product)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msproduct) |  | 10.1385/1597452750 (table 5) | Sources | Best mass | Best formula | Loss     | Loss formula | Interpreted loss | Interpreted formula     | Final      |
/// |-------|---------------------------------------------------------------------------------------------------------------------------|-----------------------------|------------------------------------------------|------------------------------------------|-----------------------------------------------------------------------|-----------------------------------------|-----------------------------|------------------------------|-------------------|--------------------------------------------------------------------------|---------------------------------------------------------------------|-------------------------|-------------------------|-------------------|--------------------------------------------------------------------------|------------------------------|---------|----------:|--------------|----------|--------------|------------------|-------------------------|------------|
/// | A     | 44                                                                                                                        | 44                          |                                                | 44                                       |                                                                       | 44                                      |                             | 44.05                        |                   |                                                                          | 44.0500                                                             |                         |                         |                   |                                                                          |                              | 6       |   44.0500 |              |          |              |                  |                         |            |
/// | R     | 129                                                                                                                       | 129                         |                                                | 129                                      |                                                                       | 129                                     |                             | 129.11                       |                   |                                                                          | 129.1140                                                            | 129.1135                | C5H13N4+                |                   |                                                                          |                              | 8       |  129.1138 | C5H13N4+     |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             | 185                                            |                                          | 185                                                                   |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 2       |       185 |              | -55.8862 |              | C-2O-2           |                         | C-2O-2     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 115.09                       |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |    115.09 |              | 14.0238  |              | C1H2             |                         | C1H2       |
/// |       | 112                                                                                                                       | 112                         | 112                                            | 112                                      | 112                                                                   | 112                                     | 112.09                      | 112.09                       |                   | 112.0869                                                                 | 112.0875                                                            |                         |                         |                   | C5H10N3+                                                                 | C5H10N3+                     | 12      |  112.0872 | C5H10N3+     | 17.0266  | H3N1         |                  |                         | H3N1       |
/// |       | 100                                                                                                                       | 100                         | 100                                            | 100                                      | 100                                                                   | 100                                     | 100.09                      |                              |                   | 100.0869                                                                 | 100.0875                                                            |                         |                         |                   | C4H10N3+                                                                 |                              | 10      |  100.0872 | C4H10N3+     | 29.0266  | C1H3N1       |                  |                         | C1H3N1     |
/// |       | 87                                                                                                                        | 87                          | 87                                             | 87                                       | 87                                                                    | 87                                      | 87.09                       | 87.09                        | 87.0922           | 87.0917                                                                  |                                                                     |                         |                         | C4H11N2+          | C4H11N2+                                                                 |                              | 12      |   87.0920 | C4H11N2+     | 42.0218  | C2H2N2       |                  |                         | C2H2N2     |
/// |       | 73                                                                                                                        | 73                          | 73                                             |                                          | 73                                                                    | 72                                      | 73.00                       |                              | 73.0640           |                                                                          |                                                                     |                         |                         | C2H7N3+           |                                                                          |                              | 8       |   73.0640 | C2H7N3+      | 56.0498  | C3H6N1       |                  |                         | C3H6N1     |
/// |       | 70                                                                                                                        | 70                          | 70                                             | 70                                       | 70                                                                    | 70                                      | 70.07                       | 70.07                        | 70.0657           | 70.0651                                                                  | 70.0657                                                             |                         |                         | C4H8N1+           | C4H8N1+                                                                  |                              | 13      |   70.0655 | C4H8N1+      | 59.0483  | C1H5N3       |                  |                         | C1H5N3     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 60.06                        |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |     60.06 |              | 69.0538  |              | C3H4N2O-1        | C2H6N1O1+               | C3H4N2O-1  |
/// |       |                                                                                                                           | 59                          |                                                |                                          |                                                                       | 59                                      |                             |                              | 59.0483           |                                                                          |                                                                     |                         |                         | CH5N3+            |                                                                          |                              | 4       |   59.0483 | CH5N3+       | 70.0655  | C4H8N1       |                  |                         | C4H8N1     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              | 43.0296           |                                                                          |                                                                     |                         |                         | C1H3N2+           |                                                                          |                              | 2       |   43.0296 | C1H3N2+      | 86.0842  | C4H10N2      |                  |                         | C4H10N2    |
/// |       | 29                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        29 |              | 100.1138 |              |                  | H1N2/C1H1O1/C1H3N1/C2H5 |            |
/// | N     | 87                                                                                                                        | 87                          | 87                                             | 87                                       | 87                                                                    | 87                                      | 87.09                       | 87.06                        |                   | 87.0553                                                                  | 87.0558                                                             |                         |                         |                   | C3H7N2O1+                                                                |                              | 11      |   87.0556 | C3H7N2O1+    |          |              |                  |                         |            |
/// |       | 70                                                                                                                        | 70                          | 70                                             | 70                                       |                                                                       | 70                                      |                             | 70.03                        |                   |                                                                          | 70.0293                                                             |                         |                         |                   |                                                                          | C3H4N1O1+                    | 8       |   70.0293 | C3H4N1O1+    | 17.0263  | H3N1         |                  |                         | H3N1       |
/// | D     | 88                                                                                                                        | 88                          | 88                                             | 88                                       | 88                                                                    | 88                                      | 88.04                       | 88.04                        | 88.0399           | 88.0393                                                                  | 88.0399                                                             |                         |                         | C3H6N1O2+         | C3H6N1O2+                                                                |                              | 13      |   88.0397 | C3H6N1O2+    |          |              |                  |                         |            |
/// |       | 70                                                                                                                        |                             | 70                                             | 70                                       |                                                                       | 70                                      |                             | 70.03                        |                   |                                                                          | 70.0293                                                             |                         |                         |                   |                                                                          | C3H4N1O1+                    | 7       |   70.0293 | C3H4N1O1+    | 18.0104  | H2O1         |                  |                         | H2O1       |
/// | C     | 76                                                                                                                        |                             |                                                | 76                                       |                                                                       | 76                                      |                             |                              |                   |                                                                          | 76.0221                                                             |                         |                         |                   |                                                                          |                              | 4       |   76.0221 |              |          |              |                  |                         |            |
/// | E     | 102                                                                                                                       | 102                         |                                                | 102                                      | 102                                                                   | 102                                     | 102.06                      | 102.05                       | 102.0555          | 102.0550                                                                 | 102.0555                                                            | 102.0550                | C4H8N1O2+               | C4H8N1O2+         | C4H8N1O2+                                                                |                              | 14      |  102.0553 | C4H8N1O2+    |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             |                                                | 91                                       |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        91 |              | 11.0553  |              |                  |                         |            |
/// |       |                                                                                                                           |                             |                                                | 84                                       |                                                                       |                                         |                             | 84.04                        |                   |                                                                          | 84.0449                                                             |                         |                         |                   |                                                                          | C4H6N1O1+                    | 4       |   84.0449 | C4H6N1O1+    | 18.0104  | H2O1         |                  |                         | H2O1       |
/// | Q     | 101                                                                                                                       | 101                         | 101                                            | 101                                      | 101                                                                   | 101                                     | 101.11                      | 101.11                       |                   | 101.0709                                                                 | 101.0715                                                            | 101.0709                | C4H9N2O1+               |                   | C4H9N2O1+                                                                |                              | 13      |  101.0711 | C4H9N2O1+    |          |              |                  |                         |            |
/// |       | 129                                                                                                                       | 129                         | 129                                            | 129                                      | 129                                                                   | 129                                     | 129.1                       | 129.11                       |                   | 129.0659                                                                 | 129.1028                                                            |                         |                         |                   | C5H9N2O2+                                                                |                              | 11      |  129.0844 | C5H9N2O2+    | -28.0133 | C-1O-1       |                  |                         | C-1O-1     |
/// |       | 84                                                                                                                        | 84                          | 84                                             | 84                                       | 84                                                                    | 84                                      | 84.08                       | 84.04                        | 84.0813           | 84.0444                                                                  | 84.0449                                                             |                         |                         | C5H10N1+          | C4H6N1O1+                                                                | C4H6N1O1+                    | 14      |   84.0569 | C5H10N1+     | 17.0142  | H3N1         |                  |                         | H3N1       |
/// |       | 56                                                                                                                        |                             |                                                | 56                                       |                                                                       | 56                                      |                             | 56.05                        |                   |                                                                          | 56.0500                                                             |                         |                         |                   |                                                                          |                              | 5       |   56.0500 |              | 45.0211  |              | C1H3N1O1         |                         | C1H3N1O1   |
/// | G     | 30                                                                                                                        | 30                          |                                                | 30                                       |                                                                       | 30                                      |                             | 30.03                        | 30.0344           |                                                                          | 30.0344                                                             |                         |                         | C1H4N1+           |                                                                          |                              | 8       |   30.0344 | C1H4N1+      |          |              |                  |                         |            |
/// | H     | 110                                                                                                                       | 110                         | 110                                            | 110                                      | 110                                                                   | 110                                     | 110.07                      | 110.07                       | 110.0718          | 110.0713                                                                 | 110.0718                                                            | 110.0713                | C5H8N3+                 | C5H8N3+           | C5H8N3+                                                                  |                              | 15      |  110.0716 | C5H8N3+      |          |              |                  |                         |            |
/// |       | 166                                                                                                                       | 166                         |                                                |                                          |                                                                       | 166                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       166 |              | -55.9284 |              | C-2O-2           |                         | C-2O-2     |
/// |       | 138                                                                                                                       | 138                         |                                                |                                          |                                                                       | 138                                     | 138.07                      |                              |                   | 138.0662                                                                 |                                                                     |                         |                         |                   | C6H8N3O1+                                                                |                              | 6       |  138.0662 | C6H8N3O1+    | -27.9946 | C-1O-1       |                  |                         | C-1O-1     |
/// |       | 123                                                                                                                       | 123                         |                                                |                                          |                                                                       | 123                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       123 |              | -12.9284 |              | H3O-1            |                         | H3O-1      |
/// |       | 121                                                                                                                       | 121                         |                                                |                                          |                                                                       | 121                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       121 |              | -10.9284 |              | H5O-1            |                         | H5O-1      |
/// |       | 82                                                                                                                        | 82                          |                                                |                                          |                                                                       | 82                                      |                             |                              | 82.0531           |                                                                          |                                                                     |                         |                         | C4H6N2+           |                                                                          |                              | 5       |   82.0531 | C4H6N2+      | 28.0185  | C1H2N1       |                  |                         | C1H2N1     |
/// | I/L/J | 86                                                                                                                        | 86                          | 86                                             | 86                                       | 86                                                                    | 86                                      | 86.1                        | 86.10                        | 86.0970           | 86.0964                                                                  | 86.0970                                                             |                         |                         | C5H12N+           | C5H12N1+                                                                 |                              | 13      |   86.0968 | C5H12N+      |          |              |                  |                         |            |
/// |       | 72                                                                                                                        | 72                          |                                                | 72                                       |                                                                       | 72                                      |                             |                              |                   |                                                                          | 72.0449                                                             |                         |                         |                   |                                                                          |                              | 5       |   72.0449 |              | 14.0519  |              | C1H2             |                         | C1H2       |
/// |       | 44                                                                                                                        |                             |                                                | 44                                       |                                                                       | 44                                      |                             |                              |                   |                                                                          | 44.0500                                                             |                         |                         |                   |                                                                          |                              | 4       |   44.0500 |              | 42.0468  |              | C3H6             |                         | C3H6       |
/// | K     | 101                                                                                                                       | 101                         | 101                                            | 101                                      | 101                                                                   | 101                                     | 101.11                      | 101.11                       |                   | 101.1073                                                                 | 101.1079                                                            |                         |                         |                   | C5H13N2+                                                                 |                              | 11      |  101.1076 | C5H13N2+     |          |              |                  |                         |            |
/// |       | 129                                                                                                                       | 129                         | 129                                            | 129                                      | 129                                                                   | 129                                     | 129.1                       | 129.11                       |                   | 129.1022                                                                 |                                                                     | 129.1022                | C6H13N2O1+              |                   | C6H13N2O1+                                                               |                              | 12      |  129.1022 | C6H13N2O1+   | -27.9946 | C-1O-1       |                  |                         | C-1O-1     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   | 126.0913                                                                 |                                                                     |                         |                         |                   | C7H12N1O1+                                                               |                              | 2       |  126.0913 | C7H12N1O1+   | -24.9837 | C-2H1N1O-1   |                  |                         | C-2H1N1O-1 |
/// |       | 112                                                                                                                       | 112                         |                                                | 112                                      |                                                                       | 112                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 4       |       112 |              | -10.8924 |              | H5O-1            |                         | H5O-1      |
/// |       | 84                                                                                                                        | 84                          | 84                                             | 84                                       | 84                                                                    | 84                                      | 84.08                       | 84.08                        | 84.0813           | 84.0808                                                                  | 84.0813                                                             |                         |                         | C5H10N1+          | C5H10N1+                                                                 | C5H10N1+                     | 14      |   84.0811 | C5H10N1+     | 17.0265  | H3N1         |                  |                         | H3N1       |
/// |       | 70                                                                                                                        | 70                          |                                                |                                          |                                                                       | 70                                      |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |        70 |              | 31.1076  |              | C1H5N1           |                         | C1H5N1     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 56.05                        |                   |                                                                          | 56.0500                                                             |                         |                         |                   |                                                                          |                              | 2       |   56.0500 |              | 45.0576  |              | C2H7N1           |                         | C2H7N1     |
/// | M     | 104                                                                                                                       | 104                         | 104                                            | 104                                      | 104                                                                   | 104                                     | 104.05                      | 104.06                       |                   | 104.0528                                                                 | 104.0534                                                            |                         |                         |                   | C4H10N1S1+                                                               |                              | 11      |  104.0531 | C4H10N1S1+   |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       | 70                                      |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              |         |        70 |              | 34.0531  |              | H2S1             |                         | H2S1       |
/// |       | 61                                                                                                                        | 61                          |                                                |                                          |                                                                       | 61                                      |                             |                              | 61.0112           |                                                                          |                                                                     |                         |                         | C2H5S1+           |                                                                          |                              | 5       |   61.0112 | C2H5S1+      | 43.0419  | C2H3N1       |                  |                         | C2H3N1     |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C3H6N1+                      | 1       |        ?? | C3H6N1+      | ??       | C1H4S1       |                  |                         | C1H4S1     |
/// | F     | 120                                                                                                                       | 120                         | 120                                            | 120                                      | 120                                                                   | 120                                     | 120.08                      | 120.08                       | 120.0813          | 120.0808                                                                 | 120.0813                                                            | 120.0808                | C8H10N+                 | C8H10N1+          | C8H10N1+                                                                 |                              | 15      |  120.0811 | C8H10N+      |          |              |                  |                         |            |
/// |       | 91                                                                                                                        | 91                          |                                                | 91                                       |                                                                       | 91                                      |                             |                              |                   |                                                                          | 91.0548                                                             |                         |                         |                   |                                                                          |                              | 5       |   91.0548 |              | 29.0263  |              | C1H3N1           | C7H7+                   | C1H3N1     |
/// | P     | 70                                                                                                                        | 70                          | 70                                             | 70                                       | 70                                                                    | 70                                      | 70.07                       | 70.07                        | 70.0657           | 70.0651                                                                  | 70.0657                                                             |                         |                         | C4H8N1+           | C4H8N1+                                                                  |                              | 13      |   70.0655 | C4H8N1+      |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             | 126                                            |                                          | 126                                                                   |                                         | 126.06                      |                              |                   | 126.055                                                                  |                                                                     |                         |                         |                   | C6H8N1O2+                                                                |                              | 5       |  126.0550 | C6H8N1O2+    | -55.9895 | C-2O-2       |                  |                         | C-2O-2     |
/// | S     | 60                                                                                                                        | 60                          | 60                                             | 60                                       | 60                                                                    | 60                                      | 60.04                       | 60.04                        |                   | 60.0444                                                                  | 60.0449                                                             |                         |                         |                   | C2H6N1O1+                                                                |                              | 11      |   60.0447 | C2H6N1O1+    |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C2H4N1+                      | 1       |        ?? | C2H4N1+      | ??       | H2O1         |                  |                         | H2O1       |
/// | T     | 74                                                                                                                        | 74                          | 74                                             | 74                                       | 74                                                                    | 74                                      |                             | 74.06                        | 74.0606           | 74.0600                                                                  | 74.0606                                                             |                         |                         | C3H8N1O1+         | C3H8N1O1+                                                                |                              | 12      |   74.0604 | C3H8N1O1+    |          |              |                  |                         |            |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C3H6N1O1+                    | 1       |        ?? | C3H6N1O1+    | ??       | H2N1         |                  |                         | H2N1       |
/// | W     | 159                                                                                                                       | 159                         | 159                                            | 159                                      | 159                                                                   | 159                                     | 159.09                      | 159.09                       |                   | 159.0917                                                                 | 159.0922                                                            | 159.0917                | C10H11N2+               |                   | C10H11N2+                                                                |                              | 13      |  159.0919 | C10H11N2+    |          |              |                  |                         |            |
/// |       |                                                                                                                           | 171                         |                                                |                                          |                                                                       | 171                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 2       |       171 |              | -11.9081 |              | H4O-1            |                         | H4O-1      |
/// |       | 170                                                                                                                       | 170                         | 170                                            |                                          | 170                                                                   | 170                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 5       |       170 |              | -10.9081 |              | H5O-1            |                         | H5O-1      |
/// |       | 132                                                                                                                       |                             |                                                | 132                                      |                                                                       | 132                                     |                             | 132.08                       |                   |                                                                          | 132.0813                                                            |                         |                         |                   |                                                                          |                              | 5       |  132.0813 |              | 27.0106  |              | C1H1N1           |                         | C1H1N1     |
/// |       | 130                                                                                                                       | 130                         | 130                                            | 130                                      | 130                                                                   | 130                                     |                             | 130.07                       |                   |                                                                          | 130.0657                                                            |                         |                         |                   |                                                                          |                              | 8       |  130.0657 |              | 29.0262  |              | C1H3N1           |                         | C1H3N1     |
/// |       | 117                                                                                                                       | 117                         | 117                                            | 117                                      | 117                                                                   | 117                                     |                             |                              |                   |                                                                          | 117.0578                                                            |                         |                         |                   |                                                                          |                              | 7       |  117.0578 |              | 42.0341  |              | C2H4N1           |                         | C2H4N1     |
/// |       | 100                                                                                                                       |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |       100 |              | 59.0919  |              | C3H9N1/C2H7N2    |                         |            |
/// |       |                                                                                                                           |                             |                                                | 77                                       |                                                                       | 77                                      |                             |                              |                   |                                                                          | 77.0391                                                             |                         |                         |                   |                                                                          |                              | 3       |   77.0391 |              | 82.0528  |              | C4H6N2           | C6H5                    | C4H6N2     |
/// |       | 11                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        11 |              | 148.0919 |              |                  |                         |            |
/// | Y     | 136                                                                                                                       | 136                         | 136                                            | 136                                      | 136                                                                   | 136                                     | 136.08                      | 136.08                       | 136.0762          | 136.0757                                                                 | 136.0762                                                            | 136.0757                | C8H10N1O1+              | C8H10N1O1+        | C8H10N1O1+                                                               |                              | 15      |  136.0760 | C8H10N1O1+   |          |              |                  |                         |            |
/// |       | 107                                                                                                                       | 107                         |                                                | 107                                      |                                                                       | 107                                     |                             |                              |                   |                                                                          | 107.0497                                                            |                         |                         |                   |                                                                          |                              | 5       |  107.0497 |              | 29.0263  |              | C1H3N1           |                         | C1H3N1     |
/// |       | 91                                                                                                                        | 91                          |                                                | 91                                       |                                                                       | 91                                      |                             |                              |                   |                                                                          | 91.0548                                                             |                         |                         |                   |                                                                          |                              | 5       |   91.0548 |              | 45.0212  |              | C1H3N1O1         |                         | C1H3N1O1   |
/// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              | 55.0184           |                                                                          |                                                                     |                         |                         | C3H3O1+           |                                                                          |                              | 2       |   55.0184 | C3H3O1+      | 81.0576  | C5H7N1       |                  |                         | C5H7N1     |
/// | V     | 72                                                                                                                        | 72                          | 72                                             | 72                                       | 72                                                                    | 72                                      | 72.08                       | 72.08                        | 72.0813           | 72.0808                                                                  | 72.0813                                                             |                         |                         | C4H10N1+          | C4H10N1+                                                                 |                              | 13      |   72.0811 | C4H10N1+     |          |              |                  |                         |            |
/// |       | 69                                                                                                                        |                             |                                                | 69                                       |                                                                       | 69                                      |                             |                              |                   |                                                                          | 69.0704                                                             |                         |                         |                   |                                                                          |                              | 4       |   69.0704 |              | 3.0107   |              | C1H1O-1          |                         | C1H1O-1    |
/// |       | 55                                                                                                                        |                             |                                                | 55                                       |                                                                       | 55                                      |                             |                              |                   |                                                                          | 55.0548                                                             |                         |                         |                   |                                                                          |                              | 4       |   55.0548 |              | 17.0263  |              | H3N1             |                         | H3N1       |
/// |       | 44                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        44 |              | 28.0811  |              | C1H2N1           |                         | C1H2N1     |
/// |       |                                                                                                                           |                             |                                                | 41                                       |                                                                       | 41                                      |                             |                              |                   |                                                                          | 41.0391                                                             |                         |                         |                   |                                                                          |                              | 3       |   41.0391 |              | 31.0420  |              | C1H5N1           |                         | C1H5N1     |
static IMMONIUM_LOSSES: LazyLock<Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>> = LazyLock::new(||
    // TODO: For B/Z there are common immonium ions, but the mass is the same (meaning the loss is different), find a way of representing that
    vec![(
        vec![AminoAcid::Arginine], vec![
            NeutralLoss::Gain(molecular_formula!(C 2 O 2)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 2)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 2 H 2 N 2)),
            NeutralLoss::Loss(molecular_formula!(C 3 H 6 N 2)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 3)),
            NeutralLoss::Loss(molecular_formula!(C 3 H 4 N 2 O -1)),
            NeutralLoss::Loss(molecular_formula!(C 4 H 8 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 4 H 10 N 2)),
        ]),
        (vec![AminoAcid::Asparagine], vec![NeutralLoss::Loss(molecular_formula!(H 3 N 1))]),
        (vec![AminoAcid::AsparticAcid, AminoAcid::GlutamicAcid, AminoAcid::Serine], vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]),
        (vec![AminoAcid::Glutamine], vec![
            NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1 O 1)),
        ]),
        (vec![AminoAcid::Histidine], vec![
            NeutralLoss::Gain(molecular_formula!(C 2 O 2)),
            NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(H 3 O -1)),
            NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 2 N 1)),
        ]),
        (vec![AminoAcid::Leucine, AminoAcid::Isoleucine, AminoAcid::AmbiguousLeucine], vec![
            NeutralLoss::Loss(molecular_formula!(C 1 H 2)),
            NeutralLoss::Loss(molecular_formula!(C 3 H 6)),
        ]),
        (vec![AminoAcid::Lysine], vec![
            NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(C -2 H 1 N 1 O -1)),
            NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 2 H 7 N 1)),
        ]),
        (vec![AminoAcid::Methionine], vec![
            NeutralLoss::Loss(molecular_formula!(H 2 S 1)),
            NeutralLoss::Loss(molecular_formula!(C 2 H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 4 S 1)),
        ]),
        (vec![AminoAcid::Phenylalanine], vec![NeutralLoss::Gain(molecular_formula!(C 2 O 2))]),
        (vec![AminoAcid::Threonine], vec![NeutralLoss::Loss(molecular_formula!(H 2 N 1))]),
        (vec![AminoAcid::Tryptophan], vec![
            NeutralLoss::Loss(molecular_formula!(H 4 O -1)),
            NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 1 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 2 H 4 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 4 H 6 N 2)),
        ]),
        (vec![AminoAcid::Tyrosine], vec![
            NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1 O 1)),
            NeutralLoss::Loss(molecular_formula!(C 5 H 7 N 1)),
        ]),
        (vec![AminoAcid::Valine], vec![
            NeutralLoss::Loss(molecular_formula!(C 1 H 1 O -1)),
            NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 2 N 1)),
            NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 1)),
        ]),
    ]);

pub fn immonium_losses() -> &'static Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)> {
    LazyLock::force(&IMMONIUM_LOSSES)
}

/// Generate all uncharged diagnostic ions for this monosaccharide.
/// According to: <https://doi.org/10.1016/j.trac.2018.09.007>.
static GLYCAN_LOSSES: LazyLock<Vec<(MonoSaccharide, bool, Vec<NeutralLoss>)>> =
    LazyLock::new(|| {
        vec![
            (
                MonoSaccharide::new(BaseSugar::Hexose(None), &[]),
                false,
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 4 O 2)),
                    NeutralLoss::Loss(molecular_formula!(C 1 H 6 O 3)),
                    NeutralLoss::Loss(molecular_formula!(C 2 H 6 O 3)),
                ],
            ),
            (
                MonoSaccharide::new(BaseSugar::Hexose(None), &[GlycanSubstituent::NAcetyl]),
                false,
                vec![
                    NeutralLoss::Loss(molecular_formula!(H 2 O 1)),
                    NeutralLoss::Loss(molecular_formula!(H 4 O 2)),
                    NeutralLoss::Loss(molecular_formula!(C 2 H 4 O 2)),
                    NeutralLoss::Loss(molecular_formula!(C 1 H 6 O 3)),
                    NeutralLoss::Loss(molecular_formula!(C 2 H 6 O 3)),
                    NeutralLoss::Loss(molecular_formula!(C 4 H 8 O 4)),
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
                ),
                false,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
            (
                MonoSaccharide::new(
                    BaseSugar::Nonose(None),
                    &[
                        GlycanSubstituent::Amino,
                        GlycanSubstituent::Glycolyl,
                        GlycanSubstituent::Acid,
                    ],
                ),
                false,
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))],
            ),
        ]
    });

pub fn glycan_losses() -> &'static Vec<(MonoSaccharide, bool, Vec<NeutralLoss>)> {
    LazyLock::force(&GLYCAN_LOSSES)
}
