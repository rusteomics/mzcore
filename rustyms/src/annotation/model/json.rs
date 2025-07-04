use std::ops::RangeInclusive;

use itertools::Itertools;
use serde_json::Value;

use crate::{
    annotation::model::{
        ChargeRange, GlycanModel, GlycanPeptideFragment, Location, PrimaryIonSeries,
        SatelliteIonSeries, SatelliteLocation,
    },
    error::{Context, CustomError},
    fragment::{FragmentKind, NeutralLoss},
    parse_json::{ParseJson, use_serde},
    prelude::{AminoAcid, FragmentationModel},
};

impl ParseJson for FragmentationModel {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::show(map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","))
            };
            Ok(Self {
                a: PrimaryIonSeries::from_json_value(map.remove("a").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'a' is missing",
                        context(&map),
                    )
                })?)?,
                b: PrimaryIonSeries::from_json_value(map.remove("b").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'b' is missing",
                        context(&map),
                    )
                })?)?,
                c: PrimaryIonSeries::from_json_value(map.remove("c").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'c' is missing",
                        context(&map),
                    )
                })?)?,
                d: SatelliteIonSeries::from_json_value(map.remove("d").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'd' is missing",
                        context(&map),
                    )
                })?)?,
                v: SatelliteIonSeries::from_json_value(map.remove("v").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'v' is missing",
                        context(&map),
                    )
                })?)?,
                w: SatelliteIonSeries::from_json_value(map.remove("w").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'w' is missing",
                        context(&map),
                    )
                })?)?,
                x: PrimaryIonSeries::from_json_value(map.remove("x").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'x' is missing",
                        context(&map),
                    )
                })?)?,
                y: PrimaryIonSeries::from_json_value(map.remove("y").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'y' is missing",
                        context(&map),
                    )
                })?)?,
                z: PrimaryIonSeries::from_json_value(map.remove("z").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'z' is missing",
                        context(&map),
                    )
                })?)?,
                precursor: <(
                        Vec<NeutralLoss>,
                        Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
                        (u8, Option<Vec<AminoAcid>>),
                        ChargeRange,
                    )>::from_json_value(map.remove("precursor").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'precursor' is missing",
                        context(&map),
                    )
                })?)?,
                immonium: Option::from_json_value(map.remove("immonium").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'immonium' is missing",
                        context(&map),
                    )
                })?)?,
                modification_specific_neutral_losses: bool::from_json_value(map.remove("modification_specific_neutral_losses").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'modification_specific_neutral_losses' is missing",
                        context(&map),
                    )
                })?)?,
                modification_specific_diagnostic_ions: Option::from_json_value(map.remove("modification_specific_diagnostic_ions").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'modification_specific_diagnostic_ions' is missing",
                        context(&map),
                    )
                })?)?,
                glycan: GlycanModel::from_json_value(map.remove("glycan").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'glycan' is missing",
                        context(&map),
                    )
                })?)?,
                allow_cross_link_cleavage: bool::from_json_value(map.remove("allow_cross_link_cleavage").ok_or_else(|| {
                    CustomError::error(
                        "Invalid FragmentationModel",
                        "The required property 'allow_cross_link_cleavage' is missing",
                        context(&map),
                    )
                })?)?,
            })
        } else {
            Err(CustomError::error(
                "Invalid Database SimpleModification",
                "The value has to be a map",
                Context::show(value),
            ))
        }
    }
}

impl ParseJson for PrimaryIonSeries {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::show(map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","))
            };
            Ok(Self {
                location: Location::from_json_value(map.remove("location").ok_or_else(|| {
                    CustomError::error(
                        "Invalid PrimaryIonSeries",
                        "The required property 'location' is missing",
                        context(&map),
                    )
                })?)?,
                neutral_losses: Vec::from_json_value(map.remove("neutral_losses").ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid PrimaryIonSeries",
                            "The required property 'neutral_losses' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                amino_acid_neutral_losses: Vec::from_json_value(
                    map.remove("amino_acid_neutral_losses").ok_or_else(|| {
                        CustomError::error(
                            "Invalid PrimaryIonSeries",
                            "The required property 'amino_acid_neutral_losses' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                amino_acid_side_chain_losses: <(u8, Option<Vec<AminoAcid>>)>::from_json_value(
                    map.remove("amino_acid_side_chain_losses").ok_or_else(|| {
                        CustomError::error(
                            "Invalid PrimaryIonSeries",
                            "The required property 'amino_acid_side_chain_losses' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                charge_range: ChargeRange::from_json_value(
                    map.remove("charge_range").ok_or_else(|| {
                        CustomError::error(
                            "Invalid PrimaryIonSeries",
                            "The required property 'charge_range' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                allowed_variants: Vec::from_json_value(
                    map.remove("allowed_variants").ok_or_else(|| {
                        CustomError::error(
                            "Invalid PrimaryIonSeries",
                            "The required property 'allowed_variants' is missing",
                            context(&map),
                        )
                    })?,
                )?,
            })
        } else {
            Err(CustomError::error(
                "Invalid Database SimpleModification",
                "The value has to be a map",
                Context::show(value),
            ))
        }
    }
}

impl ParseJson for Location {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for SatelliteLocation {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for SatelliteIonSeries {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::show(map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","))
            };
            Ok(Self {
                location: SatelliteLocation::from_json_value(map.remove("location").ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'location' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                neutral_losses: Vec::from_json_value(map.remove("neutral_losses").ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'neutral_losses' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                amino_acid_neutral_losses: Vec::from_json_value(
                    map.remove("amino_acid_neutral_losses").ok_or_else(|| {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'amino_acid_neutral_losses' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                amino_acid_side_chain_losses: <(u8, Option<Vec<AminoAcid>>)>::from_json_value(
                    map.remove("amino_acid_side_chain_losses").ok_or_else(|| {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'amino_acid_side_chain_losses' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                charge_range: ChargeRange::from_json_value(
                    map.remove("charge_range").ok_or_else(|| {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'charge_range' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                allowed_variants: Vec::from_json_value(
                    map.remove("allowed_variants").ok_or_else(|| {
                        CustomError::error(
                            "Invalid SatelliteIonSeries",
                            "The required property 'allowed_variants' is missing",
                            context(&map),
                        )
                    })?,
                )?,
            })
        } else {
            Err(CustomError::error(
                "Invalid SatelliteIonSeries",
                "The value has to be a map",
                Context::show(value),
            ))
        }
    }
}

impl ParseJson for ChargeRange {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for GlycanModel {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::show(map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","))
            };
            Ok(Self {
                allow_structural: bool::from_json_value(
                    map.remove("allow_structural").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'allow_structural' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                compositional_range: RangeInclusive::from_json_value(
                    map.remove("compositional_range").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'compositional_range' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                neutral_losses: Vec::from_json_value(map.remove("neutral_losses").ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'neutral_losses' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                specific_neutral_losses: Vec::from_json_value(
                    map.remove("specific_neutral_losses").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'specific_neutral_losses' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                default_peptide_fragment: GlycanPeptideFragment::from_json_value(
                    map.remove("default_peptide_fragment").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'default_peptide_fragment' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                peptide_fragment_rules: Vec::from_json_value(
                    map.remove("peptide_fragment_rules").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'peptide_fragment_rules' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                oxonium_charge_range: ChargeRange::from_json_value(
                    map.remove("oxonium_charge_range").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'oxonium_charge_range' is missing",
                            context(&map),
                        )
                    })?,
                )?,
                other_charge_range: ChargeRange::from_json_value(
                    map.remove("other_charge_range").ok_or_else(|| {
                        CustomError::error(
                            "Invalid GlycanModel",
                            "The required property 'other_charge_range' is missing",
                            context(&map),
                        )
                    })?,
                )?,
            })
        } else {
            Err(CustomError::error(
                "Invalid GlycanModel",
                "The value has to be a map",
                Context::show(value),
            ))
        }
    }
}

impl ParseJson for GlycanPeptideFragment {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}

impl ParseJson for FragmentKind {
    fn from_json_value(value: Value) -> Result<Self, CustomError> {
        use_serde(value)
    }
}
