use crate::{
    error::CustomError,
    helper_functions::explain_number_error,
    modification::SimpleModification,
    ontologies::CustomDatabase,
    peptide::SimpleLinear,
    system::{usize::Charge, Mass, MassOverCharge, Time},
    AminoAcid, LinearPeptide, Modification, MolecularFormula,
};
use serde::{Deserialize, Serialize};

use super::{
    common_parser::{Location, OptionalColumn, OptionalLocation},
    csv::{parse_csv, CsvLine},
    placement_rule::PlacementRule,
    BoxedIdentifiedPeptideIter, IdentifiedPeptide, IdentifiedPeptideSource, MetaData,
    SequencePosition,
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid PLGS line",
    "This column is not a number but it is required to be a number in this PLGS format",
);
static CURATION_ERROR: (&str, &str) = (
    "Invalid PLGS line",
    "This column is not a curation but it is required to be Green, Yellow, or Red",
);

format_family!(
    /// The format for any PLGS file
    PLGSFormat,
    /// The data from any PLGS file
    PLGSData,
    PLGSVersion, [&VERSION_3_0], b',', None;
    required {
        protein_id: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_entry: String, |location: Location, _| Ok(location.get_string());
        protein_accession: String, |location: Location, _| Ok(location.get_string());
        protein_description: String, |location: Location, _| Ok(location.get_string());
        protein_db_type: String, |location: Location, _| Ok(location.get_string());
        protein_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_fpr: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_average_weight: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        protein_matched_products: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptides: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_digest_peptides: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_sequence_coverage: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptide_intensity_sum: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptide_intensity_top3: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_product_intensity_sum: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_fmol_on_column: Option<f64>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
        protein_ngram_on_column: Option<f64>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
        protein_auto_curate: PLGSCuration, |location: Location, _| location.parse(CURATION_ERROR);
        peptide_rank: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_pass: String, |location: Location, _| Ok(location.get_string());
        peptide_match_type: String, |location: Location, _| Ok(location.get_string());
        peptide_modifications: Vec<(SimpleModification, AminoAcid, Option<usize>)>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.ignore("None").array(';').map(|l| {
                let plus = l.as_str().find('+').ok_or_else(|| CustomError::error("Invalid PLGS modification", "A PLGS modification should be in the format 'modification+AA(pos)' and the plus '+' is missing.", l.context()))?;
                let modification = Modification::sloppy_modification(l.full_line(), l.location.start..l.location.start+plus, None, custom_database)?;
                let aa = l.as_str()[plus+1..plus+2].parse::<AminoAcid>().map_err(|()| CustomError::error("Invalid PLGS modification", "A PLGS modification should be in the format 'modification+AA(pos)' and the amino acid is not valid", l.context()))?;
                let num = &l.as_str()[plus+2..l.len()-1];
                let index = if num == "*" {None} else {
                    Some(num.parse::<usize>().map_err(|err| CustomError::error("Invalid PLGS modification", format!("A PLGS modification should be in the format 'modification+AA(pos)' and the pos is {}", explain_number_error(&err)), l.context()))? - 1)
                };
                Ok((modification, aa, index))
            }).collect::<Result<Vec<_>,_>>();
        peptide: LinearPeptide<SimpleLinear>, |location: Location, custom_database: Option<&CustomDatabase>| LinearPeptide::pro_forma(location.as_str(), custom_database).map(|p|p.into_simple_linear().unwrap());
        peptide_start: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_pi: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_component_id: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_products: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_unique_products: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_consecutive_products: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_complementary_products: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_raw_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_x_p_bond_identified: Option<bool>, |location: Location, _| Ok(location.or_empty().map(|l| l.as_str() == "Identified"));
        peptide_matched_product_intensity: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_product_theoretical: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_product_string: String, |location: Location, _| Ok(location.get_string());
        peptide_model_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        peptide_volume: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_csa: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_model_drift: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_relative_intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_auto_curate: PLGSCuration, |location: Location, _| location.parse(CURATION_ERROR);
        precursor_le_id: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_mass: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        precursor_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_charge: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_z: Charge, |location: Location, _| location.parse(NUMBER_ERROR).map(Charge::new::<crate::system::charge::e>);
        precursor_mz: MassOverCharge, |location: Location, _| location.parse(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mass_over_charge::mz>);
        precursor_fwhm: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_lift_off_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_inf_up_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_inf_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_touch_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_rms_fwhm_delta: f64, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        fragment_mass: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        fragment_type: String, |location: Location, _| Ok(location.get_string());
        fragment_index: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        fragment_neutral_loss: MolecularFormula, |location: Location, _| location.ignore("None").map(|l| MolecularFormula::from_pro_forma(l.full_line(), l.location.clone(), false, false)).transpose();
        fragment_description: String, |location: Location, _| Ok(location.get_string());
        fragment_sequence: String, |location: Location, _| Ok(location.get_string());
        fragment_site: String, |location: Location, _| Ok(location.get_string());
        product_rank: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        product_he_id: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        product_mass: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        product_mz: MassOverCharge, |location: Location, _| location.parse(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mass_over_charge::mz>);
        product_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        product_intensity: usize, |location: Location, _| location.parse::<usize>(NUMBER_ERROR);
        product_charge: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        product_z: Charge, |location: Location, _| location.parse(NUMBER_ERROR).map(Charge::new::<crate::system::charge::e>);
        product_fwhm: f64, |location: Location, _| location.parse::<f64>(NUMBER_ERROR);
        product_lift_off_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        product_inf_up_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        product_inf_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        product_touch_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_product_delta_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, CustomError> {
        for (m, aa, index) in &parsed.peptide_modifications {
            if let Some(index) = index {
                parsed.peptide.add_simple_modification(SequencePosition::Index(*index), m.clone());
            } else if !parsed.peptide.add_unknown_position_modification(m.clone(), Some(&[PlacementRule::AminoAcid(vec![*aa], crate::placement_rule::Position::Anywhere)]), ..)
            {
                return Err(CustomError::error(
                    "Modification of unknown position cannot be placed",
                    "There is no position where this ambiguous modification can be placed based on the placement rules in the database.",
                    crate::error::Context::show(m),
                    ));
            }
        }
        Ok(parsed)
    }
);

/// PLGS curation categories
#[derive(Debug, Clone, Copy, PartialEq, Eq, Ord, Default, PartialOrd, Serialize, Deserialize)]
pub enum PLGSCuration {
    /// Good
    #[default]
    Green,
    /// Mediocre
    Yellow,
    /// Bad
    Red,
}

impl std::str::FromStr for PLGSCuration {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_ascii_lowercase().as_str() {
            "green" => Ok(Self::Green),
            "yellow" => Ok(Self::Yellow),
            "red" => Ok(Self::Red),
            _ => Err(()),
        }
    }
}

impl std::fmt::Display for PLGSCuration {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Green => "Green",
                Self::Yellow => "Yellow",
                Self::Red => "Red",
            }
        )
    }
}

impl From<PLGSData> for IdentifiedPeptide {
    fn from(value: PLGSData) -> Self {
        Self {
            score: Some(value.peptide_score.clamp(-1.0, 1.0)),
            local_confidence: None,
            metadata: MetaData::PLGS(value),
        }
    }
}

/// An older version of a PLGS export
pub const VERSION_3_0: PLGSFormat = PLGSFormat {
    version: PLGSVersion::V3_0,
    protein_id: "protein.key",
    protein_entry: "protein.entry",
    protein_accession: "protein.accession",
    protein_description: "protein.description",
    protein_db_type: "protein.databasetype",
    protein_score: "protein.score",
    protein_fpr: "protein.falsepositiverate",
    protein_average_weight: "protein.avgmass",
    protein_matched_products: "protein.matchedproducts",
    protein_matched_peptides: "protein.matchedpeptides",
    protein_digest_peptides: "protein.digestpeps",
    protein_sequence_coverage: "protein.seqcover(%)",
    protein_matched_peptide_intensity_sum: "protein.matchedpeptideintensum",
    protein_matched_peptide_intensity_top3: "protein.top3matchedpeptideintensum",
    protein_matched_product_intensity_sum: "protein.matchedproductintensum",
    protein_fmol_on_column: "protein.fmoloncolumn",
    protein_ngram_on_column: "protein.ngramoncolumn",
    protein_auto_curate: "protein.autocurate",
    peptide_rank: "peptide.rank",
    peptide_pass: "peptide.pass",
    peptide_match_type: "peptide.matchtype",
    peptide_modifications: "peptide.modification",
    peptide: "peptide.seq",
    peptide_start: "peptide.seqstart",
    peptide_pi: "peptide.pi",
    peptide_component_id: "peptide.componentid",
    peptide_matched_products: "peptide.matchedproducts",
    peptide_unique_products: "peptide.uniqueproducts",
    peptide_consecutive_products: "peptide.consectivematchedproducts",
    peptide_complementary_products: "peptide.complementarymatchedproducts",
    peptide_raw_score: "peptide.rawscore",
    peptide_score: "peptide.score",
    peptide_x_p_bond_identified: "peptide.(x)-p bond",
    peptide_matched_product_intensity: "peptide.matchedproductssuminten",
    peptide_matched_product_theoretical: "peptide.matchedproductstheoretical",
    peptide_matched_product_string: "peptide.matchedproductsstring",
    peptide_model_rt: "peptide.modelrt",
    peptide_volume: "peptide.volume",
    peptide_csa: "peptide.csa",
    peptide_model_drift: "peptide.modeldrift",
    peptide_relative_intensity: "peptide.relintensity",
    peptide_auto_curate: "peptide.autocurate",
    precursor_le_id: "precursor.leid",
    precursor_mass: "precursor.mhp",
    precursor_rt: "precursor.rett",
    precursor_intensity: "precursor.inten",
    precursor_charge: "precursor.charge",
    precursor_z: "precursor.z",
    precursor_mz: "precursor.mz",
    precursor_fwhm: "precursor.fwhm",
    precursor_lift_off_rt: "precursor.liftoffrt",
    precursor_inf_up_rt: "precursor.infuprt",
    precursor_inf_down_rt: "precursor.infdownrt",
    precursor_touch_down_rt: "precursor.touchdownrt",
    precursor_rms_fwhm_delta: "prec.rmsfwhmdelta",
    fragment_mass: OptionalColumn::Optional("fragment.mhp"),
    fragment_type: OptionalColumn::Optional("fragment.fragmenttype"),
    fragment_index: OptionalColumn::Optional("fragment.fragind"),
    fragment_neutral_loss: OptionalColumn::Optional("neutral.losstype"),
    fragment_description: OptionalColumn::Optional("fragment.str"),
    fragment_sequence: OptionalColumn::Optional("fragment.seq"),
    fragment_site: OptionalColumn::Optional("fragment.fragsite"),
    product_rank: OptionalColumn::Optional("product.rank"),
    product_he_id: OptionalColumn::Optional("product.heid"),
    product_mass: OptionalColumn::Optional("product.mhp"),
    product_mz: OptionalColumn::Optional("product.m_z"),
    product_rt: OptionalColumn::Optional("product.rett"),
    product_intensity: OptionalColumn::Optional("product.inten"),
    product_charge: OptionalColumn::Optional("product.charge"),
    product_z: OptionalColumn::Optional("product.z"),
    product_fwhm: OptionalColumn::Optional("product.fwhm"),
    product_lift_off_rt: OptionalColumn::Optional("product.liftoffrt"),
    product_inf_up_rt: OptionalColumn::Optional("product.infuprt"),
    product_inf_down_rt: OptionalColumn::Optional("product.infdownrt"),
    product_touch_down_rt: OptionalColumn::Optional("product.touchdownrt"),
    precursor_product_delta_rt: OptionalColumn::Optional("precursorproduct.deltarett"),
};

/// All possible PLGS versions
#[derive(Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize)]
#[allow(non_camel_case_types)]
pub enum PLGSVersion {
    /// Current PLGS version
    #[default]
    V3_0,
}

impl std::fmt::Display for PLGSVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(
            f,
            "{}",
            match self {
                Self::V3_0 => "v3.0",
            }
        )
    }
}
