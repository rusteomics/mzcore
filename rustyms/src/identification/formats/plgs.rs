use std::{borrow::Cow, marker::PhantomData, ops::Range};

use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    chemistry::MolecularFormula,
    fragment::NeutralLoss,
    helper_functions::explain_number_error,
    identification::{
        BoxedIdentifiedPeptideIter, FastaIdentifier, FlankingSequence, IdentifiedPeptidoform,
        IdentifiedPeptidoformData, IdentifiedPeptidoformSource, IdentifiedPeptidoformVersion,
        KnownFileFormat, MetaData, PeptidoformPresent, SpectrumId, SpectrumIds,
        common_parser::{Location, OptionalColumn, OptionalLocation},
        csv::{CsvLine, parse_csv},
    },
    ontology::CustomDatabase,
    prelude::CompoundPeptidoformIon,
    sequence::{
        AminoAcid, MUPSettings, Modification, Peptidoform, PlacementRule, Position,
        SequencePosition, SimpleLinear, SimpleModification,
    },
    system::{Mass, MassOverCharge, OrderedTime, Time, isize::Charge},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid PLGS line",
    "This column is not a number but it is required to be a number in this PLGS format",
);
static IDENTIFIER_ERROR: (&str, &str) = (
    "Invalid PLGS line",
    "This column is not a valid identifier but is required to be in this PLGS format",
);
static CURATION_ERROR: (&str, &str) = (
    "Invalid PLGS line",
    "This column is not a curation but it is required to be Green, Yellow, or Red",
);

format_family!(
    PLGS,
    SimpleLinear, PeptidoformPresent, [&VERSION_3_0], b',', None;
    required {
        protein_id: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_entry: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        protein_accession: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        protein_description: FastaIdentifier<String>, |location: Location, _| location.parse(IDENTIFIER_ERROR);
        protein_db_type: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        protein_score: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_fpr: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_average_weight: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        protein_matched_products: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptides: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_digest_peptides: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_sequence_coverage: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptide_intensity_sum: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_peptide_intensity_top3: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_matched_product_intensity_sum: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        protein_fmol_on_column: Option<f32>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
        protein_ngram_on_column: Option<f32>, |location: Location, _| location.or_empty().parse(NUMBER_ERROR);
        protein_auto_curate: PLGSCuration, |location: Location, _| location.parse(CURATION_ERROR);
        peptide_rank: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_pass: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        peptide_match_type: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        peptide_modifications: ThinVec<(SimpleModification, AminoAcid, Option<usize>)>, |location: Location, custom_database: Option<&CustomDatabase>|
            location.ignore("None").array(';').map(|l| {
                let plus = l.as_str().find('+').ok_or_else(|| BoxedError::error("Invalid PLGS modification", "A PLGS modification should be in the format 'modification+AA(pos)' and the plus '+' is missing.", l.context().to_owned()))?;
                let modification = Modification::sloppy_modification(l.full_line(), l.location.start..l.location.start+plus, None, custom_database).map_err(BoxedError::to_owned)?;
                let aa = l.as_str()[plus+1..plus+2].parse::<AminoAcid>().map_err(|()| BoxedError::error("Invalid PLGS modification", "A PLGS modification should be in the format 'modification+AA(pos)' and the amino acid is not valid", l.context().to_owned()))?;
                let num = &l.as_str()[plus+3..l.len()-1];
                let index = if num == "*" {None} else {
                    Some(num.parse::<usize>().map_err(|err| BoxedError::error("Invalid PLGS modification", format!("A PLGS modification should be in the format 'modification+AA(pos)' and the pos is {}", explain_number_error(&err)), l.context().to_owned()))? - 1)
                };
                Ok((modification, aa, index))
            }).collect::<Result<ThinVec<_>,_>>();
        peptide: Peptidoform<SimpleLinear>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::pro_forma(location.as_str(), custom_database).map(|p|p.into_simple_linear().unwrap()).map_err(BoxedError::to_owned);
        peptide_start: u16, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_pi: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_component_id: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_products: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_unique_products: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_consecutive_products: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_complementary_products: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_raw_score: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_x_p_bond_identified: Option<bool>, |location: Location, _| Ok(location.or_empty().map(|l| l.as_str() == "Identified"));
        peptide_matched_product_intensity: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_product_theoretical: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_matched_product_string: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        peptide_model_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        peptide_volume: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_csa: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_model_drift: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_relative_intensity: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide_auto_curate: PLGSCuration, |location: Location, _| location.parse(CURATION_ERROR);
        precursor_le_id: u32, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_mass: Mass, |location: Location, _| location.parse(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        precursor_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        precursor_intensity: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_charge: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_z: Charge, |location: Location, _| location.parse(NUMBER_ERROR).map(Charge::new::<crate::system::charge::e>);
        precursor_mz: MassOverCharge, |location: Location, _| location.parse(NUMBER_ERROR).map(MassOverCharge::new::<crate::system::mass_over_charge::mz>);
        precursor_fwhm: f32, |location: Location, _| location.parse(NUMBER_ERROR);
        precursor_lift_off_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::s>);
        precursor_inf_up_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::s>);
        precursor_inf_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::s>);
        precursor_touch_down_rt: Time, |location: Location, _| location.parse(NUMBER_ERROR).map(Time::new::<crate::system::time::s>);
        precursor_rms_fwhm_delta: f32, |location: Location, _| location.parse(NUMBER_ERROR);
    }
    optional {
        fragment_mass: Mass, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Mass::new::<crate::system::dalton>));
        fragment_type: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        fragment_index: u32, |location: Location, _| location.or_empty().parse::<u32>(NUMBER_ERROR);
        fragment_neutral_loss: NeutralLoss, |location: Location, _| location.or_empty().ignore("None").map(|l| MolecularFormula::from_pro_forma(l.full_line(), l.location.clone(), false, false, false, true).map(|f| NeutralLoss::Loss(1, f)).map_err(BoxedError::to_owned)).transpose();
        fragment_description: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        fragment_sequence: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        fragment_site: Box<str>, |location: Location, _| Ok(location.get_boxed_str());
        product_rank: isize, |location: Location, _| location.parse::<isize>(NUMBER_ERROR);
        product_he_id: u32, |location: Location, _| location.or_empty().parse::<u32>(NUMBER_ERROR);
        product_mass: Mass, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Mass::new::<crate::system::dalton>));
        product_mz: MassOverCharge, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(MassOverCharge::new::<crate::system::mass_over_charge::mz>));
        product_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::min>));
        product_intensity: usize, |location: Location, _| location.or_empty().parse::<usize>(NUMBER_ERROR);
        product_charge: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        product_z: Charge, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Charge::new::<crate::system::charge::e>));
        product_fwhm: f32, |location: Location, _| location.or_empty().parse::<f32>(NUMBER_ERROR);
        product_lift_off_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::s>));
        product_inf_up_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::s>));
        product_inf_down_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::s>));
        product_touch_down_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::s>));
        precursor_product_delta_rt: Time, |location: Location, _| location.or_empty().parse(NUMBER_ERROR).map(|r| r.map(Time::new::<crate::system::time::s>));
    }

    fn post_process(_source: &CsvLine, mut parsed: Self, _custom_database: Option<&CustomDatabase>) -> Result<Self, BoxedError<'static>> {
        for (m, aa, index) in &parsed.peptide_modifications {
            if let Some(index) = index {
                parsed.peptide.add_simple_modification(SequencePosition::Index(*index), m.clone());
            } else if !parsed.peptide.add_unknown_position_modification(m.clone(), .., &MUPSettings{position: Some(vec![PlacementRule::AminoAcid(vec![*aa], Position::Anywhere)]), .. Default::default()})
            {
                return Err(BoxedError::error(
                        "Modification of unknown position cannot be placed",
                        "There is no position where this ambiguous modification can be placed based on the placement rules in the database.",
                        Context::show(m.to_string()),
                    ));
            }
        }
        Ok(parsed)
    }
);

/// PLGS curation categories
#[derive(Clone, Copy, Debug, Default, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize)]
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
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum PLGSVersion {
    /// Current PLGS version
    #[default]
    V3_0,
}

impl std::fmt::Display for PLGSVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<PLGSFormat> for PLGSVersion {
    fn format(self) -> PLGSFormat {
        match self {
            Self::V3_0 => VERSION_3_0,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V3_0 => "v3.0",
        }
    }
}

impl MetaData for PLGSData {
    fn compound_peptidoform_ion(&self) -> Option<Cow<'_, CompoundPeptidoformIon>> {
        Some(Cow::Owned(self.peptide.clone().into()))
    }

    fn format(&self) -> KnownFileFormat {
        KnownFileFormat::PLGS(self.version)
    }

    fn id(&self) -> String {
        self.peptide_component_id.to_string()
    }

    fn confidence(&self) -> Option<f64> {
        Some(2.0 / (1.0 + 1.3_f64.powf(-self.peptide_score)) - 1.0)
    }

    fn local_confidence(&self) -> Option<Cow<'_, [f64]>> {
        None
    }

    fn original_confidence(&self) -> Option<f64> {
        Some(self.peptide_score)
    }

    fn original_local_confidence(&self) -> Option<&[f64]> {
        None
    }

    fn charge(&self) -> Option<Charge> {
        Some(self.precursor_z)
    }

    fn mode(&self) -> Option<&str> {
        None
    }

    fn retention_time(&self) -> Option<Time> {
        Some(self.precursor_rt)
    }

    fn scans(&self) -> SpectrumIds {
        SpectrumIds::FileNotKnown(vec![SpectrumId::RetentionTime(
            OrderedTime::from(self.precursor_lift_off_rt)
                ..=OrderedTime::from(self.precursor_touch_down_rt),
        )])
    }

    fn experimental_mz(&self) -> Option<MassOverCharge> {
        Some(self.precursor_mz)
    }

    fn experimental_mass(&self) -> Option<Mass> {
        Some(self.precursor_mass)
    }

    fn protein_name(&self) -> Option<FastaIdentifier<String>> {
        Some(self.protein_description.clone())
    }

    fn protein_id(&self) -> Option<usize> {
        Some(self.protein_id)
    }

    fn protein_location(&self) -> Option<Range<u16>> {
        Some(self.peptide_start..self.peptide_start + self.peptide.len() as u16)
    }

    fn flanking_sequences(&self) -> (&FlankingSequence, &FlankingSequence) {
        (&FlankingSequence::Unknown, &FlankingSequence::Unknown)
    }

    fn database(&self) -> Option<(&str, Option<&str>)> {
        None
    }
}
