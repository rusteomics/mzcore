use std::{
    marker::PhantomData,
    path::{Path, PathBuf},
};

use crate::{
    error::CustomError,
    identification::{PeptidoformPresent, SpectrumId},
    ontology::CustomDatabase,
    sequence::{Peptidoform, SemiAmbiguous},
    system::{Mass, Ratio, Time, isize::Charge},
};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::identification::{
    BoxedIdentifiedPeptideIter, IdentifiedPeptidoform, IdentifiedPeptidoformSource,
    IdentifiedPeptidoformVersion, IdentifiedPeptidoformData,
    common_parser::Location,
    csv::{CsvLine, parse_csv},
};

static NUMBER_ERROR: (&str, &str) = (
    "Invalid Sage line",
    "This column is not a number but it is required to be a number in this Sage format",
);

format_family!(
    /// The format for any Sage file
    SageFormat,
    /// The data from any Sage file
    SageData,
    SemiAmbiguous, PeptidoformPresent, SageVersion, [&VERSION_0_14], b'\t', None;
    required {
        aligned_rt: Ratio, |location: Location, _| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::fraction>);
        decoy: bool, |location: Location, _| location.parse::<i8>(NUMBER_ERROR).map(|v| v == -1);
        delta_best: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        delta_mobility: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        delta_next: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        delta_rt_model: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        /// Experimental mass
        mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        fragment_ppm: Ratio, |location: Location, _| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::ppm>);
        hyperscore: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        ion_mobility: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        isotope_error: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        longest_b: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        longest_y: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        matched_intensity_pct: Ratio, |location: Location, _| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::percent>);
        matched_peaks: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        missed_cleavages: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        ms2_intensity: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        scan: SpectrumId, |location: Location, _|Ok(SpectrumId::Native(location.get_string()));
        peptide_q: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        peptide: Peptidoform<SemiAmbiguous>, |location: Location, custom_database: Option<&CustomDatabase>| Peptidoform::pro_forma(location.as_str(), custom_database).map(|p|p.into_semi_ambiguous().unwrap());
        poisson: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        posterior_error: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        predicted_mobility: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        predicted_rt: Ratio, |location: Location, _| location.parse(NUMBER_ERROR).map(Ratio::new::<crate::system::ratio::fraction>);
        protein_q: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        proteins: Vec<String>, |location: Location, _| Ok(location.get_string().split(';').map(ToString::to_string).collect_vec());
        /// PSM ID
        id: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        rank: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        raw_file: PathBuf, |location: Location, _| Ok(Path::new(&location.get_string()).to_owned());
        rt: Time, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Time::new::<crate::system::time::min>);
        sage_discriminant_score: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        scored_candidates: usize, |location: Location, _| location.parse(NUMBER_ERROR);
        semi_enzymatic: bool, |location: Location, _| location.parse::<u8>(NUMBER_ERROR).map(|n| n != 0);
        spectrum_q: f64, |location: Location, _| location.parse(NUMBER_ERROR);
        theoretical_mass: Mass, |location: Location, _| location.parse::<f64>(NUMBER_ERROR).map(Mass::new::<crate::system::dalton>);
        z: Charge, |location: Location, _| location.parse::<isize>(NUMBER_ERROR).map(Charge::new::<crate::system::e>);
    }
    optional { }
);

impl From<SageData> for IdentifiedPeptidoform<SemiAmbiguous, PeptidoformPresent> {
    fn from(value: SageData) -> Self {
        Self {
            score: Some(value.sage_discriminant_score.clamp(-1.0, 1.0)),
            local_confidence: None,
            metadata: IdentifiedPeptidoformData::Sage(value),
            complexity_marker: PhantomData,
            peptidoform_availability_marker: PhantomData,
        }
    }
}

/// An older version of a Sage export
pub const VERSION_0_14: SageFormat = SageFormat {
    version: SageVersion::V0_14,
    id: "psm_id",
    peptide: "peptide",
    proteins: "proteins",
    raw_file: "filename",
    scan: "scannr",
    rank: "rank",
    decoy: "label",
    mass: "expmass",
    theoretical_mass: "calcmass",
    z: "charge",
    missed_cleavages: "missed_cleavages",
    semi_enzymatic: "semi_enzymatic",
    isotope_error: "isotope_error",
    fragment_ppm: "fragment_ppm",
    hyperscore: "hyperscore",
    delta_next: "delta_next",
    delta_best: "delta_best",
    rt: "rt",
    aligned_rt: "aligned_rt",
    predicted_rt: "predicted_rt",
    delta_rt_model: "delta_rt_model",
    ion_mobility: "ion_mobility",
    predicted_mobility: "predicted_mobility",
    delta_mobility: "delta_mobility",
    matched_peaks: "matched_peaks",
    longest_b: "longest_b",
    longest_y: "longest_y",
    matched_intensity_pct: "matched_intensity_pct",
    scored_candidates: "scored_candidates",
    poisson: "poisson",
    sage_discriminant_score: "sage_discriminant_score",
    posterior_error: "posterior_error",
    spectrum_q: "spectrum_q",
    peptide_q: "peptide_q",
    protein_q: "protein_q",
    ms2_intensity: "ms2_intensity",
};

/// All possible Sage versions
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum SageVersion {
    /// Current sage version
    #[default]
    V0_14,
}

impl std::fmt::Display for SageVersion {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.name())
    }
}

impl IdentifiedPeptidoformVersion<SageFormat> for SageVersion {
    fn format(self) -> SageFormat {
        match self {
            Self::V0_14 => VERSION_0_14,
        }
    }
    fn name(self) -> &'static str {
        match self {
            Self::V0_14 => "v0.14",
        }
    }
}
