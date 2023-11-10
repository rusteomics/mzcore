
use serde::{Deserialize, Serialize};

use crate::chemistry::constants::*;
use crate::chemistry::element::Element;
use crate::chemistry::peptide::LinearPeptide;
use crate::chemistry::table::biomolecule_atom_table;

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum ActivationType {
    CID,
    ECD,
    ETD,
    HCD,
    PSD,
}

impl std::fmt::Display for ActivationType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum MsAnalyzer {
    FTMS,
    TRAP,
}

impl std::fmt::Display for MsAnalyzer {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Enumeration that defines ion series direction
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum FragmentIonSeriesDirection {
    NTerminal,
    CTerminal,
    Unspecified,
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum FragmentIonSeries {
    //M,
    //M_H2O,
    //M_NH3,
    a,
    //#[strum(serialize = "a-H2O")]
    a_H2O,
    //#[strum(serialize = "a-NH3")]
    a_NH3,
    b,
    //#[strum(serialize = "b-H2O")]
    b_H2O,
    //#[strum(serialize = "b-NH3")]
    b_NH3,
    c,
    c·,
    c_m1,
    c_p1,
    c_p2,
    c_H2O,
    c_NH3,
    d,
    v,
    w,
    x,
    x_H2O,
    x_NH3,
    y,
    //#[strum(serialize = "y-H2O")]
    y_H2O,
   //#[strum(serialize = "y-NH3")]
    y_NH3,
    ya,
    yb,
    z,
    z_H2O,
    z_NH3,
    z·,
    z_p1,
    z_p2,
    z_p3,
    immonium,
}

impl FragmentIonSeries {

    pub fn get_ion_mono_mass_shift(&self) -> f64 {

        // TODO: define as constants?
        let atom_table = biomolecule_atom_table();
        let h_mono_mass = atom_table.atom_by_element.get(&Element::H).unwrap().monoisotopic_mass();
        let n_mono_mass = atom_table.atom_by_element.get(&Element::N).unwrap().monoisotopic_mass();
        let o_mono_mass = atom_table.atom_by_element.get(&Element::O).unwrap().monoisotopic_mass();

        match self {
            //M => 0.0,
            //M_H2O => - WATER_MONO_MASS,
            //M_NH3 => - NH3_MONO_MASS,
            Self::a =>  - (WATER_MONO_MASS + CO_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1'),
            Self::a_H2O => - (2.0 * WATER_MONO_MASS + CO_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1' + 'H-2O-1'),
            Self::a_NH3 => - (CO_MONO_MASS + NH3_MONO_MASS + WATER_MONO_MASS), //Composition(formula='H-2O-1' + 'C-1O-1' + 'N-1H-3'),
            Self::b => - WATER_MONO_MASS, //Composition(formula='H-2O-1'),
            Self::b_H2O => -2.0 * WATER_MONO_MASS, //Composition(formula='H-2O-1' + 'H-2O-1'),
            Self::b_NH3 => - (WATER_MONO_MASS + NH3_MONO_MASS), //Composition(formula='H-2O-1' + 'N-1H-3'),
            Self::c => - WATER_MONO_MASS + NH3_MONO_MASS,
            Self::c· => - WATER_MONO_MASS + NH3_MONO_MASS + PROTON_MASS,
            Self::c_m1 => - WATER_MONO_MASS + NH3_MONO_MASS - PROTON_MASS,
            Self::c_p1 => - WATER_MONO_MASS + NH3_MONO_MASS + PROTON_MASS,
            Self::c_p2 => - WATER_MONO_MASS + NH3_MONO_MASS + 2.0 * PROTON_MASS,
            Self::c_H2O => - 2.0 * WATER_MONO_MASS + NH3_MONO_MASS,
            Self::c_NH3 => - WATER_MONO_MASS,
            Self::d => 0.0, // FIXME
            Self::v => 0.0, // FIXME
            Self::w => 0.0, // FIXME
            Self::x => - WATER_MONO_MASS + CO2_MONO_MASS,
            Self::x_H2O => - 2.0 * WATER_MONO_MASS + CO2_MONO_MASS,
            Self::x_NH3 => - WATER_MONO_MASS - NH3_MONO_MASS + CO2_MONO_MASS,
            Self::y => 0.0,
            Self::y_H2O => - WATER_MONO_MASS,
            Self::y_NH3 => - NH3_MONO_MASS,
            Self::ya => 0.0, // FIXME
            Self::yb => 0.0, // FIXME
            Self::z => - WATER_MONO_MASS + o_mono_mass - n_mono_mass,
            Self::z_H2O => - 2.0 * WATER_MONO_MASS + o_mono_mass - n_mono_mass - h_mono_mass,
            Self::z_NH3 => - WATER_MONO_MASS - NH3_MONO_MASS + o_mono_mass - n_mono_mass - h_mono_mass,
            Self::z· => - WATER_MONO_MASS + o_mono_mass - n_mono_mass, // TODO: check me
            Self::z_p1 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + h_mono_mass,
            Self::z_p2 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + 2.0 * h_mono_mass,
            Self::z_p3 => - WATER_MONO_MASS + o_mono_mass - n_mono_mass + 3.0 * h_mono_mass,
            Self::immonium => 0.0 // FIXME
        }
    }

    pub fn is_n_terminal(&self) -> Option<bool> {
        use FragmentIonSeriesDirection as SeriesDir;

        match self.get_ion_series_direction() {
            SeriesDir::NTerminal => Some(true),
            SeriesDir::CTerminal => Some(false),
            SeriesDir::Unspecified => None
        }
    }

    pub fn get_ion_series_direction(&self) -> FragmentIonSeriesDirection {

        //use FragmentIonSeries::*;
        use FragmentIonSeriesDirection as SeriesDir;

        match self {
            Self::a => SeriesDir::NTerminal,
            Self::a_H2O => SeriesDir::NTerminal,
            Self::a_NH3 => SeriesDir::NTerminal,
            Self::b => SeriesDir::NTerminal,
            Self::b_H2O => SeriesDir::NTerminal,
            Self::b_NH3 => SeriesDir::NTerminal,
            Self::c => SeriesDir::NTerminal,
            Self::c· => SeriesDir::NTerminal,
            Self::c_m1 => SeriesDir::NTerminal,
            Self::c_p1 => SeriesDir::NTerminal,
            Self::c_p2 => SeriesDir::NTerminal,
            Self::c_H2O => SeriesDir::NTerminal,
            Self::c_NH3 => SeriesDir::NTerminal,
            Self::d => SeriesDir::Unspecified, // FIXME
            Self::v => SeriesDir::Unspecified, // FIXME
            Self::w => SeriesDir::Unspecified, // FIXME
            Self::x => SeriesDir::CTerminal,
            Self::x_H2O => SeriesDir::CTerminal,
            Self::x_NH3 => SeriesDir::CTerminal,
            Self::y => SeriesDir::CTerminal,
            Self::y_H2O => SeriesDir::CTerminal,
            Self::y_NH3 => SeriesDir::CTerminal,
            Self::ya => SeriesDir::CTerminal,
            Self::yb => SeriesDir::CTerminal,
            Self::z => SeriesDir::CTerminal,
            Self::z_H2O => SeriesDir::CTerminal,
            Self::z_NH3 => SeriesDir::CTerminal,
            Self::z· => SeriesDir::CTerminal,
            Self::z_p1 => SeriesDir::CTerminal,
            Self::z_p2 => SeriesDir::CTerminal,
            Self::z_p3 => SeriesDir::CTerminal,
            Self::immonium => SeriesDir::Unspecified,
        }
    }
}

impl std::fmt::Display for FragmentIonSeries {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        use FragmentIonSeries::*;

        // TODO: update me
        match self {
            a => write!(f, "a"),
            a_NH3 => write!(f, "a-NH3"),
            a_H2O => write!(f, "a-H2O"),
            b => write!(f, "b"),
            b_NH3 => write!(f, "b-NH3"),
            b_H2O => write!(f, "a-H2O"),
            c => write!(f, "c"),
            c· => write!(f, "c·"),
            c_m1 => write!(f, "c-1"),
            c_p1 => write!(f, "c+1"),
            c_p2 => write!(f, "c+1"),
            c_NH3 => write!(f, "c-NH3"),
            c_H2O => write!(f, "c-H2O"),
            d => write!(f, "d"),
            v => write!(f, "v"),
            w => write!(f, "w"),
            x => write!(f, "x"),
            x_NH3 => write!(f, "x-NH3"),
            x_H2O => write!(f, "x-H2O"),
            y => write!(f, "y"),
            y_NH3 => write!(f, "y-NH3"),
            y_H2O => write!(f, "y-H2O"),
            ya => write!(f, "ya"),
            yb => write!(f, "yb"),
            z => write!(f, "z"),
            z_NH3 => write!(f, "z-NH3"),
            z_H2O => write!(f, "z-H2O"),
            z· => write!(f, "z·"),
            z_p1 => write!(f, "z+1"),
            z_p2 => write!(f, "z+2"),
            z_p3 => write!(f, "z+3"),
            immonium => write!(f, "immonium"),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum FragmentType {
    Immonium,
    Internal,
    Satellite,
    Sequence,
}

impl std::fmt::Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {

        use FragmentType::*;

        match self {
            Immonium => write!(f, "IM"),
            Internal => write!(f, "IN"),
            Satellite => write!(f, "SAT"),
            Sequence => write!(f, "SEQ"),
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct FragmentIonType {
    pub ion_series: FragmentIonSeries,
    pub neutral_loss: Option<NeutralLoss>,
    pub is_forward_ion: bool
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
struct FragmentMatch {
    pub label: String,
    pub r#type: Option<String>, // = None,
    pub moz: f64,
    pub calculated_moz: f64,
    pub intensity: f32,
    pub neutral_loss_mass: Option<f64> // = None
}

#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
pub enum NeutralLoss {
    CH4OS,
    H2O,
    H3PO4,
    HPO3,
    NH3,
}

impl std::fmt::Display for NeutralLoss {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

/*
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
struct TheoreticalFragmentIon {
    pub position: i32,
    pub moz: f64,
    pub nl_mass: f64, // = 0.0,
    pub charge: i8, //= 1,
    pub fragment_type: FragmentType, // = FragmentType.SEQUENCE,
    pub frag_series: Option<String>, // = None,
    pub is_reverse_series: bool,
    pub is_alternative_nl: bool // = false
}*/

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
struct FragmentIonTable {
    peptide: LinearPeptide,
    charge_by_ion_type: std::collections::HashMap<char,i8>,
    //sequence: Option<Vec<char>>, // = None,
    ptm_neutral_losses: Option<std::collections::HashMap<i32,f64>> // NL mass by seq position
}
