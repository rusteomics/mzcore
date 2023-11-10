#![allow(dead_code)]

use std::borrow::Cow;
use std::cmp::Ordering;
use anyhow::*;
use serde::{Deserialize, Serialize};

use crate::chemistry::api::*;
use crate::chemistry::element::Element;
use crate::chemistry::isotope::Isotope;

// The atomic_number uniquely identifies an element
#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Atom {
    pub element: Element,
    pub name: String,
    pub isotopes: Vec<Isotope>,
}

impl Atom {
    pub fn new(
        element: Element,
        name: &str,
        isotopes: Vec<Isotope>,
    ) -> Result<Atom> {
        if name.is_empty() { bail!("name is empty") }
        if isotopes.is_empty() { bail!("isotopes is empty") }

        Ok(Atom {
            element: element,
            name: name.to_string(),
            isotopes: isotopes,
        })
    }

    pub fn atomic_number(&self) -> u16 { self.element as u16 }

    pub fn proton_number(&self) -> u16 { self.atomic_number() }

    pub fn get_neutron_number(&self, isotope_idx: usize) -> u16 {
        self.isotopes[isotope_idx].get_neutron_number(self.proton_number())
    }

    pub fn monoisotopic_mass(&self) -> f64 {
        self.isotopes.first().unwrap().mass
    }

    pub fn calc_average_mass(&self) -> f64 {
        let mut weighted_mass_sum: f64 = 0.0;
        let mut weight_sum: f64 = 0.0;

        for iso in self.isotopes.iter() {
            let ab = iso.abundance as f64;
            weighted_mass_sum += iso.mass * ab;
            weight_sum += ab;
        }

        weighted_mass_sum / weight_sum
    }

    pub fn calc_most_abundant_mass(&self) -> f64 {
        let most_abundant_isitope = self.isotopes.iter()
            .max_by(|x, y| {
                x.abundance.partial_cmp(&y.abundance).unwrap_or(Ordering::Equal)
            });

        most_abundant_isitope.unwrap().mass
    }

    pub fn to_isotopic_variants(&self) -> Vec<AtomIsotopicVariant> {
        self.isotopes.iter().map(|isotope| AtomIsotopicVariant::new(self.clone(), isotope.clone())).collect()
    }

}

impl HasMass for Atom {
    fn mono_mass(&self) -> f64 { self.isotopes.first().unwrap().mass }
    fn average_mass(&self) -> Option<f64> { Some(self.calc_average_mass()) }
}

impl HasNameAndSymbol for Atom {
    fn name(&self) -> Cow<str> { Cow::Borrowed(&self.name) }
    fn symbol(&self) -> Cow<str>  { Cow::Borrowed(self.element.to_str()) }
}

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct AtomIsotopicVariant { atom: Atom, isotope: Isotope }

impl AtomIsotopicVariant {
    pub fn new(atom: Atom, isotope: Isotope) -> AtomIsotopicVariant {
        AtomIsotopicVariant {
            atom,
            isotope,
        }
    }

    pub fn mass(&self) -> f64 {
        self.isotope.mass
    }

}

/*

impl Element {
    /// Get all available isotopes (N, mass, abundance)
    pub fn isotopes(self) -> &'static [(u16, Mass, f64)] {
        &elemental_data()[self as usize].2
    }

    /// The mass of the specified isotope of this element (if that isotope exists)
    pub fn mass(&self, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        Some(if isotope == 0 {
            elemental_data()[*self as usize - 1].0?
        } else {
            // Specific isotope do not change anything
            elemental_data()[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, _, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        })
    }

    /// The average weight of the specified isotope of this element (if that isotope exists)
    pub fn average_weight(&self, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        Some(if isotope == 0 {
            elemental_data()[*self as usize - 1].1?
        } else {
            // Specific isotope do not change anything
            elemental_data()[*self as usize - 1]
                .2
                .iter()
                .find(|(ii, _, _)| *ii == isotope)
                .map(|(_, m, _)| *m)?
        })
    }

    /// Gives the most abundant mass based on the number of this isotope
    pub fn most_abundant_mass(&self, n: i16, isotope: u16) -> Option<Mass> {
        if *self == Self::Electron {
            return Some(da(5.485_799_090_65e-4) * Ratio::new::<r>(f64::from(n)));
        }
        Some(
            if isotope == 0 {
                // (mass, chance)
                let mut max = None;
                for iso in &elemental_data()[*self as usize - 1].2 {
                    let chance = iso.2 * f64::from(n);
                    if max.map_or(true, |m: (Mass, f64)| chance > m.1) {
                        max = Some((iso.1, chance));
                    }
                }
                max?.0
            } else {
                // Specific isotope do not change anything
                elemental_data()[*self as usize - 1]
                    .2
                    .iter()
                    .find(|(ii, _, _)| *ii == isotope)
                    .map(|(_, m, _)| *m)?
            } * Ratio::new::<r>(f64::from(n)),
        )
    }
}

 */