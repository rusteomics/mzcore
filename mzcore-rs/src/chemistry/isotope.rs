#![allow(dead_code)]

use anyhow::*;
use serde::{Deserialize, Serialize};

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Isotope {
    pub mass_number: u16,
    pub mass: f64,
    pub abundance: f32
}

impl Isotope {
    pub fn new(mass_number: u16, mass: f64, abundance: f32) -> Result<Isotope> {
        if mass_number <= 0 { bail!("mass_number must be a strictly positive number") }
        if mass <= 0.0 { bail!("mass must be a strictly positive number") }
        if abundance < 0.0 { bail!("abundance must be a positive number") }

        Ok(Isotope {
            mass_number: mass_number,
            mass: mass,
            abundance: abundance,
        })
    }

    pub fn nucleon_number(&self) -> u16 { self.mass_number }

    pub fn get_neutron_number(&self, proton_number: u16) -> u16 {
        self.mass_number - proton_number
    }

}


