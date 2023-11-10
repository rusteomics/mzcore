use std::borrow::Cow;
///
/// Some parts of this file originates from [Sage](https://github.com/lazear/sage/blob/master/crates/sage/src/spectrum.rs)
/// Copyright (c) 2022 Michael Lazear
/// SPDX-License-Identifier: MIT
///
use serde::{Serialize, Deserialize};

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SpectrumData {
    pub mz_list: Vec<f64>,
    pub intensity_list: Vec<f32>,
}

pub trait HasSpectrumData {
    fn get_mz_list(&self) -> Cow<[f64]>;
    fn get_intensity_list(&self) -> Cow<[f32]>;
}

impl HasSpectrumData for SpectrumData {
    fn get_mz_list(&self) -> Cow<[f64]> {
        Cow::from(&self.mz_list)
    }

    fn get_intensity_list(&self) -> Cow<[f32]> {
        Cow::from(&self.intensity_list)
    }
}

#[allow(dead_code)]
impl SpectrumData {
    fn to_peaks(&self) -> Vec<Peak> {
        let peaks = self.mz_list
            .iter().copied()
            .zip(self.intensity_list.iter())
            .map(|(mz, &intensity)| {
                Peak { mz, intensity }
            })
            .collect::<Vec<_>>();

        peaks
    }
}

// --- Similar to sage definitions --- //

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct Peak {
    pub mz: f64,
    pub intensity: f32,
}

impl Eq for Peak {}

impl PartialOrd for Peak {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        self.intensity
            .partial_cmp(&other.intensity)
            .or_else(|| self.mz.partial_cmp(&other.mz))
    }
}

impl Ord for Peak {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.partial_cmp(other).unwrap_or(std::cmp::Ordering::Equal)
    }
}