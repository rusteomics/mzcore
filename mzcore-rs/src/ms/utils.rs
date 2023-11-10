#![allow(dead_code)]

///
/// Some parts of this file originates from [Sage](https://github.com/lazear/sage)
/// Copyright (c) 2022 Michael Lazear
/// SPDX-License-Identifier: MIT
///

use std::cmp::Ordering;
use std::ops::Mul;

use anyhow::*;
use serde::{Serialize, Deserialize};

use crate::chemistry::constants::PROTON_MASS;

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]//#[serde(rename_all = "lowercase")]
pub enum MassTolWindow {
    ppm(f64, f64),
    mmu(f64, f64),
    Da(f64, f64),
}

impl MassTolWindow {
    /// Compute the (`lower`, `upper`) window (in Da)
    /// for a monoisotopic mass and a given tolerance
    pub fn bounds(&self, center: f64) -> (f64, f64) {
        match self {
            MassTolWindow::Da(lo, hi) => (center + lo, center + hi),
            MassTolWindow::mmu(lo, hi) => (center + lo / 1000.0, center + hi / 1000.0),
            MassTolWindow::ppm(lo, hi) => {
                let delta_lo = center * lo / 1_000_000.0;
                let delta_hi = center * hi / 1_000_000.0;
                (center + delta_lo, center + delta_hi)
            }
        }
    }

    pub fn contains(&self, center: f64, rhs: f64) -> bool {
        let (lo, hi) = self.bounds(center);
        rhs >= lo && rhs <= hi
    }

    pub fn ppm_to_delta_mass(center: f64, ppm: f64) -> f64 {
        ppm * center / 1_000_000.0
    }
}

impl Mul<f64> for MassTolWindow {
    type Output = MassTolWindow;

    fn mul(self, rhs: f64) -> Self::Output {
        match self {
            MassTolWindow::Da(lo, hi) => MassTolWindow::Da(lo * rhs, hi * rhs),
            MassTolWindow::mmu(lo, hi) => MassTolWindow::mmu(lo * rhs, hi * rhs),
            MassTolWindow::ppm(lo, hi) => MassTolWindow::ppm(lo * rhs, hi * rhs),
        }
    }
}

#[allow(non_camel_case_types)]
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]//#[serde(rename_all = "lowercase")]
pub enum MassTolUnit {
    Da,
    mmu,
    ppm
}

impl MassTolUnit {
    fn new(unit: &str) -> Option<MassTolUnit> {
        match unit {
            "Da"  => Some(MassTolUnit::Da),
            "mmu" => Some(MassTolUnit::mmu),
            "ppm" => Some(MassTolUnit::ppm),
            _     => None
        }
    }
}

impl std::fmt::Display for MassTolUnit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

pub fn calc_mz_tol_in_daltons(mz: f64, mz_tol: f64, tol_unit: MassTolUnit) -> f64 {
    match tol_unit {
        MassTolUnit::Da => mz_tol,
        MassTolUnit::mmu => mz_tol / 1000.0,
        MassTolUnit::ppm => mz_tol * mz / 1_000_000.0
    }
}

pub fn calc_mz_tol_in_ppm(mz: f64, mz_tol: f64, tol_unit: MassTolUnit) -> f64 {
    match tol_unit {
        MassTolUnit::Da => mz_tol * 1e6 / mz,
        MassTolUnit::mmu => mz_tol * 1000.0 / mz,
        MassTolUnit::ppm => mz_tol
    }
}

pub fn mz_to_mass(mz: f64, charge: i32) -> f64 {
    let z = charge as f64;
    mz * z.abs() - z * PROTON_MASS
}
pub fn mass_to_mz(mass: f64, charge: i32) -> f64 {
    if charge == 1 {
        mass + PROTON_MASS
    } else {
        let z = charge as f64;
        (mass + z * PROTON_MASS) / z.abs()
    }
}

/// Return the widest `left` and `right` indices into a `slice` (sorted by the
/// function `key`) such that all values between `low` and `high` are
/// contained in `slice[left..right]`
///
/// # Invariants
///
/// * `slice[left] <= low || left == 0`
/// * `slice[right] <= high && (slice[right+1] > high || right == slice.len())`
/// * `0 <= left <= right <= slice.len()`
#[inline]
pub fn binary_search_slice<T, F, S>(slice: &[T], key: F, low: S, high: S) -> (usize, usize)
where
    F: Fn(&T, &S) -> Ordering,
{
    let left_idx = match slice.binary_search_by(|a| key(a, &low)) {
        Result::Ok(idx) | Result::Err(idx) => {
            let mut idx = idx.saturating_sub(1);
            while idx > 0 && key(&slice[idx], &low) != Ordering::Less {
                idx -= 1;
            }
            idx
        }
    };

    let right_idx = match slice[left_idx..].binary_search_by(|a| key(a, &high)) {
        Result::Ok(idx) | Err(idx) => {
            let mut idx = idx + left_idx;
            while idx < slice.len() && key(&slice[idx], &high) != Ordering::Greater {
                idx = idx.saturating_add(1);
            }
            idx.min(slice.len())
        }
    };
    (left_idx, right_idx)
}


#[cfg(test)]
mod tests {

    use super::{MassTolWindow};

    #[test]
    fn tolerances() {
        assert_eq!(
            MassTolWindow::ppm(-10.0, 20.0).bounds(1000.0),
            (999.99, 1000.02)
        );
        assert_eq!(
            MassTolWindow::ppm(-10.0, 10.0).bounds(487.0),
            (486.99513, 487.00487)
        );
        assert_eq!(
            MassTolWindow::ppm(-50.0, 50.0).bounds(1000.0),
            (999.95, 1000.05)
        );
    }
}