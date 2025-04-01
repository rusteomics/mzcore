use std::num::NonZeroU16;
use std::sync::OnceLock;

use crate::system::{da, fraction, Ratio};

include!("shared/element.rs");

impl Element {
    /// Validate this isotope to have a defined mass
    pub fn is_valid(self, isotope: Option<NonZeroU16>) -> bool {
        if self == Self::Electron {
            isotope.is_none()
        } else {
            isotope.map_or_else(
                || elemental_data()[self as usize - 1].0.is_some(),
                |isotope| {
                    elemental_data()[self as usize - 1]
                        .2
                        .iter()
                        .any(|(ii, _, _)| *ii == isotope.get())
                },
            )
        }
    }

    /// Get all available isotopes (N, mass, abundance)
    pub fn isotopes(self) -> &'static [(u16, Mass, f64)] {
        &elemental_data()[self as usize - 1].2
    }

    /// The mass of the specified isotope of this element (if that isotope exists)
    pub fn mass(self, isotope: Option<NonZeroU16>) -> Option<Mass> {
        if self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        isotope.map_or_else(
            || elemental_data()[self as usize - 1].0,
            |isotope| {
                // Specific isotope do not change anything
                elemental_data()[self as usize - 1]
                    .2
                    .iter()
                    .find(|(ii, _, _)| *ii == isotope.get())
                    .map(|(_, m, _)| *m)
            },
        )
    }

    /// The average weight of the specified isotope of this element (if that isotope exists)
    pub fn average_weight(self, isotope: Option<NonZeroU16>) -> Option<Mass> {
        if self == Self::Electron {
            return Some(da(5.485_799_090_65e-4));
        }
        isotope.map_or_else(
            || elemental_data()[self as usize - 1].1,
            |isotope| {
                // Specific isotope do not change anything
                elemental_data()[self as usize - 1]
                    .2
                    .iter()
                    .find(|(ii, _, _)| *ii == isotope.get())
                    .map(|(_, m, _)| *m)
            },
        )
    }

    /// Gives the most abundant mass based on the number of this isotope
    pub fn most_abundant_mass(self, isotope: Option<NonZeroU16>, n: i32) -> Option<Mass> {
        if self == Self::Electron {
            return Some(da(5.485_799_090_65e-4) * Ratio::new::<fraction>(f64::from(n)));
        }
        Some(
            if let Some(isotope) = isotope {
                // Specific isotope do not change anything
                elemental_data()[self as usize - 1]
                    .2
                    .iter()
                    .find(|(ii, _, _)| *ii == isotope.get())
                    .map(|(_, m, _)| *m)?
            } else {
                // (mass, chance)
                let mut max = None;
                for iso in &elemental_data()[self as usize - 1].2 {
                    let chance = iso.2 * f64::from(n);
                    if max.is_none_or(|m: (Mass, f64)| chance > m.1) {
                        max = Some((iso.1, chance));
                    }
                }
                max?.0
            } * Ratio::new::<fraction>(f64::from(n)),
        )
    }
}

/// Get the elemental data
/// # Panics
/// It panics if the elemental data that is passed at compile time is not formatted correctly.
pub fn elemental_data() -> &'static ElementalData {
    ELEMENTAL_DATA_CELL
        .get_or_init(|| bincode::deserialize(include_bytes!("databases/elements.dat")).unwrap())
}
static ELEMENTAL_DATA_CELL: OnceLock<ElementalData> = OnceLock::new();

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod test {
    #[test]
    fn hill_notation() {
        assert_eq!(
            molecular_formula!(C 6 O 5 H 10).hill_notation(),
            "C6H10O5".to_string()
        );
    }
}
