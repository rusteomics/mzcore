use std::num::NonZeroU16;

use crate::{
    chemistry::{Element, MolecularFormula},
    quantities::{Multi, Tolerance},
    system::Mass,
};

/// Find the molecular formulas that fit the mass given the tolerance using only the provided elements.
pub fn find_formulas(
    mass: Mass,
    tolerance: Tolerance<Mass>,
    elements: &[(Element, Option<NonZeroU16>)],
) -> Multi<MolecularFormula> {
    let bounds = tolerance.bounds(mass);
    let mut options: Vec<(Mass, MolecularFormula)> = Vec::new();

    for (element, isotope) in elements.iter().copied() {
        if !element.is_valid(isotope) {
            continue;
        }
        let mass = element.mass(isotope).unwrap();
        let mut new_options = Vec::with_capacity(options.len());
        if mass <= bounds.1 {
            new_options.extend((1..=(bounds.1 / mass).value.floor() as i32).map(|n| {
                (
                    mass * f64::from(n),
                    MolecularFormula::new(&[(element, isotope, n)], &[]).unwrap(),
                )
            }));
        }
        for option in &options {
            let rest = bounds.1 - option.0;
            if mass <= rest {
                new_options.extend((1..=(rest / mass).value.floor() as i32).map(|n| {
                    let mut new_formula = option.1.clone();
                    let _ = new_formula.add((element, isotope, n));
                    (option.0 + mass * f64::from(n), new_formula)
                }));
            }
        }
        options.extend_from_slice(&new_options);
    }
    options.sort_by(|a, b| a.0.value.total_cmp(&b.0.value));
    options
        .into_iter()
        .skip_while(|(mass, _)| *mass < bounds.0)
        .take_while(|(mass, _)| *mass < bounds.1)
        .map(|(_, formula)| formula)
        .collect()
}
