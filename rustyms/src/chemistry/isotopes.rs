use std::cmp::Ordering;

use itertools::Itertools;
use ndarray::{Array1, Axis, arr1, concatenate, s};
use probability::distribution::{Binomial, Discrete};

use crate::{chemistry::MolecularFormula, system::Mass};

impl MolecularFormula {
    /// Get the isotopic distribution, using the natural distribution as defined by CIAAW.
    /// All elements are considered. The return is an array with the probability per offset.
    /// The first element of the array is the base peak, every consecutive peak is 1 Dalton heavier.
    /// The probability is normalized to (approximately) 1 total area.
    ///
    /// This approximation slightly overestimates the tail end of the distribution. Especially
    /// for species with multiple higher mass isotopes as it does not take the number of already
    /// chosen atoms for lower weighed isotopes into account.
    #[expect(clippy::missing_panics_doc)]
    pub fn isotopic_distribution(&self, threshold: f64) -> Array1<f64> {
        let mut result = arr1(&[1.0]);
        for (element, isotope, amount) in self.elements() {
            if isotope.is_some() || *amount <= 0 {
                // TODO: think about negative numbers?
                continue;
            }
            let amount = usize::try_from(*amount).unwrap();
            let isotopes = element
                .isotopes()
                .iter()
                .filter(|i| i.2 != 0.0)
                .collect_vec();
            if isotopes.len() < 2 {
                // Only a single species, so no distribution is needed
                continue;
            }
            // Get the probability and base offset (weight) for all non base isotopes
            let base = isotopes[0];
            let isotopes = isotopes
                .into_iter()
                .skip(1)
                .map(|i| (i.0 - base.0, i.2))
                .collect_vec();

            for isotope in isotopes {
                // Generate distribution (take already chosen into account?)
                let binomial = Binomial::new(amount, isotope.1);

                // See how many numbers are below the threshold from the end of the distribution
                let tail = (0..=amount)
                    .rev()
                    .map(|t| binomial.mass(t))
                    .take_while(|a| *a < threshold)
                    .count();

                // Get all numbers start to the tail threshold
                let mut distribution: Array1<f64> = (0..=amount - tail)
                    .map(|t| binomial.mass(t))
                    .flat_map(|a| {
                        // Interweave the probability of this isotope with the mass difference to generate the correct distribution
                        std::iter::once(a)
                            .chain(std::iter::repeat(0.0))
                            .take(isotope.0 as usize)
                    })
                    .collect();

                // Make the lengths equal
                match result.len().cmp(&distribution.len()) {
                    Ordering::Less => {
                        result
                            .append(
                                Axis(0),
                                Array1::zeros(distribution.len() - result.len()).view(),
                            )
                            .unwrap();
                    }
                    Ordering::Greater => {
                        distribution
                            .append(
                                Axis(0),
                                Array1::zeros(result.len() - distribution.len()).view(),
                            )
                            .unwrap();
                    }
                    Ordering::Equal => (),
                }

                // Combine distribution with previous distribution
                let mut new = Array1::zeros(result.len());
                for (i, a) in distribution.into_iter().enumerate() {
                    new += &(concatenate(
                        Axis(0),
                        &[
                            Array1::zeros(i).view(),
                            result.slice(s![0..result.len() - i]),
                        ],
                    )
                    .unwrap()
                        * a);
                }

                result = new;
            }
        }
        result
    }

    // TODO: Calculated mass is incorrect when mixing elements with many isotopes (>2).
    // Enable tests below to check for the error.
    /// Get the isotopic distribution, using the natural distribution as defined by CIAAW.
    /// All elements are considered. The return is an array with the probability per offset.
    /// The first element of the array is the base peak, every consecutive peak is ~1 Da heavier.
    /// The probability is normalized to (approximately) 1 total area.
    ///
    /// This approximation slightly overestimates the tail end of the distribution. Especially
    /// for species with multiple higher mass isotopes as it does not take the number of already
    /// chosen atoms for lower weighed isotopes into account.
    ///
    /// The mass returned is the average mass of all combinations of isotopes that generate that
    /// offset. This follows the `+iA` definition of `mzPAF`.
    #[expect(clippy::missing_panics_doc)]
    #[expect(dead_code)]
    fn isotopic_distribution_with_mass(&self, threshold: f64) -> Array1<(Mass, f64)> {
        let mut full_distribution = arr1(&[(Mass::default(), 1.0_f64)]);
        for (element, isotope, amount) in self.elements() {
            if isotope.is_some() || *amount <= 0 {
                // TODO: think about negative numbers?
                continue;
            }
            let amount = usize::try_from(*amount).unwrap();
            let isotopes = element
                .isotopes()
                .iter()
                .filter(|i| i.2 != 0.0)
                .collect_vec();
            if isotopes.len() < 2 {
                // Only a single species, so no distribution is needed
                continue;
            }
            // Get (offset, mass, probability) for the isotopes
            let base = isotopes[0];
            let isotopes = isotopes
                .into_iter()
                .skip(1)
                .map(|i| (i.0 - base.0, i.1, i.2))
                .collect_vec();

            let mut all_isotopes_distribution = arr1(&[(Mass::default(), 0, 1.0_f64)]);
            for (isotope_offset, isotope_mass, isotope_probability) in isotopes {
                // Generate distribution (take already chosen into account?)
                let binomial = Binomial::new(amount, isotope_probability);

                // See how many numbers are below the threshold from the end of the distribution
                let tail = (0..=amount)
                    .rev()
                    .map(|t| binomial.mass(t))
                    .take_while(|a| *a < threshold)
                    .count();

                // Get all numbers start to the tail threshold
                let distribution: Array1<f64> = (0..=amount - tail)
                    .map(|t| binomial.mass(t))
                    .flat_map(|a| {
                        // Interweave the probability of this isotope with the mass difference to generate the correct distribution
                        std::iter::once(a)
                            .chain(std::iter::repeat(0.0))
                            .take(isotope_offset as usize)
                    })
                    .collect();

                // Make the result fit the data
                all_isotopes_distribution
                    .append(
                        Axis(0),
                        Array1::from_elem(
                            (all_isotopes_distribution.len() + distribution.len() - 1)
                                .saturating_sub(all_isotopes_distribution.len()),
                            (Mass::default(), 0, 0.0),
                        )
                        .view(),
                    )
                    .unwrap();

                // Combine distribution with previous distribution
                let mut temporary_stack =
                    Array1::from_elem(all_isotopes_distribution.len(), (Mass::default(), 0, 0.0));

                dbg!(&distribution);
                for (shift, shift_probability) in distribution.into_iter().enumerate() {
                    // The number of this element to add
                    let num = (shift / isotope_offset as usize).min(amount);

                    // The mass of this TODO: result in too high masses if multiple isotopes (> 2) exist
                    let isotope_total_mass = isotope_mass * num as f64;
                    println!(
                        "E {element} {:.4} i {isotope_offset} p {isotope_probability} s {shift} n {num} p {shift_probability} im {:.4}",
                        isotope_mass.value, isotope_total_mass.value
                    );
                    temporary_stack = temporary_stack
                        .into_iter()
                        .zip(
                            std::iter::repeat_n((Mass::default(), 0, 0.0), shift)
                                .chain(all_isotopes_distribution.iter().copied()),
                        )
                        .enumerate()
                        .map(|(index, (new, old))| {
                            (
                                new.0 + (old.0 + isotope_total_mass) * (old.2 * shift_probability),
                                new.1 + if index == shift { num } else { 0 },
                                old.2.mul_add(shift_probability, new.2),
                            )
                        })
                        .collect();
                }
                all_isotopes_distribution = temporary_stack;
            }
            all_isotopes_distribution = all_isotopes_distribution
                .into_iter()
                .filter(|v| v.2 >= threshold)
                .collect();
            dbg!(&all_isotopes_distribution);
            full_distribution
                .append(
                    Axis(0),
                    Array1::from_elem(
                        all_isotopes_distribution
                            .len()
                            .saturating_sub(full_distribution.len()),
                        (Mass::default(), 0.0),
                    )
                    .view(),
                )
                .unwrap();

            let mut temporary_stack =
                Array1::from_elem(full_distribution.len(), (Mass::default(), 0.0));
            for (shift, (mass, total_isotopes, shift_probability)) in
                all_isotopes_distribution.into_iter().enumerate()
            {
                let mass = mass / shift_probability
                    + base.1 * amount.saturating_sub(total_isotopes) as f64;
                temporary_stack = temporary_stack
                    .into_iter()
                    .zip(
                        std::iter::repeat_n((Mass::default(), 0.0), shift)
                            .chain(full_distribution.iter().copied()),
                    )
                    .enumerate()
                    .map(|(index, (new, old))| {
                        (
                            new.0
                                + if index == shift {
                                    mass + old.0
                                } else {
                                    Mass::default()
                                },
                            old.1.mul_add(shift_probability, new.1),
                        )
                    })
                    .collect();
            }
            full_distribution = temporary_stack;
            full_distribution = full_distribution
                .into_iter()
                .filter(|v| v.1 >= threshold)
                .collect();
            dbg!(&full_distribution);
        }
        full_distribution
    }
}

#[cfg(never)] // Set to test to reenable the tests
#[allow(clippy::missing_panics_doc)]
mod test {
    #[test]
    fn with_mass() {
        // let formula = molecular_formula!(C 18 H 33 S 2 O 6 N 8);
        // let formula = molecular_formula!(O 100);
        // let formula = molecular_formula!(Se 10);
        let formula = molecular_formula!(H 202 O 10);

        let distribution = formula.isotopic_distribution_with_mass(0.001);
        dbg!(distribution);

        todo!();
    }

    #[test]
    fn h_mass_and_probability() {
        let formula = molecular_formula!(H 100);
        let expected = vec![
            (molecular_formula!([1 H 100]).monoisotopic_mass(), 0.985_604),
            (
                molecular_formula!([1 H 99] [2 H 1]).monoisotopic_mass(),
                0.014_293,
            ),
        ];
        let distribution = formula.isotopic_distribution_with_mass(0.001);

        assert_eq!(distribution.len(), expected.len());
        for (gotten, expected) in distribution.iter().zip(expected) {
            assert!((gotten.0.value - expected.0.value).abs() < 1E-4);
            assert!((gotten.1 - expected.1).abs() < 1E-4);
        }
    }

    #[test]
    fn o_mass_and_probability() {
        let formula = molecular_formula!(O 25);
        let expected = vec![
            (molecular_formula!([16 O 25]).monoisotopic_mass(), 0.941_043),
            (
                molecular_formula!([16 O 24] [17 O 1]).monoisotopic_mass(),
                0.009_026,
            ),
            (
                molecular_formula!([16 O 24] [18 O 1]).monoisotopic_mass(),
                0.048_209,
            ),
            (
                molecular_formula!([16 O 23] [17 O 1] [18 O 1]).monoisotopic_mass(),
                0.001_185,
            ),
        ];
        let distribution = formula.isotopic_distribution_with_mass(0.001);
        dbg!(&distribution);
        assert_eq!(distribution.len(), expected.len());
        for (offset, (gotten, expected)) in distribution.iter().zip(expected).enumerate() {
            assert!(
                (gotten.0.value - expected.0.value).abs() < 1E-4,
                "Wrong mass at offset: {offset} Expected: {:.6} Gotten: {:.6}",
                expected.0.value,
                gotten.0.value,
            );
            assert!(
                (gotten.1 - expected.1).abs() < 1E-4,
                "Wrong probability at offset: {offset}"
            );
        }
    }

    #[test]
    fn ho_mass_and_probability() {
        let formula = molecular_formula!(H 100 O 25);
        let expected = vec![
            (
                molecular_formula!([1 H 100] [16 O 25]).monoisotopic_mass(),
                0.927_495,
            ),
            (
                molecular_formula!([1 H 100] [16 O 24] [17 O 1]).monoisotopic_mass(), // TODO: is mixture
                0.022_346,
            ),
            (
                molecular_formula!([1 H 100] [16 O 24] [18 O 1]).monoisotopic_mass(), // TODO: is mixture
                0.047_515,
            ),
            (
                molecular_formula!([1 H 100] [16 O 23] [17 O 1] [18 O 1]).monoisotopic_mass(), // TODO: is mixture
                0.001_145,
            ),
        ];
        let distribution = formula.isotopic_distribution_with_mass(0.001);
        dbg!(&distribution);
        assert_eq!(distribution.len(), expected.len());
        for (offset, (gotten, expected)) in distribution.iter().zip(expected).enumerate() {
            assert!(
                (gotten.0.value - expected.0.value).abs() < 1E-4,
                "Wrong mass at offset: {offset} Expected: {:.6} Gotten: {:.6}",
                expected.0.value,
                gotten.0.value,
            );
            assert!(
                (gotten.1 - expected.1).abs() < 1E-4,
                "Wrong probability at offset: {offset}"
            );
        }
    }
}
