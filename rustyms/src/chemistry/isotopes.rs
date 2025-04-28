use crate::chemistry::MolecularFormula;
use itertools::Itertools;
use ndarray::{arr1, concatenate, s, Array1, Axis};
use probability::distribution::{Binomial, Discrete};
use std::cmp::Ordering;

impl MolecularFormula {
    /// Get the isotopic distribution, using the natural distribution as defined by CIAAW.
    /// All elements are considered. The return is an array with the probability per offset.
    /// The first element of the array is the base peak, every consecutive peak is 1 Dalton heavier.
    /// The probability is normalized to (approximately) 1 total area.
    ///
    /// This approximation slightly overestimates the tail end of the distribution. Especially
    /// for species with multiple higher mass isotopes as it does not take the number of already
    /// chosen atom for lower weighed isotopes into account.
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
}
