use std::{cmp::Ordering, collections::HashMap, hash::Hash};

use crate::{
    annotation::model::ChargeRange,
    chemistry::{Chemical, Element, MolecularFormula},
    sequence::SequencePosition,
    system::isize::Charge,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct CachedCharge {
    charge: MolecularCharge,
    options: HashMap<Charge, Vec<MolecularCharge>>,
    number: Charge,
}

impl CachedCharge {
    pub fn charge(&self) -> Charge {
        self.number
    }

    /// Get all options resulting in this exact charge
    pub fn options(&mut self, charge: Charge) -> &[MolecularCharge] {
        self.options
            .entry(charge)
            .or_insert_with(|| self.charge.options(charge))
    }

    /// Get all options
    pub fn range(&mut self, range: ChargeRange) -> Vec<MolecularCharge> {
        let mut options = Vec::new();
        for c in range.charges_iter(self.charge()) {
            options.extend_from_slice(self.options(c));
        }
        options
    }
}

impl From<MolecularCharge> for CachedCharge {
    fn from(value: MolecularCharge) -> Self {
        let n = value.charge();
        Self {
            charge: value,
            options: HashMap::new(),
            number: n,
        }
    }
}

impl From<&MolecularCharge> for CachedCharge {
    fn from(value: &MolecularCharge) -> Self {
        Self {
            charge: value.clone(),
            options: HashMap::new(),
            number: value.charge(),
        }
    }
}

/// A selection of ions that together define the charge of a peptide
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug, Serialize, Deserialize, Hash)]
pub struct MolecularCharge {
    /// The ions that together define the charge of the peptide.
    /// The first number is the amount of times this adduct ion occurs, the molecular formula is the full formula for the adduct ion.
    /// The charge for each ion is saved as the number of electrons missing or gained in the molecular formula.
    pub charge_carriers: Vec<(isize, MolecularFormula)>,
}

impl MolecularCharge {
    /// Create a default charge state with only protons
    #[expect(clippy::missing_panics_doc)] // Cannot panic
    pub fn proton(charge: isize) -> Self {
        Self {
            charge_carriers: vec![(
                charge,
                MolecularFormula::new(&[(Element::H, None, 1), (Element::Electron, None, -1)], &[])
                    .unwrap(),
            )],
        }
    }

    /// Create a charge state with the given ions
    pub fn new(charge_carriers: &[(isize, MolecularFormula)]) -> Self {
        Self {
            charge_carriers: charge_carriers.to_vec(),
        }
    }

    /// Get all options resulting in this exact charge
    /// # Panics
    /// If the charge is not at least 1.
    pub fn options(&self, charge: Charge) -> Vec<Self> {
        assert!(charge.value > 0);
        let own_charge = self.charge();
        let remainder = charge.value.rem_euclid(own_charge.value);
        let quotient = charge.value.div_euclid(own_charge.value).max(0);

        let mut too_low_options: Vec<Vec<(isize, MolecularFormula)>> = Vec::new();
        let mut options = Vec::new();
        for carrier in &self.charge_carriers {
            let mut new_too_low_options = Vec::new();
            if too_low_options.is_empty() {
                for n in 0..=carrier.0 {
                    let charge = n * carrier.1.charge();
                    match charge.value.cmp(&remainder) {
                        Ordering::Less => new_too_low_options.push(vec![(n, carrier.1.clone())]),
                        Ordering::Equal => options.push(vec![(n, carrier.1.clone())]),
                        Ordering::Greater => (),
                    }
                }
            } else {
                for n in 0..=carrier.0 {
                    for o in &too_low_options {
                        let mut new = o.clone();
                        new.push((n, carrier.1.clone()));
                        let full_charge = new
                            .iter()
                            .fold(Charge::default(), |acc, (amount, formula)| {
                                acc + *amount * formula.charge()
                            });

                        let charge = n * carrier.1.charge() + full_charge;
                        match charge.value.cmp(&remainder) {
                            Ordering::Less => new_too_low_options.push(new),
                            Ordering::Equal => options.push(new),
                            Ordering::Greater => (),
                        }
                    }
                }
            }
            too_low_options = new_too_low_options;
        }

        options
            .into_iter()
            .map(|charge_carriers| {
                let mut charge_carriers = charge_carriers;
                charge_carriers.extend(
                    std::iter::repeat_n(self.charge_carriers.clone(), quotient as usize).flatten(),
                );
                Self { charge_carriers }.simplified()
            })
            .collect()
    }

    /// Get the total charge of these charge carriers
    pub fn charge(&self) -> Charge {
        self.charge_carriers
            .iter()
            .fold(Charge::default(), |acc, (amount, formula)| {
                acc + *amount * formula.charge()
            })
    }

    // The elements will be sorted on ion and deduplicated
    #[must_use]
    fn simplified(mut self) -> Self {
        self.charge_carriers.retain(|el| el.0 != 0);
        self.charge_carriers.sort_by(|a, b| a.1.cmp(&b.1));
        // Deduplicate
        let mut max = self.charge_carriers.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &self.charge_carriers[index];
            let next = &self.charge_carriers[index + 1];
            if this.1 == next.1 {
                self.charge_carriers[index].0 += next.0;
                self.charge_carriers.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        self.charge_carriers.retain(|el| el.0 != 0);
        self
    }
}

impl Chemical for MolecularCharge {
    fn formula_inner(
        &self,
        _sequence_index: SequencePosition,
        _peptidoform_index: usize,
    ) -> MolecularFormula {
        self.charge_carriers
            .iter()
            .map(|(n, mol)| mol.clone() * *n as i32)
            .sum::<MolecularFormula>()
    }
}

impl std::fmt::Display for MolecularCharge {
    /// Is not guaranteed to fully conform to the ProForma standard. Because the data structure accepts more than the standard.
    /// So adducts with other than +1/-1 charge states, or adducts with complex formula (not a single element) will not adhere to the standard.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            self.charge_carriers
                .iter()
                .map(|c| c.1.charge().value * c.0)
                .sum::<isize>()
        )?;
        if !self.charge_carriers.iter().all(|c| {
            c.1 == MolecularFormula::new(
                &[(Element::H, None, 1), (Element::Electron, None, -1)],
                &[],
            )
            .unwrap()
        }) {
            write!(f, "[")?;
            let mut first = true;
            for (amount, formula) in &self.charge_carriers {
                if first {
                    first = false;
                } else {
                    write!(f, ",")?;
                }
                let charge = formula.charge().value;
                write!(f, "{amount}{formula}{charge:+}")?;
            }
            write!(f, "]")?;
        }
        Ok(())
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use crate::{
        chemistry::{Chemical, MolecularCharge},
        molecular_formula,
    };

    #[test]
    fn simple_charge_options() {
        let mc = MolecularCharge::new(&[(1, molecular_formula!(H 1 Electron -1))]);
        let options = mc.options(crate::system::isize::Charge::new::<crate::system::e>(1));
        assert_eq!(options.len(), 1);
        assert_eq!(options[0].formula(), molecular_formula!(H 1 Electron -1));
    }
}
