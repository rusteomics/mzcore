use std::{cmp::Ordering, collections::HashMap, hash::Hash};

use crate::{
    chemistry::{ChargeRange, Chemical, Element, MolecularFormula},
    sequence::SequencePosition,
    system::isize::Charge,
};
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

/// A [`MolecularCharge`] that caches the options for each charge, to not calculate this every time
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct CachedCharge {
    charge: MolecularCharge,
    options: HashMap<Charge, Vec<MolecularCharge>>,
    number: Charge,
}

impl CachedCharge {
    /// Get the main charge
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
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MolecularCharge {
    /// The ions that together define the charge of the peptide.
    /// The first number is the amount of times this adduct ion occurs, the molecular formula is the full formula for the adduct ion.
    /// The charge for each ion is saved as the number of electrons missing or gained in the molecular formula.
    pub charge_carriers: ThinVec<(isize, MolecularFormula)>,
}

impl MolecularCharge {
    /// Create a default charge state with only protons
    pub fn proton(charge: Charge) -> Self {
        if charge.value == 0 {
            Self {
                charge_carriers: ThinVec::new(),
            }
        } else {
            Self {
                charge_carriers: vec![(charge.value, molecular_formula!(H 1 :z+1))].into(),
            }
        }
    }

    /// Check if this molecular charge only consists of protons
    pub fn is_proton(&self) -> bool {
        self.charge_carriers
            .iter()
            .all(|(_, m)| *m == molecular_formula!(H 1 :z+1))
    }

    /// Create a charge state with the given ions
    pub fn new(charge_carriers: &[(isize, MolecularFormula)]) -> Self {
        Self {
            charge_carriers: charge_carriers.into(),
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
                Self {
                    charge_carriers: charge_carriers.into(),
                }
                .simplified()
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
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_proton() {
            write!(f, "{}", self.charge().value)
        } else if self.charge_carriers.iter().any(|(o, _)| *o != 0) {
            write!(f, "[")?;
            let mut first = true;
            for (amount, formula) in &self.charge_carriers {
                if *amount == 0 {
                    continue;
                }
                if first {
                    first = false;
                } else {
                    write!(f, ",")?;
                }
                write!(f, "{formula}")?;
                if *amount != 1 {
                    write!(f, "^{amount}")?;
                }
            }
            write!(f, "]")
        } else {
            Ok(())
        }
    }
}

impl crate::space::Space for MolecularCharge {
    fn space(&self) -> crate::space::UsedSpace {
        self.charge_carriers.space()
    }
}

impl From<Vec<(isize, MolecularFormula)>> for MolecularCharge {
    fn from(value: Vec<(isize, MolecularFormula)>) -> Self {
        Self {
            charge_carriers: value.into(),
        }
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use crate::chemistry::{Chemical, MolecularCharge};

    #[test]
    fn simple_charge_options() {
        let mc = MolecularCharge::new(&[(1, molecular_formula!(H 1 :z+1))]);
        let options = mc.options(crate::system::isize::Charge::new::<crate::system::e>(1));
        assert_eq!(options.len(), 1);
        assert_eq!(options[0].formula(), molecular_formula!(H 1 :z+1));
    }
}
