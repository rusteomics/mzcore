use crate::Element;

/// A molecular formula, a selection of elements of specified isotopes together forming a structure
#[derive(Debug, Clone, PartialEq, Default)]
pub struct MolecularFormula {
    /// Save all constituent parts as the element in question, the isotope (or 0 for natural distribution), and the number of this part
    elements: Vec<(crate::Element, u16, i16)>,
    /// Any addition mass, defined to be monoisotopic
    additional_mass: f64,
}

/// Any item that has a clearly defined molecular formula
pub trait Chemical {
    /// Get the molecular formula
    fn formula(&self) -> MolecularFormula;
}

impl MolecularFormula {
    /// Create a new molecular formula, the elements will be sorted on element/isotope and deduplicated
    pub fn new(elements: &[(crate::Element, u16, i16)]) -> Self {
        let result = Self {
            elements: elements.to_vec(),
            additional_mass: 0.0,
        };
        result.simplify()
    }

    // The elements will be sorted on element/isotope and deduplicated
    #[must_use]
    fn simplify(mut self) -> Self {
        self.elements.retain(|e| e.2 != 0);
        self.elements.sort_by(|a, b| {
            if a.0 == b.0 {
                // If the elements are the same sort on the isotope number
                a.1.cmp(&b.1)
            } else {
                a.0.cmp(&b.0)
            }
        });
        // Deduplicate
        let mut max = self.elements.len() - 1;
        let mut index = 0;
        while index < max {
            let this = self.elements[index];
            let next = self.elements[index + 1];
            if this.0 == next.0 && this.1 == next.1 {
                self.elements[index].2 += next.2;
                self.elements.remove(index + 1);
                max -= 1;
            } else {
                index += 1;
            }
        }
        self
    }

    /// Get an empty molecular formula with only a mass of unspecified origin
    pub const fn with_additional_mass(additional_mass: f64) -> Self {
        Self {
            elements: Vec::new(),
            additional_mass,
        }
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    pub fn add(&mut self, element: (crate::Element, u16, i16)) {
        let mut index = 0;
        let mut done = false;
        let (e, i, n) = element;
        while !done {
            let base = self.elements.get(index).copied();
            if let Some((re, ri, _)) = base {
                if e > re || (e == re && i > ri) {
                    index += 1;
                } else if e == re && i == ri {
                    self.elements[index].2 += n;
                    done = true;
                } else {
                    self.elements.insert(index, (e, i, n));
                    done = true;
                }
            } else {
                self.elements.push((e, i, n));
                done = true;
            }
        }
    }

    /// Get the elements making this formula
    pub fn elements(&self) -> &[(Element, u16, i16)] {
        &self.elements
    }

    /// Create a new molecular formula with the given global isotope modifications
    #[must_use]
    pub fn with_global_isotope_modifications(&self, substitutions: &[(Element, u16)]) -> Self {
        let mut new_elements = self.elements.clone();
        for item in &mut new_elements {
            for (substitute_element, substitute_species) in substitutions {
                if item.0 == *substitute_element {
                    item.1 = *substitute_species;
                }
            }
        }
        let result = Self {
            elements: new_elements,
            additional_mass: self.additional_mass,
        };
        result.simplify()
    }
}

impl std::ops::Add<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass += rhs.additional_mass;

        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if e > re || (e == re && i > ri) {
                    index_result += 1;
                } else if e == re && i == ri {
                    result.elements[index_result].2 += n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (e, i, n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((e, i, n));
                index_rhs += 1;
            }
        }
        result
    }
}

impl std::ops::Add<Self> for MolecularFormula {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}

impl std::ops::Sub<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn sub(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass -= rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if e > re || (e == re && i > ri) {
                    index_result += 1;
                } else if e == re && i == ri {
                    result.elements[index_result].2 -= n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (e, i, -n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((e, i, -n));
                index_rhs += 1;
            }
        }
        result
    }
}

impl std::ops::Sub<Self> for MolecularFormula {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}

impl std::ops::Mul<i16> for MolecularFormula {
    type Output = Self;
    fn mul(mut self, rhs: i16) -> Self::Output {
        self.additional_mass *= f64::from(rhs);
        self.elements.iter_mut().for_each(|part| part.2 *= rhs);
        self
    }
}

impl std::ops::AddAssign<&Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: &Self) {
        let mut index_self = 0;
        let mut index_rhs = 0;
        self.additional_mass += rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (e, i, n) = rhs.elements[index_rhs];
            if index_self < self.elements.len() {
                let (re, ri, _) = self.elements[index_self];
                if e > re || (e == re && i > ri) {
                    index_self += 1;
                } else if e == re && i == ri {
                    self.elements[index_self].2 += n;
                    index_rhs += 1;
                } else {
                    self.elements.insert(index_self, (e, i, n));
                    index_rhs += 1;
                }
            } else {
                self.elements.push((e, i, n));
                index_rhs += 1;
            }
        }
    }
}

impl std::ops::AddAssign<Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl std::iter::Sum<Self> for MolecularFormula {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res += v);
        res
    }
}

macro_rules! molecular_formula {
    ($($tail:tt)*) => {
        formula_internal!([$($tail)*] -> [])
    };
}

macro_rules! formula_internal {
    ([$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        formula_internal!([$($tail)*] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:expr] -> [$($output:tt)*]) =>{
        formula_internal!([] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([] -> [$($output:tt)*]) =>{
        MolecularFormula::new(&[$($output)*])
    };
}
