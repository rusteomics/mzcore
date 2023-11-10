#![allow(dead_code)]
///
/// Some parts of this file originates from [rustyms](https://github.com/snijderlab/rustyms)
/// Copyright (c) 2023 Douwe Schulte and contributors
/// MIT License

use anyhow::*;
use std::ops::{Add, AddAssign, Mul, Sub};
use serde::{Deserialize, Serialize};

use crate::chemistry::element::*;
use crate::chemistry::api::Chemical;
use crate::chemistry::unimod::parse_unimod_composition;

#[derive(Clone, Copy, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ElementCount {
    pub element: Element,
    pub isotope_index: u8, // 0 for natural distribution
    pub count: f32, // floating numbers allows for averagine calculations
}

impl ElementCount {
    pub fn new(element: Element, isotope_index: u8, count: f32) -> Self {
        ElementCount {
            element,
            isotope_index: isotope_index,
            count: count,
        }
    }
    pub fn new_monoisotope(element: Element, count: f32) -> Self {
        ElementCount {
            element,
            isotope_index: 0, // zero for mono isotope
            count: count,
        }
    }
}

/// A molecular formula, a selection of elements of specified isotopes together forming a structure
#[derive(Clone, Debug, Default, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct ElementalComposition {
    /// Save all constituent parts as the element in question, the isotope (or 0 for natural distribution), and the number of this part
    pub element_counts: Vec<ElementCount>,
    /// Any additional mass, defined to be monoisotopic
    pub additional_mass: f64,
}

impl ElementalComposition {

    /// Create a new molecular formula, the elements will be sorted on element/isotope and deduplicated
    pub fn new(element_counts: &[ElementCount]) -> Self {
        let result = Self {
            element_counts: element_counts.to_vec(),
            additional_mass: 0.0,
        };
        result.simplify()
    }

    pub fn from_monoisotope_tuples(element_counts: &[(Element, i16)]) -> Self {
        let result = Self {
            element_counts: element_counts.iter().map(|&(element, count)| {
                ElementCount::new_monoisotope(element, count as f32)
            }).collect(),
            additional_mass: 0.0,
        };
        result.simplify()
    }

    pub fn from_tuples(element_counts: &[(Element, u8, i16)]) -> Self {
        let result = Self {
            element_counts: element_counts.iter().map(|&(element, isotope_index, count)| {
                ElementCount { element, count: count as f32, isotope_index }
            }).collect(),
            additional_mass: 0.0,
        };
        result.simplify()
    }

    pub fn parse_unimod_composition(composition: &str) -> Result<Self> {
        let (mut elem_comp, glycan_comp) = parse_unimod_composition(composition)?;
        for (glycan, glycan_count) in glycan_comp {
            elem_comp += glycan.composition() * glycan_count;
        }

        Ok(elem_comp)
    }

    // The elements will be sorted on element/isotope and deduplicated
    #[must_use]
    fn simplify(mut self) -> Self {
        self.element_counts.retain(|el| el.count != 0.0);

        self.element_counts.sort_by(|a, b| {
            if a.element == b.element {
                // If the elements are the same sort on the isotope number
                a.isotope_index.cmp(&b.isotope_index)
            } else {
                a.element.cmp(&b.element)
            }
        });

        // Deduplicate
        let mut max = self.element_counts.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = self.element_counts[index];
            let next = self.element_counts[index + 1];
            if this.element == next.element && this.isotope_index == next.isotope_index {
                self.element_counts[index].count += next.count;
                self.element_counts.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }

        self.element_counts.retain(|el| el.count != 0.0);

        self
    }

    /// Get an empty molecular formula with only a mass of unspecified origin
    pub const fn with_additional_mass(additional_mass: f64) -> Self {
        Self {
            element_counts: Vec::new(),
            additional_mass,
        }
    }

    /// Add the given element to this formula (while keeping it ordered and simplified)
    pub fn add(&mut self, element_count: ElementCount) {
        let mut index = 0;
        let mut done = false;
        let ElementCount {element: el, isotope_index: i, count: n} = element_count;
        while !done {
            let base = self.element_counts.get(index).copied();
            if let Some(ElementCount {element: re, isotope_index: ri, count: _}) = base {
                if el > re || (el == re && i > ri) {
                    index += 1;
                } else if el == re && i == ri {
                    self.element_counts[index].count += n;
                    done = true;
                } else {
                    self.element_counts.insert(index, ElementCount::new(el, i, n));
                    done = true;
                }
            } else {
                self.element_counts.push(ElementCount::new(el, i, n));
                done = true;
            }
        }
    }

    /// Get the elements making this formula
    /*pub fn element_counts(&self) -> &[ElementCount] {
        &self.element_counts
    }*/

    /// Create a new molecular formula with the given global isotope modifications
    #[must_use]
    pub fn with_global_isotope_modifications(&self, substitutions: &[(Element, u8)]) -> Self {
        let mut new_elements = self.element_counts.clone();
        for item in &mut new_elements {
            for (substitute_element, substitute_species) in substitutions {
                if item.element == *substitute_element {
                    item.isotope_index = *substitute_species;
                }
            }
        }
        let result = Self {
            element_counts: new_elements,
            additional_mass: self.additional_mass,
        };
        result.simplify()
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons).
    /// The inverse of that number is given as the charge.
    pub fn charge(&self) -> i16 {
        -self
            .element_counts
            .iter()
            .find(|el| el.element == Element::Electron)
            .map_or(0, |el| el.count as i16)
    }
}

impl Add<&ElementalComposition> for &ElementalComposition {
    type Output = ElementalComposition;
    fn add(self, rhs: &ElementalComposition) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass += rhs.additional_mass;

        while index_rhs < rhs.element_counts.len() {
            let ElementCount {element: el, isotope_index: i, count: n} = rhs.element_counts[index_rhs];
            if index_result < result.element_counts.len() {
                let ElementCount {element: re, isotope_index: ri, count: _} = result.element_counts[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.element_counts[index_result].count += n;
                    index_rhs += 1;
                } else {
                    result.element_counts.insert(index_result, ElementCount::new(el, i, n));
                    index_rhs += 1;
                }
            } else {
                result.element_counts.push(ElementCount::new(el, i, n));
                index_rhs += 1;
            }
        }
        result.element_counts.retain(|el| el.count != 0.0);
        result
    }
}

impl Sub<&ElementalComposition> for &ElementalComposition {
    type Output = ElementalComposition;
    fn sub(self, rhs: &ElementalComposition) -> Self::Output {
        let mut result = (*self).clone();
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass -= rhs.additional_mass;
        while index_rhs < rhs.element_counts.len() {
            let ElementCount {element: el, isotope_index: i, count: n} = rhs.element_counts[index_rhs];
            if index_result < result.element_counts.len() {
                let ElementCount {element: re, isotope_index: ri, count: _} = result.element_counts[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.element_counts[index_result].count -= n;
                    index_rhs += 1;
                } else {
                    result.element_counts.insert(index_result, ElementCount::new(el, i, -n));
                    index_rhs += 1;
                }
            } else {
                result.element_counts.push(ElementCount::new(el, i, -n));
                index_rhs += 1;
            }
        }
        result.element_counts.retain(|el| el.count != 0.0);
        result
    }
}

impl Mul<&i16> for &ElementalComposition {
    type Output = ElementalComposition;
    fn mul(self, rhs: &i16) -> Self::Output {
        ElementalComposition {
            additional_mass: self.additional_mass * f64::from(*rhs),
            element_counts: self
                .element_counts
                .iter()
                .copied()
                .map(|part| ElementCount::new(part.element, part.isotope_index, part.count * f32::from(*rhs)))
                .collect(),
        }
    }
}

impl_binop_ref_cases!(impl Add, add for ElementalComposition, ElementalComposition, ElementalComposition);
impl_binop_ref_cases!(impl Sub, sub for ElementalComposition, ElementalComposition, ElementalComposition);
impl_binop_ref_cases!(impl Mul, mul for ElementalComposition, i16, ElementalComposition);

impl AddAssign<&Self> for ElementalComposition {
    fn add_assign(&mut self, rhs: &Self) {
        let mut index_self = 0;
        let mut index_rhs = 0;
        self.additional_mass += rhs.additional_mass;
        while index_rhs < rhs.element_counts.len() {
            let ElementCount {element: el, isotope_index: i, count: n} = rhs.element_counts[index_rhs];
            if index_self < self.element_counts.len() {
                let ElementCount {element: re, isotope_index: ri, count: _} = self.element_counts[index_self];
                if el > re || (el == re && i > ri) {
                    index_self += 1;
                } else if el == re && i == ri {
                    self.element_counts[index_self].count += n;
                    index_rhs += 1;
                } else {
                    self.element_counts.insert(index_self, ElementCount::new(el, i, n));
                    index_rhs += 1;
                }
            } else {
                self.element_counts.push(ElementCount::new(el, i, n));
                index_rhs += 1;
            }
        }
    }
}

impl AddAssign<Self> for ElementalComposition {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl std::iter::Sum<Self> for ElementalComposition {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res += v);
        res
    }
}


/// Easily define molecular formulas using the following syntax: `<element> <num>` or `(<isotope>)<element> <num>`
/// ```
/// # use mzcore::chemistry::composition::molecular_formula;
/// molecular_formula!(C 12 (13)C 1 H 24);
/// ```
macro_rules! molecular_formula {
    ($($tail:tt)*) => {
        __formula_internal__!([$($tail)*] -> [])
    };
}

/// Internal code for the [`molecular_formula`] macro.
macro_rules! __formula_internal__ {
    ([$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        __formula_internal__!([$($tail)*] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        __formula_internal__!([$($tail)*] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([$e:ident $n:expr] -> [$($output:tt)*]) =>{
        __formula_internal__!([] -> [$($output)*(Element::$e, 0, $n),])
    };
    ([($i:literal)$e:ident $n:expr] -> [$($output:tt)*]) =>{
        __formula_internal__!([] -> [$($output)*(Element::$e, $i, $n),])
    };
    ([] -> [$($output:tt)*]) =>{
        ElementalComposition::from_tuples(&[$($output)*])
    };
}

pub(crate) use __formula_internal__;
pub(crate) use molecular_formula;

// TODO: put in a shared utility when used from different places
/// Implement a binary operator for all ref cases after the implementation for the ref-ref case (assumes deref operator works)
#[allow(unused)]
macro_rules! impl_binop_ref_cases {
    (impl $imp:ident, $method:ident for $t:ty, $u:ty, $o:ty) => {
        impl<'a> $imp<$u> for &'a $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(self, &other)
            }
        }

        impl<'a> $imp<&'a $u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: &'a $u) -> $o {
                $imp::$method(&self, other)
            }
        }

        impl $imp<$u> for $t {
            type Output = $o;

            #[inline]
            fn $method(self, other: $u) -> $o {
                $imp::$method(&self, &other)
            }
        }
    };
}

pub(in crate::chemistry::composition) use impl_binop_ref_cases;

// Additional definitions to define an atom-based composition (mainly used for mass calculation)
use crate::chemistry::atom::AtomIsotopicVariant;

#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct AtomicComposition {
    pub atoms: Vec<(AtomIsotopicVariant, f32)>,
    /// Any additional mass, defined to be monoisotopic
    pub additional_mass: f64,
}

impl AtomicComposition {
    fn calc_mass(&self) -> f64 {
        self.atoms.iter().fold(0.0,|sum, (atom, count)| {
            sum + atom.mass() * (*count as f64)
        })
    }
}