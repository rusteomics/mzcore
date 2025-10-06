use std::{
    fmt::Write,
    hash::Hash,
    num::NonZeroU16,
    ops::{Add, AddAssign, Mul, Neg, Sub, SubAssign},
};

use ordered_float::OrderedFloat;
use serde::{Deserialize, Serialize};
use thin_vec::ThinVec;

use crate::{
    chemistry::Element,
    glycan::{GlycanPosition, MonoSaccharide},
    quantities::Multi,
    sequence::{AminoAcid, CrossLinkName, SequencePosition},
};

/// A molecular formula, a selection of elements of specified isotopes together forming a structure
#[allow(clippy::unsafe_derive_deserialize)]
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct MolecularFormula {
    /// Save all constituent parts as the element in question, the isotope (or None for natural distribution), and the number of this part
    /// The elements will be sorted on element/isotope and deduplicated, guaranteed to only contain valid isotopes.
    pub(in super::super) elements: ThinVec<(Element, Option<NonZeroU16>, i32)>,
    /// Any addition mass, defined to be monoisotopic
    pub(in super::super) additional_mass: OrderedFloat<f64>,
    /// The labels of sources of ambiguity/multiplicity
    #[serde(default)]
    pub(in super::super) labels: ThinVec<AmbiguousLabel>,
}

/// Keep track of what ambiguous option is used
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum AmbiguousLabel {
    /// An ambiguous amino acid, with the actual amino acid used tracked
    AminoAcid {
        /// Which amino acid is used
        option: AminoAcid,
        /// What location in the sequence are we talking about
        sequence_index: usize,
        /// Peptidoform index
        peptidoform_index: usize,
        /// Peptidoform ion index
        peptidoform_ion_index: usize,
    },
    /// An ambiguous modification, with the actual position
    Modification {
        /// Which ambiguous modification
        id: usize,
        /// Which location
        sequence_index: SequencePosition,
        /// Peptidoform index
        peptidoform_index: usize,
        /// Peptidoform ion index
        peptidoform_ion_index: usize,
    },
    /// The actual charge used, when there are multiple charge carriers
    ChargeCarrier(MolecularFormula),
    /// An intact cross-link
    CrossLinkBound(CrossLinkName),
    /// A broken cross-link, having the name and the stub that was left in its place
    CrossLinkBroken(CrossLinkName, MolecularFormula),
    /// A glycan fragment on a peptide fragment, with the Y breakages that lead to that fragment
    GlycanFragment(Vec<GlycanPosition>),
    /// A glycan fragment on a peptide fragment, with the monosaccharides that make up the fragment
    GlycanFragmentComposition(Vec<(MonoSaccharide, isize)>),
}

/// Any item that has a clearly defined single molecular formula
pub trait Chemical {
    /// Get the molecular formula
    #[allow(dead_code)]
    fn formula(&self) -> MolecularFormula {
        self.formula_inner(SequencePosition::default(), 0)
    }

    /// Get the molecular formula while keeping track of all ambiguous labels
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula;
}

impl<T: Chemical> Chemical for &[T] {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

impl<T: Chemical> Chemical for &Vec<T> {
    fn formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> MolecularFormula {
        self.iter()
            .map(|f| f.formula_inner(sequence_index, peptidoform_index))
            .sum()
    }
}

/// Any item that has a number of potential chemical formulas
pub trait MultiChemical {
    /// Get all possible molecular formulas
    fn formulas(&self) -> Multi<MolecularFormula> {
        self.formulas_inner(SequencePosition::default(), 0, 0)
    }

    /// Get all possible molecular formulas while keeping track of all ambiguous labels
    fn formulas_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Multi<MolecularFormula>;

    /// Get the charge of this chemical, it returns None if no charge is defined.
    #[allow(dead_code)]
    fn charge(&self) -> Option<crate::system::isize::Charge> {
        self.formulas()
            .first()
            .map(MolecularFormula::charge)
            .filter(|c| c.value != 0)
    }

    /// Return a single formula if this `MultiChemical` has only one possible formula
    fn single_formula(&self) -> Option<MolecularFormula> {
        let formulas = self.formulas();
        (formulas.len() == 1).then_some(formulas.to_vec().pop().unwrap())
    }

    /// Return a single formula if this `MultiChemical` has only one possible formula while keeping track of all ambiguous labels
    #[allow(dead_code)]
    fn single_formula_inner(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
        peptidoform_ion_index: usize,
    ) -> Option<MolecularFormula> {
        let formulas =
            self.formulas_inner(sequence_index, peptidoform_index, peptidoform_ion_index);
        (formulas.len() == 1).then_some(formulas.to_vec().pop().unwrap())
    }
}

impl MolecularFormula {
    /// Create a new molecular formula, if the chosen isotopes are not valid it returns None
    pub fn new(
        elements: &[(Element, Option<NonZeroU16>, i32)],
        labels: &[AmbiguousLabel],
    ) -> Option<Self> {
        if elements.iter().any(|e| !e.0.is_valid(e.1)) {
            None
        } else {
            let result = Self {
                elements: elements.into(),
                additional_mass: 0.0.into(),
                labels: labels.into(),
            };
            Some(result.simplify())
        }
    }

    // The elements will be sorted on element/isotope and deduplicated
    #[must_use]
    fn simplify(mut self) -> Self {
        self.elements.retain(|el| el.2 != 0);
        self.elements.sort_by(|a, b| {
            if a.0 == b.0 {
                // If the elements are the same sort on the isotope number
                a.1.cmp(&b.1)
            } else {
                a.0.cmp(&b.0)
            }
        });
        // Deduplicate
        let mut max = self.elements.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = self.elements[index];
            let next = self.elements[index + 1];
            if this.0 == next.0 && this.1 == next.1 {
                self.elements[index].2 += next.2;
                self.elements.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        self.elements
            .retain(|el: &(Element, Option<std::num::NonZero<u16>>, i32)| el.2 != 0);
        self
    }

    /// Get an empty molecular formula with only a mass of unspecified origin
    #[must_use]
    pub fn with_additional_mass(additional_mass: f64) -> Self {
        Self {
            elements: ThinVec::new(),
            additional_mass: OrderedFloat(additional_mass),
            labels: ThinVec::new(),
        }
    }

    /// Add the following labels to this formula
    #[must_use]
    pub fn with_labels(mut self, labels: &[AmbiguousLabel]) -> Self {
        self.labels.extend_from_slice(labels);
        self
    }

    /// Add the following label to this formula
    #[must_use]
    pub(crate) fn with_label(mut self, label: AmbiguousLabel) -> Self {
        self.labels.push(label);
        self
    }

    /// The labels of sources of ambiguity/multiplicity
    pub fn labels(&self) -> &[AmbiguousLabel] {
        &self.labels
    }

    /// Add the given element to this formula (while keeping it ordered and simplified).
    /// If the isotope for the added element is not valid it returns `false`. It also
    /// does so if the addition of this element overflows the maximal number of atoms
    /// in a formula.
    #[must_use]
    pub fn add(&mut self, element: (Element, Option<NonZeroU16>, i32)) -> bool {
        if element.0.is_valid(element.1) {
            let mut index = 0;
            let mut done = false;
            let (el, i, n) = element;
            while !done {
                let base = self.elements.get(index).copied();
                if let Some((re, ri, _)) = base {
                    if el > re || (el == re && i > ri) {
                        index += 1;
                    } else if el == re && i == ri {
                        if let Some(n) = self.elements[index].2.checked_add(n) {
                            self.elements[index].2 = n;
                        } else {
                            return false;
                        }
                        done = true;
                    } else {
                        self.elements.insert(index, (el, i, n));
                        done = true;
                    }
                } else {
                    self.elements.push((el, i, n));
                    done = true;
                }
            }
            true
        } else {
            false
        }
    }

    /// Add the given monoisotopic weight to this formula
    pub fn add_mass(&mut self, mass: OrderedFloat<f64>) {
        self.additional_mass += mass;
    }

    /// Get the elements making this formula
    pub fn elements(&self) -> &[(Element, Option<NonZeroU16>, i32)] {
        &self.elements
    }

    /// Get the additional mass of this formula
    pub const fn additional_mass(&self) -> OrderedFloat<f64> {
        self.additional_mass
    }

    /// Create a new molecular formula with the given global isotope modifications. If the given isotope is not valid for this element it returns `None`.
    #[must_use]
    pub fn with_global_isotope_modifications(
        &self,
        substitutions: &[(Element, Option<NonZeroU16>)],
    ) -> Option<Self> {
        if substitutions.iter().all(|e| e.0.is_valid(e.1)) {
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
                labels: self.labels.clone(),
            };
            Some(result.simplify())
        } else {
            None
        }
    }

    /// Get the number of electrons (the only charged species, any ionic species is saved as that element +/- the correct number of electrons).
    /// The inverse of that number is given as the charge.
    pub fn charge(&self) -> crate::system::isize::Charge {
        -self
            .elements
            .iter()
            .find(|el| el.0 == Element::Electron)
            .map_or_else(crate::system::isize::Charge::default, |el| {
                crate::system::isize::Charge::new::<crate::system::charge::e>(el.2 as isize)
            })
    }

    /// Check if the formula is empty (no elements and no additional mass)
    pub fn is_empty(&self) -> bool {
        self.elements.is_empty() && self.additional_mass == 0.0
    }

    /// The generic backbone to do the Hill notation sorting
    #[allow(dead_code)]
    pub(in super::super) fn hill_notation_generic(
        &self,
        f: impl Fn(&(Element, Option<NonZeroU16>, i32), &mut String),
    ) -> String {
        let mut buffer = String::new();
        if let Some(carbon) = self
            .elements
            .iter()
            .find(|e| e.0 == Element::C && e.1.is_none())
        {
            if carbon.2 != 0 {
                f(carbon, &mut buffer);
            }
            if let Some(hydrogen) = self
                .elements
                .iter()
                .find(|e| e.0 == Element::H && e.1.is_none())
                && hydrogen.2 != 0
            {
                f(hydrogen, &mut buffer);
            }
            for element in self.elements.iter().filter(|e| {
                !((e.0 == Element::H || e.0 == Element::C || e.0 == Element::Electron)
                    && e.1.is_none())
            }) {
                if element.2 != 0 {
                    f(element, &mut buffer);
                }
            }
        } else {
            for element in &self.elements {
                if element.2 != 0 && element.0 != Element::Electron {
                    f(element, &mut buffer);
                }
            }
        }
        if self.additional_mass != 0.0 {
            write!(&mut buffer, "{:+}", self.additional_mass).unwrap();
        }
        if self.charge().value != 0 {
            write!(&mut buffer, ":z{:+}", self.charge().value).unwrap();
        }
        buffer
    }
}

impl Neg for &MolecularFormula {
    type Output = MolecularFormula;
    fn neg(self) -> Self::Output {
        let mut res = self.clone();
        for element in &mut res.elements {
            element.2 = -element.2;
        }
        res
    }
}

impl Neg for MolecularFormula {
    type Output = Self;
    fn neg(mut self) -> Self::Output {
        for element in &mut self.elements {
            element.2 = -element.2;
        }
        self
    }
}

impl Add<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn add(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        result.labels.extend_from_slice(&rhs.labels);
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass += rhs.additional_mass;

        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.elements[index_result].2 += n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (el, i, n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((el, i, n));
                index_rhs += 1;
            }
        }
        result.elements.retain(|el| el.2 != 0);
        result
    }
}

impl Sub<&MolecularFormula> for &MolecularFormula {
    type Output = MolecularFormula;
    fn sub(self, rhs: &MolecularFormula) -> Self::Output {
        let mut result = (*self).clone();
        result.labels.extend_from_slice(&rhs.labels);
        let mut index_result = 0;
        let mut index_rhs = 0;
        result.additional_mass -= rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_result < result.elements.len() {
                let (re, ri, _) = result.elements[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    result.elements[index_result].2 -= n;
                    index_rhs += 1;
                } else {
                    result.elements.insert(index_result, (el, i, -n));
                    index_rhs += 1;
                }
            } else {
                result.elements.push((el, i, -n));
                index_rhs += 1;
            }
        }
        result.elements.retain(|el| el.2 != 0);
        result
    }
}

impl Mul<&i32> for &MolecularFormula {
    type Output = MolecularFormula;
    fn mul(self, rhs: &i32) -> Self::Output {
        self.checked_mul_i32(*rhs)
            .expect("Overflow in multiplying MolecularFormula")
    }
}

impl Mul<&u16> for &MolecularFormula {
    type Output = MolecularFormula;
    fn mul(self, rhs: &u16) -> Self::Output {
        self.checked_mul_u16(*rhs)
            .expect("Overflow in multiplying MolecularFormula")
    }
}

impl MolecularFormula {
    /// Do checked multiplication to handle overflows gracefully
    pub fn checked_mul_i32(&self, rhs: i32) -> Option<Self> {
        Some(Self {
            additional_mass: self.additional_mass * f64::from(rhs),
            elements: self
                .elements
                .iter()
                .copied()
                .map(|part| Some((part.0, part.1, part.2.checked_mul(rhs)?)))
                .collect::<Option<_>>()?,
            labels: self.labels.clone(),
        })
    }

    /// Do checked multiplication to handle overflows gracefully
    pub fn checked_mul_u16(&self, rhs: u16) -> Option<Self> {
        Some(Self {
            additional_mass: self.additional_mass * f64::from(rhs),
            elements: self
                .elements
                .iter()
                .copied()
                .map(|part| Some((part.0, part.1, part.2.checked_mul(i32::from(rhs))?)))
                .collect::<Option<_>>()?,
            labels: self.labels.clone(),
        })
    }
}

impl Mul<&i8> for &MolecularFormula {
    type Output = MolecularFormula;
    fn mul(self, rhs: &i8) -> Self::Output {
        MolecularFormula {
            additional_mass: self.additional_mass * f64::from(*rhs),
            elements: self
                .elements
                .iter()
                .copied()
                .map(|part| (part.0, part.1, part.2 * i32::from(*rhs)))
                .collect(),
            labels: self.labels.clone(),
        }
    }
}

impl_binop_ref_cases!(impl Add, add for MolecularFormula, MolecularFormula, MolecularFormula);
impl_binop_ref_cases!(impl Sub, sub for MolecularFormula, MolecularFormula, MolecularFormula);
impl_binop_ref_cases!(impl Mul, mul for MolecularFormula, i32, MolecularFormula);
impl_binop_ref_cases!(impl Mul, mul for MolecularFormula, u16, MolecularFormula);
impl_binop_ref_cases!(impl Mul, mul for MolecularFormula, i8, MolecularFormula);

impl AddAssign<&Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: &Self) {
        let mut index_self = 0;
        let mut index_rhs = 0;
        self.additional_mass += rhs.additional_mass;
        self.labels.extend_from_slice(&rhs.labels);
        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_self < self.elements.len() {
                let (re, ri, _) = self.elements[index_self];
                if el > re || (el == re && i > ri) {
                    index_self += 1;
                } else if el == re && i == ri {
                    self.elements[index_self].2 += n;
                    index_rhs += 1;
                } else {
                    self.elements.insert(index_self, (el, i, n));
                    index_rhs += 1;
                }
            } else {
                self.elements.push((el, i, n));
                index_rhs += 1;
            }
        }
    }
}

impl SubAssign<&Self> for MolecularFormula {
    fn sub_assign(&mut self, rhs: &Self) {
        self.labels.extend_from_slice(&rhs.labels);
        let mut index_result = 0;
        let mut index_rhs = 0;
        self.additional_mass -= rhs.additional_mass;
        while index_rhs < rhs.elements.len() {
            let (el, i, n) = rhs.elements[index_rhs];
            if index_result < self.elements.len() {
                let (re, ri, _) = self.elements[index_result];
                if el > re || (el == re && i > ri) {
                    index_result += 1;
                } else if el == re && i == ri {
                    self.elements[index_result].2 -= n;
                    index_rhs += 1;
                } else {
                    self.elements.insert(index_result, (el, i, -n));
                    index_rhs += 1;
                }
            } else {
                self.elements.push((el, i, -n));
                index_rhs += 1;
            }
        }
        self.elements.retain(|el| el.2 != 0);
    }
}

impl AddAssign<Self> for MolecularFormula {
    fn add_assign(&mut self, rhs: Self) {
        *self += &rhs;
    }
}

impl SubAssign<Self> for MolecularFormula {
    fn sub_assign(&mut self, rhs: Self) {
        *self -= &rhs;
    }
}

impl std::iter::Sum<Self> for MolecularFormula {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res += v);
        res
    }
}

#[macro_export]
/// Easily define molecular formulas using the following syntax: `<element> <num>` or `[<isotope> <element> <num>]`.
/// The spaces are required by the Rust compiler.
/// ```
/// # use rustyms::*;
/// molecular_formula!(C 12 [13 C 1] H 24);
/// ```
/// # Panics
/// It panics if the defined molecular formula is not valid. A formula is not valid if not existing isotopes are used
/// or if an element is used that does not have a defined molecular weight (does not have natural abundance).
macro_rules! molecular_formula {
    ($($tail:tt)*) => {
        $crate::formula_internal!([$($tail)*] -> [])
    };
}

#[doc(hidden)]
#[macro_export]
/// Internal code for the [`molecular_formula`] macro.
macro_rules! formula_internal {
    ([$e:ident $n:literal $($tail:tt)*] -> [$($output:tt)*]) => {
        $crate::formula_internal!([$($tail)*] -> [$($output)*($crate::chemistry::Element::$e, None, $n),])
    };
    ([[$i:literal $e:ident $n:literal] $($tail:tt)*] -> [$($output:tt)*]) => {
        $crate::formula_internal!([$($tail)*] -> [$($output)*($crate::chemistry::Element::$e, Some(std::num::NonZeroU16::new($i).unwrap()), $n),])
    };
    ([$e:ident $n:expr] -> [$($output:tt)*]) =>{
        $crate::formula_internal!([] -> [$($output)*($crate::chemistry::Element::$e, None, $n),])
    };
    ([[$i:literal $e:ident] $n:expr] -> [$($output:tt)*]) =>{
        $crate::formula_internal!([] -> [$($output)*($crate::chemistry::Element::$e, Some(std::num::NonZeroU16::new($i).unwrap()), $n),])
    };
    ([] -> [$($output:tt)*]) =>{
        $crate::chemistry::MolecularFormula::new(&[$($output)*], &[]).unwrap()
    };
    ([($($l:expr),*)] -> [$($output:tt)*]) =>{
        $crate::chemistry::MolecularFormula::new(&[$($output)*], &[$($l),*]).unwrap()
    };
}

/// A label for a satellite ion, none for most amino acids but a or b for Thr and Ile
#[derive(
    Clone, Copy, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize, Default,
)]
pub enum SatelliteLabel {
    /// No label needed
    #[default]
    None,
    /// Heaviest of the two options
    A,
    /// Lightest of the two options
    B,
}

impl std::fmt::Display for SatelliteLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::None => "",
                Self::A => "a",
                Self::B => "b",
            }
        )
    }
}
