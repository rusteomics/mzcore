use std::ops::{Add, AddAssign, Deref, Mul, MulAssign, Neg, Sub};

use itertools::{Itertools, MinMaxResult};
use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{AmbiguousLabel, MolecularFormula},
    system::OrderedMass,
};

/// A collection of potentially multiple of the generic type, it is used be able to easily
/// combine multiple of this multi struct into all possible combinations.
#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct Multi<M>(Vec<M>);

impl<M> Multi<M> {
    /// Get the single value if there is only one value, None otherwise
    #[must_use]
    pub fn single(mut self) -> Option<M> {
        if self.len() == 1 { self.0.pop() } else { None }
    }
}

impl<M: Eq + std::hash::Hash + Clone> Multi<M> {
    /// Get all unique values
    #[must_use]
    pub fn unique(&self) -> Self {
        self.0.iter().unique().cloned().collect()
    }
}

impl<'a, M> Multi<M>
where
    &'a M: 'a,
    OrderedMass: From<&'a M>,
{
    /// Get the bounds for the mass
    pub fn mass_bounds(&'a self) -> MinMaxResult<&'a M> {
        self.0.iter().minmax_by_key(|f| OrderedMass::from(f))
    }
}

impl<M> Deref for Multi<M> {
    type Target = [M];
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<M: Default> Default for Multi<M> {
    // Default is one empty M to make the cartesian product with a default return useful results
    fn default() -> Self {
        Self(vec![M::default()])
    }
}

impl<'a, M> Neg for &'a Multi<M>
where
    &'a M: Neg<Output = M> + 'a,
{
    type Output = Multi<M>;
    fn neg(self) -> Self::Output {
        self.0.iter().map(|f| -f).collect()
    }
}

impl<M> Neg for Multi<M>
where
    M: Neg<Output = M> + Clone,
{
    type Output = Self;
    fn neg(self) -> Self::Output {
        self.0.iter().cloned().map(|f| -f).collect()
    }
}

impl<'a, M> Add<&'a M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Adds this M to all Ms in the multi
    fn add(self, rhs: &M) -> Self::Output {
        Multi(self.0.iter().map(|m| rhs.clone() + m).collect())
    }
}

impl<'a, M> Add<&'a M> for Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: &'a M) -> Self::Output {
        self.0.iter().cloned().map(|m| m + rhs).collect()
    }
}

impl<'a, M> Add<M> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        Multi(self.0.iter().map(|m| rhs.clone() + m).collect())
    }
}

impl<M> Add<M> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Adds this formula to all formulas in the multi formula
    fn add(self, rhs: M) -> Self::Output {
        self.0.iter().cloned().map(|m| m + rhs.clone()).collect()
    }
}

impl<M> AddAssign<M> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    /// Adds this formula to all formulas in the multi formula
    fn add_assign(&mut self, rhs: M) {
        *self = Self(self.0.iter().cloned().map(|m| m + rhs.clone()).collect());
    }
}

impl<M> AddAssign<M> for &mut Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    /// Adds this formula to all formulas in the multi formula
    fn add_assign(&mut self, rhs: M) {
        **self = Multi(self.0.iter().cloned().map(|m| m + rhs.clone()).collect());
    }
}

impl<'a, M> AddAssign<&'a M> for &mut Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    /// Adds this formula to all formulas in the multi formula
    fn add_assign(&mut self, rhs: &'a M) {
        **self = Multi(self.0.iter().cloned().map(|m| m + rhs).collect());
    }
}

impl<'a, M> Sub<&'a M> for &'a Multi<M>
where
    &'a M: Sub<M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &M) -> Self::Output {
        Multi(self.0.iter().map(|m| m - rhs.clone()).collect())
    }
}

impl<'a, M> Sub<&'a M> for Multi<M>
where
    M: Sub<&'a M, Output = M> + Clone + 'a,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: &'a M) -> Self::Output {
        self.0.iter().cloned().map(|m| m - rhs).collect()
    }
}

impl<'a, M> Sub<M> for &'a Multi<M>
where
    &'a M: Sub<M, Output = M> + 'a,
    M: Clone,
{
    type Output = Multi<M>;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        self.0.iter().map(|m| m - rhs.clone()).collect()
    }
}

impl<M> Sub<M> for Multi<M>
where
    M: Sub<M, Output = M> + Clone,
{
    type Output = Self;
    /// Subtracts this formula from all formulas in the multi formula
    fn sub(self, rhs: M) -> Self::Output {
        self.0.iter().cloned().map(|m| m - rhs.clone()).collect()
    }
}

impl<'a, M> Mul<&'a Multi<M>> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: &Multi<M>) -> Self::Output {
        Multi(
            self.0
                .iter()
                .cartesian_product(rhs.0.iter())
                .map(|(a, b)| b.clone() + a)
                .collect(),
        )
    }
}

impl<M> Mul<&Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: &Self) -> Self::Output {
        self.0
            .iter()
            .cloned()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<'a, M> Mul<Multi<M>> for &'a Multi<M>
where
    M: Add<&'a M, Output = M> + Clone,
{
    type Output = Multi<M>;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Multi<M>) -> Self::Output {
        self.0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<M> Mul<Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    type Output = Self;
    /// Cartesian product between the two multi formulas
    fn mul(self, rhs: Self) -> Self::Output {
        self.0
            .iter()
            .cloned()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a)
            .collect()
    }
}

impl<'a, M> MulAssign<&'a Self> for Multi<M>
where
    M: Add<M, Output = M> + 'a + Clone,
{
    fn mul_assign(&mut self, rhs: &Self) {
        let new_data = self
            .0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a.clone())
            .collect();
        self.0 = new_data;
    }
}

impl<M> MulAssign<Self> for Multi<M>
where
    M: Add<M, Output = M> + Clone,
{
    fn mul_assign(&mut self, rhs: Self) {
        let new_data = self
            .0
            .iter()
            .cartesian_product(rhs.0.iter())
            .map(|(a, b)| b.clone() + a.clone())
            .collect();
        self.0 = new_data;
    }
}

impl<M: Default> std::iter::Sum<Self> for Multi<M>
where
    Self: MulAssign<Self>,
{
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut res = Self::default();
        iter.for_each(|v| res *= v);
        res
    }
}

impl<M> From<M> for Multi<M> {
    fn from(value: M) -> Self {
        Self(vec![value])
    }
}

impl<M: Clone> From<&M> for Multi<M> {
    fn from(value: &M) -> Self {
        Self(vec![value.clone()])
    }
}

impl<M> From<Vec<M>> for Multi<M> {
    fn from(value: Vec<M>) -> Self {
        Self(value)
    }
}

impl<M: Clone> From<&[M]> for Multi<M> {
    fn from(value: &[M]) -> Self {
        Self(value.into())
    }
}

impl<'a, M> From<&'a [Self]> for Multi<M>
where
    Self: Mul<&'a Self, Output = Self> + 'a,
    M: Default + Clone + Add<M, Output = M>,
{
    /// Get all potential combination from a series of multi elements. If the series is empty it returns the default element.
    fn from(value: &'a [Self]) -> Self {
        value.iter().fold(Self::default(), |acc, a: &Self| {
            if a.len() == 1 {
                acc.iter().cloned().map(|ac| ac + a.0[0].clone()).collect()
            } else {
                acc * a
            }
        })
    }
}

impl<M> FromIterator<M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = M>>(iter: T) -> Self {
        Self(iter.into_iter().collect())
    }
}

impl<'a, M: Clone + 'a> FromIterator<&'a M> for Multi<M> {
    fn from_iter<T: IntoIterator<Item = &'a M>>(iter: T) -> Self {
        Self(iter.into_iter().cloned().collect())
    }
}

impl Multi<MolecularFormula> {
    pub(crate) fn with_label(self, label: &AmbiguousLabel) -> Self {
        Self(
            self.0
                .iter()
                .cloned()
                .map(|o| o.with_label(label.clone()))
                .collect(),
        )
    }
}
