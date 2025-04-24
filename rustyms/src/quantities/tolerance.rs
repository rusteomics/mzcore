use std::{fmt::Display, str::FromStr};

use itertools::Itertools;
use serde::{Deserialize, Serialize};
use uom::fmt::DisplayStyle;

use crate::{
    quantities::Multi,
    system::{da, Mass, MassOverCharge, OrderedMass, OrderedRatio, Ratio},
};

/// A tolerance around a given unit for searching purposes
#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
pub enum Tolerance<T> {
    /// A relative search tolerance
    Relative(OrderedRatio),
    /// An absolute tolerance defined by a constant offset from the unit (bounds are unit - tolerance, unit + tolerance)
    Absolute(T),
}

impl<T> Tolerance<T> {
    /// Create a new ppm value
    pub fn new_ppm(value: f64) -> Self {
        Self::Relative(Ratio::new::<crate::system::ratio::ppm>(value).into())
    }

    /// Create a new relative value
    pub fn new_relative(value: Ratio) -> Self {
        Self::Relative(value.into())
    }

    /// Create a new absolute value
    pub fn new_absolute(value: impl Into<T>) -> Self {
        Self::Absolute(value.into())
    }

    /// Convert this tolerance into another absolute type.
    pub fn convert<O: From<T>>(self) -> Tolerance<O> {
        match self {
            Self::Relative(r) => Tolerance::Relative(r),
            Self::Absolute(a) => Tolerance::Absolute(a.into()),
        }
    }
}

impl<T> Tolerance<T>
where
    T: std::ops::Mul<Ratio, Output = T>
        + std::ops::Sub<T, Output = T>
        + std::ops::Add<T, Output = T>
        + Copy,
{
    /// Find the bounds around a given value for this tolerance
    pub fn bounds(&self, value: impl Into<T>) -> (T, T) {
        let value = value.into();
        match self {
            Self::Relative(tolerance) => (
                value
                    * (Ratio::new::<crate::system::ratio::fraction>(1.0) - tolerance.into_inner()),
                value
                    * (Ratio::new::<crate::system::ratio::fraction>(1.0) + tolerance.into_inner()),
            ),
            Self::Absolute(tolerance) => (value - *tolerance, value + *tolerance),
        }
    }
}

impl<T: Display> Display for Tolerance<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Absolute(value) => format!("{value} abs"),
                Self::Relative(tolerance) => format!("{} rel", tolerance.value),
            }
        )
    }
}

impl Display for Tolerance<Mass> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Absolute(value) => format!(
                    "{}",
                    value.into_format_args(crate::system::mass::dalton, DisplayStyle::Abbreviation)
                ),
                Self::Relative(tolerance) => format!(
                    "{}",
                    tolerance
                        .into_format_args(crate::system::ratio::ppm, DisplayStyle::Abbreviation)
                ),
            }
        )
    }
}

impl FromStr for Tolerance<Mass> {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let num_str = String::from_utf8(
            s.bytes()
                .take_while(|c| {
                    c.is_ascii_digit()
                        || *c == b'.'
                        || *c == b'-'
                        || *c == b'+'
                        || *c == b'e'
                        || *c == b'E'
                })
                .collect::<Vec<_>>(),
        )
        .map_err(|_| ())?;
        let num = num_str.parse::<f64>().map_err(|_| ())?;
        match s[num_str.len()..].trim() {
            "ppm" => Ok(Self::Relative(
                Ratio::new::<crate::system::ratio::ppm>(num).into(),
            )),
            "da" => Ok(Self::Absolute(da(num))),
            _ => Err(()),
        }
    }
}

impl<T> TryFrom<&str> for Tolerance<T>
where
    Self: FromStr,
{
    type Error = <Self as std::str::FromStr>::Err;
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        value.parse()
    }
}

/// Check if two values are within the specified tolerance from each other.
pub trait WithinTolerance<A, B> {
    /// Check if two values are within the specified tolerance from each other.
    fn within(&self, a: &A, b: &B) -> bool;
}

impl WithinTolerance<MassOverCharge, MassOverCharge> for Tolerance<MassOverCharge> {
    fn within(&self, a: &MassOverCharge, b: &MassOverCharge) -> bool {
        match self {
            Self::Absolute(tol) => (a.value - b.value).abs() <= tol.value,
            Self::Relative(tolerance) => a.ppm(*b) <= tolerance.into_inner(),
        }
    }
}

impl WithinTolerance<Mass, Mass> for Tolerance<Mass> {
    fn within(&self, a: &Mass, b: &Mass) -> bool {
        match self {
            Self::Absolute(tol) => (a.value - b.value).abs() <= tol.value,
            Self::Relative(tolerance) => a.ppm(*b) <= tolerance.into_inner(),
        }
    }
}

impl WithinTolerance<Multi<Mass>, Multi<Mass>> for Tolerance<Mass> {
    fn within(&self, a: &Multi<Mass>, b: &Multi<Mass>) -> bool {
        a.iter()
            .cartesian_product(b.iter())
            .any(|(a, b)| self.within(a, b))
    }
}

impl WithinTolerance<Multi<Mass>, Mass> for Tolerance<Mass> {
    fn within(&self, a: &Multi<Mass>, b: &Mass) -> bool {
        a.iter().any(|a| self.within(a, b))
    }
}

impl WithinTolerance<Mass, Multi<Mass>> for Tolerance<Mass> {
    fn within(&self, a: &Mass, b: &Multi<Mass>) -> bool {
        b.iter().any(|b| self.within(a, b))
    }
}

impl WithinTolerance<Mass, Mass> for Tolerance<OrderedMass> {
    fn within(&self, a: &Mass, b: &Mass) -> bool {
        match self {
            Self::Absolute(tol) => (a.value - b.value).abs() <= tol.value,
            Self::Relative(tolerance) => a.ppm(*b) <= tolerance.into_inner(),
        }
    }
}

impl WithinTolerance<Multi<Mass>, Multi<Mass>> for Tolerance<OrderedMass> {
    fn within(&self, a: &Multi<Mass>, b: &Multi<Mass>) -> bool {
        a.iter()
            .cartesian_product(b.iter())
            .any(|(a, b)| self.within(a, b))
    }
}

impl WithinTolerance<Multi<Mass>, Mass> for Tolerance<OrderedMass> {
    fn within(&self, a: &Multi<Mass>, b: &Mass) -> bool {
        a.iter().any(|a| self.within(a, b))
    }
}

impl WithinTolerance<Mass, Multi<Mass>> for Tolerance<OrderedMass> {
    fn within(&self, a: &Mass, b: &Multi<Mass>) -> bool {
        b.iter().any(|b| self.within(a, b))
    }
}
