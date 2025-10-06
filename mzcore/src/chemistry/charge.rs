use std::ops::RangeInclusive;

use context_error::{BasicKind, BoxedError};
use serde::{Deserialize, Serialize};
use serde_json::Value;

use crate::{
    parse_json::{ParseJson, use_serde},
    system::{e, isize::Charge},
};

/// Control what charges are allowed for an ion series. Defined as an inclusive range.
/// Any charge above the precursor charge will result in the quotient time the precursor
/// charge carriers + all options for the remainder within the limits of the precursor
/// charge carriers.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub struct ChargeRange {
    /// Start point
    pub start: ChargePoint,
    /// End point (inclusive)
    pub end: ChargePoint,
}

impl ChargeRange {
    /// Get the number of possible charges for the given precursor charge.
    pub fn len(&self, precursor: Charge) -> usize {
        (self.end.to_absolute(precursor).value - self.start.to_absolute(precursor).value.max(1))
            .unsigned_abs()
    }

    /// Get all possible charges for the given precursor charge.
    pub fn charges(&self, precursor: Charge) -> RangeInclusive<Charge> {
        Charge::new::<e>(self.start.to_absolute(precursor).value.max(1))
            ..=self.end.to_absolute(precursor)
    }

    /// Get all possible charges for the given precursor charge.
    pub fn charges_iter(
        &self,
        precursor: Charge,
    ) -> impl DoubleEndedIterator<Item = Charge> + Clone {
        (self.start.to_absolute(precursor).value.max(1)..=self.end.to_absolute(precursor).value)
            .map(Charge::new::<e>)
    }

    /// Solely single charged
    pub const ONE: Self = Self {
        start: ChargePoint::Absolute(1),
        end: ChargePoint::Absolute(1),
    };
    /// Only the exact precursor charge
    pub const PRECURSOR: Self = Self {
        start: ChargePoint::Relative(0),
        end: ChargePoint::Relative(0),
    };
    /// Range from 1 to the precursor
    pub const ONE_TO_PRECURSOR: Self = Self {
        start: ChargePoint::Absolute(1),
        end: ChargePoint::Relative(0),
    };
}

/// A reference point for charge range definition.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
pub enum ChargePoint {
    /// Relative to the precursor, with the given offset.
    Relative(isize),
    /// Absolute charge.
    Absolute(isize),
}

impl ChargePoint {
    /// Get the absolute charge of this charge point given a precursor charge
    fn to_absolute(self, precursor: Charge) -> Charge {
        match self {
            Self::Absolute(a) => Charge::new::<e>(a),
            Self::Relative(r) => Charge::new::<e>(precursor.value + r),
        }
    }
}

impl ParseJson for ChargeRange {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        use_serde(value)
    }
}
