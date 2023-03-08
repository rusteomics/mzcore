use crate::system::{f64::MassOverCharge, mass_over_charge::mz};

pub struct Model {
    pub a: Location,
    pub b: Location,
    pub c: Location,
    pub d: Location,
    pub v: Location,
    pub w: Location,
    pub x: Location,
    pub y: Location,
    pub z: Location,
    pub ppm: MassOverCharge,
}
#[allow(clippy::struct_excessive_bools)]
pub struct PossibleIons {
    pub a: bool,
    pub b: bool,
    pub c: bool,
    pub d: bool,
    pub v: bool,
    pub w: bool,
    pub x: bool,
    pub y: bool,
    pub z: bool,
}

impl PossibleIons {
    pub fn size_upper_bound(&self) -> usize {
        usize::from(self.a)
            + usize::from(self.b)
            + usize::from(self.c)
            + usize::from(self.d) * 2
            + usize::from(self.v)
            + usize::from(self.w) * 2
            + usize::from(self.x)
            + usize::from(self.y)
            + usize::from(self.z) * 2
    }
}

impl Model {
    pub const fn ions(&self, index: usize, length: usize) -> PossibleIons {
        PossibleIons {
            a: self.a.possible(index, length),
            b: self.b.possible(index, length),
            c: self.c.possible(index, length),
            d: self.d.possible(index, length),
            v: self.v.possible(index, length),
            w: self.w.possible(index, length),
            x: self.x.possible(index, length),
            y: self.y.possible(index, length),
            z: self.z.possible(index, length),
        }
    }

    #[allow(clippy::too_many_arguments, clippy::many_single_char_names)]
    pub fn new(
        a: Location,
        b: Location,
        c: Location,
        d: Location,
        v: Location,
        w: Location,
        x: Location,
        y: Location,
        z: Location,
        ppm: MassOverCharge,
    ) -> Self {
        Self {
            a,
            b,
            c,
            d,
            v,
            w,
            x,
            y,
            z,
            ppm,
        }
    }

    pub fn all() -> Self {
        Self {
            a: Location::SkipN(1),
            b: Location::SkipN(1),
            c: Location::SkipN(1),
            d: Location::SkipN(1),
            v: Location::SkipN(1),
            w: Location::All,
            x: Location::All,
            y: Location::All,
            z: Location::All,
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn ethcd() -> Self {
        Self {
            a: Location::None,
            b: Location::SkipNC(1, 1),
            c: Location::SkipNC(1, 1),
            d: Location::None,
            v: Location::None,
            w: Location::SkipN(1),
            x: Location::None,
            y: Location::SkipN(1),
            z: Location::SkipN(1),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn cid_hcd() -> Self {
        Self {
            a: Location::TakeN { skip: 1, take: 1 },
            b: Location::SkipNC(1, 1),
            c: Location::None,
            d: Location::TakeN { skip: 1, take: 1 },
            v: Location::None,
            w: Location::None,
            x: Location::None,
            y: Location::SkipC(1),
            z: Location::None,
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }

    pub fn etcid() -> Self {
        Self {
            a: Location::None,
            b: Location::SkipNC(1, 1),
            c: Location::SkipNC(1, 1),
            d: Location::None,
            v: Location::None,
            w: Location::SkipN(1),
            x: Location::None,
            y: Location::SkipN(1),
            z: Location::SkipN(1),
            ppm: MassOverCharge::new::<mz>(20.0),
        }
    }
}

pub enum Location {
    SkipN(usize),
    SkipNC(usize, usize),
    TakeN { skip: usize, take: usize },
    SkipC(usize),
    TakeC(usize),
    All,
    None,
}

impl Location {
    pub const fn possible(&self, index: usize, length: usize) -> bool {
        match self {
            Self::SkipN(n) => index >= *n,
            Self::SkipNC(n, c) => index >= *n && length - index > *c,
            Self::TakeN { skip, take } => index >= *skip && index < *skip + *take,
            Self::SkipC(n) => length - index > *n,
            Self::TakeC(n) => length - index <= *n,
            Self::All => true,
            Self::None => false,
        }
    }
}
