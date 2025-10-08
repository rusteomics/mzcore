use crate::glycan::{
    BaseSugar, Configuration, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, MonoSaccharide,
    NonoseIsomer, PentoseIsomer,
};

impl MonoSaccharide {
    /// Get the shape, colour, inner modifications, and outer modifications for this monosaccharide.
    pub(super) fn get_shape(&self) -> (Shape, Colour, String, String) {
        // Common substitutions
        let mut nacetyl = 0;
        let mut acid = 0;
        let mut amino = 0;
        let mut deoxy = 0;
        // Additional needed substitutions
        let mut acetyl = 0;
        let mut glycolyl = 0;
        let mut nglycolyl = 0;
        let mut o_carboxy_ethyl = 0;
        let mut inner_modifications = if self.furanose {
            "f".to_string()
        } else {
            String::new()
        };
        if let Some(c) = &self.configuration {
            inner_modifications.push_str(match *c {
                Configuration::D => "D",
                Configuration::L => "L",
                Configuration::DD => "DD",
                Configuration::LL => "LL",
                Configuration::DL => "DL",
                Configuration::LD => "LD",
            });
        }
        let mut outer_modifications = String::new();
        for m in &self.substituents {
            match m {
                GlycanSubstituent::NAcetyl => nacetyl += 1,
                GlycanSubstituent::Acid => acid += 1,
                GlycanSubstituent::Amino => amino += 1,
                GlycanSubstituent::Deoxy => deoxy += 1,
                GlycanSubstituent::Acetyl => acetyl += 1,
                GlycanSubstituent::Glycolyl => glycolyl += 1,
                GlycanSubstituent::OCarboxyEthyl => o_carboxy_ethyl += 1,
                GlycanSubstituent::NGlycolyl => nglycolyl += 1,
                GlycanSubstituent::Didehydro => inner_modifications.push_str("en"),
                GlycanSubstituent::Alcohol => inner_modifications.push('o'), // Missing symbols: an for anhydro, on for lactone, am for lactam
                _ => outer_modifications.push_str(m.notation()),
            }
        }
        let outer_mods = |nacetyl: usize,
                          acid: usize,
                          amino: usize,
                          deoxy: usize,
                          acetyl: usize,
                          glycolyl: usize,
                          nglycolyl: usize,
                          o_carboxy_ethyl: usize| {
            [
                GlycanSubstituent::NAcetyl.notation().repeat(nacetyl),
                GlycanSubstituent::Acid.notation().repeat(acid),
                GlycanSubstituent::Amino.notation().repeat(amino),
                GlycanSubstituent::Deoxy.notation().repeat(deoxy),
                GlycanSubstituent::Acetyl.notation().repeat(acetyl),
                GlycanSubstituent::Glycolyl.notation().repeat(glycolyl),
                GlycanSubstituent::NGlycolyl.notation().repeat(nglycolyl),
                GlycanSubstituent::OCarboxyEthyl
                    .notation()
                    .repeat(o_carboxy_ethyl),
                outer_modifications,
            ]
            .join("")
        };
        match &self.base_sugar {
            BaseSugar::Pentose(isomer) => (
                Shape::Star,
                match isomer {
                    None | Some(PentoseIsomer::Xylulose) => Colour::Background,
                    Some(PentoseIsomer::Arabinose) => Colour::Green,
                    Some(PentoseIsomer::Lyxose) => Colour::Yellow,
                    Some(PentoseIsomer::Xylose) => Colour::Orange,
                    Some(PentoseIsomer::Ribose) => Colour::Pink,
                },
                inner_modifications,
                outer_mods(
                    nacetyl,
                    acid,
                    amino,
                    deoxy,
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Hexose(isomer) => {
                if o_carboxy_ethyl > 0 && nacetyl > 0 {
                    (
                        Shape::Hexagon,
                        Colour::Purple,
                        inner_modifications,
                        outer_mods(
                            nacetyl - 1,
                            acid,
                            amino,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl - 1,
                        ),
                    )
                } else if o_carboxy_ethyl > 0 && nglycolyl > 0 {
                    (
                        Shape::Hexagon,
                        Colour::LightBlue,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl - 1,
                            o_carboxy_ethyl - 1,
                        ),
                    )
                } else if o_carboxy_ethyl > 0 && amino > 0 {
                    (
                        Shape::Hexagon,
                        Colour::Brown,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino - 1,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl - 1,
                        ),
                    )
                } else if deoxy > 1 {
                    let c = match isomer {
                        Some(HexoseIsomer::Glucose) => Colour::Blue,
                        Some(HexoseIsomer::Mannose) => Colour::Green,
                        Some(HexoseIsomer::Galactose) => Colour::Orange,
                        Some(HexoseIsomer::Altrose) => Colour::Pink,
                        Some(HexoseIsomer::Allose) => Colour::Purple,
                        Some(HexoseIsomer::Talose) => Colour::LightBlue,
                        Some(_) | None => Colour::Background,
                    };
                    (
                        Shape::Rectangle,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino,
                            deoxy - 2,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if amino > 1 && deoxy > 0 {
                    (
                        Shape::Hexagon,
                        Colour::Blue,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino - 2,
                            deoxy - 1,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if nacetyl > 0 && deoxy > 0 {
                    let c = match isomer {
                        Some(HexoseIsomer::Glucose) => Colour::Blue,
                        Some(HexoseIsomer::Mannose) => Colour::Green,
                        Some(HexoseIsomer::Galactose) => Colour::Red,
                        Some(HexoseIsomer::Altrose) => Colour::Pink,
                        Some(HexoseIsomer::Talose) => Colour::LightBlue,
                        Some(_) | None => Colour::Background,
                    };
                    (
                        Shape::DividedTriangle,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl - 1,
                            acid,
                            amino,
                            deoxy - 1,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if deoxy > 0 {
                    let c = match isomer {
                        Some(HexoseIsomer::Glucose) => Colour::Blue,
                        Some(HexoseIsomer::Mannose) => Colour::Green,
                        Some(HexoseIsomer::Galactose) => Colour::Red,
                        Some(HexoseIsomer::Gulose) => Colour::Orange,
                        Some(HexoseIsomer::Altrose) => Colour::Pink,
                        Some(HexoseIsomer::Talose) => Colour::LightBlue,
                        Some(_) | None => Colour::Background,
                    };
                    (
                        Shape::Triangle,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino,
                            deoxy - 1,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if acid > 0 || amino > 0 || nacetyl > 0 {
                    let c = match isomer {
                        Some(HexoseIsomer::Glucose) => Colour::Blue,
                        Some(HexoseIsomer::Mannose) => Colour::Green,
                        Some(HexoseIsomer::Galactose) => Colour::Yellow,
                        Some(HexoseIsomer::Gulose) => Colour::Orange,
                        Some(HexoseIsomer::Altrose) => Colour::Pink,
                        Some(HexoseIsomer::Allose) => Colour::Purple,
                        Some(HexoseIsomer::Talose) => Colour::LightBlue,
                        Some(HexoseIsomer::Idose) => Colour::Brown,
                        Some(_) | None => Colour::Background,
                    };
                    let shape = if acid > 0 {
                        Shape::DividedDiamond
                    } else if amino > 0 {
                        Shape::CrossedSquare
                    } else {
                        Shape::Square
                    };
                    (
                        shape,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl - usize::from(shape == Shape::Square),
                            acid - usize::from(shape == Shape::DividedDiamond),
                            amino - usize::from(shape == Shape::CrossedSquare),
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else {
                    let (s, c) = match isomer {
                        None => (Shape::Circle, Colour::Background),
                        Some(HexoseIsomer::Glucose) => (Shape::Circle, Colour::Blue),
                        Some(HexoseIsomer::Mannose) => (Shape::Circle, Colour::Green),
                        Some(HexoseIsomer::Galactose) => (Shape::Circle, Colour::Yellow),
                        Some(HexoseIsomer::Gulose) => (Shape::Circle, Colour::Orange),
                        Some(HexoseIsomer::Altrose) => (Shape::Circle, Colour::Pink),
                        Some(HexoseIsomer::Allose) => (Shape::Circle, Colour::Purple),
                        Some(HexoseIsomer::Talose) => (Shape::Circle, Colour::LightBlue),
                        Some(HexoseIsomer::Idose) => (Shape::Circle, Colour::Brown),
                        Some(HexoseIsomer::Psicose) => (Shape::Pentagon, Colour::Pink),
                        Some(HexoseIsomer::Fructose) => (Shape::Pentagon, Colour::Green),
                        Some(HexoseIsomer::Sorbose) => (Shape::Pentagon, Colour::Orange),
                        Some(HexoseIsomer::Tagatose) => (Shape::Pentagon, Colour::Yellow),
                    };
                    (
                        s,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                }
            }
            BaseSugar::Heptose(Some(HeptoseIsomer::GlyceroMannoHeptopyranose)) => (
                Shape::Hexagon,
                Colour::Green,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    acid,
                    amino,
                    deoxy,
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Heptose(None) if acid > 1 && deoxy > 0 => (
                Shape::Hexagon,
                Colour::Orange,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    acid - 2,
                    amino,
                    deoxy - 1,
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Octose if acid > 0 && deoxy > 0 => (
                Shape::Hexagon,
                Colour::Yellow,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    acid - 1,
                    amino,
                    deoxy - 1,
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Nonose(isomer) if acid > 0 && amino > 0 => {
                if amino > 1 && deoxy > 1 {
                    (
                        Shape::FlatDiamond,
                        match isomer {
                            Some(NonoseIsomer::Pse) => Colour::Green,
                            Some(NonoseIsomer::Leg) => Colour::Yellow,
                            Some(NonoseIsomer::ELeg) => Colour::LightBlue,
                            Some(NonoseIsomer::Aci) => Colour::Pink,
                            _ => Colour::Background,
                        },
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid - 1,
                            amino - 2,
                            deoxy - 2,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else {
                    let colour = if deoxy > 0 {
                        if *isomer == Some(NonoseIsomer::Kdn) {
                            Colour::Green
                        } else {
                            Colour::Red
                        }
                    } else if acetyl > 0 {
                        Colour::Purple
                    } else if glycolyl > 0 {
                        Colour::LightBlue
                    } else {
                        Colour::Brown
                    };
                    (
                        Shape::Diamond,
                        colour,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid - 1,
                            amino - 1,
                            deoxy - usize::from(colour == Colour::Red || colour == Colour::Green),
                            acetyl - usize::from(colour == Colour::Purple),
                            glycolyl - usize::from(colour == Colour::LightBlue),
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                }
            }
            _ => (
                Shape::Hexagon,
                Colour::Background,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    acid,
                    amino,
                    deoxy,
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
        }
    }
}

/// All colours from Symbol Nomenclature For Glycans (SNFG)
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Colour {
    Background,
    Blue,
    Green,
    Yellow,
    Orange,
    Pink,
    Purple,
    LightBlue,
    Brown,
    Red,
}

impl Colour {
    /// Represented as bytes 0..=255
    pub(super) const fn rgb(self) -> [u8; 3] {
        match self {
            Self::Background => [255, 255, 255],
            Self::Blue => [0, 144, 188],
            Self::Green => [0, 166, 81],
            Self::Yellow => [255, 212, 0],
            Self::Orange => [244, 121, 32],
            Self::Pink => [246, 158, 161],
            Self::Purple => [165, 67, 153],
            Self::LightBlue => [143, 204, 233],
            Self::Brown => [161, 122, 77],
            Self::Red => [237, 28, 36],
        }
    }
}

/// All symbols from Symbol Nomenclature For Glycans (SNFG)
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum Shape {
    Circle,
    Square,
    CrossedSquare,
    DividedDiamond,
    Triangle,
    LeftPointingTriangle,
    RightPointingTriangle,
    DividedTriangle,
    Rectangle,
    Star,
    Diamond,
    FlatDiamond,
    Hexagon,
    Pentagon,
}

impl Shape {
    /// The height of a symbol as ratio to the width
    pub(super) const fn height(self) -> f32 {
        match self {
            Self::Rectangle | Self::FlatDiamond | Self::Hexagon => 0.5,
            _ => 1.0,
        }
    }
}
