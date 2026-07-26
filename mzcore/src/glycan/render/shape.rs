use std::{collections::HashMap, fmt::Write};

use itertools::Itertools;

use crate::glycan::{
    BaseSugar, Configuration, GlycanModificationDetails, GlycanSubstituent, HeptoseIsomer,
    HexoseIsomer, MonoSaccharide, NonoseIsomer, PentoseIsomer,
};

impl MonoSaccharide {
    /// Get the shape, colour, inner modifications, and outer modifications for this monosaccharide.
    pub(super) fn get_shape(
        &self,
        modification_details: GlycanModificationDetails,
    ) -> (Shape, Colour, String, String) {
        // Common substitutions
        let mut nacetyl = Vec::new();
        let mut acid = Vec::new();
        let mut amino = Vec::new();
        let mut deoxy = Vec::new();
        // Additional needed substitutions
        let mut acetyl = Vec::new();
        let mut glycolyl = Vec::new();
        let mut nglycolyl = Vec::new();
        let mut o_carboxy_ethyl = Vec::new();
        let mut didehydro = Vec::new();
        let mut alcohol = Vec::new();
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
        let mut outer_modifications: HashMap<&str, Vec<Option<u8>>> = HashMap::new();
        for (m, location) in &self.substituents {
            match m {
                GlycanSubstituent::NAcetyl => nacetyl.push(*location),
                GlycanSubstituent::Acid => acid.push(*location),
                GlycanSubstituent::Amino => amino.push(*location),
                GlycanSubstituent::Deoxy => deoxy.push(*location),
                GlycanSubstituent::Acetyl => acetyl.push(*location),
                GlycanSubstituent::Glycolyl => glycolyl.push(*location),
                GlycanSubstituent::OCarboxyEthyl => o_carboxy_ethyl.push(*location),
                GlycanSubstituent::NGlycolyl => nglycolyl.push(*location),
                GlycanSubstituent::Didehydro => didehydro.push(*location),
                GlycanSubstituent::Alcohol => alcohol.push(*location), /* Missing symbols: */
                // an for anhydro,
                // on for lactone,
                // am for lactam
                _ => outer_modifications.entry(m.notation()).or_default().push(*location),
            }
        }
        let show_mods = |locations: &[Option<u8>], name: &'static str| match modification_details {
            GlycanModificationDetails::Full => {
                if locations.is_empty() {
                    String::new()
                } else {
                    locations.iter().fold(String::new(), |mut acc, l| {
                        if !acc.is_empty() {
                            let _ = acc.write_char(',');
                        }
                        if let Some(l) = l {
                            let _ = write!(&mut acc, "{l}");
                        } else {
                            let _ = acc.write_char('?');
                        }
                        acc
                    }) + name
                }
            }
            GlycanModificationDetails::OnlyContent => name.repeat(locations.len()),
            GlycanModificationDetails::NeverShow => String::new(),
        };
        inner_modifications.push_str(&show_mods(&didehydro, "en"));
        inner_modifications.push_str(&show_mods(&alcohol, "o"));
        let outer_modifications = outer_modifications
            .into_iter()
            .sorted()
            .map(|(name, locations)| show_mods(&locations, name))
            .collect();
        let filter = |mut list: Vec<Option<u8>>, preferred: &[u8]| {
            for p in preferred {
                if let Some(p) = list
                    .iter()
                    .position(|e| *e == Some(*p))
                    .or_else(|| list.iter().position(Option::is_none))
                {
                    list.remove(p);
                } else {
                    let _ = list.pop();
                }
            }
            list
        };
        let outer_mods = |nacetyl: Vec<Option<u8>>,
                          acid: Vec<Option<u8>>,
                          amino: Vec<Option<u8>>,
                          deoxy: Vec<Option<u8>>,
                          acetyl: Vec<Option<u8>>,
                          glycolyl: Vec<Option<u8>>,
                          nglycolyl: Vec<Option<u8>>,
                          o_carboxy_ethyl: Vec<Option<u8>>| {
            [
                show_mods(&nacetyl, GlycanSubstituent::NAcetyl.notation()),
                show_mods(&acid, GlycanSubstituent::Acid.notation()),
                show_mods(&amino, GlycanSubstituent::Amino.notation()),
                show_mods(&deoxy, GlycanSubstituent::Deoxy.notation()),
                show_mods(&acetyl, GlycanSubstituent::Acetyl.notation()),
                show_mods(&glycolyl, GlycanSubstituent::Glycyl.notation()),
                show_mods(&nglycolyl, GlycanSubstituent::NGlycolyl.notation()),
                show_mods(
                    &o_carboxy_ethyl,
                    GlycanSubstituent::OCarboxyEthyl.notation(),
                ),
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
                if !o_carboxy_ethyl.is_empty() && !nacetyl.is_empty() {
                    (
                        Shape::Hexagon,
                        Colour::Purple,
                        inner_modifications,
                        outer_mods(
                            filter(nacetyl, &[2]),
                            acid,
                            amino,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            filter(o_carboxy_ethyl, &[3]),
                        ),
                    )
                } else if !o_carboxy_ethyl.is_empty() && !nglycolyl.is_empty() {
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
                            filter(nglycolyl, &[2]),
                            filter(o_carboxy_ethyl, &[3]),
                        ),
                    )
                } else if !o_carboxy_ethyl.is_empty() && !amino.is_empty() {
                    (
                        Shape::Hexagon,
                        Colour::Brown,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            filter(amino, &[2]),
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            filter(o_carboxy_ethyl, &[3]),
                        ),
                    )
                } else if deoxy.len() > 1 {
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
                            filter(deoxy, &[2, 6]),
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if amino.len() > 1 && !deoxy.is_empty() {
                    (
                        Shape::Hexagon,
                        Colour::Blue,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            filter(amino, &[2, 4]),
                            filter(deoxy, &[6]),
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if !nacetyl.is_empty() && !deoxy.is_empty() {
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
                            filter(nacetyl, &[2]),
                            acid,
                            amino,
                            filter(deoxy, &[6]),
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if !deoxy.is_empty() {
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
                            filter(deoxy, &[6]),
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if !acid.is_empty() || !amino.is_empty() || !nacetyl.is_empty() {
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
                    let (shape, nacetyl, acid, amino) = if !acid.is_empty() {
                        (Shape::DividedDiamond, nacetyl, filter(acid, &[6]), amino)
                    } else if !amino.is_empty() {
                        (Shape::CrossedSquare, nacetyl, acid, filter(amino, &[2]))
                    } else {
                        (Shape::Square, filter(nacetyl, &[2]), acid, amino)
                    };
                    (
                        shape,
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
            BaseSugar::Heptose(None) if acid.len() > 1 && !deoxy.is_empty() => (
                Shape::Hexagon,
                Colour::Orange,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    filter(acid, &[1, 6]),
                    amino,
                    filter(deoxy, &[3]),
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Octose if !acid.is_empty() && !deoxy.is_empty() => (
                Shape::Hexagon,
                Colour::Yellow,
                inner_modifications,
                outer_mods(
                    nacetyl,
                    filter(acid, &[2]),
                    amino,
                    filter(deoxy, &[3]),
                    acetyl,
                    glycolyl,
                    nglycolyl,
                    o_carboxy_ethyl,
                ),
            ),
            BaseSugar::Nonose(isomer) if !acid.is_empty() && !amino.is_empty() => {
                if amino.len() > 1 && deoxy.len() > 1 {
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
                            filter(acid, &[2]),
                            filter(amino, &[5, 7]),
                            filter(deoxy, &[3, 9]),
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else {
                    let colour = if !deoxy.is_empty() {
                        deoxy = filter(deoxy, &[3]);
                        if *isomer == Some(NonoseIsomer::Kdn) {
                            Colour::Green
                        } else {
                            Colour::Red
                        }
                    } else if !acetyl.is_empty() {
                        acetyl = filter(acetyl, &[5]);
                        Colour::Purple
                    } else if !glycolyl.is_empty() {
                        glycolyl = filter(glycolyl, &[5]);
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
                            filter(acid, &[2]),
                            filter(amino, &[5]),
                            deoxy,
                            acetyl,
                            glycolyl,
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
