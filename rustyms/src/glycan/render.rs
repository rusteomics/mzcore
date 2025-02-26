use std::{f32::consts::PI, fmt::Write};

use itertools::Itertools;

use crate::{
    fragment::GlycanPosition,
    glycan::{
        BaseSugar, Configuration, GlycanStructure, GlycanSubstituent, HeptoseIsomer, HexoseIsomer,
        MonoSaccharide, NonoseIsomer, PentoseIsomer,
    },
};

impl MonoSaccharide {
    /// Get the shape, colour, inner modifications, and outer modifications for this monosaccharide.
    fn get_shape(&self) -> (Shape, Colour, String, String) {
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
            if *c == Configuration::D {
                inner_modifications.push('D');
            } else {
                inner_modifications.push('L');
            }
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
                } else if acid > 0 {
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
                    (
                        Shape::DividedDiamond,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid - 1,
                            amino,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if amino > 0 {
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
                    (
                        Shape::CrossedSquare,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl,
                            acid,
                            amino - 1,
                            deoxy,
                            acetyl,
                            glycolyl,
                            nglycolyl,
                            o_carboxy_ethyl,
                        ),
                    )
                } else if nacetyl > 0 {
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
                    (
                        Shape::Square,
                        c,
                        inner_modifications,
                        outer_mods(
                            nacetyl - 1,
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
                            Some(NonoseIsomer::Leg) => {
                                if self.epi {
                                    Colour::LightBlue
                                } else {
                                    Colour::Yellow
                                }
                            }
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
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Colour {
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
    /// Represented as percentages 0..=100
    const fn cmyk(self) -> [u8; 4] {
        match self {
            Self::Background => [0, 0, 0, 0],
            Self::Blue => [100, 50, 0, 0],
            Self::Green => [100, 0, 100, 0],
            Self::Yellow => [0, 15, 100, 0],
            Self::Orange => [0, 65, 100, 0],
            Self::Pink => [0, 47, 24, 0],
            Self::Purple => [38, 88, 0, 0],
            Self::LightBlue => [41, 5, 3, 0],
            Self::Brown => [32, 48, 76, 13],
            Self::Red => [0, 100, 100, 0],
        }
    }

    /// Represented as bytes 0..=255
    const fn rgb(self) -> [u8; 3] {
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
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
enum Shape {
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
    const fn height(self) -> f32 {
        match self {
            Self::Rectangle | Self::FlatDiamond | Self::Hexagon => 0.5,
            _ => 1.0,
        }
    }
}

impl GlycanStructure {
    /// Build the rendered glycan, this can be used to create SVG images of this glycan. The footnotes list is used to gather modification texts that are too big to place in line. The caller will have to find their own way of displaying this to the user.
    pub fn render(&self, footnotes: &mut Vec<String>) -> RenderedGlycan {
        self.inner_render(0, &[], footnotes)
    }

    /// Build the rendered glycan.
    fn inner_render(
        &self,
        depth: usize,
        path: &[usize],
        footnotes: &mut Vec<String>,
    ) -> RenderedGlycan {
        let (shape, colour, inner_modifications, outer_modifications) = self.sugar.get_shape();
        // Automatically make footnotes out of long outer modification texts
        let outer_modifications = if outer_modifications.len() > 6 {
            let index = footnotes.iter().position(|e| *e == outer_modifications);
            index.map_or_else(
                || {
                    let index = footnotes.len();
                    footnotes.push(outer_modifications);
                    OuterModifications::Footnote(index)
                },
                OuterModifications::Footnote,
            )
        } else if !outer_modifications.is_empty() {
            OuterModifications::Text(outer_modifications)
        } else {
            OuterModifications::Empty
        };

        if self.branches.is_empty() {
            RenderedGlycan {
                y: 0,
                x: 0.0,
                mid_point: 0.5,
                width: 1.0,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                position: GlycanPosition {
                    inner_depth: depth,
                    series_number: depth,
                    branch: path.to_vec(),
                    attachment: None,
                },
                title: self.sugar.to_string(),
                branch_index: 0,
                branches: Vec::new(),
                sides: Vec::new(),
            }
        } else {
            let mut y_depth = 0;
            let mut branches = Vec::new();
            let mut sides = Vec::new();
            for (branch_index, branch) in self.branches.iter().enumerate() {
                let mut new_path = path.to_vec();
                new_path.push(branch_index);
                let mut rendered = branch.inner_render(
                    depth + 1,
                    if self.branches.len() > 1 {
                        &new_path
                    } else {
                        path
                    },
                    footnotes,
                );
                rendered.branch_index = branch_index;
                if rendered.is_sideways() && sides.len() < 2 {
                    if sides.is_empty() && rendered.shape == Shape::Triangle {
                        rendered.shape = Shape::LeftPointingTriangle;
                    } else if sides.len() == 1 && rendered.shape == Shape::Triangle {
                        rendered.shape = Shape::RightPointingTriangle;
                    }
                    sides.push(rendered);
                } else {
                    y_depth = y_depth.max(rendered.y);
                    branches.push(rendered);
                }
            }
            // Update all branch placements
            let mut displacement = 0.0;
            for branch in &mut branches {
                branch.transpose(y_depth - branch.y, displacement);
                displacement += branch.width;
            }
            if !branches.is_empty() {
                y_depth += 1;
            }
            // Determine the center point for this sugar
            let mut center = match branches.len() {
                0 => 0.5,
                1 => branches[0].mid_point,
                n => {
                    // Find the median midpoint of the branches
                    (branches[n / 2 - (n + 1) % 2].x
                        + branches[n / 2 - (n + 1) % 2].mid_point
                        + branches[n / 2].x
                        + branches[n / 2].mid_point)
                        / 2.0
                }
            };
            let mut width = branches.last().map_or(1.0, |b| b.x + b.width);
            if !sides.is_empty() {
                sides[0].transpose(y_depth, center + 0.5);
                width = width.max(center + 0.5 + sides[0].width);
            }
            if sides.len() == 2 {
                let mut x = center - 0.5 - sides[1].width;
                if x < 0.0 {
                    let shift = -x;
                    center += shift;
                    for branch in &mut branches {
                        branch.transpose(0, shift);
                    }
                    sides[0].transpose(0, shift);
                    width += shift;
                    x = 0.0;
                }
                sides[1].transpose(y_depth, x);
            }
            RenderedGlycan {
                y: y_depth,
                x: 0.0,
                mid_point: center,
                width,
                shape,
                colour,
                inner_modifications,
                outer_modifications,
                position: GlycanPosition {
                    inner_depth: depth,
                    series_number: depth,
                    branch: path.to_vec(),
                    attachment: None,
                },
                title: self.sugar.to_string(),
                branch_index: 0,
                branches,
                sides,
            }
        }
    }
}

/// A rendered glycan, contains all information needed to render this to svg or a bitmap.
#[derive(Debug, Clone)]
pub struct RenderedGlycan {
    /// The depth of this sugar along the main axis of the glycan, starting at 0 at the top (in the leaves)
    y: usize,
    /// The sideways placement of this whole tree starting at 0 at the leftmost monosaccharide, 1.0 is the width of one monosaccharide
    x: f32,
    /// The sideways placement of this sugar within this tree, for the absolute sideways placement of this sugar add this to `x`
    mid_point: f32,
    /// The total width of the (sub)tree with all of its branches and sides
    width: f32,
    /// The shape of the monosaccharide
    shape: Shape,
    /// The colour of the monosaccharide
    colour: Colour,
    /// Text to be shown inside the monosaccharide
    inner_modifications: String,
    /// Text to be shown outside the monosaccharide
    outer_modifications: OuterModifications,
    /// The position of this sugar
    position: GlycanPosition,
    /// Full name of the glycan
    title: String,
    /// The index into the branches of the parent monosaccharide
    branch_index: usize,
    /// All branches that go up the tree
    branches: Vec<RenderedGlycan>,
    /// All branches that go to the side (Fucoses)
    sides: Vec<RenderedGlycan>,
}

#[derive(Debug, Clone)]
/// Modifications that are to be shown outside of the monosaccharide
enum OuterModifications {
    /// Too long of a text, or it did not fit, so show as a footnote
    Footnote(usize),
    /// Text
    Text(String),
    /// No modification
    Empty,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
/// THe direction of the rendered glycan
pub enum GlycanDirection {
    /// A top down tree, with the root at the bottom
    TopDown,
    /// A left to right tree, with the root at the right
    LeftToRight,
}

/// A subtree of a rendered glycan, used to restrict the canvas for glycan fragments
struct SubTree<'a> {
    /// The root for this sub tree
    tree: &'a RenderedGlycan,
    /// Total depth of the glycans with the breaks applied
    depth: usize,
    /// The horizontal offset from the left
    left_offset: f32,
    /// The horizontal offset from the right
    right_offset: f32,
    /// If this fragment is topped by a breaking symbol, needed to calculate the correct height for the canvas
    break_top: bool,
    /// All breaking branches, standardised to the linked root
    branch_breaks: Vec<(usize, &'a [usize])>,
}

impl RenderedGlycan {
    /// Transpose this glycan and all of its branches
    fn transpose(&mut self, y: usize, x: f32) {
        self.y += y;
        self.x += x;
        for branch in &mut self.branches {
            branch.transpose(y, x);
        }
        for side in &mut self.sides {
            side.transpose(y, x);
        }
    }

    /// Check if this sugar should be rendered to the side of the parent sugar
    fn is_sideways(&self) -> bool {
        self.colour == Colour::Red
            && self.shape == Shape::Triangle
            && self.branches.is_empty()
            && self.sides.is_empty()
    }

    /// Get the subtree starting on the given position, return None if the starting position is not valid, it also indicates the depth of this subtree for the given branch breakages and if a break tops the structure
    fn get_subtree<'a>(
        &'a self,
        start: &'a GlycanPosition,
        branch_breaks: &'a [GlycanPosition],
    ) -> Option<SubTree<'a>> {
        /// Calculate the maximal depth, break top, left and right offset
        fn canvas_size(
            tree: &RenderedGlycan,
            breakages: &[(usize, &[usize])],
        ) -> (usize, bool, f32, f32) {
            let lx = (tree.x + tree.mid_point - 0.5).max(0.0);
            let rx = (tree.width - tree.mid_point - 0.5).max(0.0);
            // The tree is cut here
            if breakages.iter().any(|b| b.0 == 0) {
                return (0, true, lx, rx);
            };

            let total_branches = tree.branches.len() + tree.sides.len();
            let (depth, break_top, left_offset, right_offset) = match total_branches {
                0 => (0, false, lx, rx),
                1 => tree.branches.first().map_or((0, false, lx, rx), |branch| {
                    canvas_size(
                        branch,
                        &breakages.iter().map(|b| (b.0 - 1, b.1)).collect_vec(),
                    )
                }),
                _ => tree
                    .branches
                    .iter()
                    .enumerate()
                    .map(|(i, branch)| {
                        (
                            i,
                            canvas_size(
                                branch,
                                &breakages
                                    .iter()
                                    .filter(|b| b.1.first() == Some(&branch.branch_index))
                                    .map(|b| (b.0 - 1, &b.1[1..]))
                                    .collect_vec(),
                            ),
                        )
                    })
                    .fold((0, false, lx, rx), |acc, (i, v)| {
                        (
                            acc.0.max(v.0),
                            if v.0 >= acc.0 { v.1 } else { acc.1 },
                            if i == 0 { v.2 } else { acc.2 },
                            if i == tree.branches.len() - 1 {
                                v.3
                            } else {
                                acc.3
                            },
                        )
                    }),
            };
            (
                depth + 1,
                break_top,
                if tree.sides.len() == 2 {
                    left_offset.min(tree.x + tree.mid_point - 1.5).max(0.0)
                } else {
                    left_offset
                },
                if tree.sides.is_empty() {
                    right_offset
                } else {
                    right_offset.min(tree.width - tree.mid_point - 1.5).max(0.0)
                },
            )
        }

        let mut tree = self;
        let mut depth = 0;
        let mut branch_choices = start.branch.clone();
        while depth < start.inner_depth {
            depth += 1;

            let total_branches = tree.branches.len() + tree.sides.len();
            match total_branches {
                0 => return None,
                1 => tree = tree.branches.first().or_else(|| tree.sides.first())?,
                _ => {
                    let index = branch_choices.pop()?;
                    tree = tree
                        .branches
                        .iter()
                        .find(|b| b.branch_index == index)
                        .or_else(|| tree.sides.iter().find(|b| b.branch_index == index))?;
                }
            }
        }

        let rules = branch_breaks
            .iter()
            .filter(|b| b.inner_depth >= start.inner_depth && b.branch.starts_with(&start.branch))
            .map(|b| {
                (
                    b.inner_depth - start.inner_depth,
                    &b.branch[start.branch.len()..],
                )
            })
            .collect_vec();
        let (depth, break_top, left_offset, right_offset) = canvas_size(tree, &rules);
        Some(SubTree {
            tree,
            depth,
            left_offset,
            right_offset,
            break_top,
            branch_breaks: rules,
        })
    }

    /// Render this glycan as SVG. The SVG will be appended to the given buffer.
    ///  * `output`: the buffer to append the SVG to.
    ///  * `basis`: the text to draw at the root of the tree.
    ///  * `column_size`: the size (in pixels) of one block in the glycan, the full size with the padding and sugar size included.
    ///  * `sugar_size`: the size (in pixels) of a monosaccharide.
    ///  * `stroke_size`: the size (in pixels) of the strokes in the graphic.
    ///  * `direction`: the direction the draw the image in.
    ///  * `root_break`: the first monosaccharide to be included in the rendering, used to draw glycan fragments.
    ///  * `branch_breaks`: the list of breaks in the branches of the fragment, the fragment will include the indicated glycan position.
    ///  * `foreground`: the colour to be used for the foreground, in RGB order.
    ///  * `background`: the colour to be used for the background, in RGB order, this is used to fill 'empty' sugars if the isomeric state is unknown.
    ///  * `footnotes`: used to gather modification texts that are too big to place in line. The caller will have to find their own way of displaying this to the user.
    ///
    /// # Errors
    /// If the underlying buffer errors the error is returned. Otherwise `Ok(false)` is returned if the given `root_break` is not valid, and `Ok(true)` is returned if the rendering was fully successful.
    pub fn to_svg(
        &self,
        mut output: impl Write,
        basis: Option<String>,
        column_size: f32,
        sugar_size: f32,
        stroke_size: f32,
        direction: GlycanDirection,
        root_break: Option<GlycanPosition>,
        branch_breaks: &[GlycanPosition],
        foreground: [u8; 3],
        background: [u8; 3],
        footnotes: &mut Vec<String>,
    ) -> Result<bool, std::fmt::Error> {
        fn render_element(
            buffer: &mut impl Write,
            element: &RenderedGlycan,
            column_size: f32,
            sugar_size: f32,
            stroke_size: f32,
            direction: GlycanDirection,
            x_offset: f32,
            y_offset: f32,
            breaks: &[(usize, &[usize])],
            foreground: &str,
            background: &str,
            incoming_stroke: (f32, f32, f32, f32),
            footnotes: &mut Vec<String>,
        ) -> Result<(), std::fmt::Error> {
            let x = pick_direction("x", "y", direction);
            let y = pick_direction("y", "x", direction);
            let raw_x = element.x - x_offset;
            let raw_y = element.y as f32 - y_offset;

            let total_branches = element.branches.len() + element.sides.len();
            let mut strokes = vec![incoming_stroke];
            // First all lines to get good stacking behaviour
            for (side, branch) in element
                .branches
                .iter()
                .map(|b| (false, b))
                .chain(element.sides.iter().map(|b| (true, b)))
            {
                let origin_x = (raw_x + element.mid_point) * column_size;
                let origin_y = (raw_y + 0.5) * column_size;
                let base_x = (branch.x + branch.mid_point - x_offset) * column_size;
                if (total_branches == 1 && breaks.iter().any(|b| b.0 == 1))
                    || breaks
                        .iter()
                        .any(|b| b.0 == 1 && b.1.first() == Some(&branch.branch_index))
                {
                    let base_y =
                        (raw_y - 0.5 + f32::from(side)).mul_add(column_size, stroke_size * 0.5);
                    let angle = f32::atan2(base_y - origin_y, base_x - origin_x);
                    write!(
                        buffer,
                        "<line {x}1=\"{origin_x}\" {y}1=\"{origin_y}\" {x}2=\"{base_x}\" {y}2=\"{base_y}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                    )?;
                    let x1 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).cos(), base_x);
                    let y1 = (sugar_size / 2.0).mul_add(-0.5f32.mul_add(PI, -angle).sin(), base_y);
                    let x2 = (sugar_size / 2.0).mul_add(-0.5f32.mul_add(PI, -angle).cos(), base_x);
                    let y2 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).sin(), base_y);
                    write!(
                        buffer,
                        "<line {x}1=\"{x1}\" {y}1=\"{y1}\" {x}2=\"{x2}\" {y}2=\"{y2}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                    )?;
                    let x3 = (stroke_size * 2.0).mul_add(-angle.cos(), x1);
                    let y3 = (stroke_size * 2.0).mul_add(-angle.sin(), y1);
                    write!(
                        buffer,
                        "<line {x}1=\"{x1}\" {y}1=\"{y1}\" {x}2=\"{x3}\" {y}2=\"{y3}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                    )?;
                    let offset = 0.25f32.mul_add(column_size, stroke_size);
                    let r = 0.25f32
                        .mul_add(column_size, -stroke_size)
                        .min(sugar_size * 0.25);
                    let adjusted_x = offset.mul_add(-angle.cos(), base_x);
                    let adjusted_y = offset.mul_add(-angle.sin(), base_y);
                    write!(
                        buffer,
                        "<circle r=\"{r}\" c{x}=\"{adjusted_x}\" c{y}=\"{adjusted_y}\" fill=\"transparent\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                    )?;
                    strokes.push(pick_box(
                        (
                            origin_x.min(base_x),
                            base_y.min(origin_y),
                            origin_x.max(base_x),
                            base_y.max(origin_y),
                        ),
                        direction,
                    ));
                } else {
                    let base_y = ((branch.y as f32) - y_offset + 0.5) * column_size;
                    write!(
                        buffer,
                        "<line {x}1=\"{origin_x}\" {y}1=\"{origin_y}\" {x}2=\"{base_x}\" {y}2=\"{base_y}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                    )?;
                    strokes.push(pick_box(
                        (
                            origin_x.min(base_x),
                            base_y.min(origin_y),
                            origin_x.max(base_x),
                            base_y.max(origin_y),
                        ),
                        direction,
                    ));
                }
            }

            // Render the sugar
            let fill = if element.colour == Colour::Background {
                background.to_string()
            } else {
                let colour = element.colour.rgb();
                format!("rgb({},{},{})", colour[0], colour[1], colour[2])
            };
            let title = format!(
                " data-sugar=\"{}\" data-position=\"{}-{}\"",
                element.title,
                element.position.inner_depth,
                element.position.branch.iter().join(",")
            );
            match element.shape {
                Shape::Circle => write!(
                    buffer,
                    "<circle r=\"{}\" c{x}=\"{}\" c{y}=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                    sugar_size / 2.0,
                    (raw_x + element.mid_point) * column_size,
                    (raw_y + 0.5) * column_size
                )?,
                Shape::Square => write!(
                    buffer,
                    "<rect {x}=\"{}\" {y}=\"{}\" width=\"{}\" height=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                    (raw_x + element.mid_point).mul_add(column_size, - sugar_size / 2.0),
                    (raw_y + 0.5).mul_add(column_size, - sugar_size / 2.0),
                    sugar_size,
                    sugar_size
                )?,
                Shape::Rectangle => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    write!(
                        buffer,
                        "<rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                        base_x - sugar_size / 2.0,
                        base_y - sugar_size / 4.0,
                        sugar_size,
                        sugar_size / 2.0
                    )?;},
                Shape::Triangle => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y2}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::LeftPointingTriangle => {
                    let x1 = (raw_x + element.mid_point).mul_add(column_size,- sugar_size / 2.0);
                    let x2 = x1 + sugar_size;
                    let y1 = (raw_y + 0.5).mul_add(column_size, - sugar_size / 2.0);
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                    points(&[(x1, y2), (x2, y1), (x2, y3)], direction)
                )?;},
                Shape::RightPointingTriangle => {
                    let x1 = (raw_x + element.mid_point).mul_add(column_size,- sugar_size / 2.0);
                    let x2 = x1 + sugar_size;
                    let y1 = (raw_y + 0.5).mul_add(column_size, - sugar_size / 2.0);
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                    points(&[(x1, y1), (x2, y2), (x1, y3)], direction)
                )?;},
                Shape::Diamond => {
                    let x1 = (raw_x + element.mid_point).mul_add(column_size,- sugar_size / 2.0);
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = (raw_y + 0.5).mul_add(column_size, - sugar_size / 2.0);
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                    points(&[(x2,y1),(x3,y2),(x2,y3),(x1,y2)], direction)
                )?;},
                Shape::FlatDiamond => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 4.0;
                    let y2 = y1 + sugar_size / 4.0;
                    let y3 = y1 + sugar_size / 2.0;

                    write!(
                    buffer,
                    "<polygon points=\"{x2} {y1} {x3} {y2} {x2} {y3} {x1} {y2}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::Hexagon => {
                    let a = sugar_size / 2.0 / 3.0_f32.sqrt();
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + a;
                    let x3 = x1 + sugar_size - a;
                    let x4 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 4.0;
                    let y2 = y1 + sugar_size / 4.0;
                    let y3 = y1 + sugar_size / 2.0;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y1} {x4} {y2} {x3} {y3} {x2} {y3}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::Pentagon => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let a = (18.0/360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let b = (18.0/360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let c = (36.0/360.0*2.0*PI).cos() * sugar_size / 2.0;
                    let d = (36.0/360.0*2.0*PI).sin() * sugar_size / 2.0;
                    let x1 = base_x - a;
                    let x2 = base_x - d;
                    let x3 = base_x;
                    let x4 = base_x + d;
                    let x5 = base_x + a;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0 - b;
                    let y3 = y1 + sugar_size / 2.0 + c;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x3} {y1} {x5} {y2} {x4} {y3} {x2} {y3}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::Star => {
                    // The Phi constant, the ratio for the "golden ratio"
                    const PHI: f32 = 1.618_034_f32;
                    // Calculate sizes of parts of the pentagram
                    let a = (18.0/360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let b = (18.0/360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let c = (36.0/360.0*2.0*PI).cos() * sugar_size / 2.0;
                    let d = (36.0/360.0*2.0*PI).sin() * sugar_size / 2.0;
                    let e = 2.0 * a / (54.0/360.0*2.0*PI).sin() / (1.0 + 1.0 / PHI);
                    let f = (18.0/360.0 * 2.0 * PI).cos().mul_add(e, -(sugar_size / 2.0));
                    let g = (18.0/360.0 * 2.0 * PI).sin() * e;
                    let h = (sugar_size / 2.0 - b) * (18.0/360.0*2.0*PI).tan();
                    let j = (18.0/360.0*2.0*PI).tan() * g;
                    // Calculate the positions of the pentagram points
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - a;
                    let x2 = base_x - d;
                    let x3 = base_x - g;
                    let x4 = base_x - h;
                    let x5 = base_x;
                    let x6 = base_x + h;
                    let x7 = base_x + g;
                    let x8 = base_x + d;
                    let x9 = base_x + a;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0 - b;
                    let y3 = y1 + sugar_size / 2.0 + j;
                    let y4 = y1 + sugar_size / 2.0 + f;
                    let y5 = y1 + sugar_size / 2.0 + c;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x4} {y2} {x5} {y1} {x6} {y2} {x9} {y2} {x7} {y3} {x8} {y5} {x5} {y4} {x2} {y5} {x3} {y3}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::CrossedSquare => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let y1 = base_y - sugar_size / 2.0;
                    let x2 = x1 + sugar_size;
                    let y2 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y1} {x2} {y1} {x2} {y2}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y1} {x1} {y2} {x2} {y2}\" fill=\"{background}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y1} {x2} {y1} {x2} {y2} {x1} {y2}\" fill=\"transparent\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::DividedDiamond => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y- sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{x1} {y2} {x2} {y1} {x3} {y2}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y2} {x2} {y3} {x3} {y2}\" fill=\"{background}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x1} {y2} {x2} {y1} {x3} {y2} {x2} {y3}\" fill=\"transparent\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
                Shape::DividedTriangle => {
                    let (base_x, base_y) = pick_point((((raw_x + element.mid_point) * column_size), (raw_y + 0.5) * column_size), direction);
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size;

                    write!(
                    buffer,
                    "<polygon points=\"{x2} {y1} {x3} {y2} {x2} {y2}\" fill=\"{fill}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x2} {y1} {x1} {y2} {x2} {y2}\" fill=\"{background}\"  stroke=\"{foreground}\" stroke-width=\"{stroke_size}\" stroke-linejoin=\"bevel\"/><polygon points=\"{x2} {y1} {x3} {y2} {x1} {y2}\" fill=\"transparent\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"{title}/>",
                )?;},
            }
            if !element.inner_modifications.is_empty() {
                write!(buffer, "<text {x}=\"{}\" {y}=\"{}\" fill=\"{foreground}\" text-anchor=\"middle\" font-style=\"italic\" font-size=\"{}px\" dominant-baseline=\"middle\">{}</text>",
                    (raw_x + element.mid_point) * column_size,
                    (raw_y + 0.5) * column_size,
                    sugar_size / 2.0,
                    element.inner_modifications)?;
            }
            if let Some((pos_x, pos_y, alignment, text)) = text_location(
                &element.outer_modifications,
                element.shape,
                pick_point(
                    (
                        (raw_x + element.mid_point - 0.5) * column_size,
                        raw_y * column_size,
                    ),
                    direction,
                ),
                column_size,
                sugar_size,
                stroke_size,
                &strokes,
                footnotes,
            ) {
                write!(buffer, "<text x=\"{pos_x}\" y=\"{pos_y}\" fill=\"{foreground}\" text-anchor=\"{alignment}\" font-size=\"{}px\" dominant-baseline=\"hanging\">{text}</text>", sugar_size / 2.0)?;
            }
            // Render all connected sugars
            for (index, branch) in element
                .branches
                .iter()
                .chain(element.sides.iter())
                .enumerate()
            {
                if !((total_branches == 1 && breaks.iter().any(|b| b.0 == 1))
                    || breaks
                        .iter()
                        .any(|b| b.0 == 1 && b.1.first() == Some(&branch.branch_index)))
                {
                    render_element(
                        buffer,
                        branch,
                        column_size,
                        sugar_size,
                        stroke_size,
                        direction,
                        x_offset,
                        y_offset,
                        &breaks
                            .iter()
                            .filter(|b| {
                                total_branches > 1 && b.1.first() == Some(&branch.branch_index)
                                    || total_branches == 1
                            })
                            .map(|b| (b.0 - 1, &b.1[usize::from(total_branches > 1)..]))
                            .collect_vec(),
                        foreground,
                        background,
                        strokes[index + 1],
                        footnotes,
                    )?;
                }
            }
            Ok(())
        }

        let base = GlycanPosition {
            inner_depth: 0,
            series_number: 0,
            branch: Vec::new(),
            attachment: None,
        };
        let Some(sub_tree) = self.get_subtree(root_break.as_ref().unwrap_or(&base), branch_breaks)
        else {
            return Ok(false);
        };

        let depth = sub_tree.depth as f32 + f32::from(sub_tree.break_top);
        let width =
            (sub_tree.tree.x + sub_tree.tree.width - sub_tree.left_offset - sub_tree.right_offset)
                * column_size;
        let height = (depth + f32::from(basis.is_some() && root_break.is_none())).mul_add(
            column_size,
            if root_break.is_some() {
                3.5 * stroke_size
            } else {
                0.0
            },
        );

        let x = pick_direction("x", "y", direction);
        let y = pick_direction("y", "x", direction);
        let w = pick_direction("width", "height", direction);
        let h = pick_direction("height", "width", direction);

        write!(output, "<svg {w}=\"{}\" {h}=\"{}\">", width, height)?;

        let foreground = format!("rgb({},{},{})", foreground[0], foreground[1], foreground[2]);
        let background = format!("rgb({},{},{})", background[0], background[1], background[2]);
        let stroke = if root_break.is_some() {
            let base_x =
                (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset) * column_size;
            let base_y = depth.mul_add(column_size, stroke_size * 3.0);
            write!(
                output,
                "<line {x}1=\"{base_x}\" {y}1=\"{}\" {x}2=\"{base_x}\" {y}2=\"{}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                (depth - 0.5) * column_size,
                base_y,
            )?;
            write!(
                output,
                "<line {x}1=\"{}\" {y}1=\"{base_y}\" {x}2=\"{}\" {y}2=\"{base_y}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                base_x - sugar_size / 2.0,
                base_x + sugar_size / 2.0,
                
            )?;
            write!(
                output,
                "<line {x}1=\"{bs_x}\" {y}1=\"{}\" {x}2=\"{bs_x}\" {y}2=\"{}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                stroke_size.mul_add(-2.0, base_y),
                base_y,
                bs_x=base_x - sugar_size / 2.0,
            )
            ?;
            (base_x, (depth - 0.5) * column_size, base_x, base_y)
        } else if let Some(basis) = basis {
            let base_x =
                (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset) * column_size;
            let base_y = depth.mul_add(column_size, (column_size - sugar_size) / 2.0);
            write!(
                output,
                "<line {x}1=\"{base_x}\" {y}1=\"{}\" {x}2=\"{base_x}\" {y}2=\"{base_y}\" stroke=\"{foreground}\" stroke-width=\"{stroke_size}\"/>",
                (depth - 0.5) * column_size,
            )
            ?;
            if direction == GlycanDirection::TopDown {
                write!(output, "<text x=\"{base_x}\" y=\"{}\" fill=\"{foreground}\" text-anchor=\"middle\" font-size=\"{}px\" dominant-baseline=\"ideographic\">{basis}</text>",
                        base_y  + sugar_size,
                        sugar_size)?;
            } else {
                write!(output, "<text y=\"{base_x}\" x=\"{}\" fill=\"{foreground}\" text-anchor=\"end\" font-size=\"{}px\" dominant-baseline=\"middle\">{basis}</text>",
                (depth + 1.0) * column_size,
                sugar_size)?;
            }
            (base_x, (depth - 0.5) * column_size, base_x, base_y)
        } else {
            (
                sub_tree.tree.x + sub_tree.tree.mid_point,
                (depth + 0.5) * column_size,
                sub_tree.tree.x + sub_tree.tree.mid_point,
                (depth + 0.5) * column_size,
            )
        };
        render_element(
            &mut output,
            sub_tree.tree,
            column_size,
            sugar_size,
            stroke_size,
            direction,
            sub_tree.left_offset,
            sub_tree.tree.y as f32 - (depth - 1.0),
            &sub_tree.branch_breaks,
            &foreground,
            &background,
            pick_box(stroke, direction),
            footnotes,
        )?;

        write!(output, "</svg>")?;

        Ok(true)
    }
}

/// Determine the best location for text, returns the x, y, text-anchor, and the contents
fn text_location(
    outer_modifications: &OuterModifications,
    shape: Shape,
    position: (f32, f32),
    column_size: f32,
    sugar_size: f32,
    stroke_size: f32,
    strokes: &[(f32, f32, f32, f32)], // x1, y1, x2, y2
    footnotes: &mut Vec<String>,
) -> Option<(f32, f32, &'static str, String)> {
    let text = match outer_modifications {
        OuterModifications::Empty => return None,
        OuterModifications::Footnote(index) => (index + 1).to_string(), // Human numbering
        OuterModifications::Text(text) => text.clone(),
    };

    // Stay within box, dodge strokes, if not fitting fall back to adding to footnotes
    let text_height = sugar_size / 2.0;
    let text_width = (text.len() as f32) * text_height; // Rule of thumb, on average text is thinner than square, so this should be a good upper limit
    let vertical_padding = shape.height().mul_add(-sugar_size, column_size) / 2.0;

    if vertical_padding >= text_height {
        let mut options = vec![
            (
                (
                    (column_size - text_width).mul_add(0.5, position.0),
                    position.1,
                    (column_size + text_width).mul_add(0.5, position.0),
                    position.1 + text_height,
                ),
                "middle",
            ),
            (
                (
                    (column_size - text_width).mul_add(0.5, position.0),
                    position.1 + column_size - text_height,
                    (column_size + text_width).mul_add(0.5, position.0),
                    position.1 + column_size,
                ),
                "middle",
            ),
            (
                (
                    position.0 + column_size.mul_add(0.5, -stroke_size) - text_width,
                    position.1,
                    position.0 + column_size.mul_add(0.5, -stroke_size),
                    position.1 + text_height,
                ),
                "end",
            ),
            (
                (
                    stroke_size.mul_add(2.0, column_size.mul_add(0.5, position.0)),
                    position.1 + column_size - text_height,
                    stroke_size.mul_add(2.0, column_size.mul_add(0.5, position.0) + text_width),
                    position.1 + column_size,
                ),
                "start",
            ),
            (
                (
                    position.0,
                    position.1,
                    position.0 + text_width,
                    position.1 + text_height,
                ),
                "start",
            ),
            (
                (
                    position.0 + column_size - text_width,
                    position.1,
                    position.0 + column_size,
                    position.1 + text_height,
                ),
                "end",
            ),
            (
                (
                    position.0,
                    position.1 + column_size - text_height,
                    position.0 + text_width,
                    position.1 + column_size,
                ),
                "start",
            ),
            (
                (
                    position.0 + column_size - text_width,
                    position.1 + column_size - text_height,
                    position.0 + column_size,
                    position.1 + column_size,
                ),
                "end",
            ),
        ];

        // Remove any options that put text outside of the box
        options.retain(|(option, _)| {
            option.0 >= position.0
                && option.2 <= position.0 + column_size
                && option.1 >= position.1
                && option.3 <= position.1 + column_size
        });

        for stroke in strokes {
            options.retain(|(option, _)| !hitbox_test(*option, *stroke));
        }
        if let Some((option, anchor)) = options.first() {
            let (x, y) = match *anchor {
                "end" => (option.2, option.1),
                "middle" => ((option.0 + option.2) / 2.0, option.1),
                _ => (option.0, option.1),
            };
            return Some((x, y, anchor, text));
        }
    }

    if let OuterModifications::Text(text) = outer_modifications {
        let index = footnotes
            .iter()
            .position(|p| *p == *text)
            .unwrap_or_else(|| {
                footnotes.push(text.clone());
                footnotes.len() - 1
            });

        text_location(
            &OuterModifications::Footnote(index),
            shape,
            position,
            column_size,
            sugar_size,
            stroke_size,
            strokes,
            footnotes,
        )
    } else {
        Some((
            position.0 + column_size,
            position.1 + column_size - text_height,
            "end",
            text,
        ))
    }
}

/// Test if two boxes hit
fn hitbox_test(box1: (f32, f32, f32, f32), box2: (f32, f32, f32, f32)) -> bool {
    debug_assert!(box1.0 <= box1.2, "Invalid boxes: {box1:?} {box2:?}");
    debug_assert!(box2.0 <= box2.2, "Invalid boxes: {box1:?} {box2:?}");
    debug_assert!(box1.1 <= box1.3, "Invalid boxes: {box1:?} {box2:?}");
    debug_assert!(box2.1 <= box2.3, "Invalid boxes: {box1:?} {box2:?}");
    box1.2 > box2.0 && box1.0 < box2.2 && box1.3 > box2.1 && box1.1 < box2.3
}

fn pick_direction<T>(a: T, b: T, direction: GlycanDirection) -> T {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        b
    }
}

fn pick_point<T>(a: (T, T), direction: GlycanDirection) -> (T, T) {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        (a.1, a.0)
    }
}

fn pick_box<T>(a: (T, T, T, T), direction: GlycanDirection) -> (T, T, T, T) {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        (a.1, a.0, a.3, a.2)
    }
}

fn points(points: &[(f32, f32)], direction: GlycanDirection) -> String {
    points
        .iter()
        .map(|p| {
            let (a, b) = pick_point(*p, direction);
            format!("{a} {b}",)
        })
        .join(" ")
}

#[test]
fn test_rendering() {
    const COLUMN_SIZE: f32 = 30.0;
    const SUGAR_SIZE: f32 = 15.0;
    const STROKE_SIZE: f32 = 1.5;
    let mut html = String::new();
    let mut footnotes = Vec::new();
    write!(&mut html, "<html lang=\"en\"><head><meta charset=\"UTF-8\"><meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\"><title>Glycan render test</title></head><body>").unwrap();

    let codes = [
        ("G01670UQ", "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"),
        ("G13523IF", "Fuc(?1-?)Gal(?1-?)GalNAc(?1-"),
        ("G00613DO", "GlcN(b1-4)GlcNAc(b1-4)GlcNAc(b1-4)GlcNAc6S(?1-"),
        ("G00621IU", "Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Gal(a1-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-8)Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-"),
        ("G01464QV", "Rha2,3,4Ac3(a1-2)[Xyl(b1-3)]Ara(a1-"),
        ("G04421VO", "Fruf(b2-1a)[Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)Glc(a1-2)]Glc"),
        ("G04458LN", "Kdn(a2-3)Gal(b1-4)ManNAc(b1-2)[Kdn(a2-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Kdn(a2-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Gc(a2-3)Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-"),
        ("G69524KC", "Xyl(?1-?)Ara(?1-?)[Gal(?1-?)]GlcA"),
        ("G37707YH", "Fuc(a1-2)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-2)[Gal(a1-3)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[GlcNAc(b1-4)][Neu5Gc(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-3/6)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(?1-"),
        ("G07370RP", "Rha(a1-3)Qui(b1-4)Rha(a1-2)Glc(b1-2)[Rha(a1-6)]Glc(b1-"),
        ("G11504PZ", "Dig3CMe(b1-3)Oli(b1-3)Oli(b1-"),
        ("G64699IM", "GlcA(b1-3)GalNAc(b1-4)4eLeg?5,7Ac2(a2-"),
        ("G14402AU", "D-Araf(b1-5)Dha(?2-3)[GalA(a1-4)GalA(a1-4)]GalA(a1-4)GalA"),
        ("G08395BZ", "Glc(b1-2a)[Ido(b1-3)]Psif"),
        ("G49642ZT", "Man(?1-?)[Man(?1-?)]Man(?1-?)[Man(?1-?)]Man(?1-?)GlcNAc(?1-?)[Fuc(?1-?)][Fuc(?1-?)]GlcNAc(?1-"),
        ("G59426OB", "Hex(?1-?)HexNAc(?1-?)HexA(?1-?)Gal(?1-?)GalNAc-ol"),
        ("G75424NV", "Hex?(?1-?)Hex?NAc(?1-?)[Hex?NAc(?1-?)]Hex?(?1-?)[Hex?(?1-?)[Hex?(?1-?)]Hex?(?1-?)][Hex?NAc(?1-?)]Hex?(?1-?)Hex?NAc(?1-?)Hex?NAc(?1-"),
        ("G36128WO", "Ido(b1-3)ManNAc(?1-3)[Ido(b1-3)L-AllNAc(b1-3)Ido(b1-4)AltNAc(b1-6)]Tal(b1-4)D-Ido(?1-"),
        ("G83422GV", "L-6dTal(a1-3)[Fuc(a1-2)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)]GlcNAc(b1-3)Gal(b1-3)[Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)]GalNAc(a1-"),
        ("G09073GJ","GalNAc(?1-?)GlcA2,3NAc2(?1-?)D-FucNAc"),
        ("G00069DT","Neu(a2-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)GlcNAc(b1-3)Gal(b1-4)Glc(b1-"),
    ];

    for (_, iupac) in &codes {
        let structure = GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0).unwrap();
        if !structure
            .render(&mut footnotes)
            .to_svg(
                &mut html,
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                GlycanDirection::TopDown,
                None,
                &[],
                [66, 66, 66],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap()
        {
            panic!("Rendering failed")
        }
    }

    for (_, iupac) in &codes {
        let structure = GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0).unwrap();
        if !structure
            .render(&mut footnotes)
            .to_svg(
                &mut html,
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                GlycanDirection::LeftToRight,
                None,
                &[],
                [66, 66, 66],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap()
        {
            panic!("Rendering failed")
        }
    }

    for (index, root, breaks) in [
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            Vec::new(),
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 4,
                series_number: 4,
                branch: vec![1],
                attachment: None,
            }],
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 5,
                    series_number: 5,
                    branch: vec![0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![1],
                    attachment: None,
                },
            ],
        ),
        (
            0,
            Some(GlycanPosition {
                inner_depth: 4,
                series_number: 4,
                branch: vec![1],
                attachment: None,
            }),
            Vec::new(),
        ),
        (
            14,
            Some(GlycanPosition {
                inner_depth: 0,
                series_number: 0,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }],
        ),
        (
            14,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }),
            vec![GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: vec![0],
                attachment: None,
            }],
        ),
        (
            16,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: Vec::new(),
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 4,
                    series_number: 4,
                    branch: vec![1, 1],
                    attachment: None,
                },
            ],
        ),
        (
            18,
            Some(GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: vec![0],
                attachment: None,
            }),
            vec![
                GlycanPosition {
                    inner_depth: 3,
                    series_number: 3,
                    branch: vec![0, 0],
                    attachment: None,
                },
                GlycanPosition {
                    inner_depth: 6,
                    series_number: 6,
                    branch: vec![0, 1],
                    attachment: None,
                },
            ],
        ),
        (
            1,
            None,
            vec![GlycanPosition {
                inner_depth: 1,
                series_number: 1,
                branch: Vec::new(),
                attachment: None,
            }],
        ),
        (
            1,
            None,
            vec![GlycanPosition {
                inner_depth: 2,
                series_number: 2,
                branch: Vec::new(),
                attachment: None,
            }],
        ),
    ] {
        let structure =
            GlycanStructure::from_short_iupac(codes[index].1, 0..codes[index].1.len(), 0).unwrap();
        if !structure
            .render(&mut footnotes)
            .to_svg(
                &mut html,
                Some("pep".to_string()),
                COLUMN_SIZE,
                SUGAR_SIZE,
                STROKE_SIZE,
                GlycanDirection::TopDown,
                root,
                &breaks,
                [0, 0, 0],
                [255, 255, 255],
                &mut footnotes,
            )
            .unwrap()
        {
            panic!("Rendering failed")
        }
    }

    write!(&mut html, "<hr>").unwrap();
    if !footnotes.is_empty() {
        write!(&mut html, "<ol>").unwrap();
        for note in footnotes {
            write!(&mut html, "<li>{note}</li>").unwrap();
        }
        write!(&mut html, "</ol><hr>").unwrap();
    }
    for (code, _) in &codes {
        write!(
            &mut html,
            "<image src=\"https://image.glycosmos.org/snfg/png/{code}\"/>"
        )
        .unwrap();
    }
    write!(&mut html, "</body></html>").unwrap();
    std::fs::write("../rendered_glycans.html", html).unwrap();
}
