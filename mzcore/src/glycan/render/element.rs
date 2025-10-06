use std::f32::consts::PI;

use itertools::Itertools;

use crate::glycan::{
    GlycanBranchIndex, GlycanBranchMassIndex, GlycanDirection, GlycanPosition,
    render::{
        absolute::{AbsolutePositionedGlycan, OuterModifications},
        shape::{Colour, Shape},
    },
};

/// A rendered glycan, contains all information needed to render this to svg or a bitmap.
#[derive(Debug)]
pub struct RenderedGlycan {
    /// The size of the canvas
    pub(super) size: (f32, f32),
    /// All elements to be rendered
    pub(super) elements: Vec<Element>,
    /// The background colour
    pub(super) background: [u8; 3],
    /// Midpoint in pixels from the right for a top down glycan or in pixels from the top for a left to right glycan
    pub midpoint: f32,
}

#[derive(Clone, Debug)]
pub(super) enum Element {
    Line {
        from: (f32, f32),
        to: (f32, f32),
        stroke: [u8; 3],
        stroke_size: f32,
    },
    Circle {
        r: f32,
        center: (f32, f32),
        fill: Option<[u8; 3]>,
        stroke: [u8; 3],
        stroke_size: f32,
        svg_header: String,
    },
    Rectangle {
        top: (f32, f32),
        w: f32,
        h: f32,
        fill: [u8; 3],
        stroke: [u8; 3],
        stroke_size: f32,
        svg_header: String,
    },
    Polygon {
        points: Vec<(f32, f32)>,
        fill: [u8; 3],
        stroke: [u8; 3],
        stroke_size: f32,
        svg_header: String,
        bevel: bool,
    },
    Curve {
        start: (f32, f32),
        points: Vec<(f32, f32, f32, f32)>,
        stroke: [u8; 3],
        stroke_size: f32,
    },
    Text {
        text: String,
        position: (f32, f32),
        anchor: TextAnchor,
        baseline: TextBaseline,
        fill: [u8; 3],
        size: f32,
        italic: bool,
    },
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TextAnchor {
    Start,
    Middle,
    End,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(super) enum TextBaseline {
    Hanging,
    Middle,
    Ideographic,
}

/// The symbol or text to use at the base of a glycan.
///
#[doc = include_str!("../../../images/glycan_root.svg")]
///
/// _Glycan [G01670UQ](http://glytoucan.org/Structures/Glycans/G01670UQ) using the different root types: None, Line, Symbol, Text("pep"), Text("N"), Text("Arg")_
///
/// ```rust
/// # use rustyms::glycan::{GlycanStructure, GlycanDirection, GlycanRoot, GlycanSelection};
/// const COLUMN_SIZE: f32 = 30.0;
/// const SUGAR_SIZE: f32 = 15.0;
/// const STROKE_SIZE: f32 = 1.5;
/// let mut output = String::new();
/// let mut footnotes = Vec::new();
/// let short_iupac = "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"; // Definition for G01670UQ
/// let structure = GlycanStructure::from_short_iupac(short_iupac, 0..short_iupac.len(), 0).unwrap();
/// for root in [
///     GlycanRoot::None,
///     GlycanRoot::Line,
///     GlycanRoot::Symbol,
///     GlycanRoot::Text("pep".to_string()),
///     GlycanRoot::Text("N".to_string()),
///     GlycanRoot::Text("Arg".to_string()),
/// ] {
///     let rendered = structure
///         .render(
///             root,
///             COLUMN_SIZE,
///             SUGAR_SIZE,
///             STROKE_SIZE,
///             GlycanDirection::TopDown,
///             GlycanSelection::FULL,
///             [0, 0, 0],
///             [255, 255, 255],
///             &mut footnotes,
///         )
///         .unwrap();
///     rendered.to_svg(&mut output).unwrap();
/// }
/// ```
/// This examples shows how to generate SVGs for all the different root types as seen in the above picture.
/// Note that this writes all SVGs after each other to the variable `output`. Also note that this writes
/// all modifications that did not fit inside the image in the variable `footnotes` and this will need to
/// be dealt with by the caller, as indicated in [`GlycanStructure::render`](crate::glycan::GlycanStructure::render).
#[derive(Clone, Debug, Default, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum GlycanRoot {
    /// No symbol, this will also not draw a line from the root sugar
    #[default]
    None,
    /// No symbol, but this will draw a line from the root sugar
    Line,
    /// A tilde ('~') like symbol to indicate the full peptidoform
    Symbol,
    /// A piece of text, take care to not make this too big as it will be cut off in the image.
    /// Commonly used options are 'pep' to indicate the full peptidoform, or to indicate the
    /// attached amino acid any of 'Arg', or 'N'.
    Text(String),
}

/// The selected (part) of a glycan to render, using [`Self::FULL`] is a shortcut to get the full glycan.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GlycanSelection<'a> {
    /// A subtree of the glycan, with potentially a break of the root of the subtree and breaks in the branches.
    /// If no breaks are specified the full glycan is shown. The root is the first monosaccharide to be included
    /// in the rendering. The fragment will not include the indicated glycan positions for the branch breaks.
    Subtree(Option<&'a GlycanPosition>, &'a [GlycanPosition]),
    /// A single sugar, all it branches will be shown as broken.
    SingleSugar(&'a GlycanPosition),
}

impl GlycanSelection<'static> {
    /// A shorthand for a full glycan.
    pub const FULL: Self = Self::Subtree(None, &[]);
}

impl AbsolutePositionedGlycan {
    /// Render this glycan to the internal rendering representation, returns None if the root break contains an invalid position.
    #[expect(clippy::many_single_char_names, clippy::too_many_arguments)] // Doing geometry
    pub(super) fn render<'a>(
        &'a self,
        basis: GlycanRoot,
        column_size: f32,
        sugar_size: f32,
        stroke_size: f32,
        direction: GlycanDirection,
        selection: GlycanSelection<'a>,
        foreground: [u8; 3],
        background: [u8; 3],
        footnotes: &'a mut Vec<String>,
    ) -> Option<RenderedGlycan> {
        fn render_element(
            buffer: &mut Vec<Element>,
            element: &AbsolutePositionedGlycan,
            column_size: f32,
            sugar_size: f32,
            stroke_size: f32,
            direction: GlycanDirection,
            x_offset: f32,
            y_offset: f32,
            breaks: &[(usize, Vec<(GlycanBranchIndex, GlycanBranchMassIndex)>)],
            foreground: [u8; 3],
            background: [u8; 3],
            incoming_stroke: (f32, f32, f32, f32),
            footnotes: &mut Vec<String>,
        ) {
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
                        .any(|b| b.0 == 1 && b.1.first().map(|b| b.0) == Some(branch.branch_index))
                {
                    let base_y =
                        (raw_y - 0.5 + f32::from(side)).mul_add(column_size, stroke_size * 0.5);
                    let angle = f32::atan2(base_y - origin_y, base_x - origin_x);
                    buffer.push(Element::Line {
                        from: pick_point((origin_x, origin_y), direction),
                        to: pick_point((base_x, base_y), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    let x1 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).cos(), base_x);
                    let y1 =
                        (sugar_size / 2.0).mul_add((-0.5f32).mul_add(PI, -angle).sin(), base_y);
                    let x2 =
                        (sugar_size / 2.0).mul_add((-0.5f32).mul_add(PI, -angle).cos(), base_x);
                    let y2 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).sin(), base_y);
                    buffer.push(Element::Line {
                        from: pick_point((x1, y1), direction),
                        to: pick_point((x2, y2), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    let x3 = (stroke_size * 2.0).mul_add(-angle.cos(), x1);
                    let y3 = (stroke_size * 2.0).mul_add(-angle.sin(), y1);
                    buffer.push(Element::Line {
                        from: pick_point((x1, y1), direction),
                        to: pick_point((x3, y3), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    let offset = 0.25f32.mul_add(column_size, stroke_size);
                    let r = 0.25f32
                        .mul_add(column_size, -stroke_size)
                        .min(sugar_size * 0.25);
                    let adjusted_x = offset.mul_add(-angle.cos(), base_x);
                    let adjusted_y = offset.mul_add(-angle.sin(), base_y);
                    buffer.push(Element::Circle {
                        r,
                        center: pick_point((adjusted_x, adjusted_y), direction),
                        fill: None,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                    });
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
                    buffer.push(Element::Line {
                        from: pick_point((origin_x, origin_y), direction),
                        to: pick_point((base_x, base_y), direction),
                        stroke: foreground,
                        stroke_size,
                    });
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
                background
            } else {
                element.colour.rgb()
            };
            let title = format!(
                " data-sugar=\"{}\" data-position=\"{}-{}\"",
                element.title,
                element.position.inner_depth,
                element.position.branch.iter().map(|b| b.0).join(",")
            );
            match element.shape {
                Shape::Circle => buffer.push(Element::Circle {
                    r: sugar_size / 2.0,
                    center: pick_point(
                        (
                            (raw_x + element.mid_point) * column_size,
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    ),
                    fill: Some(fill),
                    stroke: foreground,
                    stroke_size,
                    svg_header: title,
                }),
                Shape::Square => buffer.push(Element::Rectangle {
                    top: pick_point(
                        (
                            (raw_x + element.mid_point).mul_add(column_size, -sugar_size / 2.0),
                            (raw_y + 0.5).mul_add(column_size, -sugar_size / 2.0),
                        ),
                        direction,
                    ),
                    w: sugar_size,
                    h: sugar_size,
                    fill,
                    stroke: foreground,
                    stroke_size,
                    svg_header: title,
                }),
                Shape::Rectangle => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    buffer.push(Element::Rectangle {
                        top: (base_x - sugar_size / 2.0, base_y - sugar_size / 4.0),
                        w: sugar_size,
                        h: sugar_size / 2.0,
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                    });
                }
                Shape::Triangle => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x2, y1), (x3, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::LeftPointingTriangle => {
                    let x1 = (raw_x + element.mid_point).mul_add(column_size, -sugar_size / 2.0);
                    let x2 = x1 + sugar_size;
                    let y1 = (raw_y + 0.5).mul_add(column_size, -sugar_size / 2.0);
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![
                            pick_point((x1, y2), direction),
                            pick_point((x2, y1), direction),
                            pick_point((x2, y3), direction),
                        ],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::RightPointingTriangle => {
                    let x1 = (raw_x + element.mid_point).mul_add(column_size, -sugar_size / 2.0);
                    let x2 = x1 + sugar_size;
                    let y1 = (raw_y + 0.5).mul_add(column_size, -sugar_size / 2.0);
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![
                            pick_point((x1, y1), direction),
                            pick_point((x2, y2), direction),
                            pick_point((x1, y3), direction),
                        ],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::Diamond => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![(x2, y1), (x3, y2), (x2, y3), (x1, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::FlatDiamond => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 4.0;
                    let y2 = y1 + sugar_size / 4.0;
                    let y3 = y1 + sugar_size / 2.0;

                    buffer.push(Element::Polygon {
                        points: vec![(x2, y1), (x3, y2), (x2, y3), (x1, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::Hexagon => {
                    let a = sugar_size / 2.0 / 3.0_f32.sqrt();
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + a;
                    let x3 = x1 + sugar_size - a;
                    let x4 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 4.0;
                    let y2 = y1 + sugar_size / 4.0;
                    let y3 = y1 + sugar_size / 2.0;

                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x2, y1), (x3, y1), (x4, y2), (x3, y3), (x2, y3)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::Pentagon => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let a = (18.0 / 360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let b = (18.0 / 360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let c = (36.0 / 360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let d = (36.0 / 360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let x1 = base_x - a;
                    let x2 = base_x - d;
                    let x3 = base_x;
                    let x4 = base_x + d;
                    let x5 = base_x + a;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0 - b;
                    let y3 = y1 + sugar_size / 2.0 + c;

                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x3, y1), (x5, y2), (x4, y3), (x2, y3)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::Star => {
                    // The Phi constant, the ratio for the "golden ratio"
                    const PHI: f32 = 1.618_034_f32;
                    // Calculate sizes of parts of the pentagram
                    let a = (18.0 / 360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let b = (18.0 / 360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let c = (36.0 / 360.0 * 2.0 * PI).cos() * sugar_size / 2.0;
                    let d = (36.0 / 360.0 * 2.0 * PI).sin() * sugar_size / 2.0;
                    let e = 2.0 * a / (54.0 / 360.0 * 2.0 * PI).sin() / (1.0 + 1.0 / PHI);
                    let f = (18.0 / 360.0 * 2.0 * PI)
                        .cos()
                        .mul_add(e, -(sugar_size / 2.0));
                    let g = (18.0 / 360.0 * 2.0 * PI).sin() * e;
                    let h = (sugar_size / 2.0 - b) * (18.0 / 360.0 * 2.0 * PI).tan();
                    let j = (18.0 / 360.0 * 2.0 * PI).tan() * g;
                    // Calculate the positions of the pentagram points
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
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

                    buffer.push(Element::Polygon {
                        points: vec![
                            (x1, y2),
                            (x4, y2),
                            (x5, y1),
                            (x6, y2),
                            (x9, y2),
                            (x7, y3),
                            (x8, y5),
                            (x5, y4),
                            (x2, y5),
                            (x3, y3),
                        ],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::CrossedSquare => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let y1 = base_y - sugar_size / 2.0;
                    let x2 = x1 + sugar_size;
                    let y2 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![(x1, y1), (x2, y1), (x2, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x1, y1), (x1, y2), (x2, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x1, y1), (x2, y1), (x2, y2), (x1, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::DividedDiamond => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size / 2.0;
                    let y3 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x2, y1), (x3, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x2, y3), (x3, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x1, y2), (x2, y1), (x3, y2), (x2, y3)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
                Shape::DividedTriangle => {
                    let (base_x, base_y) = pick_point(
                        (
                            ((raw_x + element.mid_point) * column_size),
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    );
                    let x1 = base_x - sugar_size / 2.0;
                    let x2 = x1 + sugar_size / 2.0;
                    let x3 = x1 + sugar_size;
                    let y1 = base_y - sugar_size / 2.0;
                    let y2 = y1 + sugar_size;

                    buffer.push(Element::Polygon {
                        points: vec![(x2, y1), (x3, y2), (x2, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x2, y1), (x1, y2), (x2, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: String::new(),
                        bevel: true,
                    });
                    buffer.push(Element::Polygon {
                        points: vec![(x2, y1), (x3, y2), (x1, y2)],
                        fill,
                        stroke: foreground,
                        stroke_size,
                        svg_header: title,
                        bevel: false,
                    });
                }
            }
            if !element.inner_modifications.is_empty() {
                buffer.push(Element::Text {
                    text: element.inner_modifications.clone(),
                    position: pick_point(
                        (
                            (raw_x + element.mid_point) * column_size,
                            (raw_y + 0.5) * column_size,
                        ),
                        direction,
                    ),
                    anchor: TextAnchor::Middle,
                    baseline: TextBaseline::Middle,
                    fill: foreground,
                    size: sugar_size / 2.0,
                    italic: true,
                });
            }
            if let Some((pos_x, pos_y, anchor, text)) = text_location(
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
                buffer.push(Element::Text {
                    text,
                    position: (pos_x, pos_y),
                    anchor,
                    baseline: TextBaseline::Hanging,
                    fill: foreground,
                    size: sugar_size / 2.0,
                    italic: false,
                });
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
                        .any(|b| b.0 == 1 && b.1.first().map(|b| b.0) == Some(branch.branch_index)))
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
                                (total_branches > 1
                                    && b.1.first().map(|b| b.0) == Some(branch.branch_index)
                                    || total_branches == 1)
                                    && b.0 > 0
                            })
                            .map(|b| (b.0 - 1, b.1[usize::from(total_branches > 1)..].to_vec()))
                            .collect_vec(),
                        foreground,
                        background,
                        strokes[index + 1],
                        footnotes,
                    );
                }
            }
        }

        let sub_tree = self.get_subtree(selection)?;

        let width =
            (sub_tree.tree.x + sub_tree.tree.width - sub_tree.left_offset - sub_tree.right_offset)
                * column_size;
        let depth = sub_tree.depth as f32 + if sub_tree.break_top { 0.75 } else { 0.0 };
        let height = depth.mul_add(
            column_size,
            if sub_tree.break_bottom {
                3.5 * stroke_size
            } else {
                (match basis {
                    GlycanRoot::None => 0.0_f32,
                    GlycanRoot::Line | GlycanRoot::Symbol => 0.5,
                    GlycanRoot::Text(_) => 1.0,
                }) * column_size
            },
        );

        let size = pick_point((width, height), direction);

        let mut buffer = Vec::new();
        let stroke = if sub_tree.break_bottom {
            let base_x =
                (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset) * column_size;
            let base_y = depth.mul_add(column_size, stroke_size * 3.0);
            buffer.push(Element::Line {
                from: pick_point((base_x, (depth - 0.5) * column_size), direction),
                to: pick_point((base_x, base_y), direction),
                stroke: foreground,
                stroke_size,
            });
            buffer.push(Element::Line {
                from: pick_point((base_x - sugar_size / 2.0, base_y), direction),
                to: pick_point((base_x + sugar_size / 2.0, base_y), direction),
                stroke: foreground,
                stroke_size,
            });
            let bs_x = base_x - sugar_size / 2.0;
            buffer.push(Element::Line {
                from: pick_point((bs_x, stroke_size.mul_add(-2.0, base_y)), direction),
                to: pick_point((bs_x, base_y), direction),
                stroke: foreground,
                stroke_size,
            });
            (base_x, (depth - 0.5) * column_size, base_x, base_y)
        } else {
            match basis {
                GlycanRoot::None => (
                    sub_tree.tree.x + sub_tree.tree.mid_point,
                    (depth + 0.5) * column_size,
                    sub_tree.tree.x + sub_tree.tree.mid_point,
                    (depth + 0.5) * column_size,
                ),
                GlycanRoot::Line => {
                    let base_x = (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset)
                        * column_size;
                    let base_y = depth.mul_add(column_size, (column_size - sugar_size) / 2.0);
                    buffer.push(Element::Line {
                        from: pick_point((base_x, (depth - 0.5) * column_size), direction),
                        to: pick_point((base_x, base_y), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    (base_x, (depth - 0.5) * column_size, base_x, base_y)
                }
                GlycanRoot::Symbol => {
                    let base_x = (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset)
                        * column_size;
                    let base_y = depth.mul_add(column_size, (column_size - sugar_size) / 2.0);
                    buffer.push(Element::Line {
                        from: pick_point((base_x, (depth - 0.5) * column_size), direction),
                        to: pick_point((base_x, base_y), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    buffer.push(Element::Curve {
                        start: pick_point(
                            (base_x - (sugar_size * 0.75).min(column_size * 0.5), base_y),
                            direction,
                        ),
                        points: vec![
                            pick_double_point(
                                (
                                    sugar_size.mul_add(-0.5, base_x),
                                    sugar_size.mul_add(0.5, base_y),
                                    base_x,
                                    base_y,
                                ),
                                direction,
                            ),
                            pick_double_point(
                                (
                                    sugar_size.mul_add(0.5, base_x),
                                    sugar_size.mul_add(-0.5, base_y),
                                    base_x + (sugar_size * 0.75).min(column_size * 0.5),
                                    base_y,
                                ),
                                direction,
                            ),
                        ],
                        stroke: foreground,
                        stroke_size,
                    });
                    (base_x, (depth - 0.5) * column_size, base_x, base_y)
                }
                GlycanRoot::Text(basis) => {
                    let base_x = (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset)
                        * column_size;
                    let base_y = depth.mul_add(column_size, (column_size - sugar_size) / 2.0);
                    buffer.push(Element::Line {
                        from: pick_point((base_x, (depth - 0.5) * column_size), direction),
                        to: pick_point((base_x, base_y), direction),
                        stroke: foreground,
                        stroke_size,
                    });
                    if direction == GlycanDirection::TopDown {
                        buffer.push(Element::Text {
                            text: basis,
                            position: (base_x, base_y + sugar_size),
                            anchor: TextAnchor::Middle,
                            baseline: TextBaseline::Ideographic,
                            fill: foreground,
                            size: sugar_size,
                            italic: false,
                        });
                    } else {
                        buffer.push(Element::Text {
                            text: basis,
                            position: ((depth + 1.0) * column_size, base_x),
                            anchor: TextAnchor::End,
                            baseline: TextBaseline::Middle,
                            fill: foreground,
                            size: sugar_size,
                            italic: false,
                        });
                    }
                    (base_x, (depth - 0.5) * column_size, base_x, base_y)
                }
            }
        };

        // If the full glycan has broken off immediately draw the break symbol
        if sub_tree.branch_breaks.iter().any(|r| r.0 == 0) {
            let origin_y = depth.mul_add(column_size, (column_size - sugar_size) / 2.0);
            let base_x =
                (sub_tree.tree.x + sub_tree.tree.mid_point - sub_tree.left_offset) * column_size;
            let base_y = (depth - 0.5).mul_add(column_size, stroke_size * 0.5);
            let angle = -PI / 2.0; // Always straight
            buffer.push(Element::Line {
                from: pick_point((base_x, origin_y), direction),
                to: pick_point((base_x, base_y), direction),
                stroke: foreground,
                stroke_size,
            });
            let x1 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).cos(), base_x);
            let y1 = (sugar_size / 2.0).mul_add((-0.5f32).mul_add(PI, -angle).sin(), base_y);
            let x2 = (sugar_size / 2.0).mul_add((-0.5f32).mul_add(PI, -angle).cos(), base_x);
            let y2 = (sugar_size / 2.0).mul_add(0.5f32.mul_add(PI, -angle).sin(), base_y);
            buffer.push(Element::Line {
                from: pick_point((x1, y1), direction),
                to: pick_point((x2, y2), direction),
                stroke: foreground,
                stroke_size,
            });
            let x3 = (stroke_size * 2.0).mul_add(-angle.cos(), x1);
            let y3 = (stroke_size * 2.0).mul_add(-angle.sin(), y1);
            buffer.push(Element::Line {
                from: pick_point((x1, y1), direction),
                to: pick_point((x3, y3), direction),
                stroke: foreground,
                stroke_size,
            });
            let offset = 0.25f32.mul_add(column_size, stroke_size);
            let r = 0.25f32
                .mul_add(column_size, -stroke_size)
                .min(sugar_size * 0.25);
            let adjusted_x = offset.mul_add(-angle.cos(), base_x);
            let adjusted_y = offset.mul_add(-angle.sin(), base_y);
            buffer.push(Element::Circle {
                r,
                center: pick_point((adjusted_x, adjusted_y), direction),
                fill: None,
                stroke: foreground,
                stroke_size,
                svg_header: String::new(),
            });
        } else {
            render_element(
                &mut buffer,
                sub_tree.tree,
                column_size,
                sugar_size,
                stroke_size,
                direction,
                sub_tree.left_offset,
                sub_tree.tree.y as f32 - (depth - 1.0),
                &sub_tree.branch_breaks,
                foreground,
                background,
                pick_box(stroke, direction),
                footnotes,
            );
        }

        Some(RenderedGlycan {
            size,
            elements: buffer,
            background,
            midpoint: (sub_tree.tree.mid_point - sub_tree.left_offset) * column_size,
        })
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
) -> Option<(f32, f32, TextAnchor, String)> {
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
                TextAnchor::Middle,
            ),
            (
                (
                    (column_size - text_width).mul_add(0.5, position.0),
                    position.1 + column_size - text_height,
                    (column_size + text_width).mul_add(0.5, position.0),
                    position.1 + column_size,
                ),
                TextAnchor::Middle,
            ),
            (
                (
                    position.0 + column_size.mul_add(0.5, -stroke_size) - text_width,
                    position.1,
                    position.0 + column_size.mul_add(0.5, -stroke_size),
                    position.1 + text_height,
                ),
                TextAnchor::End,
            ),
            (
                (
                    stroke_size.mul_add(2.0, column_size.mul_add(0.5, position.0)),
                    position.1 + column_size - text_height,
                    stroke_size.mul_add(2.0, column_size.mul_add(0.5, position.0) + text_width),
                    position.1 + column_size,
                ),
                TextAnchor::Start,
            ),
            (
                (
                    position.0,
                    position.1,
                    position.0 + text_width,
                    position.1 + text_height,
                ),
                TextAnchor::Start,
            ),
            (
                (
                    position.0 + column_size - text_width,
                    position.1,
                    position.0 + column_size,
                    position.1 + text_height,
                ),
                TextAnchor::End,
            ),
            (
                (
                    position.0,
                    position.1 + column_size - text_height,
                    position.0 + text_width,
                    position.1 + column_size,
                ),
                TextAnchor::Start,
            ),
            (
                (
                    position.0 + column_size - text_width,
                    position.1 + column_size - text_height,
                    position.0 + column_size,
                    position.1 + column_size,
                ),
                TextAnchor::End,
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
                TextAnchor::Start => (option.0, option.1),
                TextAnchor::Middle => (f32::midpoint(option.0, option.2), option.1),
                TextAnchor::End => (option.2, option.1),
            };
            return Some((x, y, *anchor, text));
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
            TextAnchor::End,
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

fn pick_point<T>(a: (T, T), direction: GlycanDirection) -> (T, T) {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        (a.1, a.0)
    }
}

fn pick_double_point<T>(a: (T, T, T, T), direction: GlycanDirection) -> (T, T, T, T) {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        (a.1, a.0, a.3, a.2)
    }
}

fn pick_box<T>(a: (T, T, T, T), direction: GlycanDirection) -> (T, T, T, T) {
    if direction == GlycanDirection::TopDown {
        a
    } else {
        (a.1, a.0, a.3, a.2)
    }
}
