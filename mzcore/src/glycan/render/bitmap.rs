use itertools::Itertools;
use swash::{
    FontRef,
    scale::{Render, ScaleContext, Source},
};
use zeno::{Fill, Format, Mask, PathBuilder, Point, Scratch, Stroke, Vector};

use crate::glycan::{RenderedGlycan, render::element::Element};

use super::element::{TextAnchor, TextBaseline};

impl RenderedGlycan {
    /// Render this glycan as an RGBA bitmap.
    ///  * `format`: the used strategy for antialiasing.
    ///  * `font`: the font for rendering text.
    ///  * `context`: the context for caching rendering text.
    /// # Panics
    /// If the glyph renderer failed. See [`swash::scale::Render::render`].
    pub fn to_bitmap(
        &self,
        format: Format,
        font: FontRef,
        context: &mut ScaleContext,
    ) -> (Vec<u8>, usize) {
        let mask_factor = if format == Format::Alpha { 1 } else { 4 };
        let image_width = self.size.0.ceil() as usize;
        let mut image = std::iter::repeat_n(
            [
                self.background[0],
                self.background[1],
                self.background[2],
                0,
            ],
            image_width * self.size.1.ceil() as usize,
        )
        .flatten()
        .collect_vec();

        let mut scratch = Scratch::new();
        let mut stroke_mask = Vec::new();
        let mut fill_mask = Vec::new();
        for element in &self.elements {
            // Draw into the mask
            let (x, y, mask_width, fill, stroke) = match element {
                Element::Line {
                    from,
                    to,
                    stroke,
                    stroke_size,
                } => {
                    let xmin = (from.0.min(to.0) - stroke_size).floor();
                    let xmax = (from.0.max(to.0) + stroke_size).ceil();
                    let ymin = (from.1.min(to.1) - stroke_size).floor();
                    let ymax = (from.1.max(to.1) + stroke_size).ceil();
                    let width = (xmax - xmin) as usize;
                    let height = (ymax - ymin) as usize;
                    let commands = vec![
                        zeno::Command::MoveTo(Vector::new(
                            from.0 - xmin + stroke_size / 2.0,
                            from.1 - ymin + stroke_size / 2.0,
                        )),
                        zeno::Command::LineTo(Vector::new(
                            to.0 - xmin + stroke_size / 2.0,
                            to.1 - ymin + stroke_size / 2.0,
                        )),
                        zeno::Command::Close,
                    ];
                    stroke_mask.fill(0);
                    stroke_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Stroke::new(*stroke_size))
                        .size(width as u32, height as u32)
                        .render_into(&mut stroke_mask, None);
                    (xmin as usize, ymin as usize, width, None, Some(*stroke))
                }
                Element::Circle {
                    r,
                    center,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header: _,
                } => {
                    let width = (center.0.fract() + r * 2.0 + stroke_size).ceil() as usize;
                    let height = (center.1.fract() + r * 2.0 + stroke_size).ceil() as usize;
                    let mut commands = Vec::new();
                    commands.add_circle(
                        (
                            center.0.fract() + r + stroke_size / 2.0,
                            center.1.fract() + r + stroke_size / 2.0,
                        ),
                        *r,
                    );
                    if fill.is_some() {
                        fill_mask.fill(0);
                        fill_mask.resize(height * width * mask_factor, 0);
                        Mask::with_scratch(&commands, &mut scratch)
                            .format(format)
                            .style(Fill::NonZero)
                            .size(width as u32, height as u32)
                            .render_into(&mut fill_mask, None);
                    }
                    stroke_mask.fill(0);
                    stroke_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Stroke::new(*stroke_size))
                        .size(width as u32, height as u32)
                        .render_into(&mut stroke_mask, None);
                    (
                        (center.0 - r) as usize,
                        (center.1 - r) as usize,
                        width,
                        *fill,
                        Some(*stroke),
                    )
                }
                Element::Rectangle {
                    top,
                    w,
                    h,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header: _,
                } => {
                    let width = (top.0.fract() + w + stroke_size).ceil() as usize;
                    let height = (top.1.fract() + h + stroke_size).ceil() as usize;
                    let mut commands = Vec::new();
                    commands.add_rect(
                        (
                            top.0.fract() + stroke_size / 2.0,
                            top.1.fract() + stroke_size / 2.0,
                        ),
                        *w,
                        *h,
                    );
                    fill_mask.fill(0);
                    fill_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Fill::NonZero)
                        .size(width as u32, height as u32)
                        .render_into(&mut fill_mask, None);
                    stroke_mask.fill(0);
                    stroke_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Stroke::new(*stroke_size))
                        .size(width as u32, height as u32)
                        .render_into(&mut stroke_mask, None);
                    (
                        (top.0 - stroke_size / 2.0) as usize,
                        (top.1 - stroke_size / 2.0) as usize,
                        width,
                        Some(*fill),
                        Some(*stroke),
                    )
                }
                Element::Polygon {
                    points,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header: _,
                    bevel,
                } => {
                    let (xmin, xmax, ymin, ymax) = points
                        .iter()
                        .fold((f32::MAX, f32::MIN, f32::MAX, f32::MIN), |acc, (x, y)| {
                            (acc.0.min(*x), acc.1.max(*x), acc.2.min(*y), acc.3.max(*y))
                        });
                    let xmin = (xmin - stroke_size).floor();
                    let xmax = (xmax + stroke_size).ceil();
                    let ymin = (ymin - stroke_size).floor();
                    let ymax = (ymax + stroke_size).ceil();
                    let width = (xmax - xmin) as usize;
                    let height = (ymax - ymin) as usize;
                    let mut commands = Vec::with_capacity(points.len() + 2);
                    commands.push(zeno::Command::MoveTo(Point::new(
                        points[0].0 - xmin + stroke_size / 2.0,
                        points[0].1 - ymin + stroke_size / 2.0,
                    )));
                    for point in points {
                        commands.push(zeno::Command::LineTo(Point::new(
                            point.0 - xmin + stroke_size / 2.0,
                            point.1 - ymin + stroke_size / 2.0,
                        )));
                    }
                    commands.push(zeno::Command::Close);
                    fill_mask.fill(0);
                    fill_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Fill::NonZero)
                        .size(width as u32, height as u32)
                        .render_into(&mut fill_mask, None);
                    stroke_mask.fill(0);
                    stroke_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Stroke::new(*stroke_size).join(if *bevel {
                            zeno::Join::Bevel
                        } else {
                            zeno::Join::Miter
                        }))
                        .size(width as u32, height as u32)
                        .render_into(&mut stroke_mask, None);
                    (
                        xmin as usize,
                        ymin as usize,
                        width,
                        Some(*fill),
                        Some(*stroke),
                    )
                }
                Element::Text {
                    text,
                    position,
                    anchor,
                    baseline,
                    fill,
                    size,
                    italic: _, // Needs a separate font
                } => {
                    let mut scaler = context.builder(font).size(*size).hint(true).build();
                    let metrics = font.metrics(&[]);
                    let normalisation_factor = size / f32::from(metrics.units_per_em);
                    let y_offset = (match baseline {
                        TextBaseline::Hanging => metrics.ascent,
                        TextBaseline::Middle => metrics.ascent - metrics.x_height / 2.0,
                        TextBaseline::Ideographic => metrics.ascent + metrics.descent,
                    })
                    .mul_add(-normalisation_factor, position.1);
                    let mut width = 0.0;
                    for c in text.chars() {
                        let id = font.charmap().map(c);
                        width += font.glyph_metrics(&[]).advance_width(id);
                    }

                    let x_offset = (match anchor {
                        TextAnchor::Start => 0.0,
                        TextAnchor::Middle => width / 2.0,
                        TextAnchor::End => width,
                    })
                    .mul_add(-normalisation_factor, position.0);

                    let mut offset = 0.0;
                    for c in text.chars() {
                        let id = font.charmap().map(c);
                        let glyph_metrics = font.glyph_metrics(&[]);
                        let mask = Render::new(&[Source::Outline])
                            .format(format)
                            .offset(Vector::new(
                                (x_offset + offset).fract(),
                                y_offset.fract() - 1.0,
                            ))
                            .render(&mut scaler, id)
                            .unwrap();
                        draw_mask(
                            (&mut image, image_width),
                            (&mask.data, mask.placement.width as usize),
                            (x_offset + offset + mask.placement.left as f32) as usize,
                            (y_offset + mask.placement.top as f32) as usize,
                            *fill,
                            format,
                        );

                        offset += glyph_metrics.advance_width(id) * normalisation_factor;
                    }
                    (0, 0, 0, None, None)
                }
                Element::Curve {
                    start,
                    points,
                    stroke,
                    stroke_size,
                } => {
                    let (xmin, xmax, ymin, ymax) = points.iter().fold(
                        (f32::MAX, f32::MIN, f32::MAX, f32::MIN),
                        |acc, (a, b, x, y)| {
                            (
                                acc.0.min(*x).min(*a),
                                acc.1.max(*x).max(*a),
                                acc.2.min(*y).min(*b),
                                acc.3.max(*y).max(*b),
                            )
                        },
                    );
                    let xmin = (xmin - stroke_size).floor();
                    let xmax = (xmax + stroke_size).ceil();
                    let ymin = (ymin - stroke_size).floor();
                    let ymax = (ymax + stroke_size).ceil();
                    let width = (xmax - xmin) as usize;
                    let height = (ymax - ymin) as usize;
                    let mut commands = Vec::with_capacity(points.len() + 1);
                    commands.push(zeno::Command::MoveTo(Point::new(
                        start.0 - xmin + stroke_size / 2.0,
                        start.1 - ymin + stroke_size / 2.0,
                    )));
                    for point in points {
                        commands.push(zeno::Command::QuadTo(
                            Point::new(
                                point.0 - xmin + stroke_size / 2.0,
                                point.1 - ymin + stroke_size / 2.0,
                            ),
                            Point::new(
                                point.2 - xmin + stroke_size / 2.0,
                                point.3 - ymin + stroke_size / 2.0,
                            ),
                        ));
                    }
                    stroke_mask.fill(0);
                    stroke_mask.resize(height * width * mask_factor, 0);
                    Mask::with_scratch(&commands, &mut scratch)
                        .format(format)
                        .style(Stroke::new(*stroke_size))
                        .size(width as u32, height as u32)
                        .render_into(&mut stroke_mask, None);
                    (xmin as usize, ymin as usize, width, None, Some(*stroke))
                }
            };
            if let Some(fill) = fill {
                draw_mask(
                    (&mut image, image_width),
                    (&fill_mask, mask_width),
                    x,
                    y,
                    fill,
                    format,
                );
            }
            if let Some(stroke) = stroke {
                draw_mask(
                    (&mut image, image_width),
                    (&stroke_mask, mask_width),
                    x,
                    y,
                    stroke,
                    format,
                );
            }
        }
        (image, image_width)
    }
}

/// Draw the specified mask onto the specified image
#[allow(clippy::identity_op, clippy::needless_pass_by_value)] // I like the + 0 in position calculations for symmetry reasons and image tuple looks makes no sense to pass by reference
fn draw_mask(
    image: (&mut [u8], usize),
    mask: (&[u8], usize),
    x: usize,
    y: usize,
    colour: [u8; 3],
    format: Format,
) {
    let mask_factor = if format == Format::Alpha { 1 } else { 4 };
    let mask_height = mask.0.len() / mask_factor / mask.1;
    for r in 0..mask_height {
        for w in 0..mask.1 {
            let image_pos = ((r + y) * image.1 + (w + x)) * 4;
            let mask_pos = (r * mask.1 + w) * mask_factor;

            if image_pos >= image.0.len() || mask_pos >= mask.0.len() {
                continue;
            }

            if format == Format::Alpha {
                image.0[image_pos + 0] = blend(mask.0[mask_pos], colour[0], image.0[image_pos + 0]);
                image.0[image_pos + 1] = blend(mask.0[mask_pos], colour[1], image.0[image_pos + 1]);
                image.0[image_pos + 2] = blend(mask.0[mask_pos], colour[2], image.0[image_pos + 2]);
            } else {
                image.0[image_pos + 0] =
                    blend(mask.0[mask_pos + 0], colour[0], image.0[image_pos + 0]);
                image.0[image_pos + 1] =
                    blend(mask.0[mask_pos + 1], colour[1], image.0[image_pos + 1]);
                image.0[image_pos + 2] =
                    blend(mask.0[mask_pos + 2], colour[2], image.0[image_pos + 2]);
            }
            image.0[image_pos + 3] = 255;
        }
    }
}

const fn blend(alpha: u8, foreground: u8, background: u8) -> u8 {
    (((alpha as u16 * foreground as u16) + (255 - alpha) as u16 * background as u16) / 255) as u8
}
