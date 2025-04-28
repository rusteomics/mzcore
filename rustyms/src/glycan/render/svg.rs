use std::fmt::Write;

use itertools::Itertools;

use crate::glycan::{
    RenderedGlycan,
    render::element::{Element, TextAnchor, TextBaseline},
};

impl RenderedGlycan {
    /// Render this glycan as SVG. The SVG will be appended to the given buffer.
    ///  * `output`: the buffer to append the SVG to.
    ///
    /// # Errors
    /// If the underlying buffer errors the error is returned.
    pub fn to_svg(&self, mut output: impl Write) -> Result<(), std::fmt::Error> {
        fn clr(clr: Option<&[u8; 3]>) -> String {
            if let Some([r, g, b]) = clr {
                format!("rgb({r},{g},{b})")
            } else {
                "transparent".to_string()
            }
        }

        write!(
            output,
            "<svg width=\"{}\" height=\"{}\" style=\"--midpoint:{}\" class=\"glycan\">",
            self.size.0, self.size.1, self.midpoint
        )?;
        for element in &self.elements {
            match element {
                Element::Line {
                    from,
                    to,
                    stroke,
                    stroke_size,
                } => write!(
                    output,
                    "<line x1=\"{}\" y1=\"{}\" x2=\"{}\" y2=\"{}\" stroke=\"{}\" stroke-width=\"{stroke_size}px\"/>",
                    from.0,
                    from.1,
                    to.0,
                    to.1,
                    clr(Some(stroke))
                )?,
                Element::Circle {
                    r,
                    center,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header,
                } => write!(
                    output,
                    "<circle r=\"{r}\" cx=\"{}\" cy=\"{}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{stroke_size}px\"{svg_header}/>",
                    center.0,
                    center.1,
                    clr(fill.as_ref()),
                    clr(Some(stroke))
                )?,
                Element::Rectangle {
                    top,
                    w,
                    h,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header,
                } => write!(
                    output,
                    "<rect x=\"{}\" y=\"{}\" width=\"{w}\" height=\"{h}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{stroke_size}px\"{svg_header}/>",
                    top.0,
                    top.1,
                    clr(Some(fill)),
                    clr(Some(stroke))
                )?,
                Element::Polygon {
                    points,
                    fill,
                    stroke,
                    stroke_size,
                    svg_header,
                    bevel,
                } => write!(
                    output,
                    "<polygon points=\"{}\" fill=\"{}\" stroke=\"{}\" stroke-width=\"{stroke_size}px\"{svg_header}{}/>",
                    points.iter().map(|(a, b)| format!("{a} {b}")).join(" "),
                    clr(Some(fill)),
                    clr(Some(stroke)),
                    if *bevel {
                        " stroke-linejoin=\"bevel\""
                    } else {
                        ""
                    }
                )?,
                Element::Text {
                    text,
                    position,
                    anchor,
                    baseline,
                    fill,
                    size,
                    italic,
                } => write!(
                    output,
                    "<text x=\"{}\" y=\"{}\" fill=\"{}\" font-size=\"{size}px\" text-anchor=\"{}\" dominant-baseline=\"{}\"{}>{text}</text>",
                    position.0,
                    position.1,
                    clr(Some(fill)),
                    match anchor {
                        TextAnchor::Start => "start",
                        TextAnchor::Middle => "middle",
                        TextAnchor::End => "End",
                    },
                    match baseline {
                        TextBaseline::Hanging => "hanging",
                        TextBaseline::Middle => "middle",
                        TextBaseline::Ideographic => "ideographic",
                    },
                    if *italic {
                        " font-style=\"italic\""
                    } else {
                        ""
                    }
                )?,
                Element::Curve {
                    start,
                    points,
                    stroke,
                    stroke_size,
                } => write!(
                    output,
                    "<path d=\"M {} {} {}\" stroke=\"{}\" stroke-width=\"{stroke_size}px\" fill=\"none\"/>",
                    start.0,
                    start.1,
                    points
                        .iter()
                        .map(|(a, b, c, d)| format!("Q {a} {b} {c} {d}"))
                        .join(" "),
                    clr(Some(stroke))
                )?,
            }
        }
        write!(output, "</svg>")?;
        Ok(())
    }
}
