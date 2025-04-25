use serde::*;
use std::{
    fmt,
    ops::{Bound, RangeBounds},
};

/// A struct to define the context of an error message
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq, Eq, Hash)]
pub enum Context {
    /// When no context can be given
    None,
    /// When only a line (e.g. in a file) can be shown
    Show {
        /// The line to be shown to the user (e.g. filename)
        line: String,
    },
    /// When a full line is faulty and no special position can be annotated
    FullLine {
        /// The line number to recognise where the error is located
        line_index: usize,
        /// The line to show the issue itself
        line: String,
    },
    /// When a special position can be annotated on a line.
    /// ```text
    ///      |
    /// 104  | ATOM      O  N   MET A   1      27.251  24.447   2.594  1.00 11.79           N
    ///      |        ^^^^
    ///        <-   -><-->
    /// ```
    /// The first space (annotated by `<-`, `->`) is the offset, in this case 7. The
    /// second space is the length, in this case 4.
    Line {
        /// The line index to recognise where the error is located.
        line_index: Option<usize>,
        /// The line to show the issue itself.
        line: String,
        /// The offset of the special position to be annotated.
        offset: usize,
        /// The length of the special position to be annotated.
        length: usize,
    },
    /// To show multiple lines where an error occurred.
    Range {
        /// The index of the first line
        start_line_index: usize,
        /// The lines to show
        lines: Vec<String>,
        /// The possible offset of the first line, will be padded with spaces
        offset: usize,
    },
    /// To show multiple lines where an error occurred.
    RangeHighlights {
        /// The index of the first line
        start_line_index: usize,
        /// The lines to show
        lines: Vec<String>,
        /// Highlights defined by the line (relative to the set of lines given), start column in that line and length of highlight
        highlights: Vec<(usize, usize, usize)>,
    },
    /// To show multiple contexts
    Multiple {
        /// The contexts to show
        contexts: Vec<(Option<String>, Context)>,
    },
}

#[expect(clippy::needless_pass_by_value)] // the impl ToString should be passed like this, otherwise &str gives errors
impl Context {
    /// Creates a new context when no context can be given
    pub const fn none() -> Self {
        Self::None
    }

    /// Creates a new context when only a line (eg filename) can be shown
    pub fn show(line: impl ToString) -> Self {
        Self::Show {
            line: line.to_string().replace('\t', " "),
        }
    }

    /// Creates a new context when a full line is faulty and no special position can be annotated
    pub fn full_line(line_index: usize, line: impl ToString) -> Self {
        Self::FullLine {
            line_index,
            line: line.to_string().replace('\t', " "),
        }
    }

    /// Creates a new context when a special position can be annotated on a line
    pub fn line(
        line_index: Option<usize>,
        line: impl ToString,
        offset: usize,
        length: usize,
    ) -> Self {
        Self::Line {
            line_index,
            line: line.to_string().replace('\t', " "),
            offset,
            length,
        }
    }

    /// Create a context highlighting a certain range on a single line
    pub fn line_range(
        line_index: Option<usize>,
        line: impl ToString,
        range: impl RangeBounds<usize>,
    ) -> Self {
        let line = line.to_string();
        match (range.start_bound(), range.end_bound()) {
            (Bound::Unbounded, Bound::Unbounded) => {
                line_index.map_or_else(|| Self::show(&line), |i| Self::full_line(i, &line))
            }
            (start, end) => {
                let start = match start {
                    Bound::Excluded(n) => n + 1,
                    Bound::Included(n) => *n,
                    Bound::Unbounded => 0,
                };
                Self::line(
                    line_index,
                    &line,
                    start,
                    match end {
                        Bound::Excluded(n) => n - 1,
                        Bound::Included(n) => *n,
                        Bound::Unbounded => line.chars().count(),
                    }
                    .saturating_sub(start),
                )
            }
        }
    }

    /// Creates a new context to highlight a certain position
    #[expect(clippy::unwrap_used, clippy::missing_panics_doc)]
    pub fn position(pos: &FilePosition<'_>) -> Self {
        if pos.text.is_empty() {
            Self::Line {
                line_index: Some(pos.line_index),
                line: String::new(),
                offset: 0,
                length: 3,
            }
        } else {
            Self::Line {
                line_index: Some(pos.line_index),
                line: pos
                    .text
                    .lines()
                    .next()
                    .unwrap()
                    .to_string()
                    .replace('\t', " "),
                offset: 0,
                length: 3,
            }
        }
    }

    /// Creates a new context from a start and end point within a single file
    pub fn range(start: &FilePosition<'_>, end: &FilePosition<'_>) -> Self {
        if start.line_index == end.line_index {
            Self::Line {
                line_index: Some(start.line_index),
                line: start.text[..(end.column - start.column)].to_string(),
                offset: start.column,
                length: end.column - start.column,
            }
        } else {
            Self::Range {
                start_line_index: start.line_index,
                lines: start
                    .text
                    .lines()
                    .take(end.line_index - start.line_index)
                    .map(ToString::to_string)
                    .collect::<Vec<String>>(),
                offset: start.column,
            }
        }
    }

    /// Overwrite the line number with the given number, if applicable
    #[must_use]
    pub fn overwrite_line_number(self, line_index: usize) -> Self {
        match self {
            Self::FullLine { line, .. } => Self::FullLine { line_index, line },
            Self::Line {
                line,
                offset,
                length,
                ..
            } => Self::Line {
                line_index: Some(line_index),
                line,
                offset,
                length,
            },
            Self::Range { lines, offset, .. } => Self::Range {
                start_line_index: line_index,
                offset,
                lines,
            },
            Self::RangeHighlights {
                lines, highlights, ..
            } => Self::RangeHighlights {
                start_line_index: line_index,
                lines,
                highlights,
            },
            Self::Multiple { contexts } => Self::Multiple {
                contexts: contexts
                    .into_iter()
                    .map(|(l, c)| (l, c.overwrite_line_number(line_index)))
                    .collect(),
            },
            n => n,
        }
    }

    /// Display this context, with an optional note after the context.
    /// # Errors
    /// If the underlying formatter errors.
    fn display(&self, f: &mut fmt::Formatter<'_>, note: Option<&str>) -> fmt::Result {
        const MAX_COLS: usize = 95;
        let mut tail = true; // End with a tailing line ╵
        #[expect(
            clippy::cast_sign_loss,
            clippy::cast_precision_loss,
            clippy::cast_possible_truncation
        )]
        let get_margin = |n| ((n + 1) as f64).log10().max(1.0).ceil() as usize;
        let margin = match self {
            Self::None | Self::Multiple { .. } => 0,
            Self::Show { .. } => 2,
            Self::FullLine { line_index: n, .. } => get_margin(*n),
            Self::Line { line_index: n, .. } => n.map_or(0, get_margin),
            Self::Range {
                start_line_index: n,
                lines: l,
                ..
            }
            | Self::RangeHighlights {
                start_line_index: n,
                lines: l,
                ..
            } => get_margin(n + l.len()),
        };
        match self {
            Self::None => {
                return Ok(());
            }
            Self::Show { line } => {
                write!(f, "\n{:pad$} ╷\n{:pad$} │ {}", "", "", line, pad = margin)?;
            }
            Self::FullLine { line_index, line } => write!(
                f,
                "\n{:pad$} ╷\n{:<pad$} │ {}",
                "",
                line_index + 1,
                line,
                pad = margin
            )?,
            Self::Line {
                line_index,
                line,
                offset,
                length,
            } => {
                let (start, end) = if line.len() > MAX_COLS {
                    let pad = MAX_COLS.saturating_sub(*length) / 2;
                    let mut start = offset.saturating_sub(pad).min(line.len());
                    while !line.is_char_boundary(start) {
                        start -= 1;
                    }
                    let mut end = (start + MAX_COLS).min(line.len());
                    while !line.is_char_boundary(end) {
                        end += 1;
                    }
                    (start, end)
                } else {
                    (0, line.len())
                };
                write!(
                    f,
                    "\n{:pad$} ╷\n{:<pad$} │ {}{}{}\n{:pad$} · {}{}",
                    "",
                    line_index.map_or(String::new(), |n| (n + 1).to_string()),
                    if start == 0 { "" } else { "…" },
                    &line[start..end],
                    if end == line.len() { "" } else { "…" },
                    "",
                    " ".repeat(*offset - start + usize::from(start != 0)),
                    if *length == 0 {
                        "└".to_string()
                    } else {
                        "‾".repeat(*length)
                    },
                    pad = margin
                )?;
            }
            Self::Range {
                start_line_index,
                lines,
                offset,
            } => {
                write!(f, "\n{:pad$} ╷", "", pad = margin)?;
                let mut number = *start_line_index + 1;
                write!(
                    f,
                    "\n{:<pad$} │ {}{}",
                    number,
                    " ".repeat(*offset),
                    lines[0],
                    pad = margin
                )?;
                for line in lines.iter().skip(1) {
                    number += 1;
                    write!(f, "\n{number:<margin$} │ {line}")?;
                }
            }
            Self::RangeHighlights {
                start_line_index: start_linenumber,
                lines,
                highlights,
            } => {
                write!(f, "\n{:pad$} ╷", "", pad = margin)?;
                let mut number = *start_linenumber + 1;
                let mut highlights_peek = highlights.iter().peekable();

                for (index, line) in lines.iter().enumerate() {
                    number += 1;
                    write!(f, "\n{number:<margin$} │ {line}")?;
                    let mut first = true;
                    let mut last_offset = 0;
                    while let Some(high) = highlights_peek.peek() {
                        if high.0 > index {
                            break;
                        }
                        if let Some(high) = highlights_peek.next() {
                            if first {
                                write!(f, "\n{:pad$} · ", "", pad = margin)?;
                                first = false;
                            }
                            if last_offset < high.1 {
                                write!(
                                    f,
                                    "{}{}",
                                    " ".repeat(high.1 - last_offset),
                                    "‾".repeat(high.2)
                                )?;
                                last_offset = high.1 + high.2;
                            } else {
                                eprintln!("A highlight in a range error message is detected to overlap with a previous highlight, it is skipped.");
                                // Panicking on error gave the following very intense error message (in test code):
                                // `thread panicked while panicking. aborting. ... (exit code: 0xc000001d, STATUS_ILLEGAL_INSTRUCTION)`
                                // To prevent other people from panicking upon seeing this error message this error is not raised currently.
                            }
                        }
                    }
                }
            }
            Self::Multiple { contexts } => {
                for (note, context) in contexts {
                    context.display(f, note.as_deref())?;
                }
                tail = false;
            }
        }
        // Last line
        if let Some(note) = note {
            write!(f, "\n{:pad$} ╰{}", "", note, pad = margin)
        } else if tail {
            write!(f, "\n{:pad$} ╵", "", pad = margin)
        } else {
            Ok(())
        }
    }
}

impl fmt::Display for Context {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        self.display(f, None)
    }
}

#[derive(Debug, Eq, PartialEq, Copy, Clone)]
/// A position in a file for use in parsing/lexing
pub struct FilePosition<'a> {
    /// The remaining text (as ref so no copies)
    pub text: &'a str,
    /// The current line index
    pub line_index: usize,
    /// The current column number
    pub column: usize,
}
