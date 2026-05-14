use std::ops::RangeBounds;

use context_error::{BasicKind, BoxedError, Context, CreateError};

use crate::{
    chemistry::ELEMENT_PARSE_LIST,
    helper_functions::{RangeExtension, explain_number_error, next_number},
    prelude::Element,
};

fn tokenise<'a>(
    value: &'a str,
    base_context: &Context<'a>,
    range: impl RangeBounds<usize>,
) -> Result<Wurcs, BoxedError<'a, BasicKind>> {
    let (mut index, end) = range.bounds(value.len());
    if value[index..end].starts_with("WURCS=2.0/") {
        index += 10;
        if !value[index..end].is_ascii() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "Defintition contains non-ascii characters",
                base_context.clone().add_highlight((0, index..end)),
            ));
        }

        let Some((counts_str, _)) = value[index..end].split_once('/') else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A WURCS 2.0 string should start with the counts followed by a slash",
                base_context.clone().add_highlight((0, index..end)),
            ));
        };
        let counts = tokenise_counts(counts_str, base_context, index, end)?;
        let mut unique_residues = Vec::with_capacity(counts.0 as usize);

        index += counts_str.len() + 1;
        while value[index..end].starts_with('[') {
            let (len, unique_res) = tokenise_unique_res(value, base_context, index, end)?;
            index += len;
            unique_residues.push(unique_res);
        }

        if !value[index..end].starts_with('/') {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "The separator '/' between the unique residues and the residue counts is missing",
                base_context.clone().add_highlight((0, index, 1)),
            ));
        }
        index += 1;

        let Some((residues_str, _)) = value[index..end].split_once('/') else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "The separator '/' after the residue counts is missing",
                base_context.clone().add_highlight((0, index..end)),
            ));
        };

        let mut residue_counts = Vec::with_capacity(unique_residues.len());
        let mut offset = 0;
        for res in residues_str.split('-') {
            let num = res.parse().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    format!("The res count {}", explain_number_error(&err)),
                    base_context
                        .clone()
                        .add_highlight((0, index + offset, res.len())),
                )
            })?;
            offset += res.len() + 1;
            residue_counts.push(num);
        }

        index += residues_str.len() + 1;

        let mut linkage = Vec::with_capacity(counts.2 as usize);
        let mut offset = 0;
        while index + offset < end {
            let mut glips = Vec::new();
            // Parse GLIP
            loop {
                let mut res_index = 0;
                let mut len = 0;
                for c in value[index + offset..end].as_bytes() {
                    if c.is_ascii_lowercase() {
                        len += 1;
                        res_index = res_index * 52 + (c - b'a');
                    } else if c.is_ascii_uppercase() {
                        len += 1;
                        res_index = res_index * 52 + (c - b'A' + 26);
                    } else {
                        break;
                    }
                }
                if len == 0 {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        "Missing residue index for GLIP",
                        base_context.clone().add_highlight((0, index + offset, 1)),
                    ));
                }
                offset += len;
                let (len, position) =
                    maybe_number(value, index + offset..end, base_context, "GLIP Position")?;
                offset += len;
                let extra = if let Some(direction) = value
                    .as_bytes()
                    .get(index + offset)
                    .and_then(|v| (*v).try_into().ok())
                {
                    offset += 1;
                    let (len, num) = next_number::<false, false, u8>(value, index + offset..end)
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The star index is missing"),
                                base_context.clone().add_highlight((0, index + offset, 1)),
                            )
                        })
                        .map(|(len, _, num)| {
                            num.map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid WURCS 2.0",
                                    format!("The star index {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((0, index + offset, len)),
                                )
                            })
                            .map(|n| (len, n))
                        })
                        .flatten()?;
                    offset += len;
                    Some((direction, num))
                } else {
                    None
                };
                glips.push(GLIP {
                    res_index,
                    position,
                    direction: extra.clone().map_or(Direction::Obvious, |(d, _)| d),
                    star_index: extra.map_or(0, |(_, i)| i),
                });
                if value.as_bytes().get(index + offset) == Some(&b'-') {
                    offset += 1;
                    continue;
                } else {
                    break;
                }
            }
            // Parse mod
            let (len, modification) = tokenise_map(value, index + offset..end, base_context)?;
            offset += len;
            linkage.push(Linkage {
                lips: glips,
                modification,
            });
            if let Some(ch) = value.as_bytes().get(index + offset) {
                if *ch == b'_' {
                    offset += 1;
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        "Invalid character after MAP",
                        base_context.clone().add_highlight((0, index + offset, 1)),
                    ));
                }
            }
        }

        // TODO: think about if it is needed to verify the counts

        Ok(Wurcs {
            residues: unique_residues,
            sequence: residue_counts,
            linkage,
        })
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "A WURCS 2.0 string should start with 'WURCS=2.0/'",
            base_context.clone().add_highlight((0, index..end)),
        ))
    }
}

fn tokenise_counts<'a>(
    value: &'a str,
    base_context: &Context<'a>,
    index: usize,
    end: usize,
) -> Result<(u8, u8, u8, bool), BoxedError<'a, BasicKind>> {
    let mut split = value.split(',');
    let mut offset = 0;

    let unique_res_count = split
        .next()
        .ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A WURCS 2.0 unique count should have three numbers",
                base_context.clone().add_highlight((0, index..end)),
            )
        })
        .map(|v| {
            offset += v.len() + 1;
            v.parse::<u8>().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    format!("The unique res count {}", explain_number_error(&err)),
                    base_context
                        .clone()
                        .add_highlight((0, index + offset - v.len() - 1, v.len())),
                )
            })
        })
        .flatten()?;
    let res_count = split
        .next()
        .ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A WURCS 2.0 unique count should have three numbers",
                base_context.clone().add_highlight((0, index..end)),
            )
        })
        .map(|v| {
            offset += v.len() + 1;
            v.parse::<u8>().map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    format!("The res count {}", explain_number_error(&err)),
                    base_context
                        .clone()
                        .add_highlight((0, index + offset - v.len() - 1, v.len())),
                )
            })
        })
        .flatten()?;
    let (lin_count, uncertain) = split
        .next()
        .ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A WURCS 2.0 unique count should have three numbers",
                base_context.clone().add_highlight((0, index..end)),
            )
        })
        .map(|v| {
            offset += v.len() + 1;
            let (uncertain, number) = if v.ends_with('+') {
                (true, &v[..v.len() - 1])
            } else {
                (false, v)
            };
            number
                .parse::<u8>()
                .map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        format!("The lin count {}", explain_number_error(&err)),
                        base_context.clone().add_highlight((
                            0,
                            index + offset - v.len() - 1 - usize::from(uncertain),
                            number.len(),
                        )),
                    )
                })
                .map(|n| (n, uncertain))
        })
        .flatten()?;

    if split.next().is_none() {
        Ok((unique_res_count, res_count, lin_count, uncertain))
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "There are too many counts provided",
            base_context.clone().add_highlight((0, index + offset..end)),
        ))
    }
}

/// Assumes to start at '['
fn tokenise_unique_res<'a>(
    value: &'a str,
    base_context: &Context<'a>,
    index: usize,
    end: usize,
) -> Result<(usize, Residue), BoxedError<'a, BasicKind>> {
    // Parse the skeleton
    let skeleton_length = value[index + 1..end]
        .as_bytes()
        .iter()
        .take_while(|c| **c != b'-' && **c != b'_' && **c != b']')
        .count();
    if skeleton_length < 2 {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "Skeleton code has to be at least two characters",
            base_context
                .clone()
                .add_highlight((0, index + 1, skeleton_length)),
        ));
    }

    let start: TerminalCarbon = value.as_bytes()[index + 1].try_into().map_err(|()| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "Invalid terminal carbon symbol at the start of the skeleton code",
            base_context.clone().add_highlight((0, index + 1, 1)),
        )
    })?;

    let backbone = value.as_bytes()[index + 2..index + skeleton_length]
        .into_iter()
        .enumerate()
        .map(|(i, s)| {
            (*s).try_into().map_err(|()| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    "Invalid carbon symbol in skeleton code",
                    base_context.clone().add_highlight((0, index + 2 + i, 1)),
                )
            })
        })
        .collect::<Result<Vec<Carbon>, _>>()?;
    let skeleton_end: TerminalCarbon = value.as_bytes()[index + skeleton_length]
        .try_into()
        .map_err(|()| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "Invalid terminal carbon symbol at the end of the skeleton code",
                base_context
                    .clone()
                    .add_highlight((0, index + skeleton_length, 1)),
            )
        })?;

    // Parse the anomeric info
    let mut offset = skeleton_length + 1;
    let anomeric = if value.as_bytes()[index + offset] == b'-' {
        offset += 1;
        // If ? then unknown
        let (len, num) = maybe_number(
            value,
            index + offset..end,
            base_context,
            "anomeric location",
        )?;
        offset += len;
        let center = value.as_bytes()[index + offset].try_into().map_err(|()| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "Invalid anomeric symbol in unique residue",
                base_context.clone().add_highlight((0, index + offset, 1)),
            )
        })?;
        offset += 1;
        Some((num, center))
    } else {
        None
    };

    // Pare any following mods (separated by '_')
    let mut mods = Vec::new();
    while value.as_bytes()[index + offset] == b'_' {
        offset += 1;
        let mut lips = Vec::new();
        // Parse LIP
        loop {
            let (len, position) =
                maybe_number(value, index + offset..end, base_context, "LIP Position")?;
            offset += len;
            let extra = if let Ok(direction) = value.as_bytes()[index + offset].try_into() {
                offset += 1;
                let (len, num) = next_number::<false, false, u8>(value, index + offset..end)
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid WURCS 2.0",
                            format!("The star index is missing"),
                            base_context.clone().add_highlight((0, index + offset, 1)),
                        )
                    })
                    .map(|(len, _, num)| {
                        num.map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The star index {}", explain_number_error(&err)),
                                base_context.clone().add_highlight((0, index + offset, len)),
                            )
                        })
                        .map(|n| (len, n))
                    })
                    .flatten()?;
                offset += len;
                Some((direction, num))
            } else {
                None
            };
            lips.push(LIP {
                position,
                direction: extra.clone().map_or(Direction::Obvious, |(d, _)| d),
                star_index: extra.map_or(0, |(_, i)| i),
            });
            if value.as_bytes()[index + offset] == b'-' {
                offset += 1;
                continue;
            } else {
                break;
            }
        }
        // Parse mod
        let (len, modification) = tokenise_map(value, index + offset..end, base_context)?;
        offset += len;
        mods.push(Mod { lips, modification });
    }

    // End at the closing ']'
    if value.as_bytes()[index + offset] == b']' {
        Ok((
            offset + 1,
            Residue {
                start,
                skeleton: backbone,
                end: skeleton_end,
                anomeric,
                mods,
            },
        ))
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "Unique residue was not closed, it should be closed with ']'",
            base_context.clone().add_highlight((0, index + offset, 1)),
        ))
    }
}

fn maybe_number<'a>(
    value: &'a str,
    range: std::ops::Range<usize>,
    base_context: &Context<'a>,
    number: &'static str,
) -> Result<(usize, Option<u8>), BoxedError<'a, BasicKind>> {
    if value[range.clone()].starts_with('?') {
        Ok((1, None))
    } else {
        next_number::<false, false, u8>(value, range.clone())
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    format!("The {number} is missing"),
                    base_context
                        .clone()
                        .add_highlight((0, range.clone().start, 1)),
                )
            })
            .map(|(len, _, num)| {
                num.map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        format!("The {number} {}", explain_number_error(&err)),
                        base_context
                            .clone()
                            .add_highlight((0, range.clone().start, len)),
                    )
                })
                .map(|n| (len, Some(n)))
            })
            .flatten()
    }
}

fn tokenise_map<'a>(
    value: &str,
    range: std::ops::Range<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, Vec<MAPSymbol>), BoxedError<'a, BasicKind>> {
    let mut tokens = Vec::new();
    let mut offset = 0;
    'outer: loop {
        if range.start + offset >= range.end {
            return Ok((offset, tokens));
        }
        for (name, el) in ELEMENT_PARSE_LIST {
            if value[range.start + offset..range.end].starts_with(name) {
                offset += name.len();
                tokens.push(MAPSymbol::Element(*el));
                continue 'outer;
            }
        }
        match value.as_bytes()[range.start + offset] {
            b'*' => {
                offset += 1;
                if let Some((len, index)) =
                    next_number::<false, false, u8>(value, range.start + offset..range.end)
                        .map(|(len, _, num)| {
                            num.map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid WURCS 2.0",
                                    format!("The star index {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((
                                        0,
                                        range.start + offset,
                                        len,
                                    )),
                                )
                            })
                            .map(|n| (len, n))
                        })
                        .transpose()?
                {
                    offset += len;
                    tokens.push(MAPSymbol::Star(Some(index)));
                } else {
                    tokens.push(MAPSymbol::Star(None));
                }
            }
            b'/' => {
                offset += 1;
                let (len, num) =
                    next_number::<false, false, u8>(value, range.start + offset..range.end)
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The branch index is missing"),
                                base_context
                                    .clone()
                                    .add_highlight((0, range.start + offset, 1)),
                            )
                        })
                        .map(|(len, _, num)| {
                            num.map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid WURCS 2.0",
                                    format!("The branch index {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((
                                        0,
                                        range.start + offset,
                                        len,
                                    )),
                                )
                            })
                            .map(|n| (len, n))
                        })
                        .flatten()?;
                offset += len;
                tokens.push(MAPSymbol::Branch(num));
            }
            b'$' => {
                offset += 1;
                let (len, num) =
                    next_number::<false, false, u8>(value, range.start + offset..range.end)
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The cyclic index is missing"),
                                base_context
                                    .clone()
                                    .add_highlight((0, range.start + offset, 1)),
                            )
                        })
                        .map(|(len, _, num)| {
                            num.map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid WURCS 2.0",
                                    format!("The cyclic index {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((
                                        0,
                                        range.start + offset,
                                        len,
                                    )),
                                )
                            })
                            .map(|n| (len, n))
                        })
                        .flatten()?;
                offset += len;
                tokens.push(MAPSymbol::Cyclic(num));
            }
            b'^' => {
                offset += 1;
                tokens.push(MAPSymbol::Chirality(
                    match value.as_bytes()[range.start + offset] {
                        b'R' => Chirality::R,
                        b'S' => Chirality::S,
                        b'X' => Chirality::Unknown,
                        b'E' => Chirality::E,
                        b'Z' => Chirality::Z,
                        _ => {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                "The chirality or gemoetrical isomerism indicator is invalid",
                                base_context
                                    .clone()
                                    .add_highlight((0, range.start + offset, 1)),
                            ));
                        }
                    },
                ));
            }
            b'=' => {
                offset += 1;
                tokens.push(MAPSymbol::DoubleBond);
            }
            b'#' => {
                offset += 1;
                tokens.push(MAPSymbol::TripleBond);
            }
            b'(' => {
                offset += 1;
                tokens.push(MAPSymbol::AromaticStart);
            }
            b')' => {
                offset += 1;
                tokens.push(MAPSymbol::AromaticEnd);
            }
            _ => return Ok((offset, tokens)),
        }
    }
}

#[derive(Debug)]
struct Wurcs {
    residues: Vec<Residue>,
    sequence: Vec<u8>,
    linkage: Vec<Linkage>,
}

#[derive(Debug)]
struct Residue {
    start: TerminalCarbon,
    skeleton: Vec<Carbon>,
    end: TerminalCarbon,
    anomeric: Option<(Option<u8>, AnomericSymbol)>,
    mods: Vec<Mod>,
}

#[derive(Debug, Copy, Clone)]
enum Carbon {
    /// 'O'
    Ketone,
    /// '1'
    HydroxyLeft,
    /// '2'
    HydroxyRight,
    /// 'd' deoxy
    Methylene,
    /// 'a' anomeric
    Hemiketal,
    /// 'U'
    KetoneOrHemiketal,
}

impl TryFrom<u8> for Carbon {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'O' => Ok(Self::Ketone),
            b'1' => Ok(Self::HydroxyLeft),
            b'2' => Ok(Self::HydroxyRight),
            b'd' => Ok(Self::Methylene),
            b'a' => Ok(Self::Hemiketal),
            b'U' => Ok(Self::KetoneOrHemiketal),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Copy, Clone)]
enum TerminalCarbon {
    /// 'o'
    Aldehyde,
    /// 'h'
    Hydroxymethyl,
    /// 'A'
    Carboxyl,
    /// 'm'
    Methyl,
    /// 'a'
    Hemiacetal,
    /// 'u'
    AldehydeOrHemiacetal,
}

impl TryFrom<u8> for TerminalCarbon {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'o' => Ok(Self::Aldehyde),
            b'h' => Ok(Self::Hydroxymethyl),
            b'A' => Ok(Self::Carboxyl),
            b'm' => Ok(Self::Methyl),
            b'a' => Ok(Self::Hemiacetal),
            b'u' => Ok(Self::AldehydeOrHemiacetal),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Copy, Clone)]
enum AnomericSymbol {
    /// 'a'
    Alpha,
    /// 'b'
    Beta,
    /// 'u'
    Up,
    /// 'd'
    Down,
    /// 'x'
    Unknown,
    /// 'o'
    None,
}

impl TryFrom<u8> for AnomericSymbol {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'a' => Ok(Self::Alpha),
            b'b' => Ok(Self::Beta),
            b'u' => Ok(Self::Up),
            b'd' => Ok(Self::Down),
            b'x' => Ok(Self::Unknown),
            b'o' => Ok(Self::None),
            _ => Err(()),
        }
    }
}

#[derive(Debug)]
struct Mod {
    lips: Vec<LIP>,
    modification: Vec<MAPSymbol>,
}

#[derive(Debug, Copy, Clone)]
struct LIP {
    position: Option<u8>,
    direction: Direction,
    star_index: u8,
}

#[derive(Debug, Copy, Clone)]
enum Direction {
    /// 'u'
    Upside,
    /// 'd'
    Downside,
    /// 't'
    Tres,
    /// 'a'
    Same,
    /// 'b'
    Opposite,
    /// 'c'
    Third,
    /// 'x'
    Unknown,
    /// 'e'
    Entgegen,
    /// 'z'
    Zusammen,
    /// 'f'
    UnknownGeometricalIsomerism,
    /// 'n'
    Obvious,
}

impl TryFrom<u8> for Direction {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'u' => Ok(Self::Upside),
            b'd' => Ok(Self::Downside),
            b't' => Ok(Self::Tres),
            b'a' => Ok(Self::Same),
            b'b' => Ok(Self::Opposite),
            b'c' => Ok(Self::Third),
            b'x' => Ok(Self::Unknown),
            b'e' => Ok(Self::Entgegen),
            b'z' => Ok(Self::Zusammen),
            b'f' => Ok(Self::UnknownGeometricalIsomerism),
            b'n' => Ok(Self::Obvious),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Copy, Clone)]
enum MAPSymbol {
    Element(Element),
    /// '*' or '*n'
    Star(Option<u8>),
    /// '/n'
    Branch(u8),
    /// '$n'
    Cyclic(u8),
    /// '='
    DoubleBond,
    /// '#'
    TripleBond,
    Chirality(Chirality),
    AromaticStart,
    AromaticEnd,
}

#[derive(Debug, Copy, Clone)]
enum Chirality {
    R,
    S,
    E,
    Z,
    Unknown,
}

#[derive(Debug)]
// TODO: misses ambiguous linkage
struct Linkage {
    lips: Vec<GLIP>,
    modification: Vec<MAPSymbol>,
}

#[derive(Debug, Copy, Clone)]
struct GLIP {
    /// Encoded using base 52
    res_index: u8,
    position: Option<u8>,
    direction: Direction,
    star_index: u8,
}

#[cfg(test)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::glycan::wurcs::{Wurcs, tokenise};

    fn test_tokenise(value: &str) -> Result<Wurcs, BoxedError<'_, BasicKind>> {
        tokenise(value, &Context::default().lines(0, value), ..)
    }

    #[test]
    fn tokenise1() {
        let t = test_tokenise("WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/");
        println!("{t:?}");
        assert!(t.is_ok());
        let t = test_tokenise("WURCS=2.0/1,1,0/[Aad22112h-2a_2-6_5*NCC/3=O]/1/");
        println!("{t:?}");
        assert!(t.is_ok());
        let t = test_tokenise(
            "WURCS=2.0/3,5,4/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3/a4-b1_b4-c1_c3-d1_c6-e1",
        );
        println!("{t:?}");
        assert!(t.is_ok());
        todo!()
    }
}
