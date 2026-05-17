#![allow(dead_code)]
use std::ops::RangeBounds;

use context_error::{BasicKind, BoxedError, Context, CreateError};

use crate::{
    chemistry::ELEMENT_PARSE_LIST,
    helper_functions::{RangeExtension, explain_number_error, next_number},
    prelude::Element,
};

pub fn tokenise_wurcs<'a>(
    value: &'a str,
    base_context: &Context<'a>,
    range: impl RangeBounds<usize>,
) -> Result<Wurcs, BoxedError<'a, BasicKind>> {
    let (mut index, end) = range.bounds(value.len() - 1);
    if value[index..=end].starts_with("WURCS=2.0/") {
        index += 10;
        if !value[index..=end].is_ascii() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "Defintition contains non-ascii characters",
                base_context.clone().add_highlight((0, index..=end)),
            ));
        }

        let Some((counts_str, _)) = value[index..=end].split_once('/') else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A WURCS 2.0 string should start with the counts followed by a slash",
                base_context.clone().add_highlight((0, index..=end)),
            ));
        };
        let counts = tokenise_counts(counts_str, base_context, index, end)?;
        let mut unique_residues = Vec::with_capacity(counts.0 as usize);

        index += counts_str.len() + 1;
        while value[index..=end].starts_with('[') {
            let (len, unique_res) = tokenise_unique_res(value, base_context, index, end)?;
            index += len;
            unique_residues.push(unique_res);
        }

        if !value[index..=end].starts_with('/') {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "The separator '/' between the unique residues and the residue counts is missing",
                base_context.clone().add_highlight((0, index, 1)),
            ));
        }
        index += 1;

        let Some((residues_str, _)) = value[index..=end].split_once('/') else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "The separator '/' after the residue counts is missing",
                base_context.clone().add_highlight((0, index..=end)),
            ));
        };

        let mut residue_sequence = Vec::with_capacity(unique_residues.len());
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
            if num as usize > unique_residues.len() || num == 0 {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid WURCS 2.0",
                    "The sequence refers to a nonexisting unique residue",
                    base_context
                        .clone()
                        .add_highlight((0, index + offset, res.len())),
                ));
            }
            offset += res.len() + 1;
            residue_sequence.push(num);
        }

        index += residues_str.len() + 1;

        let mut linkage = Vec::with_capacity(counts.2 as usize);
        let mut offset = 0;
        while index + offset <= end {
            let mut glips = Vec::new();
            // Parse GLIP
            loop {
                let mut alternate = Vec::new();
                let res_alt_open = if value.as_bytes().get(index + offset) == Some(&b'{') {
                    offset += 1;
                    true
                } else {
                    false
                };
                // Check if %
                let backbone_prob = if value.as_bytes().get(index + offset) == Some(&b'%') {
                    let (len, probability) =
                        parse_probability(value, index + offset + 1..=end, base_context)?;
                    offset += len + 1;
                    Some(probability)
                } else {
                    None
                };
                // Parse glip
                let (len, glip) = parse_glip(value, index + offset..=end, base_context)?;
                alternate.push(glip);
                offset += len;
                // Check if % (if not already parsed?)
                let opposite_prob = if value.as_bytes().get(index + offset) == Some(&b'%') {
                    let (len, probability) =
                        parse_probability(value, index + offset + 1..=end, base_context)?;
                    offset += len + 1;
                    Some(probability)
                } else {
                    None
                };
                // Check if alternate
                while index + offset <= end && value.as_bytes().get(index + offset) == Some(&b'|') {
                    offset += 1;
                    let (len, glip) = parse_glip(value, index + offset..=end, base_context)?;
                    offset += len;
                    alternate.push(glip);
                }
                let res_alt_close = if value.as_bytes().get(index + offset) == Some(&b'}') {
                    offset += 1;
                    true
                } else {
                    false
                };
                if alternate.len() == 1 {
                    if res_alt_open || res_alt_close {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid WURCS 2.0",
                            "A residue alternate GLIP must contain more than one GLIP",
                            base_context
                                .clone()
                                .add_highlight((0, index..=index + offset)),
                        ));
                    }
                    glips.push(match (backbone_prob, opposite_prob) {
                        (Some(_), Some(_)) => {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                "A probability range can only be provided for the backbone or opposite side, not both at the same time",
                                base_context
                                    .clone()
                                    .add_highlight((0, index..=index + offset)),
                            ));
                        }
                        (Some(prob), None) => GLIPOption::Statistic(true, prob, alternate[0]),
                        (None, Some(prob)) => GLIPOption::Statistic(false, prob, alternate[0]),
                        (None, None) => GLIPOption::Known(alternate[0]),
                    });
                } else {
                    if backbone_prob.is_some() || opposite_prob.is_some() {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid WURCS 2.0",
                            "A GLIP cannot be both statistic and alternate at the same time",
                            base_context
                                .clone()
                                .add_highlight((0, index..=index + offset)),
                        ));
                    }
                    glips.push(match (res_alt_open, res_alt_close) {
                        (true, true) => {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                "A residue alternate GLIP can be can only use curly braces at one side, not both at the same time",
                                base_context
                                    .clone()
                                    .add_highlight((0, index..=index + offset)),
                            ));
                        }
                        (true, false) => GLIPOption::RESAlternative(true, alternate),
                        (false, true) => GLIPOption::RESAlternative(false, alternate),
                        (false, false) => GLIPOption::Alternative(alternate),
                    });
                }
                if index + offset <= end && value.as_bytes().get(index + offset) == Some(&b'-') {
                    offset += 1;
                    continue;
                } else {
                    break;
                }
            }
            // Parse mod
            let (len, modification) = tokenise_map(value, index + offset..=end, base_context)?;
            offset += len;
            if index + offset <= end && value.as_bytes().get(index + offset) == Some(&b'~') {
                offset += 1;
                let (len, repeat) = parse_repeat(value, index + offset..=end, base_context)?;
                offset += len;
                linkage.push(Linkage::Repeated(
                    repeat,
                    LIN {
                        lips: glips,
                        modification,
                    },
                ));
            } else {
                linkage.push(Linkage::Known(LIN {
                    lips: glips,
                    modification,
                }));
            }
            if index + offset <= end
                && let Some(ch) = value.as_bytes().get(index + offset)
            {
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
            sequence: residue_sequence,
            linkage,
        })
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "A WURCS 2.0 string should start with 'WURCS=2.0/'",
            base_context.clone().add_highlight((0, index..=end)),
        ))
    }
}

fn parse_glip<'a>(
    value: &'a str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, GLIP), BoxedError<'a, BasicKind>> {
    let mut offset = 0;
    let mut res_index = 0;
    let mut len = 0;
    for c in value[range.start() + offset..=*range.end()].as_bytes() {
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
            base_context
                .clone()
                .add_highlight((0, range.start() + offset, 1)),
        ));
    }
    offset += len;
    let (len, position) = maybe_number(
        value,
        range.start() + offset..=*range.end(),
        base_context,
        "GLIP Position",
        '?',
    )?;
    offset += len;
    let extra = if let Some(direction) = value
        .as_bytes()
        .get(range.start() + offset)
        .and_then(|v| (*v).try_into().ok())
    {
        offset += 1;
        let (len, num) =
            next_number::<false, false, u8>(value, range.start() + offset..=*range.end())
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        format!("The star index is missing"),
                        base_context
                            .clone()
                            .add_highlight((0, range.start() + offset, 1)),
                    )
                })
                .map(|(len, _, num)| {
                    num.map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid WURCS 2.0",
                            format!("The star index {}", explain_number_error(&err)),
                            base_context
                                .clone()
                                .add_highlight((0, range.start() + offset, len)),
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
    Ok((
        offset,
        GLIP {
            res_index,
            position,
            direction: extra.clone().map_or(Direction::Obvious, |(d, _)| d),
            star_index: extra.map_or(0, |(_, i)| i),
        },
    ))
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
                base_context.clone().add_highlight((0, index..=end)),
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
                base_context.clone().add_highlight((0, index..=end)),
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
                base_context.clone().add_highlight((0, index..=end)),
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
            base_context
                .clone()
                .add_highlight((0, index + offset..=end)),
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
    let skeleton_length = value[index + 1..=end]
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

    // What does 'x' mean? Maybe just unknown? Also handle '<x>' as unknown length.
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
            index + offset..=end,
            base_context,
            "anomeric location",
            '?',
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
            let mut alternate = Vec::new();
            // Check if %
            let backbone_prob = if value.as_bytes().get(index + offset) == Some(&b'%') {
                let (len, probability) =
                    parse_probability(value, index + offset + 1..=end, base_context)?;
                offset += len + 1;
                Some(probability)
            } else {
                None
            };
            // Parse lip
            let (len, lip) = parse_lip(value, index + offset..=end, base_context)?;
            alternate.push(lip);
            offset += len;
            // Check if % (if not already parsed?)
            let opposite_prob = if value.as_bytes().get(index + offset) == Some(&b'%') {
                let (len, probability) =
                    parse_probability(value, index + offset + 1..=end, base_context)?;
                offset += len + 1;
                Some(probability)
            } else {
                None
            };
            // Check if alternate
            while index + offset <= end && value.as_bytes().get(index + offset) == Some(&b'|') {
                offset += 1;
                let (len, lip) = parse_lip(value, index + offset..=end, base_context)?;
                offset += len;
                alternate.push(lip);
            }
            if alternate.len() == 1 {
                lips.push(match (backbone_prob, opposite_prob) {
                        (Some(_), Some(_)) => {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                "A probability range can only be provided for the backbone or opposite side, not both at the same time",
                                base_context
                                    .clone()
                                    .add_highlight((0, index..=index + offset)),
                            ));
                        }
                        (Some(prob), None) => LIPOption::Statistic(true, prob, alternate[0]),
                        (None, Some(prob)) => LIPOption::Statistic(false, prob, alternate[0]),
                        (None, None) => LIPOption::Known(alternate[0]),
                    });
            } else {
                if backbone_prob.is_some() || opposite_prob.is_some() {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        "A LIP cannot be both statistic and alternate at the same time",
                        base_context
                            .clone()
                            .add_highlight((0, index..=index + offset)),
                    ));
                }
                lips.push(LIPOption::Alternative(alternate));
            }
            if index + offset <= end && value.as_bytes().get(index + offset) == Some(&b'-') {
                offset += 1;
                continue;
            } else {
                break;
            }
        }
        // Parse mod
        let (len, modification) = tokenise_map(value, index + offset..=end, base_context)?;
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

fn parse_lip<'a>(
    value: &'a str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, LIP), BoxedError<'a, BasicKind>> {
    let mut offset = 0;
    let (len, position) = maybe_number(
        value,
        range.start() + offset..=*range.end(),
        base_context,
        "GLIP Position",
        '?',
    )?;
    offset += len;
    let extra = if let Some(direction) = value
        .as_bytes()
        .get(range.start() + offset)
        .and_then(|v| (*v).try_into().ok())
    {
        offset += 1;
        let (len, num) =
            next_number::<false, false, u8>(value, range.start() + offset..=*range.end())
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        format!("The star index is missing"),
                        base_context
                            .clone()
                            .add_highlight((0, range.start() + offset, 1)),
                    )
                })
                .map(|(len, _, num)| {
                    num.map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid WURCS 2.0",
                            format!("The star index {}", explain_number_error(&err)),
                            base_context
                                .clone()
                                .add_highlight((0, range.start() + offset, len)),
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
    Ok((
        offset,
        LIP {
            position,
            direction: extra.clone().map_or(Direction::Obvious, |(d, _)| d),
            star_index: extra.map_or(0, |(_, i)| i),
        },
    ))
}

fn maybe_number<'a>(
    value: &str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
    number: &'static str,
    unknown_symbol: char,
) -> Result<(usize, Option<u8>), BoxedError<'a, BasicKind>> {
    if value[range.clone()].starts_with(unknown_symbol) {
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
                        .add_highlight((0, *range.clone().start(), 1)),
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
                            .add_highlight((0, *range.clone().start(), len)),
                    )
                })
                .map(|n| (len, Some(n)))
            })
            .flatten()
    }
}

fn tokenise_map<'a>(
    value: &str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, Vec<MAPSymbol>), BoxedError<'a, BasicKind>> {
    let mut tokens = Vec::new();
    let mut offset = 0;
    'outer: loop {
        if range.start() + offset >= *range.end() {
            return Ok((offset, tokens));
        }
        for (name, el) in ELEMENT_PARSE_LIST {
            if value[range.start() + offset..=*range.end()].starts_with(name) {
                offset += name.len();
                tokens.push(MAPSymbol::Element(*el));
                continue 'outer;
            }
        }
        match value.as_bytes()[range.start() + offset] {
            b'*' => {
                offset += 1;
                if let Some((len, index)) =
                    next_number::<false, false, u8>(value, range.start() + offset..=*range.end())
                        .map(|(len, _, num)| {
                            num.map_err(|err| {
                                BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid WURCS 2.0",
                                    format!("The star index {}", explain_number_error(&err)),
                                    base_context.clone().add_highlight((
                                        0,
                                        range.start() + offset,
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
                    next_number::<false, false, u8>(value, range.start() + offset..*range.end())
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The branch index is missing"),
                                base_context
                                    .clone()
                                    .add_highlight((0, range.start() + offset, 1)),
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
                                        range.start() + offset,
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
                    next_number::<false, false, u8>(value, range.start() + offset..*range.end())
                        .ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid WURCS 2.0",
                                format!("The cyclic index is missing"),
                                base_context
                                    .clone()
                                    .add_highlight((0, range.start() + offset, 1)),
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
                                        range.start() + offset,
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
                offset += 2;
                tokens.push(MAPSymbol::Chirality(
                    match value.as_bytes()[range.start() + offset - 1] {
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
                                    .add_highlight((0, range.start() + offset, 1)),
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

fn parse_probability<'a>(
    value: &str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, Probability), BoxedError<'a, BasicKind>> {
    let parse_num = |offset: usize| {
        match value.as_bytes().get(*range.start() + offset) {
            Some(b'?') => Ok((1, None)),
            Some(b'.') => {
                let mut o = 1;
                let mut len = 0;
                let mut num = 0.0;

                while let Some(ch) = value.as_bytes().get(range.start() + offset + o)
                    && ch.is_ascii_digit()
                {
                    len += 1;
                    o += 1;
                    let v = *ch - b'0';
                    num += (v as f32).powi(-len);
                }

                if len == 0 {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid WURCS 2.0",
                        "A probability range cannot contain an empty number",
                        base_context
                            .clone()
                            .add_highlight((0, *range.start() + offset, 1)),
                    ));
                }

                Ok((o, Some(num)))
            } //parse
            _ => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid WURCS 2.0",
                "A probability range number should start with a '.' or a '?'",
                base_context
                    .clone()
                    .add_highlight((0, *range.start() + offset, 1)),
            )),
        }
    };

    let (offset, first) = parse_num(0)?;
    let (offset, v) = if value.as_bytes().get(range.start() + offset) == Some(&b'-') {
        let (len, second) = parse_num(offset + 1)?;
        (offset + len + 1, Probability::Range(first, second))
    } else {
        (offset, Probability::Single(first))
    };

    if value.as_bytes().get(range.start() + offset) == Some(&b'%') {
        Ok((offset + 1, v))
    } else {
        Err(BoxedError::new(
            BasicKind::Error,
            "Invalid WURCS 2.0",
            "A probability range should end with '%'",
            base_context
                .clone()
                .add_highlight((0, *range.start() + offset, 1)),
        ))
    }
}

fn parse_repeat<'a>(
    value: &str,
    range: std::ops::RangeInclusive<usize>,
    base_context: &Context<'a>,
) -> Result<(usize, Repeat), BoxedError<'a, BasicKind>> {
    let (offset, first) = maybe_number(value, range.clone(), base_context, "repeat start", 'n')?;
    Ok(
        if value.as_bytes().get(range.start() + offset) == Some(&b':') {
            let (len, second) = maybe_number(
                value,
                range.start() + offset + 1..=*range.end(),
                base_context,
                "repeat start",
                'n',
            )?;
            (offset + len + 1, Repeat::Range(first, second))
        } else {
            (offset, Repeat::Single(first))
        },
    )
}

#[derive(Debug)]
pub struct Wurcs {
    residues: Vec<Residue>,
    sequence: Vec<u8>,
    linkage: Vec<Linkage>,
}

#[derive(Debug)]
struct Residue {
    start: TerminalCarbon,
    skeleton: Vec<Carbon>, // TODO: handle the '<x>' case of unknown length all the same carbon
    end: TerminalCarbon,
    anomeric: Option<(Option<u8>, AnomericSymbol)>,
    mods: Vec<Mod>,
}

/// The carbon descriptors extended with the options in: https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.6b00650/suppl_file/ci6b00650_si_001.pdf section 2.8.
#[derive(Debug, Copy, Clone)]
enum Carbon {
    /// 'd' deoxy `H-C-H`
    Methylene,
    /// 'C' `X-C-X`
    Dual,
    /// '1' `X-C-H`
    HydroxyLeft,
    /// '2' `H-C-X`
    HydroxyRight,
    /// '3' `C(X)(H)`
    HydroxyOpposite,
    /// '4' `C(H)(X)`
    HydroxySame,
    /// 'x' one of 1 or 2
    HydroxyUnknown,
    /// '5' `X-C-Y`
    DualLeft,
    /// '6' `Y-C-X`
    DualRight,
    /// '7' `C(X)(Y)`
    DualOpposite,
    /// '8' `C(Y)(X)`
    DualSame,
    /// 'X' one of 5 or 6
    DualUnknown,
    /// 'O' `C=O`
    Ketone,
    /// `e`
    DoubleHydroxyEntgegen,
    /// `z`
    DoubleHydroxyZusammen,
    /// `n`
    DoubleHydroxyNoIsomer,
    /// `f`
    DoubleHydroxyUnknown,
    /// `E`
    DoubleEntgegen,
    /// `Z`
    DoubleZusammen,
    /// `N`
    DoubleNoIsomer,
    /// `F`
    DoubleUnknown,
    /// `K`
    DoubleBonded,
    /// `T`
    TripleBonded,
    /// 'a' anomeric
    Hemiketal,
    /// 'U'
    KetoneOrHemiketal,
    /// 'Q' can be any of the other carbon descriptors
    Unknown,
}

impl TryFrom<u8> for Carbon {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'd' => Ok(Self::Methylene),
            b'C' => Ok(Self::Dual),
            b'1' => Ok(Self::HydroxyLeft),
            b'2' => Ok(Self::HydroxyRight),
            b'3' => Ok(Self::HydroxyOpposite),
            b'4' => Ok(Self::HydroxySame),
            b'x' => Ok(Self::HydroxyUnknown),
            b'5' => Ok(Self::DualLeft),
            b'6' => Ok(Self::DualRight),
            b'7' => Ok(Self::DualOpposite),
            b'8' => Ok(Self::DualSame),
            b'X' => Ok(Self::DualUnknown),
            b'O' => Ok(Self::Ketone),
            b'e' => Ok(Self::DoubleHydroxyEntgegen),
            b'z' => Ok(Self::DoubleHydroxyZusammen),
            b'n' => Ok(Self::DoubleHydroxyNoIsomer),
            b'f' => Ok(Self::DoubleHydroxyUnknown),
            b'E' => Ok(Self::DoubleEntgegen),
            b'Z' => Ok(Self::DoubleZusammen),
            b'N' => Ok(Self::DoubleNoIsomer),
            b'F' => Ok(Self::DoubleUnknown),
            b'K' => Ok(Self::DoubleBonded),
            b'T' => Ok(Self::TripleBonded),
            b'a' => Ok(Self::Hemiketal),
            b'U' => Ok(Self::KetoneOrHemiketal),
            b'Q' => Ok(Self::Unknown),
            _ => Err(()),
        }
    }
}

#[derive(Debug, Copy, Clone)]
#[allow(non_camel_case_types)]
enum TerminalCarbon {
    /// 'm'
    CHHH,
    /// 'M'
    CXXX,
    /// 'h'
    CHHX,
    /// 'c'
    CXXH,
    /// 'C'
    CXXY,
    /// '1'
    CXYH,
    /// '2'
    CYXH,
    /// '3'
    CXYH_opposite,
    /// '4'
    CYXH_same,
    /// 'x'
    CXYH_unknown,
    /// '5'
    CXYZ,
    /// '6'
    CYXZ,
    /// '7'
    CXYZ_opposite,
    /// '8'
    CYXZ_same,
    /// 'X'
    CXYZ_unknown,
    /// 'o' -C=XH
    CXH,
    /// 'A' -C=XY
    CXY,
    /// 'n' =CHH
    CHH,
    /// 'N' =CXX
    CXX,
    /// 'e'
    Double_CXH_entgegen,
    /// 'z'
    Double_CXH_zusammen,
    /// 'f'
    Double_CXH_unknown,
    /// 'E'
    Double_CXY_entgegen,
    /// 'Z'
    Double_CXY_zusammen,
    /// 'F'
    Double_CXY_unknown,
    /// 'T' ≡C-X or -C≡X
    Triple_CX,
    /// 'K'
    Double_CX,
    /// 't'
    Triple_CH,
    /// 'a'
    Hemiacetal,
    /// 'u'
    AldehydeOrHemiacetal,
    /// 'Q'
    Unknown,
}

impl TryFrom<u8> for TerminalCarbon {
    type Error = ();
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'm' => Ok(Self::CHHH),
            b'M' => Ok(Self::CXXX),
            b'h' => Ok(Self::CHHX),
            b'c' => Ok(Self::CXXH),
            b'C' => Ok(Self::CXXY),
            b'1' => Ok(Self::CXYH),
            b'2' => Ok(Self::CYXH),
            b'3' => Ok(Self::CXYH_opposite),
            b'4' => Ok(Self::CYXH_same),
            b'x' => Ok(Self::CXYH_unknown),
            b'5' => Ok(Self::CXYZ),
            b'6' => Ok(Self::CYXZ),
            b'7' => Ok(Self::CXYZ_opposite),
            b'8' => Ok(Self::CYXZ_same),
            b'X' => Ok(Self::CXYZ_unknown),
            b'o' => Ok(Self::CXH),
            b'A' => Ok(Self::CXY),
            b'n' => Ok(Self::CHH),
            b'N' => Ok(Self::CXX),
            b'e' => Ok(Self::Double_CXH_entgegen),
            b'z' => Ok(Self::Double_CXH_zusammen),
            b'f' => Ok(Self::Double_CXH_unknown),
            b'E' => Ok(Self::Double_CXY_entgegen),
            b'Z' => Ok(Self::Double_CXY_zusammen),
            b'F' => Ok(Self::Double_CXY_unknown),
            b'T' => Ok(Self::Triple_CX),
            b'K' => Ok(Self::Double_CX),
            b't' => Ok(Self::Triple_CH),
            b'a' => Ok(Self::Hemiacetal),
            b'u' => Ok(Self::AldehydeOrHemiacetal),
            b'Q' => Ok(Self::Unknown),
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
    lips: Vec<LIPOption>,
    modification: Vec<MAPSymbol>,
}

#[derive(Debug, Clone)]
enum LIPOption {
    Known(LIP),
    Statistic(bool, Probability, LIP),
    Alternative(Vec<LIP>),
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
enum Linkage {
    Known(LIN),
    Repeated(Repeat, LIN),
}

#[derive(Debug)]
// TODO: misses ambiguous linkage
struct LIN {
    lips: Vec<GLIPOption>,
    modification: Vec<MAPSymbol>,
}

#[derive(Debug, Clone)]
enum GLIPOption {
    Known(GLIP),
    Statistic(bool, Probability, GLIP),
    Alternative(Vec<GLIP>),
    RESAlternative(bool, Vec<GLIP>),
}

#[derive(Debug, Copy, Clone)]
struct GLIP {
    /// Encoded using base 52
    res_index: u8,
    /// Position or unknown
    position: Option<u8>,
    direction: Direction,
    star_index: u8,
}

#[derive(Debug, Copy, Clone)]
enum Probability {
    Single(Option<f32>),
    Range(Option<f32>, Option<f32>),
}

#[derive(Debug, Copy, Clone)]
enum Repeat {
    Single(Option<u8>),
    /// Max, min
    Range(Option<u8>, Option<u8>),
}

#[cfg(test)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::glycan::wurcs::{Wurcs, tokenise_wurcs};

    fn test_tokenise(value: &str) -> Result<Wurcs, BoxedError<'_, BasicKind>> {
        tokenise_wurcs(value, &Context::default().lines(0, value), ..)
    }

    macro_rules! test {
        ($case:literal, $name:ident) => {
            #[test]
            #[allow(non_snake_case)]
            fn $name() {
                let t = test_tokenise($case);
                println!("{}\n=>{t:?}", $case);
                assert!(t.is_ok());
            }
        };
    }

    test!(
        "WURCS=2.0/2,2,1/[a2112h-1a_1-5_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1",
        case1
    );
    test!("WURCS=2.0/1,1,0/[a2122h-1b_1-5_2*NCC/3=O]/1/", case2);
    test!("WURCS=2.0/1,1,0/[Aad22112h-2a_2-6_5*NCC/3=O]/1/", case3);
    test!(
        "WURCS=2.0/3,5,4/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5]/1-1-2-3-3/a4-b1_b4-c1_c3-d1_c6-e1",
        case4
    );
    test!(
        "WURCS=2.0/2,2,2/[u2211m][u12h]/1-2/a?|b?}*OC_a?|b?}-{a?|b?",
        case_G45070JE
    );
    test!(
        "WURCS=2.0/3,3,2/[a2112h-1a_1-5_2*NCC/3=O][a2112h-1x_1-5][Aad21122h-2x_2-6_5*NCC/3=O]/1-2-3/a3-b1_c2-a?|b?}",
        case_G46252BF
    );
}
