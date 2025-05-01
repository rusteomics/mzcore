//! WIP: mzPAF parser
#![allow(dead_code)]
use std::{ops::Range, sync::LazyLock};

use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{MolecularCharge, MolecularFormula},
    error::{Context, CustomError},
    fragment::NeutralLoss,
    helper_functions::{Characters, RangeExtension, RangeMaths, explain_number_error, next_number},
    molecular_formula,
    ontology::Ontology,
    quantities::Tolerance,
    sequence::{AminoAcid, SimpleModification},
    system::{MassOverCharge, e, isize::Charge, mz},
};
// TODO: custom errors are off as they assume the input to be Bytes but now get Character offsets.

/// Parse a mzPAF peak annotation line (can contain multiple annotations).
/// # Errors
/// When the annotation does not follow the format.
pub fn parse_mzpaf(line: &str) -> Result<Vec<PeakAnnotation>, CustomError> {
    let mut annotations = Vec::new();

    // Parse first
    let (mut range, a) = parse_annotation(line, 0..line.chars().count())?;
    annotations.push(a);

    // Parse any following
    while !range.is_empty() {
        if line.chars().nth(range.start_index()) == Some(',') {
            range.add_start(1_usize);
        } else {
            return Err(CustomError::error(
                "Invalid mzPAF annotation delimiter",
                "Different mzPAF annotations should be separated with commas ','.",
                Context::line(None, line, range.start_index(), 1),
            ));
        }
        let (r, a) = parse_annotation(line, range)?;
        range = r;
        annotations.push(a);
    }

    Ok(annotations)
}

/// Parse a single mzPAF peak annotation.
/// # Errors
/// When the annotation does not follow the format.
fn parse_annotation(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, PeakAnnotation), CustomError> {
    let (left_range, auxiliary) = if line.chars().nth(range.start_index()) == Some('&') {
        (range.add_start(1_usize), true)
    } else {
        (range.clone(), false)
    };
    let (left_range, analyte_number) = parse_analyte_number(line, left_range)?;
    let (left_range, ion) = parse_ion(line, left_range)?;
    let (left_range, neutral_losses) = parse_neutral_loss(line, left_range)?;
    // TODO: Parse isotopes (currently considered ambiguous with Iodine loss/gain)
    let (left_range, adduct_type) = parse_adduct_type(line, left_range)?;
    let (left_range, charge) = parse_charge(line, left_range)?;
    let (left_range, deviation) = parse_deviation(line, left_range)?;
    let (left_range, confidence) = parse_confidence(line, left_range)?;
    if let Some(adduct) = &adduct_type {
        if adduct.charge() != charge {
            return Err(CustomError::error(
                "Invalid mzPAF annotation",
                "The defined charge should be identical to the total charge as defined in the adduct ions",
                Context::line_range(None, line, range),
            ));
        }
    }
    Ok((
        left_range,
        PeakAnnotation {
            auxiliary,
            analyte_number,
            ion,
            neutral_losses,
            charge: adduct_type.unwrap_or_else(|| MolecularCharge::proton(charge.value)),
            deviation,
            confidence,
        },
    ))
}

/// An mzPAF single peak annotation.
#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct PeakAnnotation {
    auxiliary: bool,
    analyte_number: Option<usize>,
    ion: IonType,
    neutral_losses: Vec<NeutralLoss>,
    charge: MolecularCharge,
    deviation: Option<Tolerance<MassOverCharge>>,
    confidence: Option<f64>,
}

#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
enum IonType {
    Unknown(Option<usize>),
    MainSeries(char, usize),
    Immonium(AminoAcid, Option<SimpleModification>),
    Internal(usize, usize),
    Named(String),
    Precursor,
    Reporter(MolecularFormula),
    Formula(MolecularFormula),
}

/// Parse a mzPAF analyte number. '1@...'
/// # Errors
/// When the ion is not formatted correctly.
fn parse_analyte_number(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<usize>), CustomError> {
    next_number::<false, false, usize>(line, range.clone()).map_or_else(
        || Ok((range.clone(), None)),
        |num| {
            if line.chars().nth(num.0) != Some('@') {
                return Err(CustomError::error(
                    "Invalid mzPAF analyte number",
                    "The analyte number should be followed by an at sign '@'",
                    Context::line(None, line, num.0, 1),
                ));
            }
            Ok((
                range.add_start(num.0 + 1),
                Some(num.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF analyte number",
                        format!("The analyte number number {}", explain_number_error(&err)),
                        Context::line(None, line, 0, num.0),
                    )
                })?),
            ))
        },
    )
}

/// Parse a mzPAF ion.
/// # Errors
/// When the ion is not formatted correctly.
fn parse_ion(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, IonType), CustomError> {
    match line.chars().nth(range.start_index()) {
        Some('?') => {
            if let Some(ordinal) =
                next_number::<false, false, usize>(line, range.add_start(1_usize))
            {
                Ok((
                    range.add_start(1 + ordinal.0),
                    IonType::Unknown(Some(ordinal.2.map_err(|err| {
                        CustomError::error(
                            "Invalid mzPAF unknown ion ordinal",
                            format!("The ordinal number {}", explain_number_error(&err)),
                            Context::line(None, line, range.start_index() + 1, ordinal.0),
                        )
                    })?)),
                ))
            } else {
                Ok((range.add_start(1_usize), IonType::Unknown(None)))
            }
        }
        Some(c @ ('a' | 'b' | 'c' | 'x' | 'y' | 'z')) => {
            if let Some(ordinal) =
                next_number::<false, false, usize>(line, range.add_start(1_usize))
            {
                Ok((
                    range.add_start(1 + ordinal.0),
                    IonType::MainSeries(
                        c,
                        ordinal.2.map_err(|err| {
                            CustomError::error(
                                "Invalid mzPAF unknown ion ordinal",
                                format!("The ordinal number {}", explain_number_error(&err)),
                                Context::line(None, line, range.start_index() + 1, ordinal.0),
                            )
                        })?,
                    ),
                ))
                // TODO: potentially followed by a ProForma sequence in {}
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF main series ion ordinal",
                    "For a main series ion the ordinal should be provided, like 'a12'",
                    Context::line(None, line, range.start_index(), 1),
                ))
            }
        }
        Some('I') => {
            let amino_acid = line[range.clone()].chars().nth(1).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF immonium",
                    "The source amino acid for this immonium ion should be present like 'IA'",
                    Context::line(None, line, range.start_index(), 1),
                )
            })?;
            let modification = if line[range.clone()].chars().nth(2) == Some('[') {
                let first = line[range.clone()].char_indices().nth(3).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(3)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Some((
                    last.0 + last.1.len_utf8() - first,
                    Ontology::Unimod
                        .find_name(
                            &line[range.clone()][first..last.0 + last.1.len_utf8()],
                            None,
                        )
                        .ok_or_else(|| {
                            Ontology::Unimod.find_closest(
                                &line[range.clone()][first..last.0 + last.1.len_utf8()],
                                None,
                            )
                        })?,
                ))
            } else {
                None
            };
            Ok((
                range.add_start(2 + modification.as_ref().map_or(0, |m| m.0)),
                IonType::Immonium(
                    AminoAcid::try_from(amino_acid).map_err(|()| {
                        CustomError::error(
                            "Invalid mzPAF immonium ion",
                            "The provided amino acid is not a known amino acid",
                            Context::line(None, line, range.start_index() + 1, 1),
                        )
                    })?,
                    modification.map(|m| m.1),
                ),
            ))
        }
        Some('m') => {
            let first_ordinal = next_number::<false, false, usize>(line, range.add_start(1_usize))
                .ok_or_else(|| {
                    CustomError::error(
                        "Invalid mzPAF internal ion first ordinal",
                        "The first ordinal for an internal ion should be present",
                        Context::line(None, line, range.start_index(), 1),
                    )
                })?;
            if line[range.clone()].chars().nth(first_ordinal.0) != Some(':') {
                return Err(CustomError::error(
                    "Invalid mzPAF internal ion ordinal separator",
                    "The internal ion ordinal separator should be a colon ':', like 'm4:6'",
                    Context::line(None, line, range.start_index() + 1 + first_ordinal.0, 1),
                ));
            }
            assert!(
                line[range.clone()].chars().nth(first_ordinal.0) == Some(':'),
                "Needs to be separated by colon"
            );
            let second_ordinal = next_number::<false, false, usize>(
                line,
                range.add_start(2 + first_ordinal.0 as isize),
            )
            .ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF internal ion second ordinal",
                    "The second ordinal for an internal ion should be present",
                    Context::line(None, line, range.start_index() + 1 + first_ordinal.0, 1),
                )
            })?;
            let first_location = first_ordinal.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF internal ion first ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    Context::line(None, line, range.start_index() + 1, first_ordinal.0),
                )
            })?;
            let second_location = second_ordinal.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF internal ion second ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    Context::line(
                        None,
                        line,
                        range.start_index() + 2 + first_ordinal.0,
                        second_ordinal.0,
                    ),
                )
            })?;
            Ok((
                range.add_start(2 + first_ordinal.0 + second_ordinal.0),
                IonType::Internal(first_location, second_location),
            ))
        }
        Some('_') => {
            // Format less strings
            // TODO: Potentially recognise the following as known contaminants:
            // 0@_{y1(R)}
            // 0@_{a2(LP)}
            // 0@_{b2(LP)}

            let (len, name) = if line[range.clone()].chars().nth(1) == Some('{') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != '}')
                    .last()
                    .unwrap();
                Ok((
                    last.0 + last.1.len_utf8() - first,
                    &line[range.clone()][first..last.0 + last.1.len_utf8()],
                ))
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF named compound",
                    "A named compound must be named with curly braces '{}' after the '_'",
                    Context::line(None, line, range.start_index(), 1),
                ))
            }?;
            Ok((range.add_start(3 + len), IonType::Named(name.to_string())))
        }
        Some('p') => Ok((range.add_start(1_usize), IonType::Precursor)),
        Some('r') => {
            // Same name as neutral losses
            let (len, name) = if line[range.clone()].chars().nth(1) == Some('[') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != ']')
                    .last()
                    .unwrap();
                Ok((
                    last.0 + last.1.len_utf8() - first,
                    &line[range.clone()][first..last.0 + last.1.len_utf8()],
                ))
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF reporter ion",
                    "A reporter ion must be named with square braces '[]' after the 'r'",
                    Context::line(None, line, range.start_index(), 1),
                ))
            }?;
            MZPAF_NAMED_MOLECULES
                .iter()
                .find_map(|n| (n.0 == name).then_some(n.1.clone()))
                .map_or_else(
                    || {
                        Err(CustomError::error(
                            "Unknown mzPAF named reporter ion",
                            "Unknown name",
                            Context::line(None, line, range.start_index() + 1, len),
                        ))
                    },
                    |formula| Ok((range.add_start(3 + len), IonType::Reporter(formula))),
                )
        }
        Some('f') => {
            // Simple formula
            let formula_range = if line[range.clone()].chars().nth(1) == Some('{') {
                let first = line[range.clone()].char_indices().nth(2).unwrap().0;
                let last = line[range.clone()]
                    .char_indices()
                    .skip(2)
                    .take_while(|(_, c)| *c != '}')
                    .last()
                    .unwrap();
                Ok(range.start_index() + first..range.start_index() + last.0 + last.1.len_utf8())
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF formula",
                    "A formula must have the formula defined with curly braces '{}' after the 'f'",
                    Context::line(None, line, range.start_index(), 1),
                ))
            }?;
            let formula =
                MolecularFormula::from_pro_forma(line, formula_range.clone(), false, false, true)?;

            Ok((
                range.add_start(3 + formula_range.len()),
                IonType::Formula(formula),
            ))
        }
        Some('s') => todo!(), // TODO: return as Formula
        Some(_) => Err(CustomError::error(
            "Invalid ion",
            "An ion cannot start with this character",
            Context::line(None, line, range.start, 1),
        )),
        None => Err(CustomError::error(
            "Invalid ion",
            "An ion cannot be an empty string",
            Context::line_range(None, line, range),
        )),
    }
}

fn parse_neutral_loss(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, Vec<NeutralLoss>), CustomError> {
    let mut offset = 0;
    let mut neutral_losses = Vec::new();
    while let Some(c @ ('-' | '+')) = line.chars().nth(range.start_index() + offset) {
        if line.chars().nth(range.start_index() + 1) == Some('[') {
            let first = line.char_indices().nth(range.start_index() + 2).unwrap().0;
            let last = line
                .char_indices()
                .skip(range.start_index() + 2)
                .take_while(|(_, c)| *c != ']')
                .last()
                .unwrap();
            let name = line[first..last.0 + last.1.len_utf8()].to_ascii_lowercase();

            offset += 1 + last.0 + last.1.len_utf8() - first;

            if let Some(formula) = MZPAF_NAMED_MOLECULES
                .iter()
                .find_map(|n| (n.0 == name).then_some(n.1.clone()))
            {
                neutral_losses.push(match c {
                    '+' => NeutralLoss::Gain(formula),
                    '-' => NeutralLoss::Loss(formula),
                    _ => unreachable!(),
                });
            } else {
                return Err(CustomError::error(
                    "Unknown mzPAF named neutral loss",
                    "Unknown name",
                    Context::line(None, line, offset - name.len() - 1, name.len()),
                ));
            }
        } else {
            let first = line.char_indices().nth(range.start_index() + 1).unwrap().0;
            let last = line
                .char_indices()
                .skip(range.start_index() + 2)
                .take_while(|(_, c)| c.is_ascii_alphanumeric())
                .last()
                .unwrap();
            let formula = MolecularFormula::from_pro_forma(
                line,
                first..last.0 + last.1.len_utf8(),
                false,
                false,
                true,
            )?;
            neutral_losses.push(match c {
                '+' => NeutralLoss::Gain(formula),
                '-' => NeutralLoss::Loss(formula),
                _ => unreachable!(),
            });
            offset += 1 + last.0 + last.1.len_utf8() - first;
        }
    }
    Ok((range.add_start(offset), neutral_losses))
}

fn parse_adduct_type(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, Option<MolecularCharge>), CustomError> {
    if line.chars().nth(range.start_index()) == Some('[') {
        if line
            .chars()
            .nth(range.start_index() + 1)
            .map(|c| c.to_ascii_lowercase())
            != Some('m')
        {
            return Err(CustomError::error(
                "Invalid mzPAF adduct type",
                "The adduct type should start with 'M', as in '[M+nA]'",
                Context::line(None, line, range.start_index() + 1, 1),
            ));
        }
        let mut carriers = Vec::new();
        let mut offset = 2;
        while let Some(number) = next_number::<true, false, isize>(line, range.add_start(offset)) {
            let first = line
                .char_indices()
                .nth(range.start_index() + offset + number.0)
                .unwrap()
                .0;
            let last = line
                .char_indices()
                .skip(range.start_index() + offset + number.0)
                .take_while(|(_, c)| c.is_ascii_alphanumeric())
                .fold((0, 0, ' '), |acc, (index, c)| (acc.0 + 1, index, c));
            let formula = MolecularFormula::from_pro_forma(
                line,
                first..last.1 + last.2.len_utf8(),
                false,
                false,
                true,
            )?;
            carriers.push((
                number.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF adduct ordinal",
                        format!("The ordinal number {}", explain_number_error(&err)),
                        Context::line(None, line, range.start_index() + offset, number.0),
                    )
                })?,
                formula,
            ));
            offset += number.0 + last.0;
        }
        if line.chars().nth(range.start_index() + offset) != Some(']') {
            return Err(CustomError::error(
                "Invalid mzPAF adduct type",
                "The adduct type should be closed with ']'",
                Context::line(None, line, range.start_index() + offset, 1),
            ));
        }
        Ok((
            range.add_start(offset + 1),
            Some(MolecularCharge::new(&carriers)),
        ))
    } else {
        Ok((range, None))
    }
}

fn parse_charge(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, Charge), CustomError> {
    if line.chars().nth(range.start_index()) == Some('^') {
        let charge =
            next_number::<false, false, u32>(line, range.add_start(1_usize)).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF charge",
                    "The number after the charge symbol should be present, eg '^2'.",
                    Context::line(None, line, range.start_index(), 1),
                )
            })?;
        Ok((
            range.add_start(charge.0 + 1),
            Charge::new::<e>(charge.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF charge",
                    format!("The charge number {}", explain_number_error(&err)),
                    Context::line(None, line, range.start_index() + 1, charge.0),
                )
            })? as isize),
        ))
    } else {
        Ok((range, Charge::new::<e>(1)))
    }
}

// fn parse_adduct(
//     line: &str,
//     range: Range<usize>,
// ) -> Result<(Characters, MolecularFormula), CustomError> {
//     if line[range.clone()].chars().next() == Some('^') {
//         let charge =
//             next_number::<false, false, u32>(line, range.add_start(1)).ok_or_else(|| {
//                 CustomError::error(
//                     "Invalid mzPAF charge",
//                     "The number after the charge symbol should be present, eg '^2'.",
//                     Context::line(None, line, range.start_index(), 1),
//                 )
//             })?;
//         Ok((
//             charge.0 + 1,
//             Charge::new::<e>(charge.2.map_err(|err| {
//                 CustomError::error(
//                     "Invalid mzPAF charge",
//                     format!("The charge number {}", explain_number_error(&err)),
//                     Context::line(None, line, range.start_index() + 1, charge.0),
//                 )
//             })? as isize),
//         ))
//     } else {
//         Ok((0, Charge::new::<e>(1)))
//     }
// }

/// Parse a mzPAF deviation, either a ppm or mz deviation.
/// # Errors
/// When the deviation is not '<number>' or '<number>ppm'.
fn parse_deviation(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, Option<Tolerance<MassOverCharge>>), CustomError> {
    if line.chars().nth(range.start_index()) == Some('/') {
        let number =
            next_number::<true, true, f64>(line, range.add_start(1_usize)).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF deviation",
                    "A deviation should be a number",
                    Context::line_range(None, line, range.start..=range.start + 1),
                )
            })?;
        let deviation = number.2.map_err(|err| {
            CustomError::error(
                "Invalid mzPAF deviation",
                format!("The deviation number {err}",),
                Context::line_range(None, line, range.start + 1..range.start + 1 + number.0),
            )
        })?;
        if line
            .chars()
            .skip(range.start_index() + 1)
            .take(3)
            .collect::<String>()
            .to_ascii_lowercase()
            == "ppm"
        {
            Ok((
                range.add_start(1 + number.0 + 3),
                Some(Tolerance::new_ppm(deviation)),
            ))
        } else {
            Ok((
                range.add_start(1 + number.0),
                Some(Tolerance::new_absolute(MassOverCharge::new::<mz>(
                    deviation,
                ))),
            ))
        }
    } else {
        Ok((range, None))
    }
}

/// Parse a mzPAF confidence.
/// # Errors
/// When the deviation is not '*<number>'.
fn parse_confidence(
    line: &str,
    range: Range<Characters>,
) -> Result<(Range<Characters>, Option<f64>), CustomError> {
    if line.chars().nth(range.start_index()) == Some('*') {
        let number =
            next_number::<true, true, f64>(line, range.add_start(1_usize)).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF confidence",
                    "A confidence should be a number",
                    Context::line_range(None, line, range.start..=range.start + 1),
                )
            })?;
        let confidence = number.2.map_err(|err| {
            CustomError::error(
                "Invalid mzPAF confidence",
                format!("The confidence number {err}",),
                Context::line_range(None, line, range.start + 1..range.start + 1 + number.0),
            )
        })?;
        Ok((range.add_start(number.0 + 1), Some(confidence)))
    } else {
        Ok((range, None))
    }
}

static MZPAF_NAMED_MOLECULES: LazyLock<Vec<(&str, MolecularFormula)>> = LazyLock::new(|| {
    vec![
        ("hex", molecular_formula!(C 6 H 10 O 5)),
        ("hexnac", molecular_formula!(C 8 H 13 N 1 O 5)),
        ("dhex", molecular_formula!(C 6 H 10 O 4)),
        ("neuac", molecular_formula!(C 11 H 17 N 1 O 8)),
        ("neugc", molecular_formula!(C 11 H 17 N 1 O 9)),
        ("tmt126", molecular_formula!(C 8 N 1 H 15)),
        ("tmt127n", molecular_formula!(C 8 [15 N 1] H 15)),
        ("tmt127c", molecular_formula!(C 7 [13 C 1] N 1 H 15)),
        ("tmt128n", molecular_formula!(C 7 [13 C 1] [15 N 1] H 15)),
        ("tmt128c", molecular_formula!(C 6 [13 C 2] N 1 H 15)),
        ("tmt129n", molecular_formula!(C 6 [13 C 2] [15 N 1] H 15)),
        ("tmt129c", molecular_formula!(C 5 [13 C 3] N 1 H 15)),
        ("tmt130n", molecular_formula!(C 5 [13 C 3] [15 N 1] H 15)),
        ("tmt130c", molecular_formula!(C 4 [13 C 4] N 1 H 15)),
        ("tmt131n", molecular_formula!(C 4 [13 C 4] [15 N 1] H 15)),
        ("tmt131c", molecular_formula!(C 3 [13 C 5] N 1 H 15)),
        ("tmt132n", molecular_formula!(C 3 [13 C 5] [15 N 1] H 15)),
        ("tmt132c", molecular_formula!(C 2 [13 C 6] N 1 H 15)),
        ("tmt133n", molecular_formula!(C 2 [13 C 6] [15 N 1] H 15)),
        ("tmt133c", molecular_formula!(C 1 [13 C 7] N 1 H 15)),
        ("tmt134n", molecular_formula!(C 1 [13 C 7] [15 N 1] H 15)),
        ("tmt134c", molecular_formula!(C 0 [13 C 8] N 1 H 15)),
        ("tmt135n", molecular_formula!(C 0 [13 C 8] [15 N 1] H 15)),
        ("tmtzero", molecular_formula!(C 12 H 20 N 2 O 2)),
        ("tmtpro_zero", molecular_formula!(C 15 H 25 N 3 O 3)),
        ("tmt2plex", molecular_formula!(C 11 [ 13 C 1] H 20 N 2 O 2)),
        (
            "tmt6plex",
            molecular_formula!(C 8 [13 C 5] H 20 N 1 [ 15 N 1] O 2),
        ),
        (
            "tmtpro",
            molecular_formula!(C 8 [13 C 7] H 25 [15 N 2] N 1 O 3),
        ),
        ("itraq113", molecular_formula!(C 6 N 2 H 12)),
        ("itraq114", molecular_formula!(C 5 [13 C 1] N 2 H 12)),
        (
            "itraq115",
            molecular_formula!(C 5 [13 C 1] N 1 [15 N 1] H 12),
        ),
        (
            "itraq116",
            molecular_formula!(C 4 [13 C 2] N 1 [15 N 1] H 12),
        ),
        (
            "itraq117",
            molecular_formula!(C 3 [13 C 3] N 1 [15 N 1] H 12),
        ),
        ("itraq118", molecular_formula!(C 3 [13 C 3] [15 N 2] H 12)),
        ("itraq119", molecular_formula!(C 4 [13 C 2] [15 N 2] H 12)),
        ("itraq121", molecular_formula!([13 C 6] [15 N 2] H 12)),
        (
            "itraq4plex",
            molecular_formula!(C 4 [13 C 3] H 12 N 1 [15 N 1] O 1),
        ),
        (
            "itraq8plex",
            molecular_formula!(C 7 [13 C 7] H 24 N 3 [15 N 1] O 3),
        ),
        ("tmt126-etd", molecular_formula!(C 7 N 1 H 15)),
        ("tmt127n-etd", molecular_formula!(C 7 [15 N 1] H 15)),
        ("tmt127c-etd", molecular_formula!(C 6 [13 C 1] N 1 H 15)),
        (
            "tmt128n-etd",
            molecular_formula!(C 6 [13 C 1] [15 N 1] H 15),
        ),
        ("tmt128c-etd", molecular_formula!(C 5 [13 C 2] N 1 H 15)),
        (
            "tmt129n-etd",
            molecular_formula!(C 5 [13 C 2] [15 N 1] H 15),
        ),
        ("tmt129c-etd", molecular_formula!(C 4 [13 C 3] N 1 H 15)),
        (
            "tmt130n-etd",
            molecular_formula!(C 4 [13 C 3] [15 N 1] H 15),
        ),
        ("tmt130c-etd", molecular_formula!(C 3 [13 C 4] N 1 H 15)),
        (
            "tmt131n-etd",
            molecular_formula!(C 3 [13 C 4] [15 N 1] H 15),
        ),
        ("tmt131c-etd", molecular_formula!(C 2 [13 C 5] N 1 H 15)),
    ]
});
