//! WIP: mzPAF parser
#![allow(dead_code)]
use std::{ops::Range, sync::LazyLock};

use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{MolecularCharge, MolecularFormula},
    error::{Context, CustomError},
    fragment::NeutralLoss,
    helper_functions::{
        Characters, RangeExtension, RangeMaths, end_of_enclosure, explain_number_error, next_number,
    },
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    quantities::Tolerance,
    sequence::{AminoAcid, CompoundPeptidoformIon, Peptidoform, SemiAmbiguous, SimpleModification},
    system::{MassOverCharge, e, isize::Charge, mz},
};
// TODO: ask about mzPAF spec:
// * neutral losses are either names encompassed by [] or formulas, but since isotopes are allowed in formulas seeing [] is ambiguous between a formula starting with an isotope or a name
// * are isotopes supported in adduct ion additions? (I would assume so, but not clearly stated)
// * page 12 styling for examples not right
// * Is 'Adenosine' supposed to be a reporter ion? It is not in the list but it is used as an example

/// Parse a mzPAF peak annotation line (can contain multiple annotations).
/// # Errors
/// When the annotation does not follow the format.
pub fn parse_mzpaf(
    line: &str,
    custom_database: Option<&CustomDatabase>,
) -> Result<Vec<PeakAnnotation>, CustomError> {
    let mut annotations = Vec::new();

    // Parse first
    let (mut range, a) = parse_annotation(line, 0..line.len(), custom_database)?;
    annotations.push(a);

    // Parse any following
    while !range.is_empty() {
        if line.as_bytes().get(range.start_index()).copied() == Some(b',') {
            range = range.add_start(1_usize);
        } else {
            return Err(CustomError::error(
                "Invalid mzPAF annotation delimiter",
                "Different mzPAF annotations should be separated with commas ','.",
                Context::line(None, line, range.start_index(), 1),
            ));
        }
        let (r, a) = parse_annotation(line, range, custom_database)?;
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
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<(Range<usize>, PeakAnnotation), CustomError> {
    let (left_range, auxiliary) = if line.as_bytes().get(range.start_index()).copied() == Some(b'&')
    {
        (range.add_start(1_usize), true)
    } else {
        (range.clone(), false)
    };
    let (left_range, analyte_number) = parse_analyte_number(line, left_range)?;
    let (left_range, ion) = parse_ion(line, left_range, custom_database)?;
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
    MainSeries(u8, Option<u8>, usize, Option<Peptidoform<SemiAmbiguous>>),
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
            if line.as_bytes().get(num.0 + 1).copied() != Some(b'@') {
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
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<(Range<usize>, IonType), CustomError> {
    match line.as_bytes().get(range.start_index()).copied() {
        Some(b'?') => {
            if let Some(ordinal) =
                next_number::<false, false, usize>(line, range.add_start(1_usize))
            {
                Ok((
                    range.add_start(1 + ordinal.0),
                    IonType::Unknown(Some(ordinal.2.map_err(|err| {
                        CustomError::error(
                            "Invalid mzPAF unknown ion ('?') ordinal",
                            format!("The ordinal number {}", explain_number_error(&err)),
                            Context::line(None, line, range.start_index() + 1, ordinal.0),
                        )
                    })?)),
                ))
            } else {
                Ok((range.add_start(1_usize), IonType::Unknown(None)))
            }
        }
        Some(c @ (b'a' | b'b' | b'c' | b'd' | b'v' | b'w' | b'x' | b'y' | b'z')) => {
            let (range, sub) = if let Some(sub @ (b'a' | b'b')) =
                line.as_bytes().get(range.start_index() + 1).copied()
            {
                if c != b'd' && c != b'w' {
                    return Err(CustomError::error(
                        "Invalid mzPAF main series ion ordinal",
                        "Only for the satellite ions 'd' and 'w' does a subtype exist, like 'wa12'",
                        Context::line(None, line, range.start_index(), 1),
                    ));
                }
                (range.add_start(2_usize), Some(sub))
            } else {
                (range.add_start(1_usize), None)
            };
            if let Some(ordinal) =
                next_number::<false, false, usize>(line, range.clone())
            {
                let (range, interpretation) = if line
                    .as_bytes()
                    .get(range.start_index() + ordinal.0)
                    .copied()
                    == Some(b'{')
                {
                    if let Some(location) =
                        end_of_enclosure(line, range.start_index() + 1 + ordinal.0, b'{', b'}')
                    {
                        let interpretation = Peptidoform::pro_forma(
                            &line[range.start_index() + ordinal.0 + 1
                                ..location],
                            custom_database,
                        );
                        (
                            range.add_start(ordinal.0 + location),
                            Some(interpretation.unwrap().into_semi_ambiguous().unwrap()),
                        )
                        // TODO: proper error handling and add checks to the length of the sequence
                    } else {
                        return Err(CustomError::error(
                            "Invalid mzPAF main series ion ordinal",
                            "The asserted interpretation should have a closed curly bracket, like '0@b2{LL}'",
                            Context::line(None, line, range.start_index() + ordinal.0, 1),
                        ));
                    }
                } else {
                    (range, None)
                };
                Ok((
                    range.add_start( ordinal.0),
                    IonType::MainSeries(
                        c,
                        sub,
                        ordinal.2.map_err(|err| {
                            CustomError::error(
                                "Invalid mzPAF ion ordinal",
                                format!("The ordinal number {}", explain_number_error(&err)),
                                Context::line(None, line, range.start_index(), ordinal.0),
                            )
                        })?,
                        interpretation,
                    ),
                ))
            } else {
                Err(CustomError::error(
                    "Invalid mzPAF main series ion ordinal",
                    "For a main series ion the ordinal should be provided, like 'a12'",
                    Context::line(None, line, range.start_index(), 1),
                ))
            }
        }
        Some(b'I') => {
            let amino_acid = line[range.clone()].chars().next().ok_or_else(|| {
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
                    last.0 + last.1.len_utf8() - first + 2, // Length of mod + [ + ]
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
        Some(b'm') => {
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
        Some(b'_') => {
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
        Some(b'p') => Ok((range.add_start(1_usize), IonType::Precursor)),
        Some(b'r') => {
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
                    &line[range.clone()][first + range.start_index()..range.start_index() + last.0 + last.1.len_utf8()],
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
                .find_map(|n| (n.0.eq_ignore_ascii_case(name)).then_some(n.1.clone()))
                .map_or_else(
                    || {
                        Err(CustomError::error(
                            "Unknown mzPAF named reporter ion",
                            "Unknown name",
                            Context::line(None, line, range.start_index() + 2, len),
                        ))
                    },
                    |formula| Ok((range.add_start(3 + len), IonType::Reporter(formula))),
                )
        }
        Some(b'f') => {
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
                MolecularFormula::from_pro_forma(line, formula_range.clone(), false, false, true, false)?;

            Ok((
                range.add_start(3 + formula_range.len()),
                IonType::Formula(formula),
            ))
        }
        Some(b's') => todo!("SIMLES not (yet) supported"), // TODO: return as Formula
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

// TODO needs to backtrack once it detects an isotope
fn parse_neutral_loss(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Vec<NeutralLoss>), CustomError> {
    let mut offset = 0;
    let mut neutral_losses = Vec::new();
    while let Some(c @ (b'-' | b'+')) = line.as_bytes().get(range.start_index() + offset).copied() {
        let mut amount = 1;
        // Parse leading number to detect how many times this loss occured
        if let Some(num) = next_number::<false, false, u16>(line, range.add_start(1 + offset)) {
            amount = i32::from(num.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF neutral loss leading amount",
                    format!("The neutral loss amount number {}", explain_number_error(&err)),
                    Context::line(None, line, range.start_index() + 1 + offset, range.start_index() + 1 + offset + num.0),
                )
            })?);
            offset += num.0;
        }

        if line.as_bytes().get(range.start_index() + 1 + offset).copied() == Some(b'[') {
            let last =
                end_of_enclosure(line, range.start_index() + 2 + offset, b'[', b']').ok_or_else(|| {
                    CustomError::error(
                        "Unknown mzPAF named neutral loss",
                        "Opening bracket for neutral loss name was not closed",
                        Context::line(None, line, range.start_index() + 1 + offset, 1),
                    )
                })?;
            let first = range.start_index() + 2 + offset;
            let name = &line[first..last];

            offset += 3 + last - first; // Sign, brackets, and the name

            if let Some(formula) = MZPAF_NAMED_MOLECULES
                .iter()
                .find_map(|n| (n.0.eq_ignore_ascii_case(name)).then_some(n.1.clone()))
            {
                println!("Found: ['{}'] at {}..{}", &formula,  range.start_index() + offset - name.len() - 1,  range.start_index() +  offset);
                neutral_losses.push(match c {
                    b'+' => NeutralLoss::Gain(formula * amount),
                    b'-' => NeutralLoss::Loss(formula * amount),
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
            let first = range.start_index() + 1 + offset;
            let last = line[first..]
                .char_indices()
                .take_while(|(_, c)| c.is_ascii_alphanumeric()) // TODO check if isotopes are allowed here, theoretically yes, but not clearly specified in the spec
                .last()
                .unwrap();
            let formula = MolecularFormula::from_pro_forma( // TODO 'i' is seen as iodine while the thing should casing specific
                line,
                first..first + last.0 + last.1.len_utf8(),
                false,
                false,
                true,
                false,
            )?;
            println!("Found: '{}' at {}..{}", &formula,  range.start_index() + offset, range.start_index() + offset + 1 + last.0 + last.1.len_utf8());
            neutral_losses.push(match c {
                b'+' => NeutralLoss::Gain(formula * amount),
                b'-' => NeutralLoss::Loss(formula * amount),
                _ => unreachable!(),
            });
            offset += 1 + last.0 + last.1.len_utf8();
        }
    }
    Ok((range.add_start(offset), neutral_losses))
}

fn parse_adduct_type(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<MolecularCharge>), CustomError> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'[') {
        if line
            .as_bytes()
            .get(range.start_index() + 1).copied()
            != Some(b'M')
        {
            return Err(CustomError::error(
                "Invalid mzPAF adduct type",
                "The adduct type should start with 'M', as in '[M+nA]'",
                Context::line(None, line, range.start_index() + 1, 1),
            ));
        }
        let mut carriers = Vec::new();
        let mut offset = 2;
        println!("{}", &line[range.add_start(offset)]);
         while let Some(c @ (b'-' | b'+')) = line.as_bytes().get(range.start_index() + offset).copied() {
            offset += 1; // The sign
            let mut amount = 1;
            // Parse leading number to detect how many times this adduct occured
            if let Some(num) = next_number::<false, false, u16>(line, range.add_start(1 + offset)) {
                amount = i32::from(num.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF adduct leading amount",
                        format!("The adduct amount number {}", explain_number_error(&err)),
                        Context::line(None, line, range.start_index() + 1 + offset, range.start_index() + 1 + offset + num.0),
                    )
                })?);
                offset += num.0;
            }
            if c == b'-' {
                amount *= -1;
            }

            let first = range.start_index() + offset;
            let last = line[first..]
                .char_indices()
                .take_while(|(_, c)| c.is_ascii_alphanumeric()) // TODO are isotopes allowed here?
                .last().map_or(0, |last| last.0 + last.1.len_utf8());
            let formula = MolecularFormula::from_pro_forma(
                line,
                first..first + last,
                false,
                false,
                true,
                false,
            )?;
            carriers.push((
                amount as isize,
                formula,
            ));
            offset += last;
        }
        if line.as_bytes().get(range.start_index() + offset).copied() != Some(b']') {
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
    range: Range<usize>,
) -> Result<(Range<usize>, Charge), CustomError> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'^') {
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
    range: Range<usize>,
) -> Result<(Range<usize>, Option<Tolerance<MassOverCharge>>), CustomError> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'/') {
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
        if line[range.start_index() + 1 + number.0..].starts_with("ppm")
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

/// Create a parse test based on a given case and its name.
#[macro_export]
macro_rules! mzpaf_test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::fragment::parse_mzpaf($case, None);
            println!("{}", $case);
            assert!(res.is_ok(), "{}", res.err().unwrap());
            // let back = res.as_ref().unwrap().to_string();
            // let res_back = $crate::fragment::parse_mzpaf(&back, None);
            // assert_eq!(res, res_back, "{} != {back}", $case);
        }
    };
    (ne $case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let res = $crate::fragment::parse_mzpaf($case, None);
            println!("{}\n{:?}", $case, res);
            assert!(res.is_err());
        }
    };
}

mzpaf_test!("b2-H2O/3.2ppm,b4-H2O^2/3.2ppm", positive_1);
mzpaf_test!("b2-H2O/3.2ppm*0.75,b4-H2O^2/3.2ppm*0.25", positive_2);
mzpaf_test!("1@y12/0.13,2@b9-NH3/0.23", positive_3);
mzpaf_test!("0@y1{K}", positive_4);
mzpaf_test!("0@y1{K}-NH3", positive_5);
mzpaf_test!("y1/-1.4ppm", positive_6);
mzpaf_test!("y1/-0.0002", positive_7);
mzpaf_test!("y4-H2O+2i[M+H+Na]^2", positive_8);
mzpaf_test!("?", positive_9);
mzpaf_test!("?^3", positive_10);
mzpaf_test!("?+2i^4", positive_11);
mzpaf_test!("?17", positive_12);
mzpaf_test!("?17+i/1.45ppm", positive_13);
mzpaf_test!("?17-H2O/-0.87ppm", positive_14);
mzpaf_test!("0@b2{LL}", positive_15);
mzpaf_test!("0@y1{K}", positive_16);
mzpaf_test!("0@b2{LC[Carbamidomethyl]}", positive_17);
mzpaf_test!("0@b1{[Acetyl]-M}", positive_18);
mzpaf_test!("0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", positive_19);
mzpaf_test!("m3:6", positive_20);
mzpaf_test!("b3-C2H3NO", positive_21);
mzpaf_test!("m3:6-CO", positive_22);
mzpaf_test!("m3:6-CO-H2O^2", positive_23);
mzpaf_test!("m3:4/1.1ppm,m4:5/1.1ppm", positive_24);
mzpaf_test!("m3:5", positive_25);
mzpaf_test!("IY", positive_26);
mzpaf_test!("IH", positive_27);
mzpaf_test!("IL-CH2", positive_28);
mzpaf_test!("IC[Carbamidomethyl]", positive_29);
mzpaf_test!("IY[Phospho]", positive_30);
mzpaf_test!("p^2", positive_31);
mzpaf_test!("p-H3PO4^2", positive_32);
mzpaf_test!("p^4", positive_33);
mzpaf_test!("p+H^3", positive_34);
mzpaf_test!("p^3 ", positive_35);
mzpaf_test!("p+2H^2", positive_36);
mzpaf_test!("p^2 ", positive_37);
mzpaf_test!("p+H^2", positive_38);
mzpaf_test!("p+3H", positive_39);
mzpaf_test!("p+2H", positive_40);
mzpaf_test!("p+H", positive_41);
mzpaf_test!("p", positive_42);
mzpaf_test!("r[TMT127N]", positive_43);
mzpaf_test!("r[iTRAQ114]", positive_44);
mzpaf_test!("r[TMT6plex]", positive_45);
mzpaf_test!("r[Hex]", positive_46);
mzpaf_test!("r[Adenosine]", positive_47);
mzpaf_test!("0@_{Urocanic Acid}", positive_48);
mzpaf_test!("f{C13H9}/-0.55ppm", positive_49);
mzpaf_test!("f{C12H9N}/0.06ppm", positive_50);
mzpaf_test!("f{C13H9N}/-2.01ppm", positive_51);
mzpaf_test!("f{C13H10N}/-0.11ppm", positive_52);
mzpaf_test!("f{C13H11N}/-0.09ppm", positive_53);
mzpaf_test!("f{C13H12N}/0.26ppm", positive_54);
mzpaf_test!("f{C14H10N}/0.19ppm", positive_55);
mzpaf_test!("f{C14H11N}/0.45ppm", positive_56);
mzpaf_test!("f{C14H10NO}/0.03ppm", positive_57);
mzpaf_test!("f{C16H22O}+i^3", positive_58);
mzpaf_test!("f{C15[13C1]H22O}^3", positive_59);
mzpaf_test!("s{CN=C=O}[M+H]/-0.55ppm", positive_60);
mzpaf_test!("s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm", positive_61);
mzpaf_test!("p-[Hex]", positive_62);
mzpaf_test!("y2+CO-H2O", positive_63);
mzpaf_test!("y2-H2O-NH3", positive_64);
mzpaf_test!("p-[Hex]", positive_65);
mzpaf_test!("p-[TMT6plex]-2H2O-HPO3", positive_66);
mzpaf_test!("p-2[iTRAQ115]", positive_67);
mzpaf_test!("p-[iTRAQ116]-CO-H2O-HPO3", positive_68);
mzpaf_test!("y4[M+Na]", positive_69);
mzpaf_test!("y5-H2O[M+H+Na]^2", positive_70);
mzpaf_test!("&1@y7/-0.002", positive_71);
mzpaf_test!("&y7/-0.001", positive_72);
mzpaf_test!("y7/0.000*0.95", positive_73);
mzpaf_test!("&y7/0.001", positive_74);
mzpaf_test!("&y7/0.002", positive_75);
mzpaf_test!("b6-H2O/-0.005,&y7/0.003", positive_76);
mzpaf_test!("y12-H2O^2/7.4ppm*0.70", positive_77);
mzpaf_test!("y12/3.4ppm*0.85,b9-NH3/5.2ppm*0.05", positive_78);
