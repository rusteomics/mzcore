//! WIP: mzPAF parser
#![allow(dead_code)]
use std::{num::NonZeroU16, ops::Range, sync::LazyLock};

use serde::{Deserialize, Serialize};

use crate::{
    chemistry::{ELEMENT_PARSE_LIST, Element, MolecularCharge, MolecularFormula},
    error::{Context, CustomError},
    fragment::{
        DiagnosticPosition, Fragment, FragmentType, NeutralLoss, PeptidePosition, SatelliteLabel,
    },
    helper_functions::{
        Characters, RangeExtension, RangeMaths, end_of_enclosure, explain_number_error, next_number,
    },
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    prelude::{Chemical, CompoundPeptidoformIon, MultiChemical, SequenceElement},
    quantities::Tolerance,
    sequence::{
        AminoAcid, Peptidoform, SemiAmbiguous, SimpleModification, SimpleModificationInner,
    },
    system::{Mass, MassOverCharge, OrderedMassOverCharge, e, isize::Charge, mz},
};

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
    let (left_range, isotopes) = parse_isotopes(line, left_range)?;
    let (left_range, adduct_type) = parse_adduct_type(line, left_range)?;
    let (left_range, charge) = parse_charge(line, left_range)?;
    let (left_range, deviation) = parse_deviation(line, left_range)?;
    let (left_range, confidence) = parse_confidence(line, left_range)?;
    Ok((
        left_range,
        PeakAnnotation {
            auxiliary,
            analyte_number,
            ion,
            neutral_losses,
            isotopes,
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
    isotopes: Vec<(i32, Isotope)>,
    charge: MolecularCharge,
    deviation: Option<Tolerance<OrderedMassOverCharge>>,
    confidence: Option<f64>,
}

impl PeakAnnotation {
    fn to_fragment(self, interpretation: CompoundPeptidoformIon) -> Fragment {
        // Get the peptidoform (assume no cross-linkers)
        let peptidoform = self.analyte_number.filter(|n| *n > 0).and_then(|n| {
            interpretation
                .peptidoform_ions()
                .get(n - 1)
                .and_then(|p| p.peptidoforms().first())
        });
        let (formula, ion) = match self.ion {
            IonType::Unknown(series) => (None, FragmentType::Unknown(series)),
            IonType::MainSeries(series, sub, ordinal, specific_interpretation) => {
                let specific_interpretation: Option<Peptidoform<crate::sequence::Linked>> =
                    specific_interpretation.map(Into::into);
                let peptidoform = specific_interpretation.as_ref().or(peptidoform);
                let sequence_length = peptidoform.map_or(0, Peptidoform::len);
                (
                    None,
                    match series {
                        b'a' => FragmentType::a(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    ordinal - 1,
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'b' => FragmentType::b(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    ordinal - 1,
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'c' => FragmentType::c(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    ordinal - 1,
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'd' => FragmentType::d(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    ordinal - 1,
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            peptidoform.map_or(AminoAcid::Unknown, |p| {
                                p.sequence()[ordinal - 1].aminoacid.aminoacid()
                            }),
                            0,
                            0,
                            match sub {
                                None => SatelliteLabel::None,
                                Some(b'a') => SatelliteLabel::A,
                                Some(b'b') => SatelliteLabel::B,
                                _ => unreachable!(),
                            },
                        ),
                        b'v' => FragmentType::v(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            peptidoform.map_or(AminoAcid::Unknown, |p| {
                                p.sequence()[ordinal - 1].aminoacid.aminoacid()
                            }),
                            0,
                            0,
                        ),
                        b'w' => FragmentType::w(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            peptidoform.map_or(AminoAcid::Unknown, |p| {
                                p.sequence()[ordinal - 1].aminoacid.aminoacid()
                            }),
                            0,
                            0,
                            match sub {
                                None => SatelliteLabel::None,
                                Some(b'a') => SatelliteLabel::A,
                                Some(b'b') => SatelliteLabel::B,
                                _ => unreachable!(),
                            },
                        ),
                        b'x' => FragmentType::x(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'y' => FragmentType::y(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'z' => FragmentType::z(
                            PeptidePosition {
                                sequence_index: crate::prelude::SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        _ => unreachable!(),
                    },
                )
            }
            IonType::Immonium(aa, m) => (
                Some(
                    aa.formulas().first().unwrap()
                        + m.as_ref().map(|m| m.formula()).unwrap_or_default()
                        - molecular_formula!(C 1 O 1),
                ),
                FragmentType::Immonium(
                    None, // TODO allow empty
                    m.map_or_else(
                        || SequenceElement::new(aa.into(), None),
                        |m| SequenceElement::new(aa.into(), None).with_simple_modification(m),
                    ),
                ),
            ),
            IonType::Internal(start, end) => {
                let sequence_length = self
                    .analyte_number
                    .filter(|n| *n > 0)
                    .and_then(|n| {
                        interpretation
                            .peptidoform_ions()
                            .get(n - 1)
                            .and_then(|p| p.peptidoforms().first())
                    })
                    .map_or(0, Peptidoform::len);
                (
                    None,
                    FragmentType::Internal(
                        None,
                        PeptidePosition::n(
                            crate::prelude::SequencePosition::Index(start - 1),
                            sequence_length,
                        ),
                        PeptidePosition::n(
                            crate::prelude::SequencePosition::Index(end - 1),
                            sequence_length,
                        ),
                    ),
                )
            }
            IonType::Named(_) => (None, FragmentType::Unknown(None)),
            IonType::Formula(formula) => (Some(formula), FragmentType::Unknown(None)),
            IonType::Precursor => (None, FragmentType::Precursor),
            IonType::Reporter(formula) => (
                Some(formula),
                FragmentType::Diagnostic(DiagnosticPosition::Reporter),
            ),
        };
        Fragment {
            formula,
            charge: self.charge.charge(),
            ion,
            peptidoform_ion_index: self.analyte_number.filter(|n| *n > 0).map(|n| n - 1),
            peptidoform_index: Some(0), // TODO what to do with cross-linked stuff
            neutral_loss: self.neutral_losses,
            deviation: self.deviation,
            confidence: self.confidence.map(Into::into),
            auxiliary: self.auxiliary,
        }
    }
}

impl std::fmt::Display for PeakAnnotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.auxiliary {
            write!(f, "&")?;
        }
        if let Some(number) = self.analyte_number {
            write!(f, "{number}@")?;
        }
        write!(f, "{}", self.ion)?;
        for loss in &self.neutral_losses {
            write!(f, "{loss}")?;
        }
        let charge = self.charge.charge().value;
        if self.charge != MolecularCharge::proton(charge) {
            write!(f, "[M")?;
            for (amount, option) in &self.charge.charge_carriers {
                write!(f, "{amount}{option}")?;
            }
            write!(f, "]")?;
        }
        if charge != 1 {
            write!(f, "^{charge}")?;
        }
        write!(f, "{:?}{:?}", self.deviation, self.confidence)
    }
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

impl std::fmt::Display for IonType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Unknown(n) => write!(f, "?{}", n.map_or_else(String::new, |n| n.to_string())),
            Self::MainSeries(series, sub, ordinal, interpretation) => {
                write!(
                    f,
                    "{}{}{ordinal}{}",
                    *series as char,
                    sub.map_or_else(String::new, |s| (s as char).to_string()),
                    interpretation
                        .as_ref()
                        .map_or_else(String::new, |i| format!("{{{i}}}"))
                )
            }
            Self::Immonium(aa, m) => {
                write!(
                    f,
                    "I{aa}{}",
                    m.clone().map_or_else(String::new, |m| format!("[{m}]"))
                )
            }
            Self::Internal(start, end) => write!(f, "m{start}:{end}"),
            Self::Named(name) => write!(f, "r[{name}]"),
            Self::Precursor => write!(f, "p"),
            Self::Reporter(formula) | Self::Formula(formula) => write!(f, "f{{{formula}}}"),
        }
    }
}

#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
enum Isotope {
    General,
    Average,
    Specific(Element, NonZeroU16),
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
            if line.as_bytes().get(num.0 + range.start).copied() != Some(b'@') {
                return Err(CustomError::error(
                    "Invalid mzPAF analyte number",
                    "The analyte number should be followed by an at sign '@'",
                    Context::line(None, line, num.0 + range.start, 1),
                ));
            }
            Ok((
                range.add_start(num.0 + 1),
                Some(num.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF analyte number",
                        format!("The analyte number number {}", explain_number_error(&err)),
                        Context::line(None, line, range.start, num.0),
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
            if let Some(ordinal) = next_number::<false, false, usize>(line, range.clone()) {
                let range = range.add_start(ordinal.0);
                let (end, interpretation) = if line.as_bytes().get(range.start_index()).copied()
                    == Some(b'{')
                {
                    if let Some(location) =
                        end_of_enclosure(line, range.start_index() + 1, b'{', b'}')
                    {
                        let interpretation = Peptidoform::pro_forma(
                            &line[range.start_index() + 1..location],
                            custom_database,
                        );
                        interpretation.and_then(|i| {
                            i.into_semi_ambiguous()
                                .ok_or_else(|| CustomError::error(
                                    "Invalid mzPAF interpretation", 
                                    "An mzPAF interpretation should be limited to `base-ProForma compliant` without any labile modifications", 
                                    Context::line_range(None, line, range.start_index()..location)))
                                .map(|i| (location + 1, Some(i)))
                        })?
                        // TODO: proper error handling and add checks to the length of the sequence
                    } else {
                        return Err(CustomError::error(
                            "Invalid mzPAF main series ion ordinal",
                            "The asserted interpretation should have a closed curly bracket, like '0@b2{LL}'",
                            Context::line(None, line, range.start_index(), 1),
                        ));
                    }
                } else {
                    (range.start_index(), None)
                };
                Ok((
                    end..range.end,
                    IonType::MainSeries(
                        c,
                        sub,
                        ordinal.2.map_err(|err| {
                            CustomError::error(
                                "Invalid mzPAF ion ordinal",
                                format!("The ordinal number {}", explain_number_error(&err)),
                                Context::line(
                                    None,
                                    line,
                                    range.start_index() - ordinal.0, // Maybe also offset for interpretation?
                                    ordinal.0,
                                ),
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
                let end = end_of_enclosure(line, range.start_index() + 3, b'[', b']').ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid mzPAF immonium modification",
                            "The square brackets are not closed",
                            Context::line(None, line, range.start_index(), 1),
                        )
                    },
                )?;
                let modification = &line[range.start_index() + 3..end];
                Some((
                    end - range.start_index(), // Length of mod + [ + ]
                    Ontology::Unimod
                        .find_name(modification, None)
                        .or_else(|| {
                            modification.parse::<f64>().ok().map(|n| {
                                std::sync::Arc::new(SimpleModificationInner::Mass(
                                    Mass::new::<crate::system::dalton>(n).into(),
                                ))
                            })
                        })
                        .ok_or_else(|| Ontology::Unimod.find_closest(modification, None))?,
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
            if line[range.clone()].chars().nth(first_ordinal.0 + 1) != Some(':') {
                return Err(CustomError::error(
                    "Invalid mzPAF internal ion ordinal separator",
                    "The internal ion ordinal separator should be a colon ':', like 'm4:6'",
                    Context::line(None, line, range.start_index() + 1 + first_ordinal.0, 1),
                ));
            }
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

            let (len, name) = if line[range.start_index() + 1..].starts_with('{') {
                let end = end_of_enclosure(line, range.start_index() + 2, b'{', b'}').ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid mzPAF named compound",
                            "The curly braces are not closed",
                            Context::line(None, line, range.start_index() + 1, 1),
                        )
                    },
                )?;
                Ok((
                    end - range.start_index(),
                    &line[range.start_index() + 2..end],
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
                    &line[range.clone()][first + range.start_index()
                        ..range.start_index() + last.0 + last.1.len_utf8()],
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
            let formula = MolecularFormula::from_pro_forma(
                line,
                formula_range.clone(),
                false,
                false,
                true,
                false,
            )?;

            Ok((
                range.add_start(3 + formula_range.len()),
                IonType::Formula(formula),
            ))
        }
        Some(b's') => Err(CustomError::error(
            "Unsupported feature",
            "SMILES strings are currently not supported in mzPAF definitions",
            Context::line(None, line, range.start, 1),
        )), // TODO: return as Formula
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
        let mut num_offset = 0;
        // Parse leading number to detect how many times this loss occured
        if let Some(num) = next_number::<false, false, u16>(line, range.add_start(1 + offset)) {
            amount = i32::from(num.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF neutral loss leading amount",
                    format!(
                        "The neutral loss amount number {}",
                        explain_number_error(&err)
                    ),
                    Context::line(
                        None,
                        line,
                        range.start_index() + 1 + offset,
                        range.start_index() + 1 + offset + num.0,
                    ),
                )
            })?);
            offset += num.0;
            num_offset = num.0;
        }

        println!("{}", &line[range.start_index() + 1 + offset..]);
        if line[range.start_index() + 1 + offset..].starts_with('i')
            || line[range.start_index() + 1 + offset..].starts_with("[M+")
        {
            return Ok((range.add_start(offset - num_offset), neutral_losses));
        }

        if line
            .as_bytes()
            .get(range.start_index() + 1 + offset)
            .copied()
            == Some(b'[')
        {
            let last = end_of_enclosure(line, range.start_index() + 2 + offset, b'[', b']')
                .ok_or_else(|| {
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
                neutral_losses.push(match c {
                    b'+' => NeutralLoss::Gain(formula * amount),
                    b'-' => NeutralLoss::Loss(formula * amount),
                    _ => unreachable!(),
                });
            } else if let Ok(formula) =
                MolecularFormula::from_pro_forma(line, first - 1..=last, false, false, true, false)
            {
                // Catches the case of a single isotope as formula
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
                .take_while(|(_, c)| c.is_ascii_alphanumeric() || *c == '[' || *c == ']')
                .last()
                .ok_or_else(|| {
                    CustomError::error(
                        "Invalid mzPAF",
                        "Empty neutral loss",
                        Context::line_range(None, line, first..),
                    )
                })?;
            let mut last = last.0 + last.1.len_utf8();
            if line[first..first + last].ends_with("[M") {
                last -= 2; // Detect any adduct types which might otherwise sneak in
            }
            let formula = MolecularFormula::from_pro_forma(
                line,
                first..first + last,
                false,
                false,
                true,
                false,
            )?;
            neutral_losses.push(match c {
                b'+' => NeutralLoss::Gain(formula * amount),
                b'-' => NeutralLoss::Loss(formula * amount),
                _ => unreachable!(),
            });
            offset += 1 + last;
        }
    }
    Ok((range.add_start(offset), neutral_losses))
}

fn parse_isotopes(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Vec<(i32, Isotope)>), CustomError> {
    let mut offset = 0;
    let mut isotopes = Vec::new();
    while let Some(c @ (b'-' | b'+')) = line.as_bytes().get(range.start_index() + offset).copied() {
        offset += 1;
        let mut amount = 1;
        // Parse leading number to detect how many times this isotope occurred
        if let Some(num) = next_number::<false, false, u16>(line, range.add_start(offset)) {
            amount = i32::from(num.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF isotope leading amount",
                    format!("The isotope amount number {}", explain_number_error(&err)),
                    Context::line(None, line, range.start_index() + offset, num.0),
                )
            })?);
            offset += num.0;
        }
        if c == b'-' {
            amount *= -1;
        }

        // Check if i
        if line.as_bytes().get(range.start_index() + offset).copied() != Some(b'i') {
            return Err(CustomError::error(
                "Invalid mzPAF isotope",
                "An isotope should be indicated with a lowercase 'i', eg '+i', '+5i', '+2iA', '+i13C'",
                Context::line(None, line, range.start_index() + offset, 1),
            ));
        }
        offset += 1;

        // Check if a specific isotope
        if let Some(num) = next_number::<false, false, NonZeroU16>(line, range.add_start(offset)) {
            let nucleon = NonZeroU16::from(num.2.map_err(|err| {
                CustomError::error(
                    "Invalid mzPAF isotope nucleon number",
                    format!("The nucleon number {}", explain_number_error(&err)),
                    Context::line(None, line, range.start_index() + offset, num.0),
                )
            })?);
            offset += num.0;

            let mut element = None;
            for (code, el) in ELEMENT_PARSE_LIST {
                if line[range.start + offset..].starts_with(code) {
                    element = Some(*el);
                    offset += code.len();
                    break;
                }
            }
            let element = element.ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF isotope element",
                    "No recognised element symbol was found",
                    Context::line(None, line, range.start_index() + offset, 1),
                )
            })?;
            if !element.is_valid(Some(nucleon)) {
                let ln = element.symbol().len();
                return Err(CustomError::error(
                    "Invalid mzPAF isotope",
                    format!(
                        "The nucleon number {nucleon} does not have a defined mass for {element}",
                    ),
                    Context::line(None, line, range.start_index() + offset - ln, ln),
                ));
            }
            isotopes.push((amount, Isotope::Specific(element, nucleon)));
        } else {
            // A or nothing
            if line.as_bytes().get(range.start_index() + offset).copied() == Some(b'A') {
                offset += 1;
                isotopes.push((amount, Isotope::Average));
            } else {
                isotopes.push((amount, Isotope::General));
            }
        }
    }
    Ok((range.add_start(offset), isotopes))
}

fn parse_adduct_type(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<MolecularCharge>), CustomError> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'[') {
        let closing =
            end_of_enclosure(line, range.start_index() + 1, b'[', b']').ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF adduct type",
                    "No closing bracket found for opening bracket of adduct type",
                    Context::line(None, line, range.start_index(), 1),
                )
            })?; // Excluding the ']' closing bracket
        println!("at: {}", &line[range.start..closing]);
        if line.as_bytes().get(range.start_index() + 1).copied() != Some(b'M') {
            return Err(CustomError::error(
                "Invalid mzPAF adduct type",
                "The adduct type should start with 'M', as in '[M+nA]'",
                Context::line(None, line, range.start_index() + 1, 1),
            ));
        }
        let mut carriers = Vec::new();
        let mut offset = 2; // '[M'
        while let Some(c @ (b'-' | b'+')) =
            line.as_bytes().get(range.start_index() + offset).copied()
        {
            offset += 1; // The sign
            if range.start_index() + offset >= closing {
                return Ok((
                    range.add_start(offset + 2),
                    Some(MolecularCharge::new(&carriers)),
                ));
            }
            let mut amount = 1;
            // Parse leading number to detect how many times this adduct occurred
            if let Some(num) = next_number::<false, false, u16>(line, range.add_start(offset)) {
                amount = i32::from(num.2.map_err(|err| {
                    CustomError::error(
                        "Invalid mzPAF adduct leading amount",
                        format!("The adduct amount number {}", explain_number_error(&err)),
                        Context::line(
                            None,
                            line,
                            range.start_index() + offset,
                            range.start_index() + offset + num.0,
                        ),
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
                .take_while(|(_, c)| c.is_ascii_alphanumeric() || *c == '[' || *c == ']')
                .last()
                .map_or(0, |last| last.0 + last.1.len_utf8())
                .min(closing - first); // Prevent the closing bracket from being used in an isotope
            println!("fo: {}", &line[first..first + last]);
            let formula = MolecularFormula::from_pro_forma(
                line,
                first..first + last,
                false,
                false,
                true,
                false,
            )?;
            carriers.push((amount as isize, formula));
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

/// Parse mzPAF charge, eg `^2` `^-1`
/// # Errors
/// If there is nu number after the caret, or if the number is invalid (outside of range and the like).
fn parse_charge(line: &str, range: Range<usize>) -> Result<(Range<usize>, Charge), CustomError> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'^') {
        let charge =
            next_number::<true, false, isize>(line, range.add_start(1_usize)).ok_or_else(|| {
                CustomError::error(
                    "Invalid mzPAF charge",
                    "The number after the charge symbol should be present, eg '^2'.",
                    Context::line(None, line, range.start_index(), 1),
                )
            })?;
        Ok((
            range.add_start(charge.0 + 1),
            Charge::new::<e>(
                if charge.1 { -1 } else { 1 }
                    * charge.2.map_err(|err| {
                        CustomError::error(
                            "Invalid mzPAF charge",
                            format!("The charge number {}", explain_number_error(&err)),
                            Context::line(None, line, range.start_index() + 1, charge.0),
                        )
                    })?,
            ),
        ))
    } else {
        Ok((range, Charge::new::<e>(1)))
    }
}

/// Parse a mzPAF deviation, either a ppm or mz deviation.
/// # Errors
/// When the deviation is not '<number>' or '<number>ppm'.
fn parse_deviation(
    line: &str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<Tolerance<OrderedMassOverCharge>>), CustomError> {
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
        if line[range.start_index() + 1 + number.0..].starts_with("ppm") {
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

// TODO: update list
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
            use itertools::Itertools;
            let res = $crate::fragment::parse_mzpaf($case, None);
            match res {
                Err(err) => {
                    println!("Failed: '{}'", $case);
                    println!("{err}");
                    panic!("Failed test")
                }
                Ok(res) => {
                    let back = res
                        .into_iter()
                        .map(|a| a.to_fragment(CompoundPeptidoformIon::default()).to_mzPAF())
                        .join(",");
                    let res = $crate::fragment::parse_mzpaf(&back, None);
                    if let Err(err) = res {
                        println!("Failed: '{}' was exported as '{back}'", $case);
                        println!("{err}");
                        panic!("Failed test")
                    }
                    // let back = res.as_ref().unwrap().to_string();
                    // let res_back = $crate::fragment::parse_mzpaf(&back, None);
                    // assert_eq!(res, res_back, "{} != {back}", $case);
                }
            };
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

mzpaf_test!("b2-H2O/3.2ppm,b4-H2O^2/3.2ppm", spec_positive_1);
mzpaf_test!("b2-H2O/3.2ppm*0.75,b4-H2O^2/3.2ppm*0.25", spec_positive_2);
mzpaf_test!("1@y12/0.13,2@b9-NH3/0.23", spec_positive_3);
mzpaf_test!("0@y1{K}", spec_positive_4);
mzpaf_test!("0@y1{K}-NH3", spec_positive_5);
mzpaf_test!("y1/-1.4ppm", spec_positive_6);
mzpaf_test!("y1/-0.0002", spec_positive_7);
mzpaf_test!("y4-H2O+2i[M+H+Na]^2", spec_positive_8);
mzpaf_test!("?", spec_positive_9);
mzpaf_test!("?^3", spec_positive_10);
mzpaf_test!("?+2i^4", spec_positive_11);
mzpaf_test!("?17", spec_positive_12);
mzpaf_test!("?17+i/1.45ppm", spec_positive_13);
mzpaf_test!("?17-H2O/-0.87ppm", spec_positive_14);
mzpaf_test!("0@b2{LL}", spec_positive_15);
mzpaf_test!("0@y1{K}", spec_positive_16);
mzpaf_test!("0@b2{LC[Carbamidomethyl]}", spec_positive_17);
mzpaf_test!("0@b1{[Acetyl]-M}", spec_positive_18);
mzpaf_test!("0@y4{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19a);
mzpaf_test!("0@y44{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19b);
mzpaf_test!("0@y444{M[Oxidation]ACK}-CH4OS[M+H+Na]^2", spec_positive_19c);
mzpaf_test!("m3:6", spec_positive_20);
mzpaf_test!("b3-C2H3NO", spec_positive_21);
mzpaf_test!("m3:6-CO", spec_positive_22);
mzpaf_test!("m3:6-CO-H2O^2", spec_positive_23);
mzpaf_test!("m3:4/1.1ppm,m4:5/1.1ppm", spec_positive24);
mzpaf_test!("m3:5", spec_positive_25);
mzpaf_test!("IY", spec_positive_26);
mzpaf_test!("IH", spec_positive_27);
mzpaf_test!("IL-CH2", spec_positive_28);
mzpaf_test!("IC[Carbamidomethyl]", spec_positive_29);
mzpaf_test!("IY[Phospho]", spec_positive_30);
mzpaf_test!("IC[+58.005]", spec_positive_31);
mzpaf_test!("p^2", spec_positive_32a);
mzpaf_test!("p^-2", spec_positive_32b);
mzpaf_test!("p-H3PO4^2", spec_positive_33);
mzpaf_test!("p^4", spec_positive_34);
mzpaf_test!("p+H^3", spec_positive_35);
mzpaf_test!("p^3", spec_positive_36);
mzpaf_test!("p+2H^2", spec_positive_37);
mzpaf_test!("p^2", spec_positive_38);
mzpaf_test!("p+H^2", spec_positive_39);
mzpaf_test!("p+3H", spec_positive_40);
mzpaf_test!("p+2H", spec_positive_41);
mzpaf_test!("p+H", spec_positive_42);
mzpaf_test!("p", spec_positive_43);
mzpaf_test!("r[TMT127N]", spec_positive_44);
mzpaf_test!("r[iTRAQ114]", spec_positive_45);
mzpaf_test!("r[TMT6plex]", spec_positive_46);
mzpaf_test!("r[Hex]", spec_positive_47);
mzpaf_test!("r[Adenosine]", spec_positive_48);
mzpaf_test!("0@_{Urocanic Acid}", spec_positive_49);
mzpaf_test!("f{C13H9}/-0.55ppm", spec_positive_50);
mzpaf_test!("f{C12H9N}/0.06ppm", spec_positive_51);
mzpaf_test!("f{C13H9N}/-2.01ppm", spec_positive_52);
mzpaf_test!("f{C13H10N}/-0.11ppm", spec_positive_53);
mzpaf_test!("f{C13H11N}/-0.09ppm", spec_positive_54);
mzpaf_test!("f{C13H12N}/0.26ppm", spec_positive_55);
mzpaf_test!("f{C14H10N}/0.19ppm", spec_positive_56);
mzpaf_test!("f{C14H11N}/0.45ppm", spec_positive_57);
mzpaf_test!("f{C14H10NO}/0.03ppm", spec_positive_58);
mzpaf_test!("f{C16H22O}+i^3", spec_positive_59);
mzpaf_test!("f{C15[13C1]H22O}^3", spec_positive_60);
mzpaf_test!("s{CN=C=O}[M+H]/-0.55ppm", spec_positive_61);
mzpaf_test!("s{COc(c1)cccc1C#N}[M+H+Na]^2/1.29ppm", spec_positive_62);
mzpaf_test!("p-[Hex]", spec_positive_63);
mzpaf_test!("y2+CO-H2O", spec_positive_64);
mzpaf_test!("y2-H2O-NH3", spec_positive_65);
mzpaf_test!("p-[TMT6plex]-2H2O-HPO3", spec_positive_66);
mzpaf_test!("p-2[iTRAQ115]", spec_positive_67);
mzpaf_test!("p-[iTRAQ116]-CO-H2O-HPO3", spec_positive_68);
mzpaf_test!("y2-[2H1]-NH3", spec_positive_69);
mzpaf_test!("y5-H2[18O1][M+Na]", spec_positive_70);
mzpaf_test!("y12+i", spec_positive_71);
mzpaf_test!("y12+i13C", spec_positive_72);
mzpaf_test!("y12+i15N", spec_positive_73);
mzpaf_test!("y12+2i13C+i15N", spec_positive_74);
mzpaf_test!("y12+iA", spec_positive_75);
mzpaf_test!("y12+2iA", spec_positive_76);
mzpaf_test!("y4[M+Na]", spec_positive_77);
mzpaf_test!("y5-H2O[M+H+Na]^2", spec_positive_78);
mzpaf_test!("y6[M+[2H2]]", spec_positive_79);
mzpaf_test!("y5[M+[15N1]H4]", spec_positive_80);
mzpaf_test!("&1@y7/-0.002", spec_positive_81);
mzpaf_test!("&y7/-0.001", spec_positive_82);
mzpaf_test!("y7/0.000*0.95", spec_positive_83);
mzpaf_test!("&y7/0.001", spec_positive_84);
mzpaf_test!("&y7/0.002", spec_positive_85);
mzpaf_test!("b6-H2O/-0.005,&y7/0.003", spec_positive_86);
mzpaf_test!("y12-H2O^2/7.4ppm*0.70", spec_positive_87);
mzpaf_test!("y12/3.4ppm*0.85,b9-NH3/5.2ppm*0.05", spec_positive_88);

mzpaf_test!(ne r"0@y4{Mar^2dation]ACK}-CH4OS", fuzz_0);
mzpaf_test!(ne r"y4-4", fuzz_1);
mzpaf_test!(ne r"0@y4{Mar^2dation]ACK", fuzz_2);
mzpaf_test!(ne r"0@y4{", fuzz_3);
mzpaf_test!(ne r"IM[", fuzz_4);
mzpaf_test!(ne r"IM[Carboxymethyl", fuzz_5);

mzpaf_test!(ne r"r[", test_0);
mzpaf_test!(ne r"y5-4552H2O555557O-42[17O1]SSSSSSSSSSSSSSSSJSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS555VSS42[17O1LSSSSSO55555V553[5V5555555SSS422H2355355V555g555555M+ ]N", test_1);
mzpaf_test!(ne r"x5-92[17O-144444443.2pp/3.2O/3.2MpmO*0.75,b4-H2O^2/3.2ppm*0 on]I-92HInOT1+][", test_2);
mzpaf_test!(ne r"x5-9B92[17O-1414444444[ion]HII,_4[17O-141444444Q[[ion]IImO]Iion]+][M", test_3);
mzpaf_test!(ne r"y5-O2[2H2O55O-4O[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS55555V55555V55555552H2O5555VSS42[17O1]SSSSSO55555V55555V55555552H2O55355V555555555M+I", test_4);
mzpaf_test!(ne r"x5-9292[17O-1414444    ion]IIICH4[17O- 14144444MUUUUUUU]IImO]Iion]+][M", test_5);
mzpaf_test!(ne r"f{}M", test_6);
mzpaf_test!(ne r"f{", test_7);
mzpaf_test!(ne r"y5-O2[2H2O55O-4O[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS55555V55555V555555552H2O5555VSS42[17O1]SSSSSO55555V55555V55555552H2O5*355V5555555555M+I", test_8);
mzpaf_test!(ne r"r[]RIa", test_9);
mzpaf_test!(ne r"x5-92B22B2[17O-1414444444[[17O-1414444444[ion]III    17O0@y4{M-1414444444[[ionon]+][M", test_10);
mzpaf_test!(ne r"x5-9292[17O41414444Pidomethyl],p-H3PO4^2,0@2[/3.0.75,b4-H2O^2/3]+][M", test_11);
mzpaf_test!(ne r"x5-9292[17O-1414464444]+][M", test_12);
mzpaf_test!(ne r"x5-929IIICH4[17O-1414444444[[ion]IImO]1414444444[[ion]IImO][889OO][M", test_13);
mzpaf_test!(ne r"y5-O2[2H2O55O-4O[17O1]SSSSSSSSSSS4SSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS55555V-5555V555552H2O55O-4O[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS    4242[17O1]SS555V42[17O1]SS95555V55555V55555552HSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS5J555V-5555V555552H2O55O-4O[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS15H1]H475bdatm*0.2O/3.2ppm*0.7O^2/3.2ppm*0.2O/3.2ppq*0.]-1=]|InjdcyhbarbbV-5555V555552H2O55O-4O[1 7O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS    4242[17O1]SS555V42[17O1]SS95555V55555V55555552HSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS55555V-5555V555552H2O55O-4O[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS15H1]H475bdatm*0.2O/3.2ppm*0.7O^2/3.2ppm*0.2O/3.2@y4{M/0010000000[-000,00b4-,-,-.dpmxppq*0.]-1=]|Injdcyhbarbbmidomethyl],pM7O1]SSSSS42[1O2[2H2O55O-4O[17O1]SSSSSS2O5555VSS42[17O1]SSSSSC55555V55555V55555552H2O5535552H2O5555VSS42[17O1]SSSSSO55555V55555V55555552H2O55355V5555555555Mmidomethyl],pM7O1]SSSSS42[1O2[2H2O55O-4O[17O1]SSSSSS2O5555VSS42[17O1]SSSSSO55555V55555V55555552H2O5535552H2O5555VSS42[17O1]SSSSSO55555V55555V55555552H2O55355V5555555555M+I", test_14);
mzpaf_test!(ne r"1@r[daAAKuy8s<745aG>Ca``vmw[L}.", test_15);
mzpaf_test!(ne r"b2-8[15N1191111111]H\15N115555555555555", test_16);
mzpaf_test!(ne r"x5-927O[17O-14444444[ion]I-92II}-T1+][", test_17);
mzpaf_test!(ne r"x5-9S[17O444444444[ion]III7:31+][M", test_18);
mzpaf_test!(ne r"x5-92[17O-141444444  ion]IIICH4OS[M+]^}-T1+][444444[imO]III}-T1+][_", test_19);
mzpaf_test!(ne r"x5-92I-92I1111[17O114444444[ion]I-92I111111111111117O-14444444[ion]II}][", test_20);
mzpaf_test!(ne r"x5-92[17O-1414444444Zion]IIIIIIIIIIIIIIIIIIICf4OS[M+]^}-T1k]S444444[imN]III}-T1+][M", test_21);
mzpaf_test!(ne r"x5-92[17O-111111111ffftfftfIIII}-T1+][M", test_22);
mzpaf_test!(ne r"0@r[Oxidati   bami3:4/1.1py4{ ", test_23);
mzpaf_test!(ne r"y5-617O111112[17O11111111][M+Na],y5[M?", test_24);
mzpaf_test!(ne r"x5-929B[17O-141ICH4[17O-1414444944[[O-1414444944[[ion]IImS]Iion]+][M", test_25);
mzpaf_test!(ne r"y5-9217O1111111+][M[17O1111111+][M", test_26);
mzpaf_test!(ne r"x5-92[17O-1414444444Iiooonn ׾A,dgmethyl],ooooarbamdgmehoooooooooooooooooooon]IIICf4OS[M+]^}-T1+][444444[imO]III}-T1+][M", test_27);
mzpaf_test!(ne r"x5-9292[17O-1414444]ioni\n]IIICH4[17O-14144;4444[[ion]IImO]Iion]+][M", test_28);
mzpaf_test!(ne r"y5-72[17O111111111#11111111111111111111111]|M/Na],y5[M+", test_29);
mzpaf_test!(ne r"x5-922IIII[17O-14444444[ion]I-92IIIIIIIIIIIIIIIIIIII0@y4{<333@>Ac_]ion]VDK}-IIII]@y4{<[dxic_]i}}}}}}}}on]ADK[", test_30);
mzpaf_test!(ne r"x5-9N[17O-1414444444[ion]IIICH4OS[M+]E}-T1+][54T1+][M", test_31);
mzpaf_test!(ne r"b22-H2O-33333O3333-H2O-33333O3333333333-Hpm*0.2", test_32);
mzpaf_test!(ne r"f{}fioft   -[M", test_33);
mzpaf_test!(ne r"x5-9S[17O-1414444444[ion]IIICV4OS[M+]^}-T1+][444444[imO]III}-T1+][M", test_34);
mzpaf_test!(ne r"x5-9292[17O-1444444ion]III}-T1+][M", test_35);
mzpaf_test!(ne r"r[]xidation]AAAAAAAAAAAAAAAAAACK}", test_36);
mzpaf_test!(ne r"f{", test_37);
mzpaf_test!(ne r"r[", test_38);
mzpaf_test!(ne r"f{}-T1-@.1----[", test_39);
mzpaf_test!(ne r"x5-9IIIC2[17O-1414444444[ion]IIICV4OS[M+]^@-T1+][444444[imO]III}-T1+][M", test_40);
mzpaf_test!(ne r"f{", test_41);
mzpaf_test!(ne r"x5-92[17O4444444 11+]-992[17O4444444 114I7:31+][M", test_42);
mzpaf_test!(ne r"r[", test_43);
mzpaf_test!(ne r"y5-927O1111119278+1111i17O-92i17O11M", test_44);
mzpaf_test!(ne r"x5-9292[17O-1414444bamiddmethyl],_{444[ion]IIICH4[17O-1414444444[[ionRIImO]Iion]+][M", test_45);
mzpaf_test!(ne r"b2-55555V5555555!555555555V555555555555", test_46);
mzpaf_test!(ne r"x5-92[17O444444444 11+]-92 7O44[ion]II71+][M", test_47);
mzpaf_test!(ne r"x5-92[17O-1414442938_ops]IIICV4OS[M+]^}-T1+][444447\pkT`SUV|.X0/__H", test_48);
mzpaf_test!(ne r"r[].iij thyl].`l", test_49);
mzpaf_test!(ne r"0@y4^1,0@f{}-[cm,m4MQOxidatio-n] K}", test_50);
mzpaf_test!(ne r"x5-9292[17O-1414444mO]Iion]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_51);
mzpaf_test!(ne r"x5-92[17O-64444444dion]I-9MMMMMMMMMMMMMMǐ-T1][", test_52);
mzpaf_test!(ne r"x5-929IIICH4[17O-141444442[17O-14144444R4[ion]IIICH4[17O-1414444444[[ion---mO]Ihon]+][M", test_53);
mzpaf_test!(ne r"0@y4^1,00@f{}-cO", test_54);
mzpaf_test!(ne r"m3:4/1.1ppm,m4:5/-1.1ppm,r[eȆr[eȆ,w,m4555-.1opm", test_55);
mzpaf_test!(ne r"r[", test_56);
mzpaf_test!(ne r"r[]mUUUUUmidomVthyl.2pp.2ppm*0m*0.2", test_57);
mzpaf_test!(ne r"x5-9217O4444444 11+]-92[27O44[ion]IIII7:31+][M", test_58);
mzpaf_test!(ne r"0@r[TM36Y", test_59);
mzpaf_test!(ne r"x5-9U92[17O-1414444444[ion]IIICH4[17O-141444443-i3Larbamidomeon92[17O-1414444444[ion]I[M", test_60);
mzpaf_test!(ne r"0@y4^1,0@f{", test_61);
mzpaf_test!(ne r"x5-92[17O44444440@y4{MBOxidationyA@ }-CH231+][M", test_62);
mzpaf_test!(ne r"x5-92[17O4444444 1+]-9292[17O4444464 1J+]-92[1[1744Nion]III7:31+][[0@y4{GG[Oxid[[[[[[[[[[[[[M", test_63);
mzpaf_test!(ne r"x5-92U555555555[17O4444444 11]-92[17O44[ion]III7:31+][M", test_64);
mzpaf_test!(ne r"0@y4{M[Oxidation]ACK}-CH4OSH3PO4^2,0@f{", test_65);
mzpaf_test!(ne r"m3:511+88[15N1191111111]H[15N11111]H[15N1191111.11]H4]H4]14]H4]", test_66);
mzpaf_test!(ne r"b2-H2O/3.2ppm*0.75,b4-OOOOO4-8[15N1191111111]H[1NNH2O^2H2", test_67);
mzpaf_test!(ne r"0@r[eȆr[eȆ,w,_{Urrdr4OSH3PȆ@_{", test_68);
mzpaf_test!(ne r"b2-0000000W02-000000000000000001000W020000000-000000000W00000W00000000000%0000000000000000000000000004000000000000000000000000004000000000", test_69);
mzpaf_test!(ne r"x5-93I333333333U333I333333333U2[17O-14144PPPPPPPPPPPPPPPP:mpPPP[i4444P4on]IIIC醆lUrrV4][M", test_70);
mzpaf_test!(ne r"r[],_{Ue~ȆlrrdS", test_71);
mzpaf_test!(ne r"m3:111+88[15N1191111111]H[15N11911115[M]H4]H4]14]H4]", test_72);
mzpaf_test!(ne r"r[", test_73);
mzpaf_test!(ne r"f{}K", test_75);
mzpaf_test!(ne r"x5-92[17O-144444444pion]I-92IǏǐ-T1+][", test_76);
mzpaf_test!(ne r"c3-iA-iA,r[NiA-iA,c3-iA,22-iA,22", test_77);
mzpaf_test!(ne r"y5-1211H11111111[17O-111H11111111111][M,}-C[M+", test_78);
mzpaf_test!(ne r"f{}", test_79);
mzpaf_test!(ne r"r[", test_80);
mzpaf_test!(ne r"y5-4289U26326223B397O-42Z172[2H2O55555O1]SSSSSSSSSSSSSSSSSV55555+555557O-42[17O1]SSSSSSSidation]ACX}31]StioCarbamidomQtiyl],_{Urocanid Ac555M5555M+Im]N", test_81);
mzpaf_test!(ne r"b2-55555P5555555H555555 55555555552-55555I 55555G55555", test_82);
mzpaf_test!(ne r"y5-425V5555555[1H2O555555557O-42[17O1]SSSSSSSSSSSSSSSSSSSSSS[M+Im]N+Na],", test_83);
mzpaf_test!(ne r"x5-9292[17O-1414444]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_84);
mzpaf_test!(ne r"x5-9297[17O-1414444444on]+]IIICH4[17O-14]III  4[17O-141444@44[[ion]IImO]Ii[ion]lM", test_85);
mzpaf_test!(ne r"y5-42S1H2O52[17O1]SSSSSSSSSSS-45555V55555555   SSfSSSnS[M+Im]N+Narrrrr],", test_86);
mzpaf_test!(ne r"y5-42[1H255555557]5557]-42[155557]-42[17O1]S-42[15557]-42[17O1]SS", test_88);
mzpaf_test!(ne r"y5-H2[17O-111H111111114{<R[M,}-C[MMMMMMMidation]ACK}-CMMMMMMM+", test_89);
mzpaf_test!(ne r"y5-2[17O1111111111?11111][M,}-C[M++", test_90);
mzpaf_test!(ne r"y5-41[1H2O55555V55557]-42[17O555555557]-42[17OQ]SS", test_91);
mzpaf_test!(ne r"m3:4/1.1ppm,p-H3PO4^2,0@f{}-CH3PO4^2,0@f{U}[R:9/0{GGGGGGGGGGGGGGGGGGGGGGGGGWOxidathooo-MMMMid}", test_93);
mzpaf_test!(ne r"y5-4834F764424832=678726983693229733333333333333331[17O-4HHHHPHHHHHHHHHHM,yH5-4834F764424832=6787M,y5", test_94);
mzpaf_test!(ne r"x5-9U92[17O-1414444444[ion]IIICH4[17OI1414444444[[ion]IImO]Iion]+][M", test_95);
mzpaf_test!(ne r"y5-42S42[17O1]SSSS1O55555V555555555[2H5555V555555555557O-42[1777777777777777SSSSSSSS42[17O1]SSSSSOCAC5555V25C-[NUjwe]gyhd	Ebid}55SSSSOCAC5555V555555M+Im]N", test_96);
mzpaf_test!(ne r"y5-42S-45555V5555555[1H2O52[17O1]SSSSSSSSSSS-45555V5555555[1H2O52[17O1]SSSSSSSSSSS-45555V555555555557SSfSSSSS[M+Im]N+Na],", test_97);
mzpaf_test!(ne r"y5-42[  1H222222222 ,y}-00 hyl4 A5ethyl],o-c3PO4^2,0@_{UroM+Nb],", test_99);
mzpaf_test!(ne r"f{", test_100);
mzpaf_test!(ne r"y5-H22[17O+2[17O+1111101111)1111][M,", test_101);
mzpaf_test!(ne r"x5-9292[17O-1414425299^hhn]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_102);
mzpaf_test!(ne r"y5-17O11O1111111111@11111][M,}-C[M+", test_103);
mzpaf_test!(ne r"y5-42[17O-422222244\\`S/Qd]M+Na],y5", test_104);
mzpaf_test!(ne r"x5-9292[17O-141444444eeeeeeeeeen]IImO]Iion]+][M", test_105);
mzpaf_test!(ne r"x5-92[17O-44414444y2{<in>_Xǐ-T1+][", test_106);
mzpaf_test!(ne r"m3:4/1.1ppm,p-H2222O^2/1.1ppm,p-H3PO4^2,0@f{", test_107);
mzpaf_test!(ne r"0@r[eȆr[e((((}2", test_108);
mzpaf_test!(ne r"r[", test_109);
mzpaf_test!(ne r"0@r[idat^3:_BS|FKP", test_110);
mzpaf_test!(ne r"m8:9/0,r[i !724*./111111111111111111111.1/0.0ppm", test_111);
mzpaf_test!(ne r"m3:4/1.1ppm,p-H3PO4^2,0@f{}", test_112);
mzpaf_test!(ne r"x5-92[17O4444444 ]-92[17O44114444 11+]-92m17O44[ion]	II7:31+]-92[17O44[ion]III7:31+][M", test_113);
mzpaf_test!(ne r"x5-929F[17O-1414444444[ion]IIICH4[17O-1uir]A[#xladru ]A[#xleLovir]A", test_114);
mzpaf_test!(ne r"r[]A]AC[alrbCadpmel],_{Uromethyl]K^]>JiBat}", test_115);
mzpaf_test!(ne r"x5-9292[17O-1414444yyyyyyyyyyyyyyyyyyyyyyyidomethyl]._~Rxuyyyyyyyyion]IImO]Iion]+][M", test_116);
mzpaf_test!(ne r"x5-9292[17O-1414476368\rhn]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_117);
mzpaf_test!(ne r"x5-9217O4444444171+]-93[ 1O44[ion]1+III7:31+][M", test_118);
mzpaf_test!(ne r"0@y4-2223O224222222223O22422222[222223O224222222223O2242222222222222ـG6Gatinssatins]AKK", test_119);
mzpaf_test!(ne r"b2-2[18O1111111111][M[ 18O1][M2[/34-H2O^2/3.2ppm*0.25", test_120);
mzpaf_test!(ne r"x5-9292[17O-1417444]]]]]444[ionBICICH4[+][M", test_121);
mzpaf_test!(ne r"0@r[e~ᆏ\,_{Ur	k  n   y ,A", test_122);
mzpaf_test!(ne r"b2-H2O555555+55W5555555555555d555", test_123);
mzpaf_test!(ne r"r[", test_124);
mzpaf_test!(ne r"b2-5555U5555555W5555/55555555d555", test_125);
mzpaf_test!(ne r"r[]", test_126);
mzpaf_test!(ne r"b2-42O555555555W5555555555555d555", test_127);
mzpaf_test!(ne r"x5-92[17O-1414442698_us/]IIICH4YS[M+]^}-T1+][444446^rlK`THQ~.N0/__V", test_128);
mzpaf_test!(ne r"x5-9292[17O-1414444n]IIICH4[17O-1414444444[[ion]II444[ion]*IICH4[17O-1414444444[[ion]IImM\Iion]+][M", test_129);
mzpaf_test!(ne r"x5-92N2[17O-114444444[[ion]IImO]IionO-1414444444[[ion]IImO]Iion]+][M", test_130);
mzpaf_test!(ne r"x5-927O27O-144O[17O-1444444444[ion]I", test_131);
mzpaf_test!(ne r"x5-92[17O-144442222$22222IǏǐ-T1+][", test_132);
mzpaf_test!(ne r"x5-9O92[17O-1414444444[kon]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_133);
mzpaf_test!(ne r"f{", test_134);
mzpaf_test!(ne r"x5-9C92[17O-1414444444[ion]IIICH4[17O-1414444444[[ionZIImO]Iion]+][M", test_135);
mzpaf_test!(ne r"x5-92[17O714444444\}on]I17O-14444444ji@ @ JIT1+]YM", test_136);
mzpaf_test!(ne r"f{} ,y3333333cdqidome/+++++++++++++oM[Ox7ppppm,m2:4/0.0ppm", test_137);
mzpaf_test!(ne r"r[", test_139);
mzpaf_test!(ne r"0@y44-1111C1111111111(1", test_140);
mzpaf_test!(ne r"x5-92[17O4444444 11+]-92[17O444444444 11+]-92[17O44[ionLII[ion]III7:31+][M", test_141);
mzpaf_test!(ne r"f{}{BK}czscx&`S", test_142);
mzpaf_test!(ne r"x9-47B3[17O-1706397377[izp]IIICH5[17O-1702895828[[tnp]YWrG]Suwp]/]^R", test_143);
mzpaf_test!(ne r"x5-927O4434444 11+]-4[i[17O4434444 11+]-4[i ]III7:31+][M", test_144);
mzpaf_test!(ne r"x5-9292ICH4[17O-1414444444[[ion]IImO]Ii4944444[M", test_145);
mzpaf_test!(ne r"x5-929O7O4444444 4444 11+11+]-92[17O44[ion]III7:C1+][M", test_146);
mzpaf_test!(ne r"x5-9292[17O-1414497377^lzn]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_148);
mzpaf_test!(ne r"0@r[eȆ,w,y4{VVVVK}-i0Da}-C", test_149);
mzpaf_test!(ne r"0@r[Tm782TB>a_^xqhKL 	 }-", test_150);
mzpaf_test!(ne r"x5-92[17O44444555H2O/3.2ppm*0.'5,b]III7:31+][M", test_151);
mzpaf_test!(ne r"r[", test_152);
mzpaf_test!(ne r"r[]-_}VpubeyhdCarbair}", test_153);
mzpaf_test!(ne r"x5-92[17O-24444444lllllllionP 2HISOT1+][", test_154);
mzpaf_test!(ne r"1@f{},D#7777777777777777777atioYYAhvgM[O", test_155);
mzpaf_test!(ne r"11@f{", test_156);
mzpaf_test!(ne r"0@r[C  ardddddd iC", test_157);
mzpaf_test!(ne r"b2-0000055I55000055I555-977777777-e45pmM0.5-977777777-e45pmM0.5", test_158);
mzpaf_test!(ne r"x5-92[17O-111111111Ion]I-9-H1+][", test_159);
mzpaf_test!(ne r"x5-9292ICH4[17O-1414444444[[innH4[17O-14  idomethyl]._~Qktee14444444[[inn]IImO]I  n +][M", test_160);
mzpaf_test!(ne r"r[", test_161);
mzpaf_test!(ne r"x5-92[17O514444444~ion]I-92IIn-T1+][", test_162);
mzpaf_test!(ne r"f{}KB", test_163);
mzpaf_test!(ne r"f{", test_164);
mzpaf_test!(ne r"y5-927O1111111+11127O1111111+1111i17O-91i17O-92i17O11M", test_165);
mzpaf_test!(ne r"x5-9292[17O-1414444c_]i   d IICH4[ O-141444A44@ {M[Oxidation]ACK}-", test_166);
mzpaf_test!(ne r"x5-917O44444442[17O4444444 11+]-92[17 114 11+]-92[17 11+O44[ion]III7:+O44[ion]III7:31+][M", test_167);
mzpaf_test!(ne r"b2-H2O/3.2ppm*0.75,b4-H2O^2/2ppm*0.75,r[tb4-H2O^2/33.2ppm*0", test_168);
mzpaf_test!(ne r"r[]+H", test_169);
mzpaf_test!(ne r"x5-92[17O-1444444444   ]I9217O-144444N[i];IT}]II}][M", test_170);
mzpaf_test!(ne r"r[", test_171);
mzpaf_test!(ne r"b2-H2O-333S333333333 3132ppm*0.25", test_172);
mzpaf_test!(ne r"x5-9292[17O-1414444? ion]IIICH4[17O-141O444444[[ion]IImO]Iion]+][M", test_173);
mzpaf_test!(ne r"1@r[M+r[M +H+RaH+RBS8", test_174);
mzpaf_test!(ne r"x5-929NO0000P00000000000000000O04[17O-1414444444[[ion]IImO]Iion]+][M", test_175);
mzpaf_test!(ne r"x5-929O[17O-1414544444[ion]14444444[ion]IIICH4[17O-IIICH4[17O-1414444%44[[ion]I ]B-014-[MO[M", test_176);
mzpaf_test!(ne r"x5-9292[17O-1414444444ion]IIICH4[17O-1414444444[[ion]IImO]Iion\+][M", test_177);
mzpaf_test!(ne r"x5-92U2[17O-1414444444<4444444444444444444444444[[io~]IImO]Iion]+][M", test_178);
mzpaf_test!(ne r"x5-9O92[27O-1414444445[ion]IIICH4[1444444[[ion]IImO4]I][Mion]+][Mon]+][M", test_179);
mzpaf_test!(ne r"x5-9292[17O-1414444]IImO]444[ion]IIICH4[17O-1414444444[[ion]IImO]Iion]+][M", test_180);
mzpaf_test!(ne r"f{}", test_181);
mzpaf_test!(ne r"f{", test_182);
mzpaf_test!(ne r"f{", test_183);
mzpaf_test!(ne r"x5-9292[17O5141444bWWCiddde4444tion]IIICH4[1WO17O5414444444[[ioation]ACCCCCCCCCCCCCCCCCCCCCCCCM", test_184);
mzpaf_test!(ne r"x5-92[17O-1414449985_lwg]IIICV4OS[M+]^}-T1+][444445_wxS^YKQ~,Y0/[]Q", test_185);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N-3333H3333333333:4,y544+[15N1]yp]", test_186);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N-33315N1333333333:4,y544+[15N1]yp]", test_187);
mzpaf_test!(ne r"m3:4,y544+[15N1]+CC[15N-333O33333333[15N1[Oxidation]FOOOOOOOO-H2]y", test_188);
mzpaf_test!(ne r"y5-9O701111111+1111i17O-92 7O11M", test_189);
mzpaf_test!(ne r"m3:4,y544+CC[15N-33333CCC[15N-333333333m333333:4,y544+[15N1]yp]", test_190);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N-33333CC[15Nm3333333UHHGHHHN4,y544+[15N1]yp]", test_191);
mzpaf_test!(ne r"f{}Kf{BKB}", test_192);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N-333[15N133333333:4,y544+33333]333[15N133333333:4yp]", test_193);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N73333333y544+15N733333333$33KS[M+]", test_194);
mzpaf_test!(ne r"r[] _{UsrdSsrdS", test_195);
mzpaf_test!(ne r"b2-H2O-33333V33333333", test_196);
mzpaf_test!(ne r"r[", test_197);
mzpaf_test!(ne r"r[]ACKOUUUiU[OUUUiUUUmamUiOwidation]ACK2,Em*0.7", test_198);
mzpaf_test!(ne r"f{", test_199);
mzpaf_test!(ne r"b2-H2O/3.2ppm*0.75,r[OUUUiUUUmamido/3.2ppm*0.7b4-H2O^-0000p.25", test_200);
mzpaf_test!(ne r"x5-9N92[17O-1414444444[ion]IIICH4[1O-1414444444[[ion]IImO]Iion]+][M", test_201);
mzpaf_test!(ne r"y4{MCI}-3333I3333334{M[3|_#8]H33333333F6/4m*{F9I3333333332`j ACN}-H4]A", test_202);
mzpaf_test!(ne r"y5-62[17O111111111][11][OO11-3333miO],H2B17OS1111111][O-3333my5[M+", test_203);
mzpaf_test!(ne r"r[", test_204);
mzpaf_test!(ne r"m3:4/1.1ppm,r[OUUUUUUUxid43:4/1.1ppm,_{}4!5/m", test_205);
mzpaf_test!(ne r"f{UWK}-1111F111111117", test_206);
mzpaf_test!(ne r"x5-9292[17O-1414496855]nhi[PQTAK7\17O-1414444444[[ion]IImO]Iion]+][M", test_207);
mzpaf_test!(ne r"_{UU ũK}-i-11111K11111111", test_208);
mzpaf_test!(ne r"f{}-", test_209);
mzpaf_test!(ne r"000000000@r[d aFCwxxon]ACK", test_210);
mzpaf_test!(ne r"y5-45V555555552[2H2O555557O-42[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4342[17O1]SS555VSS42[17O1]SSS42[17O1]SS555VSS42[17O1]SSSSSO55555V55555V55555552H2O55355V5555555555M+Im]N", test_211);
mzpaf_test!(ne r"x5-9297[17O-1414444OOOOOon]IIICH4[17O714]IIICH4[7O-141444@44[[ion]IImO]Iion]+][M", test_212);
mzpaf_test!(ne r"0@r[]a]ICK", test_213);
mzpaf_test!(ne r"x5-92[17O44444411]92[17O44[ion]III7:31+][M", test_214);
mzpaf_test!(ne r"f{}K}", test_215);
mzpaf_test!(ne r"0@r[", test_216);
mzpaf_test!(ne r"0@r[eȆ,w,_{Urrdr[eȆpm*1.23,b6H", test_217);
mzpaf_test!(ne r"0@r[OxiMdaACwxxxxxxxxqglrhhzlgjuxcfkwsq_FBW", test_218);
mzpaf_test!(ne r"0@r[K", test_219);
mzpaf_test!(ne r"0@r[dK", test_220);
mzpaf_test!(ne r"00000000000000000@r[daACwxxon]ACK", test_221);
mzpaf_test!(ne r"f{", test_222);
mzpaf_test!(ne r"x5-9HHHHHHHHHHHH292[17O-1414444444[Gon]IIICSHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH3+][M", test_223);
mzpaf_test!(ne r"x5-92IIIII[17O-1414444444[io+]IIICH4YS[M+]^}-T1+][444444[imO]III}-T1+][M", test_224);
mzpaf_test!(ne r"b2-900F000000W60000000W0@0900F000000W60000000W0@00000000000000  0", test_225);
mzpaf_test!(ne r"y5-4282H2O555557O-42[17O1]SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS4242[17O1]SS555VSS42[17O1]SSS42[17O1]SS555VSS42[17O1]SSSSSO55555V55555V55555552H2O55355V555555555 5M+Im]N", test_226);
mzpaf_test!(ne r"b2-000011100W+0000000110W0000000000000000000000001100000000", test_227);
mzpaf_test!(ne r"x5-9292[17O-1414435279]gji\TZJFG7^17O-1414444444[[ion]IImO]Iion]+][M", test_228);
mzpaf_test!(ne r"x5-92[17O4444444 11+]-92[17O44444444 11+]-92[4[ion]III7:31+][M", test_229);
mzpaf_test!(ne r"m3:4,y544+CCCC[15N73333333y544+CCC[15N7333333y544+3C[15N733333333y544+CCCC[15N73y544+CCC[15N73333333y544+3C[15N733333333033KS[M+]", test_230);
mzpaf_test!(ne r"b2-000001111W0000000000000110000000", test_231);
mzpaf_test!(ne r"b2-000001111W0000000000000010000000", test_232);
mzpaf_test!(ne r"x3-9292[17O-1414444bamicomethyl]-_}Vgkbblxa	444[ion]IIE䆛l_-]}Noqn~2,8rhy/054-@   IionR+][M", test_233);
mzpaf_test!(ne r"m3:4/1.1ppm,p-H33PO4^2,0@f{}-CH3PO,^2,0@f{U}-MMMMid}", test_234);
mzpaf_test!(ne r"x5-92[17O4444444 11]-917O4444444 11]-92[17O44[2[17O44[ion]III7:31+][M", test_235);
mzpaf_test!(ne r"r[", test_236);
mzpaf_test!(ne r"f{}bion]idatACatwon]idatACK+NNNQORWU", test_237);
mzpaf_test!(ne r"r[],_{Vgkbblxa	EfjCarbamidom%thylnVgkbblxbzu", test_238);
mzpaf_test!(ne r"x5-92[17O-044444443:4/]DII", test_239);
mzpaf_test!(ne r"x5-9292[17O-1414465576]mgm`JKVAP2\17O-1414444444[[ion]IImO]Iion]+][M", test_240);
mzpaf_test!(ne r"x5-9C92[17O-1414444444[ion2Pon]IIICon]4444[[17O-141444aIMO]Iion]+][M", test_241);
mzpaf_test!(ne r"0@y0-4WWWWWW11111111114i5n]i", test_242);
mzpaf_test!(ne r"z5-72[17O1111111111][M+NNNNNNNN+N-72+4-7+1.1p+H+Na]1]", test_243);
mzpaf_test!(ne r"0@r[OGWOxidathooooo[X:5#5CK}-C/0.qrs", test_244);
mzpaf_test!(ne r"x5-9297[17O-1414444]IIICH4[17O-141444@44[[ion]IIm444[ion]IIICH4[17O-14]IIICH4[17O-141444@44[[ion]IImO]Iion]+][M", test_245);
mzpaf_test!(ne r"x5-9O97[17O-1414444444[ion]IIICH4[17O-14]IIICH4[17O-1414i3Tlm*1.23,IImO]Iion]+][M", test_246);
mzpaf_test!(ne r"f{", test_247);
mzpaf_test!(ne r"x5-92[17O444444444 11+]-92[17O44[i~n]III7:31+][M", test_248);
mzpaf_test!(ne r"x3-9292[17O-1414444444䆛l_-]}Noqn~2,8rhy[ion]dIE䆛l_-]}Noqn~2,8rhy/054-a5,SIion]+][M", test_249);
