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
            let (end, name) = if line[range.start_index() + 1..].starts_with('[') {
                let end = end_of_enclosure(line, range.start_index() + 2, b'[', b']').ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid mzPAF reference compound",
                            "The square brackets are not closed",
                            Context::line(None, line, range.start_index() + 1, 1),
                        )
                    },
                )?;
                Ok((end, &line[range.start_index() + 2..end]))
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
                            Context::line_range(None, line, range.start_index() + 2..end),
                        ))
                    },
                    |formula| Ok((end + 1..range.end, IonType::Reporter(formula))),
                )
        }
        Some(b'f') => {
            // Simple formula
            let formula_range = if line[range.start_index() + 1..].starts_with('{') {
                let end = end_of_enclosure(line, range.start_index() + 2, b'{', b'}').ok_or_else(
                    || {
                        CustomError::error(
                            "Invalid mzPAF formula fragment",
                            "The curly braces are not closed",
                            Context::line(None, line, range.start_index() + 1, 1),
                        )
                    },
                )?;
                Ok(range.start_index() + 2..end)
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
            amount = num.2.map_err(|err| {
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
            })?;
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
                    b'+' => NeutralLoss::Gain(amount, formula),
                    b'-' => NeutralLoss::Loss(amount, formula),
                    _ => unreachable!(),
                });
            } else if let Ok(formula) =
                MolecularFormula::from_pro_forma(line, first - 1..=last, false, false, true, false)
            {
                // Catches the case of a single isotope as formula
                neutral_losses.push(match c {
                    b'+' => NeutralLoss::Gain(amount, formula),
                    b'-' => NeutralLoss::Loss(amount, formula),
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
                b'+' => NeutralLoss::Gain(amount, formula),
                b'-' => NeutralLoss::Loss(amount, formula),
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
mzpaf_test!(ne r"r[]", fuzz_6);
mzpaf_test!(ne r"r[Adensosine", fuzz_7);
mzpaf_test!(ne r"f{}", fuzz_8);
mzpaf_test!(ne r"f{C2H6", fuzz_9);
