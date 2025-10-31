//! mzPAF parser
use std::{num::NonZeroU16, ops::Range, sync::LazyLock};

use context_error::*;
use serde::{Deserialize, Serialize};

use mzcore::{
    chemistry::{
        Chemical, ELEMENT_PARSE_LIST, Element, MolecularCharge, MolecularFormula, MultiChemical,
        NeutralLoss, SatelliteLabel,
    },
    molecular_formula,
    ontology::{CustomDatabase, Ontology},
    quantities::Tolerance,
    sequence::{
        AminoAcid, Linked, PeptidePosition, Peptidoform, SemiAmbiguous, SequenceElement,
        SequencePosition, SimpleModification, SimpleModificationInner,
    },
    system::{Mass, MassOverCharge, OrderedMassOverCharge, e, isize::Charge, thomson},
};

use crate::{
    fragment::*,
    helper_functions::{
        RangeExtension, RangeMaths, end_of_enclosure, explain_number_error, next_number,
    },
    mzspeclib::AnalyteTarget,
};

impl Fragment {
    /// Parse a [mzPAF](https://www.psidev.info/mzPAF) peak annotation line (can contain multiple annotations).
    /// mzPAF version 1.0 is supported. Except for the SMILES constructs.
    ///
    /// # Errors
    /// When the annotation does not follow the format.
    pub fn mz_paf<'a>(
        line: &'a str,
        custom_database: Option<&CustomDatabase>,
        interpretation: &[(u32, AnalyteTarget)],
    ) -> Result<Vec<Self>, BoxedError<'a, BasicKind>> {
        Self::mz_paf_inner(
            &Context::none().lines(0, line),
            line,
            0..line.len(),
            custom_database,
            interpretation,
        )
    }

    /// This parses a substring of the given string as an mzPAF definition. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// When the annotation does not follow the format.
    pub(crate) fn mz_paf_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        custom_database: Option<&CustomDatabase>,
        interpretation: &[(u32, AnalyteTarget)],
    ) -> Result<Vec<Self>, BoxedError<'a, BasicKind>> {
        parse_intermediate_representation(base_context, line, range, custom_database).and_then(
            |annotations| {
                annotations
                    .into_iter()
                    .map(|a| a.into_fragment(interpretation, base_context))
                    .collect()
            },
        )
    }
}

/// Parse a mzPAF line into the internal representation, see [`parse_mz_paf`] for more information.
/// # Errors
/// When the annotation does not follow the format.
fn parse_intermediate_representation<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<Vec<PeakAnnotation>, BoxedError<'a, BasicKind>> {
    let mut annotations = Vec::new();

    // Parse first
    let (mut range, a) = parse_annotation(base_context, line, range, custom_database)?;
    annotations.push(a);

    // Parse any following
    while !range.is_empty() {
        if line.as_bytes().get(range.start_index()).copied() == Some(b',') {
            range = range.add_start(1_usize);
        } else {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF annotation delimiter",
                "Different mzPAF annotations should be separated with commas ','.",
                Context::line(None, line, range.start_index(), 1),
            ));
        }
        let (r, a) = parse_annotation(base_context, line, range, custom_database)?;
        range = r;
        annotations.push(a);
    }

    Ok(annotations)
}

/// Parse a single mzPAF peak annotation.
/// # Errors
/// When the annotation does not follow the format.
fn parse_annotation<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<(Range<usize>, PeakAnnotation), BoxedError<'a, BasicKind>> {
    let (left_range, auxiliary) = if line.as_bytes().get(range.start_index()).copied() == Some(b'&')
    {
        (range.add_start(1_usize), true)
    } else {
        (range, false)
    };
    let (left_range, analyte_number) = parse_analyte_number(base_context, line, left_range)?;
    let (left_range, ion) = parse_ion(base_context, line, left_range, custom_database)?;
    let (left_range, neutral_losses) = parse_neutral_loss(base_context, line, left_range)?;
    let (left_range, isotopes) = parse_isotopes(base_context, line, left_range)?;
    let (left_range, adduct_type) = parse_adduct_type(base_context, line, left_range)?;
    let (left_range, charge) = parse_charge(base_context, line, left_range)?;
    let (left_range, deviation) = parse_deviation(base_context, line, left_range)?;
    let (left_range, confidence) = parse_confidence(base_context, line, left_range)?;
    Ok((
        left_range,
        PeakAnnotation {
            auxiliary,
            analyte_number,
            ion,
            neutral_losses,
            isotopes,
            charge: adduct_type
                .unwrap_or_else(|| MolecularCharge::proton(Charge::new::<e>(charge.value))),
            deviation,
            confidence,
        },
    ))
}

/// An mzPAF single peak annotation.
#[derive(Clone, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct PeakAnnotation {
    pub(in crate::fragment) auxiliary: bool,
    pub(in crate::fragment) analyte_number: u32,
    pub(in crate::fragment) ion: IonType,
    pub(in crate::fragment) neutral_losses: Vec<NeutralLoss>,
    pub(in crate::fragment) isotopes: Vec<(i32, Isotope)>,
    pub(in crate::fragment) charge: MolecularCharge,
    pub(in crate::fragment) deviation: Option<Tolerance<OrderedMassOverCharge>>,
    pub(in crate::fragment) confidence: Option<f64>,
}

impl PeakAnnotation {
    /// Convert a peak annotation into a fragment.
    /// # Errors
    /// If the fragment refers to a non existing analyte, or if the referenced analyte is not a
    /// peptidoform ion when that is necessary, or when the fragment is outside of bounds for a
    /// peptidoform.
    pub(crate) fn into_fragment<'a>(
        self,
        interpretation: &[(u32, AnalyteTarget)],
        context: &Context<'a>,
    ) -> Result<Fragment, BoxedError<'a, BasicKind>> {
        // Get the peptidoform (assume no cross-linkers)
        let target = if self.analyte_number == 0 {
            None
        } else {
            Some(
                &interpretation
                    .iter()
                    .find(|(n, _)| *n == self.analyte_number)
                    .ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF analyte number",
                            "The analyte that is referenced does not exist",
                            context.clone(),
                        )
                    })?
                    .1,
            )
        };

        let (formula, ion) = match self.ion {
            IonType::Unknown(series) => (None, FragmentType::Unknown(series)),
            IonType::MainSeries(series, sub, ordinal, specific_interpretation) => {
                let (sequence_length, n_aa, c_aa) = if let Some(target) = target {
                    let specific_interpretation: Option<&Peptidoform<Linked>> =
                        specific_interpretation.as_ref().map(AsRef::as_ref);
                    let peptidoform= specific_interpretation.unwrap_or(match target {
                        AnalyteTarget::PeptidoformIon(pep) => {
                            &pep.peptidoforms()[0] // TODO: this assumes no cross-linkers
                        }
                        AnalyteTarget::MolecularFormula(_) => {
                            return Err(BoxedError::new(BasicKind::Error, "Invalid mzPAF peak annotation", "The peak annotation is a peptidoform main series fragment but the referenced analyte is a molecular formula.", context.clone()));
                        }
                        AnalyteTarget::Unknown(_) => {
                            return Err(BoxedError::new(BasicKind::Error, "Invalid mzPAF peak annotation", "The peak annotation is a peptidoform main series fragment but the referenced analyte does not contain a recognised analyte definition.", context.clone()));
                        }
                    });
                    let sequence_length = peptidoform.len();
                    let Some(n_aa) = peptidoform
                        .sequence()
                        .get(ordinal - 1)
                        .map(|a| a.aminoacid.aminoacid())
                    else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF peak annotation",
                            "The referenced position does not exist in the analyte peptidoform",
                            context.clone(),
                        ));
                    };
                    let Some(c_aa) = sequence_length
                        .checked_sub(ordinal)
                        .and_then(|i| peptidoform.sequence().get(i))
                        .map(|a| a.aminoacid.aminoacid())
                    else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF peak annotation",
                            "The referenced position does not exist in the analyte peptidoform",
                            context.clone(),
                        ));
                    };
                    (sequence_length, n_aa, c_aa)
                } else {
                    (0, AminoAcid::Unknown, AminoAcid::Unknown)
                };
                (
                    None,
                    match series {
                        b'a' => FragmentType::a(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(ordinal - 1),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'b' => FragmentType::b(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(ordinal - 1),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'c' => FragmentType::c(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(ordinal - 1),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'd' => FragmentType::d(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(ordinal - 1),
                                series_number: ordinal,
                                sequence_length,
                            },
                            n_aa,
                            0,
                            0,
                            sub,
                        ),
                        b'v' => FragmentType::v(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            c_aa,
                            0,
                            0,
                        ),
                        b'w' => FragmentType::w(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            c_aa,
                            0,
                            0,
                            sub,
                        ),
                        b'x' => FragmentType::x(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'y' => FragmentType::y(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(
                                    sequence_length.saturating_sub(ordinal),
                                ),
                                series_number: ordinal,
                                sequence_length,
                            },
                            0,
                        ),
                        b'z' => FragmentType::z(
                            PeptidePosition {
                                sequence_index: SequencePosition::Index(
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
                aa.formulas().first().map(|f| {
                    f + m.as_ref().map(|m| m.formula()).unwrap_or_default()
                        - molecular_formula!(C 1 O 1)
                }),
                FragmentType::Immonium(
                    None,
                    m.map_or_else(
                        || SequenceElement::new(aa.into(), None),
                        |m| SequenceElement::new(aa.into(), None).with_simple_modification(m),
                    ),
                ),
            ),
            IonType::Internal(start, end) => {
                let sequence_length = match target {
                    Some(AnalyteTarget::PeptidoformIon(pep)) => {
                        pep.peptidoforms()[0].len() // TODO: this assumes no cross-linkers
                    }
                    Some(AnalyteTarget::MolecularFormula(_)) => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF peak annotation",
                            "The peak annotation is a peptidoform main series fragment but the referenced analyte is a molecular formula.",
                            context.clone(),
                        ));
                    }
                    Some(AnalyteTarget::Unknown(_)) => {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF peak annotation",
                            "The peak annotation is a peptidoform main series fragment but the referenced analyte does not contain a recognised analyte definition.",
                            context.clone(),
                        ));
                    }
                    None => 0,
                };

                (
                    None,
                    FragmentType::Internal(
                        None,
                        PeptidePosition::n(SequencePosition::Index(start - 1), sequence_length),
                        PeptidePosition::n(SequencePosition::Index(end - 1), sequence_length),
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
        Ok(Fragment {
            formula,
            charge: self.charge.charge(),
            ion,
            peptidoform_ion_index: self.analyte_number.checked_sub(1).map(|a| a as usize), // TODO: figure out how to handle this properly
            peptidoform_index: Some(0),
            neutral_loss: self.neutral_losses.into(),
            deviation: self.deviation,
            confidence: self.confidence.map(Into::into),
            auxiliary: self.auxiliary,
        })
    }
}

impl std::fmt::Display for PeakAnnotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.auxiliary {
            write!(f, "&")?;
        }
        write!(f, "{}@{}", self.analyte_number, self.ion)?;
        for loss in &self.neutral_losses {
            write!(f, "{loss}")?;
        }
        let charge = self.charge.charge().value;
        if self.charge != MolecularCharge::proton(Charge::new::<e>(charge)) {
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
pub(crate) enum IonType {
    Unknown(Option<usize>),
    /// Main series, char identifier for the series, optional sattelite label, the series number, and possibly the interpretation
    MainSeries(
        u8,
        SatelliteLabel,
        usize,
        Option<Peptidoform<SemiAmbiguous>>,
    ),
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
                    sub,
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
pub(in crate::fragment) enum Isotope {
    General,
    Average,
    Specific(Element, NonZeroU16),
}

/// Parse a mzPAF analyte number. '1@...'
/// # Errors
/// When the ion is not formatted correctly.
fn parse_analyte_number<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<(Range<usize>, u32), BoxedError<'a, BasicKind>> {
    next_number::<false, false, u32>(line, range.clone()).map_or_else(
        || Ok((range.clone(), 1)),
        |num| {
            if line.as_bytes().get(num.0 + range.start).copied() != Some(b'@') {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF analyte number",
                    "The analyte number should be followed by an at sign '@'",
                    base_context
                        .clone()
                        .add_highlight((0, num.0 + range.start, 1)),
                ));
            }
            Ok((
                range.add_start(num.0 + 1),
                num.2.map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF analyte number",
                        format!("The analyte number number {}", explain_number_error(&err)),
                        base_context.clone().add_highlight((0, range.start, num.0)),
                    )
                })?,
            ))
        },
    )
}

/// Parse a mzPAF ion.
/// # Errors
/// When the ion is not formatted correctly.
fn parse_ion<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<(Range<usize>, IonType), BoxedError<'a, BasicKind>> {
    match line.as_bytes().get(range.start_index()).copied() {
        Some(b'?') => {
            if let Some(ordinal) =
                next_number::<false, false, usize>(line, range.add_start(1_usize))
            {
                Ok((
                    range.add_start(1 + ordinal.0),
                    IonType::Unknown(Some(ordinal.2.map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF unknown ion ('?') ordinal",
                            format!("The ordinal number {}", explain_number_error(&err)),
                            base_context.clone().add_highlight((
                                0,
                                range.start_index() + 1,
                                ordinal.0,
                            )),
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
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF main series ion ordinal",
                        "Only for the satellite ions 'd' and 'w' does a subtype exist, like 'wa12'",
                        base_context
                            .clone()
                            .add_highlight((0, range.start_index(), 1)),
                    ));
                }
                (
                    range.add_start(2_usize),
                    if sub == b'a' {
                        SatelliteLabel::A
                    } else {
                        SatelliteLabel::B
                    },
                )
            } else {
                (range.add_start(1_usize), SatelliteLabel::None)
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
                        )
                        .map_err(|errs| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid ProForma definition",
                                "The string could not be parsed as a ProForma definition",
                                Context::line_range(None, line, range.clone()),
                            )
                            .add_underlying_errors(errs)
                        });
                        interpretation.and_then(|(i, _)| {
                            i.into_semi_ambiguous()
                                .ok_or_else(|| BoxedError::new(BasicKind::Error,
                                    "Invalid mzPAF interpretation",
                                    "An mzPAF interpretation should be limited to `base-ProForma compliant` without any labile modifications",
                                    base_context.clone().add_highlight((0,range.start_index()..location))))
                                .map(|i| (location + 1, Some(i)))
                        })?
                        // TODO: proper error handling and add checks to the length of the sequence
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF main series ion ordinal",
                            "The asserted interpretation should have a closed curly bracket, like '0@b2{LL}'",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index(), 1)),
                        ));
                    }
                } else {
                    (range.start_index(), None)
                };
                let ordinal_num = ordinal.2.map_err(|err| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF ion ordinal",
                        format!("The ordinal number {}", explain_number_error(&err)),
                        base_context.clone().add_highlight((
                            0,
                            range.start_index() - ordinal.0, // Maybe also offset for interpretation?
                            ordinal.0,
                        )),
                    )
                })?;
                if ordinal_num == 0 {
                    Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF main series ion ordinal",
                        "The ion ordinal cannot be 0",
                        base_context.clone().add_highlight((
                            0,
                            range.start_index() - ordinal.0, // Maybe also offset for interpretation?
                            ordinal.0,
                        )),
                    ))
                } else {
                    Ok((
                        end..range.end,
                        IonType::MainSeries(c, sub, ordinal_num, interpretation),
                    ))
                }
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF main series ion ordinal",
                    "For a main series ion the ordinal should be provided, like 'a12'",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                ))
            }
        }
        Some(b'I') => {
            let amino_acid = line[range.clone()].chars().nth(1).ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF immonium",
                    "The source amino acid for this immonium ion should be present like 'IA'",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                )
            })?;
            let index = range.start_index() + 1 + amino_acid.len_utf8();
            let modification = if line[index..range.end].starts_with('[') {
                let end = end_of_enclosure(line, index + 1, b'[', b']').ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF immonium modification",
                        "The square brackets are not closed",
                        base_context.clone().add_highlight((0, index, 1)),
                    )
                })?;
                let modification = &line[index + 1..end];
                Some((
                    end - range.start_index() - 1, // Length of mod + [ + ]
                    Ontology::Unimod
                        .find_name(modification, None)
                        .or_else(|| {
                            modification.parse::<f64>().ok().map(|n| {
                                std::sync::Arc::new(SimpleModificationInner::Mass(
                                    Mass::new::<mzcore::system::dalton>(n).into(),
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
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF immonium ion",
                            "The provided amino acid is not a known amino acid",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index() + 1, 1)),
                        )
                    })?,
                    modification.map(|m| m.1),
                ),
            ))
        }
        Some(b'm') => {
            let first_ordinal = next_number::<false, false, usize>(line, range.add_start(1_usize))
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF internal ion first ordinal",
                        "The first ordinal for an internal ion should be present",
                        base_context
                            .clone()
                            .add_highlight((0, range.start_index(), 1)),
                    )
                })?;
            if line[range.clone()].chars().nth(first_ordinal.0 + 1) != Some(':') {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF internal ion ordinal separator",
                    "The internal ion ordinal separator should be a colon ':', like 'm4:6'",
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 1 + first_ordinal.0,
                        1,
                    )),
                ));
            }
            let second_ordinal = next_number::<false, false, usize>(
                line,
                range.add_start(first_ordinal.0.saturating_add(2)),
            )
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF internal ion second ordinal",
                    "The second ordinal for an internal ion should be present",
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 1 + first_ordinal.0,
                        1,
                    )),
                )
            })?;
            let first_location = first_ordinal.2.map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF internal ion first ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 1,
                        first_ordinal.0,
                    )),
                )
            })?;
            let second_location = second_ordinal.2.map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF internal ion second ordinal",
                    format!("The ordinal number {}", explain_number_error(&err)),
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 2 + first_ordinal.0,
                        second_ordinal.0,
                    )),
                )
            })?;
            if first_location == 0 {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF main series ion ordinal",
                    "The ion ordinal cannot be 0",
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 1,
                        first_ordinal.0,
                    )),
                ))
            } else if second_location == 0 {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF main series ion ordinal",
                    "The ion ordinal cannot be 0",
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + first_ordinal.0 + 2,
                        second_ordinal.0,
                    )),
                ))
            } else {
                Ok((
                    range.add_start(2 + first_ordinal.0 + second_ordinal.0),
                    IonType::Internal(first_location, second_location),
                ))
            }
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
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF named compound",
                            "The curly braces are not closed",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index() + 1, 1)),
                        )
                    },
                )?;
                Ok((
                    end - range.start_index(),
                    &line[range.start_index() + 2..end],
                ))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF named compound",
                    "A named compound must be named with curly braces '{}' after the '_'",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
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
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF reference compound",
                            "The square brackets are not closed",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index() + 1, 1)),
                        )
                    },
                )?;
                Ok((end, &line[range.start_index() + 2..end]))
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF reporter ion",
                    "A reporter ion must be named with square braces '[]' after the 'r'",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                ))
            }?;
            MZPAF_NAMED_MOLECULES
                .iter()
                .find_map(|n| (n.0.eq_ignore_ascii_case(name)).then_some(n.1.clone()))
                .map_or_else(
                    || {
                        Err(BoxedError::new(
                            BasicKind::Error,
                            "Unknown mzPAF named reporter ion",
                            "Unknown name",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index() + 2..end)),
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
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF formula fragment",
                            "The curly braces are not closed",
                            base_context
                                .clone()
                                .add_highlight((0, range.start_index() + 1, 1)),
                        )
                    },
                )?;
                Ok(range.start_index() + 2..end)
            } else {
                Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF formula",
                    "A formula must have the formula defined with curly braces '{}' after the 'f'",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                ))
            }?;
            let formula =
                MolecularFormula::from_pro_forma::<false, false>(line, formula_range.clone())?;

            Ok((
                range.add_start(3 + formula_range.len()),
                IonType::Formula(formula),
            ))
        }
        Some(b's') => Err(BoxedError::new(
            BasicKind::Error,
            "Unsupported feature",
            "SMILES strings are currently not supported in mzPAF definitions",
            base_context.clone().add_highlight((0, range.start, 1)),
        )), // TODO: return as Formula
        Some(_) => Err(BoxedError::new(
            BasicKind::Error,
            "Invalid ion",
            "An ion cannot start with this character",
            base_context.clone().add_highlight((0, range.start, 1)),
        )),
        None => Err(BoxedError::new(
            BasicKind::Error,
            "Invalid ion",
            "An ion cannot be an empty string",
            base_context.clone().add_highlight((0, range.start, 1)),
        )),
    }
}

/// Parse a neutral loss from the string.
/// # Errors
/// If the a neutral loss is detected but is invalid.
pub(in crate::fragment) fn parse_neutral_loss<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<(Range<usize>, Vec<NeutralLoss>), BoxedError<'a, BasicKind>> {
    let mut offset = 0;
    let mut neutral_losses = Vec::new();
    while let Some(c @ (b'-' | b'+')) = line.as_bytes().get(range.start_index() + offset).copied() {
        let mut amount = 1;
        let num_offset =
        // Parse leading number to detect how many times this loss occured
        if let Some(num) = next_number::<false, false, u16>(line, range.add_start(1 + offset)) {
            amount = num.2.map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF neutral loss leading amount",
                    format!(
                        "The neutral loss amount number {}",
                        explain_number_error(&err)
                    ),
                    base_context.clone().add_highlight((
                        0,
                        range.start_index() + 1 + offset,
                        range.start_index() + 1 + offset + num.0,
                    )),
                )
            })?;
            offset += num.0;
            num.0
        } else {
            0
        };

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
                    BoxedError::new(
                        BasicKind::Error,
                        "Unknown mzPAF named neutral loss",
                        "Opening bracket for neutral loss name was not closed",
                        base_context.clone().add_highlight((
                            0,
                            range.start_index() + 1 + offset,
                            1,
                        )),
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
                MolecularFormula::from_pro_forma::<false, false>(line, first - 1..=last)
            {
                // Catches the case of a single isotope as formula
                neutral_losses.push(match c {
                    b'+' => NeutralLoss::Gain(amount, formula),
                    b'-' => NeutralLoss::Loss(amount, formula),
                    _ => unreachable!(),
                });
            } else {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Unknown mzPAF named neutral loss",
                    "Unknown name",
                    base_context
                        .clone()
                        .add_highlight((0, offset - name.len() - 1, name.len())),
                ));
            }
        } else {
            let first = range.start_index() + 1 + offset;
            let last = line[first..]
                .char_indices()
                .take_while(|(_, c)| c.is_ascii_alphanumeric() || *c == '[' || *c == ']')
                .last()
                .ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF",
                        "Empty neutral loss",
                        base_context.clone().add_highlight((0, first..)),
                    )
                })?;
            let mut last = last.0 + last.1.len_utf8();
            if line[first..first + last].ends_with("[M") {
                last -= 2; // Detect any adduct types which might otherwise sneak in
            }
            let formula =
                MolecularFormula::from_pro_forma::<false, false>(line, first..first + last)?;
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

/// The parsed isotopes. First the left range, then all found isotopes as the multiplier and isotope type.
type Isotopes = (Range<usize>, Vec<(i32, Isotope)>);

/// Parse isotopes definition from the string.
/// # Errors
/// If the detected isotopes are detected but invalid.
fn parse_isotopes<'a>(
    base_context: &Context<'a>,
    line: &str,
    range: Range<usize>,
) -> Result<Isotopes, BoxedError<'a, BasicKind>> {
    let mut offset = 0;
    let mut isotopes = Vec::new();
    while let Some(c @ (b'-' | b'+')) = line.as_bytes().get(range.start_index() + offset).copied() {
        offset += 1;
        let mut amount = 1;
        // Parse leading number to detect how many times this isotope occurred
        if let Some(num) = next_number::<false, false, u16>(line, range.add_start(offset)) {
            amount = i32::from(num.2.map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF isotope leading amount",
                    format!("The isotope amount number {}", explain_number_error(&err)),
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index() + offset, num.0)),
                )
            })?);
            offset += num.0;
        }
        if c == b'-' {
            amount *= -1;
        }

        // Check if i
        if line.as_bytes().get(range.start_index() + offset).copied() != Some(b'i') {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF isotope",
                "An isotope should be indicated with a lowercase 'i', eg '+i', '+5i', '+2iA', '+i13C'",
                base_context
                    .clone()
                    .add_highlight((0, range.start_index() + offset, 1)),
            ));
        }
        offset += 1;

        // Check if a specific isotope
        if let Some(num) = next_number::<false, false, NonZeroU16>(line, range.add_start(offset)) {
            let nucleon = NonZeroU16::from(num.2.map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF isotope nucleon number",
                    format!("The nucleon number {}", explain_number_error(&err)),
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index() + offset, num.0)),
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
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF isotope element",
                    "No recognised element symbol was found",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index() + offset, 1)),
                )
            })?;
            if !element.is_valid(Some(nucleon)) {
                let ln = element.symbol().len();
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF isotope",
                    format!(
                        "The nucleon number {nucleon} does not have a defined mass for {element}",
                    ),
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index() + offset - ln, ln)),
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

/// Parse adduct types from the string.
/// # Errors
/// If and adduct type is detected but is invalid.
fn parse_adduct_type<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<MolecularCharge>), BoxedError<'a, BasicKind>> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'[') {
        let closing =
            end_of_enclosure(line, range.start_index() + 1, b'[', b']').ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF adduct type",
                    "No closing bracket found for opening bracket of adduct type",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                )
            })?; // Excluding the ']' closing bracket
        if line.as_bytes().get(range.start_index() + 1).copied() != Some(b'M') {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF adduct type",
                "The adduct type should start with 'M', as in '[M+nA]'",
                base_context
                    .clone()
                    .add_highlight((0, range.start_index() + 1, 1)),
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
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid mzPAF adduct leading amount",
                        format!("The adduct amount number {}", explain_number_error(&err)),
                        base_context.clone().add_highlight((
                            0,
                            range.start_index() + offset,
                            range.start_index() + offset + num.0,
                        )),
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
            let formula =
                MolecularFormula::from_pro_forma::<false, false>(line, first..first + last)?;
            carriers.push((amount as isize, formula));
            offset += last;
        }
        if line.as_bytes().get(range.start_index() + offset).copied() != Some(b']') {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF adduct type",
                "The adduct type should be closed with ']'",
                base_context
                    .clone()
                    .add_highlight((0, range.start_index() + offset, 1)),
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
/// If there is no number after the caret, or if the number is invalid (outside of range and the like).
fn parse_charge<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<(Range<usize>, Charge), BoxedError<'a, BasicKind>> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'^') {
        let charge =
            next_number::<true, false, isize>(line, range.add_start(1_usize)).ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF charge",
                    "The number after the charge symbol should be present, eg '^2'.",
                    base_context
                        .clone()
                        .add_highlight((0, range.start_index(), 1)),
                )
            })?;
        Ok((
            range.add_start(charge.0 + 1),
            Charge::new::<e>(
                if charge.1 { 1 } else { -1 }
                    * charge.2.map_err(|err| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid mzPAF charge",
                            format!("The charge number {}", explain_number_error(&err)),
                            base_context.clone().add_highlight((
                                0,
                                range.start_index() + 1,
                                charge.0,
                            )),
                        )
                    })?,
            ),
        ))
    } else {
        Ok((range, Charge::new::<e>(1)))
    }
}

/// Parsed deviation, the left range and the deviation if present.
type Deviation = (Range<usize>, Option<Tolerance<OrderedMassOverCharge>>);

/// Parse a mzPAF deviation, either a ppm or m/z deviation.
/// # Errors
/// When the deviation is not '<number>' or '<number>ppm'.
fn parse_deviation<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<Deviation, BoxedError<'a, BasicKind>> {
    if line.as_bytes().get(range.start_index()).copied() == Some(b'/') {
        let number =
            next_number::<true, true, f64>(line, range.add_start(1_usize)).ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF deviation",
                    "A deviation should be a number",
                    base_context
                        .clone()
                        .add_highlight((0, range.start..=range.start + 1)),
                )
            })?;
        let deviation = number.2.map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF deviation",
                format!("The deviation number {err}",),
                base_context
                    .clone()
                    .add_highlight((0, range.start + 1..range.start + 1 + number.0)),
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
                Some(Tolerance::new_absolute(MassOverCharge::new::<thomson>(
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
fn parse_confidence<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<(Range<usize>, Option<f64>), BoxedError<'a, BasicKind>> {
    // TODO: the range is in characters
    if line.chars().nth(range.start_index()) == Some('*') {
        let number =
            next_number::<true, true, f64>(line, range.add_start(1_usize)).ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid mzPAF confidence",
                    "A confidence should be a number",
                    base_context
                        .clone()
                        .add_highlight((0, range.start..=range.start + 1)),
                )
            })?;
        let confidence = number.2.map_err(|err| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid mzPAF confidence",
                format!("The confidence number {err}",),
                base_context
                    .clone()
                    .add_highlight((0, range.start + 1..range.start + 1 + number.0)),
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
        ("Hex", molecular_formula!(C 6 H 10 O 5)),
        ("HexNAc", molecular_formula!(C 8 H 13 N 1 O 5)),
        ("dHex", molecular_formula!(C 6 H 10 O 4)),
        ("NeuAc", molecular_formula!(C 11 H 17 N 1 O 8)),
        ("NeuGc", molecular_formula!(C 11 H 17 N 1 O 9)),
        ("TMT126", molecular_formula!(C 8 N 1 H 15)),
        ("TMT127N", molecular_formula!(C 8 [15 N 1] H 15)),
        ("TMT127C", molecular_formula!(C 7 [13 C 1] N 1 H 15)),
        ("TMT128N", molecular_formula!(C 7 [13 C 1] [15 N 1] H 15)),
        ("TMT128C", molecular_formula!(C 6 [13 C 2] N 1 H 15)),
        ("TMT129N", molecular_formula!(C 6 [13 C 2] [15 N 1] H 15)),
        ("TMT129C", molecular_formula!(C 5 [13 C 3] N 1 H 15)),
        ("TMT130N", molecular_formula!(C 5 [13 C 3] [15 N 1] H 15)),
        ("TMT130C", molecular_formula!(C 4 [13 C 4] N 1 H 15)),
        ("TMT131N", molecular_formula!(C 4 [13 C 4] [15 N 1] H 15)),
        ("TMT131C", molecular_formula!(C 3 [13 C 5] N 1 H 15)),
        ("TMT132N", molecular_formula!(C 3 [13 C 5] [15 N 1] H 15)),
        ("TMT132C", molecular_formula!(C 2 [13 C 6] N 1 H 15)),
        ("TMT133N", molecular_formula!(C 2 [13 C 6] [15 N 1] H 15)),
        ("TMT133C", molecular_formula!(C 1 [13 C 7] N 1 H 15)),
        ("TMT134N", molecular_formula!(C 1 [13 C 7] [15 N 1] H 15)),
        ("TMT134C", molecular_formula!(C 0 [13 C 8] N 1 H 15)),
        ("TMT135N", molecular_formula!(C 0 [13 C 8] [15 N 1] H 15)),
        ("TMTzero", molecular_formula!(C 12 H 20 N 2 O 2)),
        ("TMTpro_zero", molecular_formula!(C 15 H 25 N 3 O 3)),
        ("TMT2plex", molecular_formula!(C 11 [ 13 C 1] H 20 N 2 O 2)),
        (
            "TMT6plex",
            molecular_formula!(C 8 [13 C 5] H 20 N 1 [ 15 N 1] O 2),
        ),
        (
            "TMTpro",
            molecular_formula!(C 8 [13 C 7] H 25 [15 N 2] N 1 O 3),
        ),
        ("iTRAQ113", molecular_formula!(C 6 N 2 H 12)),
        ("iTRAQ114", molecular_formula!(C 5 [13 C 1] N 2 H 12)),
        (
            "iTRAQ115",
            molecular_formula!(C 5 [13 C 1] N 1 [15 N 1] H 12),
        ),
        (
            "iTRAQ116",
            molecular_formula!(C 4 [13 C 2] N 1 [15 N 1] H 12),
        ),
        (
            "iTRAQ117",
            molecular_formula!(C 3 [13 C 3] N 1 [15 N 1] H 12),
        ),
        ("iTRAQ118", molecular_formula!(C 3 [13 C 3] [15 N 2] H 12)),
        ("iTRAQ119", molecular_formula!(C 4 [13 C 2] [15 N 2] H 12)),
        ("iTRAQ121", molecular_formula!([13 C 6] [15 N 2] H 12)),
        (
            "iTRAQ4plex",
            molecular_formula!(C 4 [13 C 3] H 12 N 1 [15 N 1] O 1),
        ),
        (
            "iTRAQ8plex",
            molecular_formula!(C 7 [13 C 7] H 24 N 3 [15 N 1] O 3),
        ),
        ("TMT126-ETD", molecular_formula!(C 7 N 1 H 15)),
        ("TMT127N-ETD", molecular_formula!(C 7 [15 N 1] H 15)),
        ("TMT127C-ETD", molecular_formula!(C 6 [13 C 1] N 1 H 15)),
        (
            "TMT128N-ETD",
            molecular_formula!(C 6 [13 C 1] [15 N 1] H 15),
        ),
        ("TMT128C-ETD", molecular_formula!(C 5 [13 C 2] N 1 H 15)),
        (
            "TMT129N-ETD",
            molecular_formula!(C 5 [13 C 2] [15 N 1] H 15),
        ),
        ("TMT129C-ETD", molecular_formula!(C 4 [13 C 3] N 1 H 15)),
        (
            "TMT130N-ETD",
            molecular_formula!(C 4 [13 C 3] [15 N 1] H 15),
        ),
        ("TMT130C-ETD", molecular_formula!(C 3 [13 C 4] N 1 H 15)),
        (
            "TMT131N-ETD",
            molecular_formula!(C 3 [13 C 4] [15 N 1] H 15),
        ),
        ("TMT131C-ETD", molecular_formula!(C 2 [13 C 5] N 1 H 15)),
        ("sidechain_A", molecular_formula!(C 1 H 3 )),
        ("sidechain_C", molecular_formula!(C 1 H 3 S 1)),
        ("sidechain_D", molecular_formula!(C 2 H 2 O 2)),
        ("sidechain_E", molecular_formula!(C 3 H 4 O 2)),
        ("sidechain_F", molecular_formula!(C 7 H 7 )),
        ("sidechain_G", molecular_formula!(H 1)),
        ("sidechain_H", molecular_formula!(C 4 H 5 N 2)),
        ("sidechain_I", molecular_formula!(C 4 H 9 )),
        ("sidechain_J", molecular_formula!(C 4 H 9 )),
        ("sidechain_K", molecular_formula!(C 4 H 10 N 1)),
        ("sidechain_L", molecular_formula!(C 4 H 9 )),
        ("sidechain_M", molecular_formula!(C 3 H 7 S 1)),
        ("sidechain_N", molecular_formula!(C 2 H 4 N 1 O 1)),
        ("sidechain_O", molecular_formula!(C 9 H 17 N 2 O 1)),
        ("sidechain_Q", molecular_formula!(C 3 H 6 N 1 O 1)),
        ("sidechain_R", molecular_formula!(C 4 H 10 N 3)),
        ("sidechain_S", molecular_formula!(C 1 H 3 O 1)),
        ("sidechain_T", molecular_formula!(C 2 H 5 O 1)),
        ("sidechain_U", molecular_formula!(C 1 H 3 Se 1)),
        ("sidechain_V", molecular_formula!(C 3 H 7 )),
        ("sidechain_W", molecular_formula!(C 9 H 8 N 1)),
        ("sidechain_Y", molecular_formula!(C 7 H 7 O 1)),
        ("Cytosine", molecular_formula!(C 4 H 5 N 3 O 1)),
        ("Adenine", molecular_formula!(C 5 H 5 N 5)),
        ("Guanine", molecular_formula!(C 5 H 5 N 5 O 1)),
        ("Uracil", molecular_formula!(C 4 H 4 N 2 O 2)),
        ("Thymine", molecular_formula!(C 5 H 6 N 2 O 2)),
    ]
});

#[test]
#[allow(clippy::missing_panics_doc)]
fn neutral_loss() {
    assert_eq!(
        parse_neutral_loss(&Context::none(), "-H2O", 0..4),
        Ok((
            4..4,
            vec![NeutralLoss::Loss(1, molecular_formula!(H 2 O 1))]
        ))
    );
    assert_eq!(
        parse_neutral_loss(&Context::none(), "+H2O", 0..4),
        Ok((
            4..4,
            vec![NeutralLoss::Gain(1, molecular_formula!(H 2 O 1))]
        ))
    );
    assert_eq!(
        parse_neutral_loss(&Context::none(), "+NH3", 0..4),
        Ok((
            4..4,
            vec![NeutralLoss::Gain(1, molecular_formula!(N 1 H 3))]
        ))
    );
    assert_eq!(
        parse_neutral_loss(&Context::none(), "/-0.0008", 0..8),
        Ok((0..8, vec![]))
    );
}

#[test]
#[allow(clippy::missing_panics_doc)]
fn parse_correctly() {
    let pep = [(
        1_u32,
        AnalyteTarget::PeptidoformIon(
            mzcore::sequence::PeptidoformIon::pro_forma("AAAAAAAAAA", None)
                .unwrap()
                .0,
        ),
    )];
    let a = "y8^2/-0.0017";
    let (_, parse_a) = parse_annotation(&Context::none(), a, 0..a.len(), None).unwrap();
    assert!(!parse_a.auxiliary);
    assert_eq!(parse_a.analyte_number, 1);
    assert_eq!(
        parse_a.ion,
        IonType::MainSeries(b'y', SatelliteLabel::None, 8, None)
    );
    assert_eq!(parse_a.neutral_losses, Vec::new());
    assert_eq!(parse_a.isotopes, Vec::new());
    assert_eq!(parse_a.charge, MolecularCharge::proton(Charge::new::<e>(2)));
    assert_eq!(
        parse_a.deviation,
        Some(Tolerance::Absolute(
            MassOverCharge::new::<thomson>(-0.0017).into()
        )),
    );
    assert_eq!(parse_a.confidence, None);
    let frag_a = parse_a.into_fragment(&pep, &Context::none()).unwrap();
    assert_eq!(
        frag_a.ion.position(),
        Some(&PeptidePosition::c(SequencePosition::Index(2), 10))
    );

    let b = "y8+i^2/0.0002";
    let (_, parse_b) = parse_annotation(&Context::none(), b, 0..b.len(), None).unwrap();
    assert!(!parse_b.auxiliary);
    assert_eq!(parse_b.analyte_number, 1);
    assert_eq!(
        parse_b.ion,
        IonType::MainSeries(b'y', SatelliteLabel::None, 8, None)
    );
    assert_eq!(parse_b.neutral_losses, Vec::new());
    assert_eq!(parse_b.isotopes, vec![(1, Isotope::General)]);
    assert_eq!(parse_b.charge, MolecularCharge::proton(Charge::new::<e>(2)));
    assert_eq!(
        parse_b.deviation,
        Some(Tolerance::Absolute(
            MassOverCharge::new::<thomson>(0.0002).into()
        )),
    );
    assert_eq!(parse_b.confidence, None);
}
