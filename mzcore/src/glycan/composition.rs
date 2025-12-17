use std::ops::Range;

use context_error::*;

use crate::{
    ParserResult,
    chemistry::MolecularFormula,
    glycan::{GLYCAN_PARSE_LIST, MonoSaccharide},
    helper_functions::{explain_number_error, next_number, str_starts_with},
    prelude::Chemical,
};

/// All possible errors when parsing a ProForma glycan composition
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum GlycanCompositionError {
    /// A missing closing brackets for a formula '}'
    NoClosingBracket,
    /// A name was used that is not defined
    UnrecognisedName,
    /// An invalid occurrence number was specified
    InvalidNumber,
    /// An invalid formula was given for a custom monosaccharide
    InvalidFormula,
    /// Some monosaccharide occurrence was outside of bounds for [`isize`]
    OutOfRangeOccurance,
    /// An empty composition was given, this also includes cases where monosaccharides cancel each other out like `Hex2Hex-2`
    #[default]
    Empty,
    /// A synonym was given (will only be checked in strict mode)
    ImproperName,
    /// A name was given with improper capitalisation (will only be checked in strict mode)
    ImproperCapitalisation,
    /// The occurrence was missing for a glycan, this is always generated because it might mean that the composition was parsed as something not intended by the user
    ImproperMissingOccurance,
}

impl ErrorKind for GlycanCompositionError {
    type Settings = ();
    fn descriptor(&self) -> &'static str {
        match self {
            Self::ImproperCapitalisation | Self::ImproperMissingOccurance | Self::ImproperName => {
                "warning"
            }
            _ => "error",
        }
    }
    fn ignored(&self, _settings: Self::Settings) -> bool {
        false
    }
    fn is_error(&self, _settings: Self::Settings) -> bool {
        match self {
            Self::ImproperCapitalisation | Self::ImproperMissingOccurance | Self::ImproperName => {
                false
            }
            _ => true,
        }
    }
}

impl From<GlycanCompositionError> for BasicKind {
    fn from(value: GlycanCompositionError) -> Self {
        if value.is_error(()) {
            BasicKind::Error
        } else {
            BasicKind::Warning
        }
    }
}

impl MonoSaccharide {
    /// Parse the given text as a ProForma glycan composition.
    /// ```rust
    /// use mzcore::glycan::MonoSaccharide;
    /// assert!(MonoSaccharide::pro_forma_composition::<false>("HexNAc5Hex3Fuc1").is_ok());
    /// assert!(MonoSaccharide::pro_forma_composition::<false>("HexNAc4Hex5Fuc1NeuAc1{C8H13N1O5Na1:z+1}1").is_ok());
    /// ```
    /// # Errors
    /// When the composition could not be read. Or when any of the glycans numbers outside of the valid numerical range.
    /// Warns when the amount is missing.
    pub fn pro_forma_composition<const STRICT: bool>(
        value: &str,
    ) -> ParserResult<'_, Vec<(Self, isize)>, GlycanCompositionError> {
        Self::pro_forma_composition_inner::<STRICT>(
            &Context::none().lines(0, value),
            value,
            0..value.len(),
        )
    }

    /// Parse the given text as a ProForma glycan composition.
    /// # Errors
    /// When the composition could not be read. Or when any of the glycans numbers outside of the valid numerical range.
    /// Warns when the amount is missing.
    pub fn pro_forma_composition_inner<'a, const STRICT: bool>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
    ) -> ParserResult<'a, Vec<(Self, isize)>, GlycanCompositionError> {
        let mut index = range.start;
        let mut errors = Vec::new();
        let end = line.len().min(range.end);
        let mut output = Vec::new();
        while index < end {
            let start_glycan = index;
            let sugar = if line[index..end].starts_with(' ') {
                index += 1;
                continue;
            } else if line[index..end].starts_with('{') {
                let end_formula = handle!(single errors,  line[index + 1..end].find('}').ok_or_else(||BoxedError::new(
                    GlycanCompositionError::NoClosingBracket,
                    "Invalid ProForma glycan",
                    "The custom formula is not closed. No closing bracket '}' could be found.",
                    base_context.clone().add_highlight((0, index..end)),
                )));
                println!("{}", &line[index + 1..index + 1 + end_formula]);
                let formula = handle!(single errors,
                    MolecularFormula::pro_forma_inner::<true, false>(base_context, line, index + 1..index + 1+end_formula).map_err(|e| e.convert(|_| GlycanCompositionError::InvalidFormula)));
                index += end_formula + 2;
                Self::new(crate::glycan::BaseSugar::Custom(formula.into()), &[])
            } else {
                let mut found = None;
                'find_glycan: for (names, sugar) in GLYCAN_PARSE_LIST.iter() {
                    for name in names {
                        if str_starts_with::<true>(&line[index..end], name) {
                            found = Some(sugar.clone());
                            if STRICT {
                                let pro_forma_name = sugar.pro_forma_name();
                                if **name != pro_forma_name {
                                    combine_error(
                                        &mut errors,
                                        BoxedError::new(
                                            GlycanCompositionError::ImproperName,
                                            "Improper ProForma glycan",
                                            format!(
                                                "While `{name}` can be unambiguously parsed the proper name in ProForma is `{pro_forma_name}`."
                                            ),
                                            base_context.clone().add_highlight((
                                                0,
                                                index,
                                                name.len(),
                                            )),
                                        ),
                                        (),
                                    );
                                } else if !line[index..end].starts_with(&**name) {
                                    combine_error(
                                        &mut errors,
                                        BoxedError::new(
                                            GlycanCompositionError::ImproperCapitalisation,
                                            "Improper ProForma glycan",
                                            "This glycan was not written with the proper capitalisation.",
                                            base_context.clone().add_highlight((
                                                0,
                                                index,
                                                name.len(),
                                            )),
                                        ),
                                        (),
                                    );
                                }
                            }
                            index += name.len();
                            break 'find_glycan;
                        }
                    }
                }
                handle!(single errors, found.ok_or_else(|| BoxedError::new(
                    GlycanCompositionError::UnrecognisedName,
                    "Invalid ProForma glycan",
                    "No valid glycan name could be recognised.",
                    base_context.clone().add_highlight((0, index..end)),
                )))
            };
            let num = if let Some((offset, _, num)) =
                next_number::<true, false, isize>(line, index..end)
            {
                index += offset;
                handle!(single errors, num.map_err(|e| BoxedError::new(
                                   GlycanCompositionError::InvalidNumber,
                    "Invalid ProForma glycan",
                    format!("The monosaccharide occurance number {}", explain_number_error(&e)),
                    base_context.clone().add_highlight((0, index-offset..index)),
                )))
            } else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        GlycanCompositionError::ImproperMissingOccurance,
                        "Improper ProForma glycan",
                        "No amount for this glycan was specified, it is assumed to occur once.",
                        base_context.clone().add_highlight((0, start_glycan..index)),
                    ),
                    (),
                );
                1
            };

            output.push((sugar, num));
        }

        let composition = handle!(single errors, Self::simplify_composition(output).ok_or_else(|| {
            BoxedError::new(
                GlycanCompositionError::OutOfRangeOccurance,
                "Invalid ProForma glycan composition",
                format!(
                    "The occurrence of a monosaccharide species is outside of the range {} to {}",
                    isize::MIN,
                    isize::MAX
                ),
                base_context.clone().add_highlight((0, range.clone())),
            )
        }));
        if let Some(f) = composition
            .iter()
            .try_fold(MolecularFormula::default(), |acc, (s, n)| {
                s.formula()
                    .checked_mul_isize(*n)
                    .and_then(|f| acc.checked_add(&f))
            })
        {
            if composition.is_empty() || f.is_empty() {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        GlycanCompositionError::Empty,
                        "Invalid ProForma glycan composition",
                        "The glycan composition is empty",
                        base_context.clone().add_highlight((0, range)),
                    ),
                    (),
                );
                Err(errors)
            } else {
                Ok((composition, errors))
            }
        } else {
            combine_error(
                &mut errors,
                BoxedError::new(
                    GlycanCompositionError::OutOfRangeOccurance,
                    "Invalid ProForma glycan composition",
                    "The occurance of an element overflowed during calculation of the full molecular formula",
                    base_context.clone().add_highlight((0, range.clone())),
                ),
                (),
            );
            Err(errors)
        }
    }

    /// Simplify a glycan composition to be sorted and deduplicated.
    /// Returns None if overflow occurred, meaning that there where more than `isize::MAX` or less then `isize::MIN` monosaccharides for one species.
    pub(crate) fn simplify_composition(
        mut composition: Vec<(Self, isize)>,
    ) -> Option<Vec<(Self, isize)>> {
        // Sort on monosaccharide
        composition.retain(|el| el.1 != 0 && !el.0.formula().is_empty());
        composition.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate
        let mut max = composition.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &composition[index];
            let next = &composition[index + 1];
            if this.0 == next.0 {
                composition[index].1 = composition[index].1.checked_add(next.1)?;
                composition.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        composition.retain(|el| el.1 != 0);
        Some(composition)
    }

    /// Parse the given text as a MSFragger glycan composition. Examples:
    /// * HexNAc(5)Hex(3)Fuc(1)
    /// * HexNAc(4)Hex(5)Fuc(1)NeuAc(1)
    /// # Errors
    /// When the composition could not be read. Or when any of the glycans occurs outside of the valid range
    pub fn byonic_composition(text: &str) -> Result<Vec<(Self, isize)>, BoxedError<'_, BasicKind>> {
        let mut index = 0;
        let mut output = Vec::new();
        while index < text.len() {
            if text[index..].starts_with(' ') {
                index += 1;
            } else if let Some(next_open_bracket) = text[index..].find('(') {
                if let Some(next_close_bracket) = text[index + next_open_bracket..].find(')') {
                    let name = text[index..index + next_open_bracket].trim();
                    let mut sugar = None;
                    for option in GLYCAN_PARSE_LIST.as_slice() {
                        for o in &option.0 {
                            if o.eq_ignore_ascii_case(name) {
                                sugar = Some(option.1.clone());
                                break;
                            }
                        }
                    }
                    let number = text[index + next_open_bracket + 1
                        ..index + next_open_bracket + next_close_bracket]
                        .trim()
                        .parse::<isize>();
                    output.push((
                        sugar.ok_or_else(|| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid MSFragger glycan composition",
                                "The sugar name could not be recognised",
                                Context::line(None, text, index, name.len()),
                            )
                        })?,
                        number.map_err(|err| {
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid MSFragger glycan composition",
                                format!("Sugar count number is {}", explain_number_error(&err)),
                                Context::line(
                                    None,
                                    text,
                                    index + next_open_bracket + 1,
                                    next_close_bracket - 1,
                                ),
                            )
                        })?,
                    ));
                    index += next_open_bracket + next_close_bracket + 1;
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid MSFragger glycan composition",
                        "No closing bracket found ')'",
                        Context::line(None, text, index + next_open_bracket, 1),
                    ));
                }
            } else if text[index..].chars().all(|c| c.is_ascii_whitespace()) {
                break; // Allow trailing whitespace
            } else {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid MSFragger glycan composition",
                    "No opening bracket found but there is text left, the format expected is 'Sugar(Number)'",
                    Context::line(None, text, index, 1),
                ));
            }
        }

        Self::simplify_composition(output).ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid MSFragger glycan composition",
                format!(
                    "The occurrence of one monosaccharide species is outside of the range {} to {}",
                    isize::MIN,
                    isize::MAX
                ),
                Context::show(text),
            )
        })
    }

    /// Display a composition as a ProForma glycan composition
    pub fn display_composition(composition: &[(Self, isize)]) -> String {
        composition
            .iter()
            .fold(String::new(), |acc, m| acc + &format!("{}{}", m.0, m.1))
    }
}
