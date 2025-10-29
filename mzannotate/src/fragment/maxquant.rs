use std::ops::Range;

use context_error::{BasicKind, BoxedError, Context, CreateError};
use mzcore::ontology::CustomDatabase;

use crate::{
    fragment::{Fragment, PeakAnnotation},
    helper_functions::{explain_number_error, next_number},
    mzspeclib::AnalyteTarget,
};

impl Fragment {
    /// Parse a maxquant ion annotation. Some examples:
    /// ```plain
    /// y2;y5;y6;y13;y8-NH3;y30(2+);a2;b25;b2-H2O;b3-H2O;b5-H2O;b14(2+);z°10;z°12;z'12;z'13;cm5;cp27;cp31;b27
    /// ```
    ///
    /// # Errors
    /// When the annotation does not follow the format.
    pub fn maxquant<'a>(
        line: &'a str,
        custom_database: Option<&CustomDatabase>,
        interpretation: &[(u32, AnalyteTarget)],
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        Self::maxquant_inner(
            &Context::none().lines(0, line),
            line,
            0..line.len(),
            custom_database,
            interpretation,
        )
    }

    /// This parses a substring of the given string as a maxquant ion annotation. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// When the annotation does not follow the format.
    pub(crate) fn maxquant_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        custom_database: Option<&CustomDatabase>,
        interpretation: &[(u32, AnalyteTarget)],
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        parse_intermediate_representation(base_context, line, range, custom_database)
            .and_then(|annotations| annotations.into_fragment(interpretation, base_context))
    }
}

fn parse_intermediate_representation<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
    custom_database: Option<&CustomDatabase>,
) -> Result<PeakAnnotation, BoxedError<'a, BasicKind>> {
    // y2;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14;y15;y8-NH3;y30(2+);a2;b2;b3;b4;b5;b6;b7;b18;b25;b2-H2O;b3-H2O;b5-H2O;b14(2+);b22(2+);b26(2+)
    // y2;y3;y4;y7;y10;y12;y13;c2;c5;c6;c10;c14;c32;c39;c40;c41;z°10;z°12;z°14;z°19;z°34;z°36;z°37;z'12;z'13;cm5;cp27;cp31;cp33;cp34;cp35;b2;b3;b4;b5;b6;b7;b10;b11;b27
    // ion type
    let ion = line[range.clone()].chars().next().ok_or_else(|| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment annotation",
            "The fragment annotation is empty",
            base_context.clone().add_highlight((0, range.clone())),
        )
    })?;
    // modifier
    let modifier = line[range.clone()].chars().next().ok_or_else(|| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment annotation",
            "The fragment annotation onoly contains an ion type but needs to contain at least also the series number",
            base_context.clone().add_highlight((0, range.clone())),
        )
    })?;
    let modifier = if modifier.is_ascii_digit() {
        None
    } else {
        Some(modifier)
    };
    // series number
    let num_range = range.start + ion.len_utf8() + modifier.map_or(0, char::len_utf8)..range.end;
    let (length, _, num) =
        next_number::<false, false, usize>(line, num_range.clone()).ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid MaxQuant fragment annotation",
                "There is no series number",
                base_context.clone().add_highlight((0, num_range.clone())),
            )
        })?;
    let num = num.map_err(|err| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment series number",
            explain_number_error(&err),
            base_context.clone().add_highlight((0, num_range.clone())),
        )
    })?;
    // neutral losses
    let (left, losses) = crate::fragment::mzpaf::parse_neutral_loss(
        base_context,
        line,
        num_range.start + length..num_range.end,
    )?;

    // charge
    if line[left.clone()].starts_with('(') && line[left.clone()].ends_with(')') {
        let positive = match line[left.start + 1..left.end - 1].chars().last() {
            Some('+') => true,
            Some('-') => false,
            Some(_) => (),
            None => (),
        };
    }

    todo!()
}
