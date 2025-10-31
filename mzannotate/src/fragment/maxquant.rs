use std::ops::Range;

use context_error::{BasicKind, BoxedError, Context, CreateError};
use mzcore::{chemistry::SatelliteLabel, prelude::*, sequence::SimpleLinear};

use crate::{
    fragment::{Fragment, IonType, PeakAnnotation},
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
        interpretation: &Peptidoform<SimpleLinear>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        Self::maxquant_inner(
            &Context::none().lines(0, line),
            line,
            0..line.len(),
            interpretation,
        )
    }

    /// This parses a substring of the given string as a maxquant ion annotation. Additionally, this allows
    /// passing a base context to allow to set the line index and source and other properties. Note
    /// that the base context is assumed to contain the full line at line index 0.
    ///
    /// # Errors
    /// When the annotation does not follow the format.
    pub fn maxquant_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
        interpretation: &Peptidoform<SimpleLinear>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        parse_intermediate_representation(base_context, line, range).and_then(|annotations| {
            annotations.into_fragment(
                &[(
                    1,
                    AnalyteTarget::PeptidoformIon(interpretation.clone().into()),
                )],
                base_context,
            )
        })
    }
}

fn parse_intermediate_representation<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<PeakAnnotation, BoxedError<'a, BasicKind>> {
    // ion type
    let ion = line[range.clone()].chars().next().ok_or_else(|| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment annotation",
            "The fragment annotation is empty",
            base_context.clone().add_highlight((0, range.clone())),
        )
    })?;

    let ion = match u8::try_from(ion) {
        Ok(ion @ (b'a' | b'b' | b'c' | b'x' | b'y' | b'z')) => Ok(ion),
        _ => Err(BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment annotation",
            "Only a/b/c/x/y/z ion serires are expected",
            base_context
                .clone()
                .add_highlight((0, range.start..range.start + ion.len_utf8())),
        )),
    }?;

    // modifier
    // TODO: use!
    let modifier = line[range.start+1..range.end].chars().next().ok_or_else(|| {
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
    let num_range = range.start + 1 + modifier.map_or(0, char::len_utf8)..range.end;
    let (length, _, num) =
        next_number::<false, false, usize>(line, num_range.clone()).ok_or_else(|| {
            BoxedError::new(
                BasicKind::Error,
                "Invalid MaxQuant fragment annotation",
                "There is no series number",
                base_context.clone().add_highlight((0, num_range.clone())),
            )
        })?;
    let series_number = num.map_err(|err| {
        BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment series number",
            explain_number_error(&err),
            base_context.clone().add_highlight((0, num_range.clone())),
        )
    })?;

    // neutral losses
    let (left, neutral_losses) = crate::fragment::mzpaf::parse_neutral_loss(
        base_context,
        line,
        num_range.start + length..num_range.end,
    )?;

    // charge
    let charge = if line[left.clone()].starts_with('(') && line[left.clone()].ends_with(')') {
        let positive = match line[left.start + 1..left.end - 1].chars().last() {
            Some('+') => Ok(true),
            Some('-') => Ok(false),
            Some(_) => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid MaxQuant fragment annotation",
                "The charge should be tailed by the sign, for example `(2+)`",
                base_context
                    .clone()
                    .add_highlight((0, left.start + 1..left.end - 1)),
            )),
            None => Err(BoxedError::new(
                BasicKind::Error,
                "Invalid MaxQuant fragment annotation",
                "The charge is empty, it should look like `(2+)`",
                base_context
                    .clone()
                    .add_highlight((0, left.start + 1..left.end - 1)),
            )),
        }?;

        let num = line[left.start + 1..left.end - 2]
            .parse::<u8>()
            .map_err(|err| {
                BoxedError::new(
                    BasicKind::Error,
                    "Invalid MaxQuant fragment charge",
                    explain_number_error(&err),
                    base_context
                        .clone()
                        .add_highlight((0, left.start + 1..left.end - 2)),
                )
            })?;

        if positive {
            MolecularCharge::proton(mzcore::system::isize::Charge::new::<mzcore::system::e>(
                num as isize,
            ))
        } else {
            MolecularCharge::proton(mzcore::system::isize::Charge::new::<mzcore::system::e>(
                num as isize * -1,
            ))
        }
    } else if left.is_empty() {
        MolecularCharge::proton(mzcore::system::isize::Charge::new::<mzcore::system::e>(1))
    } else {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid MaxQuant fragment annotation",
            "Could not parse charge, it should look like `(2+)` but the brackets are not found",
            base_context
                .clone()
                .add_highlight((0, left.start + 1..left.end - 1)),
        ));
    };

    Ok(PeakAnnotation {
        auxiliary: false,
        analyte_number: 1,
        ion: IonType::MainSeries(ion, SatelliteLabel::None, series_number, None),
        neutral_losses,
        isotopes: Vec::new(),
        charge,
        deviation: None,
        confidence: None,
    })
}

#[test]
fn maxquant_annotations() {
    let peptide = Peptidoform::pro_forma(
        "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
        None,
    )
    .unwrap()
    .0
    .into_simple_linear()
    .unwrap();
    for annotation in "y2;y5;y6;y7;y8;y9;y10;y11;y12;y13;y14;y15;y8-NH3;y30(2+);a2;b2;b3;b4;b5;b6;b7;b18;b25;b2-H2O;b3-H2O;b5-H2O;b14(2+);b22(2+);b26(2+);y2;y3;y4;y7;y10;y12;y13;c2;c5;c6;c10;c14;c32;c39;c40;c41;z°10;z°12;z°14;z°19;z°34;z°36;z°37;z'12;z'13;cm5;cp27;cp31;cp33;cp34;cp35;b2;b3;b4;b5;b6;b7;b10;b11;b27".split(';') {
        Fragment::maxquant(annotation, &peptide).unwrap();
    }
}
