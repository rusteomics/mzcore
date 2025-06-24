use itertools::Itertools;

use crate::{identification::*, sequence::Peptidoform};

/// Test a dataset for common errors in identified peptide parsing
/// # Errors
/// * If the peptide was not identified as the correct version of the format (see parameters).
/// * See errors at [`test_identified_peptide`]
#[expect(clippy::missing_panics_doc)]
pub(super) fn test_format<
    T: IdentifiedPeptidoformSource
        + Into<IdentifiedPeptidoform<T::Complexity, T::PeptidoformAvailability>>,
>(
    reader: impl std::io::Read,
    custom_database: Option<&crate::ontology::CustomDatabase>,
    allow_mass_mods: bool,
    expect_lc: bool,
    format: Option<T::Version>,
) -> Result<usize, String>
where
    T::Format: 'static,
    T::Version: std::fmt::Display,
{
    let mut number = 0;
    for peptide in
        T::parse_reader(reader, custom_database, false, None).map_err(|e| e.to_string())?
    {
        let peptide: IdentifiedPeptidoform<T::Complexity, T::PeptidoformAvailability> =
            peptide.map_err(|e| e.to_string())?.into();
        number += 1;

        test_identified_peptidoform(&peptide, allow_mass_mods, expect_lc)?;

        if format.as_ref().is_some_and(|f| {
            peptide
                .format()
                .version()
                .is_some_and(|version| f.to_string() != version)
        }) {
            return Err(format!(
                "Peptide {} was detected as the wrong version ({} instead of {})",
                peptide.id(),
                peptide
                    .format()
                    .version()
                    .unwrap_or_else(|| "-".to_string()),
                format.unwrap(),
            ));
        }
    }
    Ok(number)
}

/// Test a peptide for common errors in identified peptide parsing
/// # Errors
/// * If the local confidence has to be there and is not there (see parameter).
/// * If the local confidence is not the same length as the peptide.
/// * If the score of the peptide is outside of the range -1.0..=1.0.
/// * If any of the local scores is outside of range -1.0..=1.0.
/// * If the peptide contains mass modifications (see parameters).
#[expect(clippy::missing_panics_doc)]
pub(super) fn test_identified_peptidoform<Complexity, PeptidoformAvailability>(
    peptidoform: &IdentifiedPeptidoform<Complexity, PeptidoformAvailability>,
    allow_mass_mods: bool,
    expect_lc: bool,
) -> Result<(), String> {
    let found_length = peptidoform
        .compound_peptidoform_ion()
        .and_then(|p| p.singular_peptidoform_ref().map(Peptidoform::len));
    if found_length != peptidoform.local_confidence().map(<[f64]>::len) {
        if expect_lc && peptidoform.local_confidence.is_none() {
            return Err(format!(
                "No local confidence was provided for peptidoform {}",
                peptidoform.id()
            ));
        } else if peptidoform.local_confidence.is_some() {
            return Err(format!(
                "The local confidence ({}) does not have the same number of elements as the peptidoform ({}) for peptidoform {}",
                peptidoform.local_confidence().map_or(0, <[f64]>::len),
                found_length.unwrap_or_default(),
                peptidoform.id()
            ));
        }
    }
    if peptidoform
        .score
        .is_some_and(|s| !(-1.0..=1.0).contains(&s))
    {
        return Err(format!(
            "The score {} for peptidoform {} is outside of range",
            peptidoform.score.unwrap(),
            peptidoform.id()
        ));
    }
    if peptidoform
        .local_confidence
        .as_ref()
        .is_some_and(|s| s.iter().any(|s| !(-1.0..=1.0).contains(s)))
    {
        return Err(format!(
            "The local score {} for peptidoform {} is outside of range",
            peptidoform.local_confidence().unwrap().iter().join(","),
            peptidoform.id()
        ));
    }
    if !allow_mass_mods
        && peptidoform.compound_peptidoform_ion().is_some_and(|p| {
            p.peptidoforms().any(|p| {
                p.sequence().iter().any(|s| {
                    s.modifications.iter().any(|m| {
                        m.simple().is_some_and(|m| {
                            matches!(**m, crate::sequence::SimpleModificationInner::Mass(_))
                        })
                    })
                })
            })
        })
    {
        return Err(format!(
            "Peptidoform {} contains mass modifications, sequence {}",
            peptidoform.id(),
            peptidoform.compound_peptidoform_ion().unwrap(),
        ));
    }
    if let Err(err) = peptidoform.compound_peptidoform_ion().map_or(Ok(()), |p| {
        p.peptidoforms()
            .try_for_each(Peptidoform::enforce_modification_rules)
    }) {
        return Err(format!(
            "Peptide {} contains misplaced modifications, sequence {}\n{}",
            peptidoform.id(),
            peptidoform.compound_peptidoform_ion().unwrap(),
            err
        ));
    }
    Ok(())
}
