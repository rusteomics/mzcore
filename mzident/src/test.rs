use itertools::Itertools;

use crate::*;
use mzcore::sequence::Peptidoform;

/// Test a dataset for common errors in PSM parsing
/// # Errors
/// * If the peptide was not identified as the correct version of the format (see parameters).
/// * See errors at [`test_psm`]
#[expect(clippy::missing_panics_doc)]
pub(super) fn test_format<
    'a,
    T: PSMSource + Into<PSM<T::Complexity, T::PeptidoformAvailability>> + 'a,
>(
    reader: impl std::io::Read + 'a,
    ontologies: &'a mzcore::ontology::Ontologies,
    allow_mass_mods: bool,
    expect_lc: bool,
    format: Option<T::Version>,
) -> Result<usize, String>
where
    T::Format: 'static,
    T::Version: std::fmt::Display,
{
    let mut number = 0;
    for peptide in T::parse_reader(reader, ontologies, false, None).map_err(|e| e.to_string())? {
        let peptide: PSM<T::Complexity, T::PeptidoformAvailability> =
            peptide.map_err(|e| e.to_string())?.into();
        number += 1;

        test_psm(&peptide, allow_mass_mods, expect_lc)?;

        if format.as_ref().is_some_and(|f| {
            peptide
                .format()
                .version()
                .is_some_and(|version| f.to_string() != version)
        }) {
            return Err(format!(
                "PSM {} was detected as the wrong version ({} instead of {})",
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

/// Test a peptide for common errors in PSM parsing
/// # Errors
/// * If the local confidence has to be there and is not there (see parameter).
/// * If the local confidence is not the same length as the peptide.
/// * If the score of the peptide is outside of the range -1.0..=1.0.
/// * If any of the local scores is outside of range -1.0..=1.0.
/// * If the peptide contains mass modifications (see parameters).
#[expect(clippy::missing_panics_doc)]
pub(super) fn test_psm<Complexity, PeptidoformAvailability>(
    peptidoform: &PSM<Complexity, PeptidoformAvailability>,
    allow_mass_mods: bool,
    expect_lc: bool,
) -> Result<(), String> {
    let found_length = peptidoform
        .compound_peptidoform_ion()
        .and_then(|p| p.singular_peptidoform_ref().map(Peptidoform::len));
    if found_length != peptidoform.local_confidence().map(|lc| lc.len()) {
        if expect_lc && peptidoform.local_confidence.is_none() {
            return Err(format!(
                "No local confidence was provided for peptidoform {}",
                peptidoform.id()
            ));
        } else if peptidoform.local_confidence.is_some() {
            return Err(format!(
                "The local confidence ({}) does not have the same number of elements as the peptidoform ({}) for peptidoform {}",
                peptidoform.local_confidence().map_or(0, |lc| lc.len()),
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
                            matches!(
                                **m,
                                mzcore::sequence::SimpleModificationInner::Mass(_, _, _)
                            )
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
    if peptidoform.compound_peptidoform_ion().is_some_and(|p| {
        p.peptidoforms()
            .any(|p| !p.enforce_modification_rules().is_empty())
    }) {
        return Err(format!(
            "Peptide {} contains misplaced modifications, sequence {}",
            peptidoform.id(),
            peptidoform.compound_peptidoform_ion().unwrap(),
        ));
    }
    #[cfg(feature = "mzannotate")]
    if let Some(mode) = peptidoform.mode()
        && let Some(model) = peptidoform.fragmentation_model()
    {
        use mzannotate::annotation::model::BuiltInFragmentationModel;

        if (model == BuiltInFragmentationModel::All || model == BuiltInFragmentationModel::None)
            && mode != "MIX"
            && mode != "HCD/ETHCD"
        {
            return Err(format!(
                "Peptide {} fragmentation model '{mode}' was not recognised",
                peptidoform.id()
            ));
        }
    } else if let Some(mode) = peptidoform.mode() {
        return Err(format!(
            "Peptide {} has no fragmentation model but mode is '{mode}'",
            peptidoform.id()
        ));
    } else if let Some(model) = peptidoform.fragmentation_model() {
        return Err(format!(
            "Peptide {} has no mode but fragmentation model is '{model}'",
            peptidoform.id()
        ));
    }
    Ok(())
}

#[ignore = "only run when interested in the sizes of the data"]
#[test]
fn test_size_plgs() {
    let mut file = open_psm_file(
        "src/test_files/plgs_v3.0_peptide.csv",
        &mzcore::ontology::STATIC_ONTOLOGIES,
        false,
    )
    .unwrap()
    .collect::<Result<Vec<_>, _>>()
    .unwrap();
    file.shrink_to_fit();
    let total = mzcore::space::Space::space(&file);
    let first = mzcore::space::Space::space(&file[0]);
    let protein = mzcore::space::Space::space(&file[0].proteins().as_ref()[0]);
    dbg!(&file[0].proteins()[0]);
    println!(
        "= PLGS Total:\n{total}PSMs: {}\n\n= PLGS First:\n{first}\n= Protein First:\n{protein}",
        file.len()
    );
    todo!();
}

#[ignore = "only run when interested in the sizes of the data"]
#[test]
fn test_size_plink() {
    let mut file = open_psm_file(
        "src/test_files/plink_v2.3_small.csv",
        &mzcore::ontology::STATIC_ONTOLOGIES,
        false,
    )
    .unwrap()
    .collect::<Result<Vec<_>, _>>()
    .unwrap();
    file.shrink_to_fit();
    let total = mzcore::space::Space::space(&file);
    let first = mzcore::space::Space::space(&file[0]);
    let protein = mzcore::space::Space::space(&file[0].proteins().as_ref()[0]);
    dbg!(&file[0].proteins()[0]);
    println!(
        "= pLink Total:\n{total}PSMs: {}\n\n= pLink First:\n{first}\n= Protein First:\n{protein}",
        file.len()
    );
    todo!();
}
