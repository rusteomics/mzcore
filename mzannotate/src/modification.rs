use mzcore::{
    chemistry::CachedCharge,
    prelude::{AminoAcid, MolecularFormula, SequencePosition},
    quantities::Multi,
    sequence::{GnoComposition, Modification, SimpleModificationInner},
};

use crate::{
    glycan::GlycanFragmention,
    prelude::{Fragment, FragmentationModel},
};

/// Generate theoretical fragments for side chains (glycans)
pub(crate) fn generate_theoretical_fragments(
    modification: &Modification,
    model: &FragmentationModel,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    charge_carriers: &mut CachedCharge,
    full_formula: &Multi<MolecularFormula>,
    attachment: Option<(AminoAcid, SequencePosition)>,
) -> Vec<Fragment> {
    match modification {
        Modification::Simple(modification) | Modification::Ambiguous { modification, .. } => {
            simple_modification_fragments(
                modification,
                model,
                peptidoform_ion_index,
                peptidoform_index,
                charge_carriers,
                full_formula,
                attachment,
            )
        }
        Modification::CrossLink { .. } => Vec::new(),
    }
}

/// Generate theoretical fragments for side chains (glycans)
pub(crate) fn simple_modification_fragments(
    modification: &SimpleModificationInner,
    model: &FragmentationModel,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    charge_carriers: &mut CachedCharge,
    full_formula: &Multi<MolecularFormula>,
    attachment: Option<(AminoAcid, SequencePosition)>,
) -> Vec<Fragment> {
    match modification {
        SimpleModificationInner::GlycanStructure(glycan)
        | SimpleModificationInner::Gno {
            composition: GnoComposition::Topology(glycan),
            ..
        } => glycan
            .clone()
            .determine_positions()
            .generate_theoretical_fragments(
                model,
                peptidoform_ion_index,
                peptidoform_index,
                charge_carriers,
                full_formula,
                attachment,
            ),
        SimpleModificationInner::Glycan(composition)
        | SimpleModificationInner::Gno {
            composition: GnoComposition::Composition(composition),
            ..
        } => crate::monosaccharide::theoretical_fragments(
            composition,
            model,
            peptidoform_ion_index,
            peptidoform_index,
            charge_carriers,
            full_formula,
            attachment,
        ),
        _ => Vec::new(),
    }
}
