use mzcore::{
    chemistry::CachedCharge,
    glycan::MonoSaccharide,
    prelude::{AminoAcid, Chemical, MolecularFormula, SequencePosition},
    quantities::Multi,
    system::isize::Charge,
};

use crate::{
    annotation::model::GlycanModel,
    fragment::{DiagnosticPosition, FragmentType},
    monosaccharide,
    prelude::{Fragment, FragmentationModel},
};

/// Generate the theoretical fragments, if any monosaccharide is present a negative number of times no fragments are generated.
pub(crate) fn theoretical_fragments(
    composition: &[(MonoSaccharide, isize)],
    model: &FragmentationModel,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    charge_carriers: &mut CachedCharge,
    full_formula: &Multi<MolecularFormula>,
    attachment: Option<(AminoAcid, SequencePosition)>,
) -> Vec<Fragment> {
    if composition.iter().any(|(_, a)| u16::try_from(*a).is_err()) {
        // u16: negative + also ensure it fits within the bounds of the molecular formula structure
        return Vec::new();
    }
    let mut fragments = Vec::new();
    let compositions =
        MonoSaccharide::composition_options(composition, model.glycan.compositional_range);

    // Generate compositional B and Y ions
    let charges_other = charge_carriers.range(model.glycan.other_charge_range);
    let charges_oxonium = charge_carriers.range(model.glycan.oxonium_charge_range);
    for fragment_composition in compositions {
        let formula: MolecularFormula = fragment_composition
            .iter()
            .map(|s| {
                s.0.formula_inner(SequencePosition::default(), peptidoform_index) * s.1 as i32
            })
            .sum();
        fragments.extend(
            Fragment::new(
                formula.clone(),
                Charge::default(),
                peptidoform_ion_index,
                peptidoform_index,
                FragmentType::BComposition(fragment_composition.clone(), attachment),
            )
            .with_charge_range_slice(&charges_oxonium)
            .flat_map(|o| o.with_neutral_losses(&model.glycan.neutral_losses)),
        );

        fragments.extend(full_formula.to_vec().iter().flat_map(|base| {
            Fragment::new(
                base - &formula,
                Charge::default(),
                peptidoform_ion_index,
                peptidoform_index,
                FragmentType::YComposition(
                    MonoSaccharide::composition_left_over(composition, &fragment_composition),
                    attachment,
                ),
            )
            .with_charge_range_slice(&charges_other)
            .flat_map(|o| o.with_neutral_losses(&model.glycan.neutral_losses))
        }));
    }

    // Generate compositional diagnostic ions
    for (sugar, _) in composition {
        fragments.extend(
            diagnostic_ions(
                sugar,
                peptidoform_ion_index,
                peptidoform_index,
                DiagnosticPosition::GlycanCompositional(sugar.clone(), attachment),
                false,
                &model.glycan,
            )
            .into_iter()
            .flat_map(|d| d.with_charge_range_slice(&charges_oxonium)),
        );
    }

    fragments
}

/// Generate all uncharged diagnostic ions for this monosaccharide.
pub(crate) fn diagnostic_ions(
    monosaccharide: &MonoSaccharide,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    position: DiagnosticPosition,
    add_base: bool,
    model: &GlycanModel,
) -> Vec<Fragment> {
    let base = Fragment::new(
        monosaccharide.formula(),
        Charge::default(),
        peptidoform_ion_index,
        peptidoform_index,
        FragmentType::Diagnostic(position),
    );
    model
        .specific_neutral_losses
        .iter()
        .filter(|(ms, precise, _)| ms.equivalent(monosaccharide, *precise))
        .flat_map(|(_, _, losses)| losses)
        .map(|loss| base.with_neutral_loss(loss))
        .chain(std::iter::repeat_n(base.clone(), usize::from(add_base)))
        .collect()
}
