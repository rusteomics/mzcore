//! Handle positioned glycan structures

use itertools::Itertools;
use mzcore::{
    chemistry::CachedCharge,
    glycan::{GlycanBreakPos, PositionedGlycanStructure},
    prelude::*,
    quantities::Multi,
    system::isize::Charge,
};

use crate::{annotation::model::GlycanModel, fragment::FragmentType, prelude::*};

use uom::num_traits::Zero;

/// Helper trait to be able to define fragmentation on glycan structures.
pub trait GlycanFragmention {
    /// Generate theoretical fragments based on the model.
    fn generate_theoretical_fragments(
        &self,
        model: &FragmentationModel,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment>;
}

/// Get uncharged diagnostic ions from all positions
fn diagnostic_ions(
    glycan: &PositionedGlycanStructure,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    attachment: Option<(AminoAcid, SequencePosition)>,
    model: &GlycanModel,
) -> Vec<Fragment> {
    let mut output = crate::monosaccharide::diagnostic_ions(
        &glycan.sugar,
        peptidoform_ion_index,
        peptidoform_index,
        crate::fragment::DiagnosticPosition::Glycan(
            glycan.position(attachment, false),
            glycan.sugar.clone(),
        ),
        true,
        model,
    );
    output.extend(glycan.branches.iter().flat_map(|b| {
        diagnostic_ions(
            b,
            peptidoform_ion_index,
            peptidoform_index,
            attachment,
            model,
        )
    }));

    output
}

impl GlycanFragmention for PositionedGlycanStructure {
    /// Generate all theoretical fragments for this glycan
    /// * `full_formula` the total formula of the whole peptide + glycan
    fn generate_theoretical_fragments(
        &self,
        model: &FragmentationModel,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        charge_carriers: &mut CachedCharge,
        full_formula: &Multi<MolecularFormula>,
        attachment: Option<(AminoAcid, SequencePosition)>,
    ) -> Vec<Fragment> {
        let charges_other = charge_carriers.range(model.glycan.other_charge_range);
        let charges_oxonium = charge_carriers.range(model.glycan.oxonium_charge_range);
        if model.glycan.allow_structural {
            {
                // Get all B fragments from this node and all its children
                let mut base_fragments =
                    leaf_fragments(self, peptidoform_ion_index, peptidoform_index, attachment)
                        .into_iter()
                        .flat_map(|f| f.with_charge_range_slice(&charges_oxonium))
                        .flat_map(|f| f.with_neutral_losses(&model.glycan.neutral_losses))
                        .collect_vec();
                // Generate all Y fragments
                base_fragments.extend(
                    self.internal_break_points(0, peptidoform_index, attachment)
                        .iter()
                        .filter(|(_, bonds, _)| {
                            bonds.iter().all(|b| !matches!(b, GlycanBreakPos::B(_)))
                                && !bonds.iter().all(|b| matches!(b, GlycanBreakPos::End(_)))
                        })
                        .flat_map(move |(f, bonds, _)| {
                            full_formula.iter().map(move |full| {
                                Fragment::new(
                                    full - self.formula_inner(
                                        SequencePosition::default(),
                                        peptidoform_index,
                                    ) + f,
                                    Charge::zero(),
                                    peptidoform_ion_index,
                                    peptidoform_index,
                                    FragmentType::Y(
                                        bonds
                                            .iter()
                                            .filter(|b| !matches!(b, GlycanBreakPos::End(_)))
                                            .map(GlycanBreakPos::position)
                                            .cloned()
                                            .collect(),
                                    ),
                                )
                            })
                        })
                        .flat_map(|f| f.with_charge_range_slice(&charges_other))
                        .flat_map(|f| f.with_neutral_losses(&model.glycan.neutral_losses)),
                );
                // Generate all diagnostic ions
                base_fragments.extend(
                    diagnostic_ions(
                        self,
                        peptidoform_ion_index,
                        peptidoform_index,
                        attachment,
                        &model.glycan,
                    )
                    .into_iter()
                    .flat_map(|f| f.with_charge_range_slice(&charges_oxonium)),
                );
                base_fragments
            }
        } else {
            Vec::new()
        }
    }
}

/// Generate all fragments without charge and neutral loss options
fn leaf_fragments(
    glycan: &PositionedGlycanStructure,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    attachment: Option<(AminoAcid, SequencePosition)>,
) -> Vec<Fragment> {
    // Find all B type fragments (with and without Y breakage)
    let mut base_fragments = glycan
        .internal_break_points(0, peptidoform_index, attachment)
        .into_iter()
        .filter(|(m, _, _)| *m != MolecularFormula::default())
        .map(|(formula, breakages, _)| {
            Fragment::new(
                formula,
                Charge::zero(),
                peptidoform_ion_index,
                peptidoform_index,
                FragmentType::B {
                    b: glycan.position(attachment, true),
                    y: breakages
                        .iter()
                        .filter(|b| matches!(b, GlycanBreakPos::Y(_)))
                        .map(GlycanBreakPos::position)
                        .cloned()
                        .collect(),
                    end: breakages
                        .iter()
                        .filter(|b| matches!(b, GlycanBreakPos::End(_)))
                        .map(GlycanBreakPos::position)
                        .cloned()
                        .collect(),
                },
            )
        })
        .collect_vec();

    // Extend with the theoretical fragments for all branches of this position
    base_fragments.extend(
        glycan
            .branches
            .iter()
            .flat_map(|b| leaf_fragments(b, peptidoform_ion_index, peptidoform_index, attachment)),
    );
    base_fragments
}
