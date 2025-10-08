use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use mzcore::{
    chemistry::{CachedCharge, DiagnosticIon},
    prelude::{
        CompoundPeptidoformIon, MolecularCharge, Peptidoform, PeptidoformIon, SequencePosition,
    },
    quantities::Multi,
    sequence::{
        AtMax, GnoComposition, HiddenInternalMethods, Linear, Linked, LinkerSpecificity,
        PeptidePosition, SimpleModificationInner,
    },
    system::isize::Charge,
};

use crate::{
    annotation::model::get_all_sidechain_losses,
    fragment::{DiagnosticPosition, FragmentType},
    glycan::GlycanFragmention,
    helper_functions::merge_hashmap,
    modification,
    prelude::{Fragment, FragmentationModel},
};

pub trait PeptidoformFragmentation {
    fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &FragmentationModel,
    ) -> Vec<Fragment>;
}

impl PeptidoformFragmentation for CompoundPeptidoformIon {
    /// Generate the theoretical fragments for this compound peptidoform.
    fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &FragmentationModel,
    ) -> Vec<Fragment> {
        let mut base = Vec::new();
        for (index, peptidoform) in self.peptidoform_ions().iter().enumerate() {
            base.extend(peptidoform_ion_inner(peptidoform, max_charge, model, index));
        }
        base
    }
}

impl PeptidoformFragmentation for PeptidoformIon {
    /// Generate the theoretical fragments for this peptidoform.
    fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &FragmentationModel,
    ) -> Vec<Fragment> {
        peptidoform_ion_inner(self, max_charge, model, 0)
    }
}

/// Generate the theoretical fragments for this peptidoform.
fn peptidoform_ion_inner(
    peptidoform_ion: &PeptidoformIon,
    max_charge: Charge,
    model: &FragmentationModel,
    peptidoform_ion_index: usize,
) -> Vec<Fragment> {
    let mut base = Vec::new();
    for (index, peptide) in peptidoform_ion.peptidoforms().iter().enumerate() {
        base.extend(generate_theoretical_fragments_inner(
            peptide,
            max_charge,
            model,
            peptidoform_ion_index,
            index,
            peptidoform_ion.peptidoforms(),
        ));
    }
    base
}

/// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
/// With the global isotope modifications applied.
/// # Panics
/// Panics if the `max_charge` is bigger than [`isize::MAX`].
pub(crate) fn generate_theoretical_fragments_inner<Complexity>(
    peptidoform: &Peptidoform<Complexity>,
    max_charge: Charge,
    model: &FragmentationModel,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    all_peptides: &[Peptidoform<Linked>],
) -> Vec<Fragment> {
    let default_charge = MolecularCharge::proton(max_charge.value);
    let mut charge_carriers: CachedCharge = peptidoform
        .get_charge_carriers()
        .unwrap_or(&default_charge)
        .into();

    let mut output = Vec::with_capacity(20 * peptidoform.sequence().len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
    for sequence_index in 0..peptidoform.sequence().len() {
        let position =
            PeptidePosition::n(SequencePosition::Index(sequence_index), peptidoform.len());
        let mut cross_links = Vec::new();
        let visited_peptides = vec![peptidoform_index];
        let (n_term, n_term_specific, n_term_seen, n_term_losses) = peptidoform.all_masses(
            ..=sequence_index,
            ..sequence_index,
            &peptidoform.get_n_term_mass(
                all_peptides,
                &visited_peptides,
                &mut cross_links,
                model.allow_cross_link_cleavage,
                peptidoform_index,
                peptidoform_ion_index,
                &model.glycan,
            ),
            all_peptides,
            &visited_peptides,
            &mut cross_links,
            model.allow_cross_link_cleavage,
            peptidoform_index,
            peptidoform_ion_index,
            &model.glycan,
        );
        let (c_term, c_term_specific, c_term_seen, c_term_losses) = peptidoform.all_masses(
            sequence_index..,
            sequence_index + 1..,
            &peptidoform.get_c_term_mass(
                all_peptides,
                &visited_peptides,
                &mut cross_links,
                model.allow_cross_link_cleavage,
                peptidoform_index,
                peptidoform_ion_index,
                &model.glycan,
            ),
            all_peptides,
            &visited_peptides,
            &mut cross_links,
            model.allow_cross_link_cleavage,
            peptidoform_index,
            peptidoform_ion_index,
            &model.glycan,
        );
        if !n_term_seen.is_disjoint(&c_term_seen) {
            continue; // There is a link reachable from both sides so there is a loop
        }
        let (modifications_total, modifications_specific, modifications_cross_links) = peptidoform
            .sequence()[sequence_index]
            .modifications
            .iter()
            .fold(
                (Multi::default(), HashMap::new(), HashSet::new()),
                |acc, m| {
                    let (f, specific, s) = m.formula_inner(
                        all_peptides,
                        &[peptidoform_index],
                        &mut cross_links,
                        model.allow_cross_link_cleavage,
                        SequencePosition::Index(sequence_index),
                        peptidoform_index,
                        peptidoform_ion_index,
                        &model.glycan,
                        Some(peptidoform.sequence()[sequence_index].aminoacid.aminoacid()),
                    );
                    (
                        acc.0 * f,
                        merge_hashmap(acc.1, specific),
                        acc.2.union(&s).cloned().collect(),
                    )
                },
            );

        output.append(&mut crate::aminoacid::fragments(
            peptidoform.sequence()[sequence_index].aminoacid.aminoacid(),
            &(n_term, n_term_specific, n_term_losses),
            &(c_term, c_term_specific, c_term_losses),
            &(modifications_total, modifications_specific),
            &mut charge_carriers,
            SequencePosition::Index(sequence_index),
            peptidoform.sequence().len(),
            &model.ions(position, peptidoform),
            peptidoform_ion_index,
            peptidoform_index,
            (
                // Allow any N terminal fragment if there is no cross-link to the C terminal side
                c_term_seen.is_disjoint(&modifications_cross_links),
                n_term_seen.is_disjoint(&modifications_cross_links),
            ),
        ));
    }
    for fragment in &mut output {
        fragment.formula = fragment.formula.as_ref().map(|f| {
            f.with_global_isotope_modifications(peptidoform.get_global())
                .expect("Invalid global isotope modification")
        });
    }

    // Generate precursor peak
    let (full_precursor, _precursor_specific, _all_cross_links) = peptidoform.formulas_inner(
        peptidoform_index,
        peptidoform_ion_index,
        all_peptides,
        &[],
        &mut Vec::new(),
        model.allow_cross_link_cleavage,
        &model.glycan,
    );
    // Allow neutral losses from modifications for the precursor
    let mut precursor_neutral_losses = if model.modification_specific_neutral_losses {
        peptidoform
            .potential_neutral_losses(.., all_peptides, peptidoform_index, &mut Vec::new())
            .into_iter()
            .map(|(n, _, _)| vec![n])
            .collect()
    } else {
        Vec::new()
    };
    // Add amino acid specific neutral losses
    precursor_neutral_losses.extend(
        model
            .precursor
            .1
            .iter()
            .filter_map(|(rule, losses)| {
                rule.iter()
                    .any(|aa| {
                        peptidoform
                            .sequence()
                            .iter()
                            .any(|seq| seq.aminoacid.aminoacid() == *aa)
                    })
                    .then_some(losses)
            })
            .flatten()
            .map(|l| vec![l.clone()]),
    );
    // Add amino acid side chain losses
    precursor_neutral_losses.extend(get_all_sidechain_losses(
        peptidoform.sequence(),
        &model.precursor.2,
    ));
    // Add all normal neutral losses
    precursor_neutral_losses.extend(model.precursor.0.iter().map(|l| vec![l.clone()]));

    output.extend(Fragment::generate_all(
        &full_precursor,
        peptidoform_ion_index,
        peptidoform_index,
        &FragmentType::Precursor,
        &Multi::default(),
        &precursor_neutral_losses,
        &mut charge_carriers,
        model.precursor.3,
    ));

    // Add glycan fragmentation to all peptide fragments
    // Assuming that only one glycan can ever fragment at the same time.
    let full_formula = peptidoform
        .formulas_inner(
            peptidoform_index,
            peptidoform_ion_index,
            all_peptides,
            &[],
            &mut Vec::new(),
            model.allow_cross_link_cleavage,
            &model.glycan,
        )
        .0;
    for (sequence_index, position) in peptidoform.sequence().iter().enumerate() {
        let attachment = (
            position.aminoacid.aminoacid(),
            SequencePosition::Index(sequence_index),
        );
        for modification in &position.modifications {
            output.extend(modification::generate_theoretical_fragments(
                modification,
                model,
                peptidoform_ion_index,
                peptidoform_index,
                &mut charge_carriers,
                &full_formula,
                Some(attachment),
            ));
        }
    }

    if let Some(charge) = model.modification_specific_diagnostic_ions {
        // Add all modification diagnostic ions
        for (dia, pos) in diagnostic_ions(peptidoform) {
            output.extend(
                Fragment {
                    formula: Some(dia.0),
                    charge: Charge::default(),
                    ion: FragmentType::Diagnostic(pos),
                    peptidoform_ion_index: Some(peptidoform_ion_index),
                    peptidoform_index: Some(peptidoform_index),
                    neutral_loss: Vec::new(),
                    deviation: None,
                    confidence: None,
                    auxiliary: false,
                }
                .with_charge_range(&mut charge_carriers, charge),
            );
        }
    }

    // Add labile glycan fragments
    for modification in peptidoform.get_labile() {
        match &**modification {
            SimpleModificationInner::Glycan(composition) => {
                output.extend(crate::monosaccharide::theoretical_fragments(
                    composition,
                    model,
                    peptidoform_ion_index,
                    peptidoform_index,
                    &mut charge_carriers,
                    &full_formula,
                    None,
                ));
            }
            SimpleModificationInner::GlycanStructure(structure)
            | SimpleModificationInner::Gno {
                composition: GnoComposition::Topology(structure),
                ..
            } => {
                output.extend(
                    structure
                        .clone()
                        .determine_positions()
                        .generate_theoretical_fragments(
                            model,
                            peptidoform_ion_index,
                            peptidoform_index,
                            &mut charge_carriers,
                            &full_formula,
                            None,
                        ),
                );
            }
            _ => (),
        }
    }

    output
}

/// Find all diagnostic ions for this full peptide
fn diagnostic_ions<Complexity>(
    peptidoform: &Peptidoform<Complexity>,
) -> Vec<(DiagnosticIon, DiagnosticPosition)> {
    peptidoform
        .iter(..)
        .flat_map(|(pos, aa)| {
            aa.diagnostic_ions(
                pos.sequence_index,
                peptidoform.get_n_term(),
                peptidoform.get_c_term(),
            )
            .into_iter()
            .map(move |diagnostic| {
                (
                    diagnostic,
                    DiagnosticPosition::Peptide(pos, aa.aminoacid.aminoacid()),
                )
            })
        })
        .chain(
            peptidoform
                .get_labile()
                .iter()
                .flat_map(move |modification| match &**modification {
                    SimpleModificationInner::Database { specificities, .. } => specificities
                        .iter()
                        .flat_map(|(_, _, diagnostic)| diagnostic)
                        .map(|diagnostic| {
                            (
                                diagnostic.clone(),
                                DiagnosticPosition::Labile(modification.clone().into()),
                            )
                        })
                        .collect::<Vec<_>>(),
                    SimpleModificationInner::Linker { specificities, .. } => specificities
                        .iter()
                        .flat_map(|rule| match rule {
                            LinkerSpecificity::Symmetric { diagnostic, .. }
                            | LinkerSpecificity::Asymmetric { diagnostic, .. } => diagnostic,
                        })
                        .map(|diagnostic| {
                            (
                                diagnostic.clone(),
                                DiagnosticPosition::Labile(modification.clone().into()),
                            )
                        })
                        .collect::<Vec<_>>(),
                    _ => Vec::new(),
                }),
        )
        .unique()
        .collect()
}

impl<Complexity: AtMax<Linear>> PeptidoformFragmentation for Peptidoform<Complexity> {
    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    /// With the global isotope modifications applied.
    ///
    /// # Panics
    /// If `max_charge` outside the range `1..=u64::MAX`.
    fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &FragmentationModel,
    ) -> Vec<Fragment> {
        generate_theoretical_fragments_inner(self, max_charge, model, 0, 0, &[])
    }
}
