use std::{collections::HashMap, sync::LazyLock};

use mzcore::{
    chemistry::{CachedCharge, NeutralLoss},
    glycan::BackboneFragmentKind,
    molecular_formula,
    prelude::{AminoAcid, IsAminoAcid, MolecularFormula, MultiChemical, SequencePosition},
    quantities::Multi,
    sequence::{BACKBONE, PeptidePosition},
};

use crate::{annotation::model::PossibleIons, fragment::FragmentType, prelude::Fragment};

/// The information from the N and C terminal to properly generate fragments
type TerminalInfo = (
    Multi<MolecularFormula>,
    HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
    Vec<Vec<NeutralLoss>>,
);

// TODO: generalise over used storage type, so using molecularformula, monoisotopic mass, or average mass, also make sure that AAs can return these numbers in a const fashion
#[expect(clippy::too_many_lines, clippy::too_many_arguments)]
pub(crate) fn fragments(
    aminoacid: AminoAcid,
    n_term: &TerminalInfo,
    c_term: &TerminalInfo,
    modifications: &(
        Multi<MolecularFormula>,
        HashMap<BackboneFragmentKind, Multi<MolecularFormula>>,
    ),
    charge_carriers: &mut CachedCharge,
    sequence_index: SequencePosition,
    sequence_length: usize,
    ions: &PossibleIons,
    peptidoform_ion_index: usize,
    peptidoform_index: usize,
    allow_terminal: (bool, bool),
) -> Vec<Fragment> {
    let mut base_fragments = Vec::with_capacity(ions.size_upper_bound());
    let n_pos = PeptidePosition::n(sequence_index, sequence_length);
    let c_pos = PeptidePosition::c(sequence_index, sequence_length);

    if allow_terminal.0 {
        if let Some(settings) = &ions.a {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::a)
                    .unwrap_or(&modifications.0)
                    - molecular_formula!(H 1 C 1 O 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::a(n_pos, 0),
                n_term.1.get(&BackboneFragmentKind::a).unwrap_or(&n_term.0),
                &n_term.2,
                charge_carriers,
                settings,
            ));
        }
        if let Some(settings) = &ions.b {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::b)
                    .unwrap_or(&modifications.0)
                    - molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::b(n_pos, 0),
                n_term.1.get(&BackboneFragmentKind::b).unwrap_or(&n_term.0),
                &n_term.2,
                charge_carriers,
                settings,
            ));
        }
        if let Some(settings) = &ions.c {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::c)
                    .unwrap_or(&modifications.0)
                    + molecular_formula!(H 2 N 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::c(n_pos, 0),
                n_term.1.get(&BackboneFragmentKind::c).unwrap_or(&n_term.0),
                &n_term.2,
                charge_carriers,
                settings,
            ));
        }
        for (aa, distance) in &ions.d.0 {
            if let Some(satellite_fragments) = aa.satellite_ion_fragments(
                sequence_index - *distance,
                peptidoform_index,
                peptidoform_ion_index,
            ) {
                for (label, formula) in satellite_fragments.iter() {
                    base_fragments.extend(Fragment::generate_series(
                        &(modifications
                            .1
                            .get(&BackboneFragmentKind::d)
                            .unwrap_or(&modifications.0)
                            * aminoacid.formulas_inner(
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            )
                            + molecular_formula!(H 1 C 1 O 1)
                            - formula),
                        peptidoform_ion_index,
                        peptidoform_index,
                        &FragmentType::d(n_pos, *aa, *distance, 0, *label),
                        n_term.1.get(&BackboneFragmentKind::d).unwrap_or(&n_term.0),
                        &n_term.2,
                        charge_carriers,
                        &ions.d.1,
                    ));
                }
            }
        }
    }
    if allow_terminal.1 {
        for (aa, distance) in &ions.v.0 {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * -aa.formulas_inner(
                    sequence_index + *distance,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) + LazyLock::force(&BACKBONE)),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::v(c_pos, *aa, *distance, 0),
                c_term.1.get(&BackboneFragmentKind::v).unwrap_or(&c_term.0),
                &c_term.2,
                charge_carriers,
                &ions.v.1,
            ));
        }
        for (aa, distance) in &ions.w.0 {
            if let Some(satellite_fragments) = aa.satellite_ion_fragments(
                sequence_index - *distance,
                peptidoform_index,
                peptidoform_ion_index,
            ) {
                for (label, formula) in satellite_fragments.iter() {
                    base_fragments.extend(Fragment::generate_series(
                        &(modifications
                            .1
                            .get(&BackboneFragmentKind::w)
                            .unwrap_or(&modifications.0)
                            * aminoacid.formulas_inner(
                                sequence_index,
                                peptidoform_index,
                                peptidoform_ion_index,
                            )
                            + molecular_formula!(H 2 N 1)
                            - formula),
                        peptidoform_ion_index,
                        peptidoform_index,
                        &FragmentType::w(c_pos, *aa, *distance, 0, *label),
                        c_term.1.get(&BackboneFragmentKind::w).unwrap_or(&c_term.0),
                        &c_term.2,
                        charge_carriers,
                        &ions.w.1,
                    ));
                }
            }
        }
        if let Some(settings) = &ions.x {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::x)
                    .unwrap_or(&modifications.0)
                    + molecular_formula!(C 1 O 1)
                    - molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::x(c_pos, 0),
                c_term.1.get(&BackboneFragmentKind::x).unwrap_or(&c_term.0),
                &c_term.2,
                charge_carriers,
                settings,
            ));
        }
        if let Some(settings) = &ions.y {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::y)
                    .unwrap_or(&modifications.0)
                    + molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::y(c_pos, 0),
                c_term.1.get(&BackboneFragmentKind::y).unwrap_or(&c_term.0),
                &c_term.2,
                charge_carriers,
                settings,
            ));
        }
        if let Some(settings) = &ions.z {
            base_fragments.extend(Fragment::generate_series(
                &(aminoacid.formulas_inner(
                    sequence_index,
                    peptidoform_index,
                    peptidoform_ion_index,
                ) * (modifications
                    .1
                    .get(&BackboneFragmentKind::z)
                    .unwrap_or(&modifications.0)
                    - molecular_formula!(H 2 N 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::z(c_pos, 0),
                c_term.1.get(&BackboneFragmentKind::z).unwrap_or(&c_term.0),
                &c_term.2,
                charge_carriers,
                settings,
            ));
        }
    }

    if allow_terminal.0
        && allow_terminal.1
        && let Some((charge, losses)) = &ions.immonium
    {
        base_fragments.extend(Fragment::generate_all(
            &(aminoacid.formulas_inner(sequence_index, peptidoform_index, peptidoform_ion_index)
                * &modifications.0
                - molecular_formula!(C 1 O 1)),
            peptidoform_ion_index,
            peptidoform_index,
            &FragmentType::Immonium(Some(n_pos), aminoacid.into()), // TODO: get the actual sequence element here
            &Multi::default(),
            &losses
                .iter()
                .filter(|(aa, _)| aa.contains(&aminoacid))
                .flat_map(|(_, l)| l.iter())
                .map(|l| vec![l.clone()])
                .collect::<Vec<_>>(),
            charge_carriers,
            *charge,
        ));
    }
    base_fragments
}
