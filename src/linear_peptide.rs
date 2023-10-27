#![warn(dead_code)]

use std::{fmt::Display, ops::RangeBounds};

use crate::{
    error::{Context, CustomError},
    modification::{AmbiguousModification, GlobalModification, ReturnModification},
    molecular_charge::MolecularCharge,
    Element, GnoComposition, MolecularFormula,
};
use itertools::Itertools;
use uom::num_traits::Zero;

use crate::{
    aminoacids::AminoAcid, modification::Modification, system::f64::*, Chemical, Fragment,
    FragmentType, Model,
};

/// A peptide with all data as provided by pro forma. Preferably generated by using the [`PeptideCollection::pro_forma`] function.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct LinearPeptide {
    /// Global isotope modifications, saved as the element and the species that
    /// all occurrence of that element will consist of. Eg (N, 15) will make
    /// all occurring nitrogens be isotope 15.
    pub global: Vec<(Element, u16)>,
    /// Labile modifications, which will not be found in the actual spectrum.
    pub labile: Vec<Modification>,
    /// N terminal modification
    pub n_term: Option<Modification>,
    /// C terminal modification
    pub c_term: Option<Modification>,
    /// The sequence of this peptide (includes local modifications)
    pub sequence: Vec<SequenceElement>,
    /// For each ambiguous modification list all possible positions it can be placed on.
    /// Indexed by the ambiguous modification id.
    pub ambiguous_modifications: Vec<Vec<usize>>,
    /// The adduct ions, if specified
    pub charge_carriers: Option<MolecularCharge>,
}

impl LinearPeptide {
    /// Get the number of amino acids making up this peptide
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if there are any amino acids in this peptide
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// The mass of the N terminal modifications. The global isotope modifications are NOT applied.
    pub fn n_term(&self) -> MolecularFormula {
        self.n_term.as_ref().map_or_else(
            || molecular_formula!(H 1),
            |m| molecular_formula!(H 1) + m.formula(),
        )
    }

    /// The mass of the C terminal modifications. The global isotope modifications are NOT applied.
    pub fn c_term(&self) -> MolecularFormula {
        self.c_term.as_ref().map_or_else(
            || molecular_formula!(H 1 O 1),
            |m| molecular_formula!(H 1 O 1) + m.formula(),
        )
    }

    pub(crate) fn enforce_modification_rules(&self) -> Result<(), CustomError> {
        for (index, element) in self.sequence.iter().enumerate() {
            element.enforce_modification_rules(index, self.sequence.len())?;
        }
        Ok(())
    }

    /// Generate all possible patterns for the ambiguous positions (Mass, String:Label).
    /// It always contains at least one pattern (being (base mass, "")).
    /// The global isotope modifications are NOT applied.
    fn ambiguous_patterns(
        &self,
        aa_range: impl RangeBounds<usize>,
        aa: &[SequenceElement],
        index: usize,
        base: MolecularFormula,
    ) -> Option<Vec<(MolecularFormula, String)>> {
        let result = self
            .ambiguous_modifications
            .iter()
            .enumerate()
            .fold(vec![Vec::new()], |acc, (id, possibilities)| {
                acc.into_iter()
                    .flat_map(|path| {
                        let mut path_clone = path.clone();
                        let options = possibilities
                            .iter()
                            .filter(|pos| aa_range.contains(pos))
                            .map(move |pos| {
                                let mut new = path.clone();
                                new.push((id, *pos));
                                new
                            });
                        options.chain(
                            possibilities
                                .iter()
                                .find(|pos| !aa_range.contains(pos))
                                .map(move |pos| {
                                    path_clone.push((id, *pos));
                                    path_clone
                                }),
                        )
                    })
                    .collect()
            })
            .into_iter()
            .map(|pattern| {
                let ambiguous_local = pattern
                    .iter()
                    .filter_map(|(id, pos)| (*pos == index).then_some(id))
                    .collect::<Vec<_>>();
                aa.iter()
                    .enumerate()
                    .try_fold(MolecularFormula::default(), |acc, (index, aa)| {
                        aa.formula(
                            &pattern
                                .iter()
                                .copied()
                                .filter_map(|(id, pos)| (pos == index).then_some(id))
                                .collect_vec(),
                        )
                        .map(|m| acc + m)
                    })
                    .map(|m| {
                        &base
                            + &m
                            + self.sequence[index]
                                .possible_modifications
                                .iter()
                                .filter(|&am| ambiguous_local.contains(&&am.id))
                                .map(|am| am.modification.formula())
                                .sum::<MolecularFormula>()
                    })
                    .map(|m| {
                        (
                            m,
                            pattern.iter().fold(String::new(), |acc, (id, pos)| {
                                format!(
                                    "{acc}{}{}@{}",
                                    if acc.is_empty() { "" } else { "," },
                                    &self.sequence[index]
                                        .possible_modifications
                                        .iter()
                                        .find(|am| am.id == *id)
                                        .map_or(String::new(), |v| v
                                            .group
                                            .as_ref()
                                            .map_or(id.to_string(), |g| g.0.clone())),
                                    pos + 1
                                )
                            }),
                        )
                    })
            })
            .collect::<Option<Vec<(MolecularFormula, String)>>>()?;
        if result.is_empty() {
            Some(vec![(base, String::new())])
        } else {
            Some(result)
        }
    }

    /// Gives the formula for the whole peptide. With the global isotope modifications applied.
    pub fn formula(&self) -> Option<MolecularFormula> {
        let mut formula = self.n_term() + self.c_term();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        for (_, pos) in self.sequence.iter().enumerate() {
            formula += pos.formula_greedy(&mut placed)?;
        }

        Some(formula.with_global_isotope_modifications(&self.global))
    }

    /// Gives the formula for the whole peptide with no C and N terminal modifications. With the global isotope modifications applied.
    pub fn bare_formula(&self) -> Option<MolecularFormula> {
        let mut formula = MolecularFormula::default();
        let mut placed = vec![false; self.ambiguous_modifications.len()];
        for (_, pos) in self.sequence.iter().enumerate() {
            formula += pos.formula_greedy(&mut placed)?;
        }

        Some(formula.with_global_isotope_modifications(&self.global))
    }

    /// Generate the theoretical fragments for this peptide, with the given maximal charge of the fragments, and the given model.
    /// With the global isotope modifications applied.
    ///
    /// # Panics
    /// If `max_charge` outside the range `1..=u64::MAX`.
    pub fn generate_theoretical_fragments(
        &self,
        max_charge: Charge,
        model: &Model,
        peptide_index: usize,
    ) -> Option<Vec<Fragment>> {
        assert!(max_charge.value >= 1.0);
        assert!(max_charge.value <= u64::MAX as f64);

        let default_charge = MolecularCharge::proton(max_charge.value as isize);
        let charge_carriers = self.charge_carriers.as_ref().unwrap_or(&default_charge);

        let mut output = Vec::with_capacity(20 * self.sequence.len() + 75); // Empirically derived required size of the buffer (Derived from Hecklib)
        for index in 0..self.sequence.len() {
            let n_term =
                self.ambiguous_patterns(0..=index, &self.sequence[0..index], index, self.n_term())?;

            let c_term = self.ambiguous_patterns(
                index..self.sequence.len(),
                &self.sequence[index + 1..self.sequence.len()],
                index,
                self.c_term(),
            )?;

            output.append(
                &mut self.sequence[index].aminoacid.fragments(
                    &n_term,
                    &c_term,
                    &self.sequence[index]
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum(),
                    charge_carriers,
                    index,
                    self.sequence.len(),
                    &model.ions(index, self.sequence.len()),
                    peptide_index,
                ),
            );
        }
        for fragment in &mut output {
            fragment.formula = fragment
                .formula
                .with_global_isotope_modifications(&self.global);
        }

        // Generate precursor peak
        output.extend(
            Fragment::new(
                self.formula()?,
                Charge::zero(),
                peptide_index,
                FragmentType::precursor,
                String::new(),
            )
            .with_charge(charge_carriers)
            .with_neutral_losses(&model.precursor),
        );

        // Add glycan fragmentation to all peptide fragments
        // Assuming that only one glycan can ever fragment at the same time,
        // and that no peptide fragmentation occurs during glycan fragmentation
        for (sequence_index, position) in self.sequence.iter().enumerate() {
            for modification in &position.modifications {
                if let Modification::GlycanStructure(glycan) = modification {
                    output.extend(
                        glycan
                            .clone()
                            .determine_positions()
                            .generate_theoretical_fragments(
                                model,
                                peptide_index,
                                charge_carriers,
                                &self.formula()?,
                                (position.aminoacid, sequence_index),
                            ),
                    );
                } else if let Modification::Gno(GnoComposition::Structure(glycan), _) = modification
                {
                    output.extend(
                        glycan
                            .clone()
                            .determine_positions()
                            .generate_theoretical_fragments(
                                model,
                                peptide_index,
                                charge_carriers,
                                &self.formula()?,
                                (position.aminoacid, sequence_index),
                            ),
                    );
                }
            }
        }

        Some(output)
    }

    pub(crate) fn apply_global_modifications(
        &mut self,
        global_modifications: &[GlobalModification],
    ) {
        let length = self.len();
        for modification in global_modifications {
            match modification {
                GlobalModification::Fixed(aa, modification) => {
                    for (_, seq) in self.sequence.iter_mut().enumerate().filter(|(index, seq)| {
                        seq.aminoacid == *aa && modification.is_possible(seq, *index, length)
                    }) {
                        seq.modifications.push(modification.clone());
                    }
                }
                GlobalModification::Free(modification) => {
                    for (_, seq) in self
                        .sequence
                        .iter_mut()
                        .enumerate()
                        .filter(|(index, seq)| modification.is_possible(seq, *index, length))
                    {
                        seq.modifications.push(modification.clone());
                    }
                }
                GlobalModification::Isotope(el, isotope) => self.global.push((*el, *isotope)),
            }
        }
    }

    /// Place all global unknown positions at all possible locations as ambiguous modifications
    pub(crate) fn apply_unknown_position_modification(
        &mut self,
        unknown_position_modifications: &[Modification],
    ) {
        for modification in unknown_position_modifications {
            let id = self.ambiguous_modifications.len();
            let length = self.len();
            #[allow(clippy::unnecessary_filter_map)]
            // Side effects so the lint does not apply here
            self.ambiguous_modifications.push(
                (0..length)
                    .filter_map(|i| {
                        if modification.is_possible(&self.sequence[i], i, length) {
                            self.sequence[i]
                                .possible_modifications
                                .push(AmbiguousModification {
                                    id,
                                    modification: modification.clone(),
                                    localisation_score: None,
                                    group: None,
                                });
                            Some(i)
                        } else {
                            None
                        }
                    })
                    .collect(),
            );
        }
    }
    /// Place all ranged unknown positions at all possible locations as ambiguous modifications
    pub(crate) fn apply_ranged_unknown_position_modification(
        &mut self,
        ranged_unknown_position_modifications: &[(usize, usize, ReturnModification)],
        ambiguous_lookup: &[(Option<String>, Option<Modification>)],
    ) {
        for (start, end, ret_modification) in ranged_unknown_position_modifications {
            let (id, modification, score, group) = match ret_modification {
                ReturnModification::Defined(def) => {
                    self.ambiguous_modifications.push(Vec::new());
                    (
                        self.ambiguous_modifications.len() - 1,
                        def.clone(),
                        None,
                        None,
                    )
                }
                ReturnModification::Preferred(i, score) => {
                    if *i >= self.ambiguous_modifications.len() {
                        self.ambiguous_modifications.push(Vec::new());
                    }
                    (
                        *i,
                        ambiguous_lookup[*i].1.clone().unwrap(),
                        *score,
                        Some((ambiguous_lookup[*i].0.clone().unwrap(), true)), // TODO: now all possible location in the range are listed as preferred
                    )
                }
                ReturnModification::Referenced(i, score) => {
                    if *i >= self.ambiguous_modifications.len() {
                        self.ambiguous_modifications.push(Vec::new());
                    }
                    (
                        *i,
                        ambiguous_lookup[*i].1.clone().unwrap(),
                        *score,
                        Some((ambiguous_lookup[*i].0.clone().unwrap(), false)),
                    )
                }
            };
            let length = self.len();
            #[allow(clippy::unnecessary_filter_map)]
            // Side effects so the lint does not apply here
            let positions = (*start..=*end)
                .filter_map(|i| {
                    if modification.is_possible(&self.sequence[i], i, length) {
                        self.sequence[i]
                            .possible_modifications
                            .push(AmbiguousModification {
                                id,
                                modification: modification.clone(),
                                localisation_score: None,
                                group: group.clone(),
                            });
                        Some(i)
                    } else {
                        None
                    }
                })
                .collect_vec();
            if let Some(score) = score {
                let individual_score = score / positions.len() as f64;
                for pos in &positions {
                    self.sequence[*pos]
                        .possible_modifications
                        .last_mut()
                        .unwrap()
                        .localisation_score = Some(individual_score);
                }
            }
            self.ambiguous_modifications[id].extend(positions);
        }
    }
}

impl Display for LinearPeptide {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut placed = Vec::new();
        if let Some(m) = &self.n_term {
            write!(f, "[{m}]-")?;
        }
        let mut last_ambiguous = None;
        for position in &self.sequence {
            placed.extend(position.display(f, &placed, last_ambiguous)?);
            last_ambiguous = position.ambiguous;
        }
        if last_ambiguous.is_some() {
            write!(f, ")")?;
        }
        if let Some(m) = &self.c_term {
            write!(f, "-[{m}]")?;
        }
        Ok(())
    }
}

/// One block in a sequence meaning an aminoacid and its accompanying modifications
#[derive(Debug, Clone, PartialEq)]
pub struct SequenceElement {
    /// The aminoacid
    pub aminoacid: AminoAcid,
    /// All present modifications
    pub modifications: Vec<Modification>,
    /// All ambiguous modifications (could be placed here or on another position)
    pub possible_modifications: Vec<AmbiguousModification>,
    /// If this aminoacid is part of an ambiguous sequence group `(QA)?` in pro forma
    pub ambiguous: Option<usize>,
}

impl SequenceElement {
    /// Create a new aminoacid without any modifications
    pub const fn new(aminoacid: AminoAcid, ambiguous: Option<usize>) -> Self {
        Self {
            aminoacid,
            modifications: Vec::new(),
            possible_modifications: Vec::new(),
            ambiguous,
        }
    }

    fn display(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        placed: &[usize],
        last_ambiguous: Option<usize>,
    ) -> Result<Vec<usize>, std::fmt::Error> {
        let mut extra_placed = Vec::new();
        if last_ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, ")")?;
        }
        if self.ambiguous.is_some() && last_ambiguous != self.ambiguous {
            write!(f, "(?")?;
        }
        write!(f, "{}", self.aminoacid.char())?;
        for m in &self.modifications {
            write!(f, "[{m}]")?;
        }
        for m in &self.possible_modifications {
            write!(
                f,
                "[{}#{}{}]",
                m.group.as_ref().map_or(
                    if placed.contains(&m.id) {
                        String::new()
                    } else {
                        extra_placed.push(m.id);
                        m.modification.to_string()
                    },
                    |group| if group.1 {
                        m.modification.to_string()
                    } else {
                        String::new()
                    }
                ),
                m.group
                    .as_ref()
                    .map_or(m.id.to_string(), |g| g.0.to_string()),
                m.localisation_score
                    .map(|v| format!("({v})"))
                    .unwrap_or_default()
            )?;
        }
        Ok(extra_placed)
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the selected ambiguous modifications, without any global isotype modifications
    pub fn formula(&self, selected_ambiguous: &[usize]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .filter(|&m| selected_ambiguous.contains(&m.id))
                        .map(|m| m.modification.formula())
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Get the molecular formula for this position (unless it is B/Z) with the ambiguous modifications placed on the very first placed (and updating this in `placed`), without any global isotype modifications
    pub fn formula_greedy(&self, placed: &mut [bool]) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            #[allow(clippy::filter_map_bool_then)] // otherwise crashes
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .filter_map(|m| {
                            (!placed[m.id]).then(|| {
                                placed[m.id] = true;
                                m.modification.formula()
                            })
                        })
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Get the molecular formula for this position (unless it is B/Z) with all ambiguous modifications, without any global isotype modifications
    pub fn formula_all(&self) -> Option<MolecularFormula> {
        if self.aminoacid == AminoAcid::B || self.aminoacid == AminoAcid::Z {
            None
        } else {
            Some(
                self.aminoacid.formula()
                    + self
                        .modifications
                        .iter()
                        .map(Chemical::formula)
                        .sum::<MolecularFormula>()
                    + self
                        .possible_modifications
                        .iter()
                        .map(|m| m.modification.formula())
                        .sum::<MolecularFormula>(),
            )
        }
    }

    /// Enforce the placement rules of predefined modifications.
    fn enforce_modification_rules(&self, index: usize, length: usize) -> Result<(), CustomError> {
        for modification in &self.modifications {
            if !modification.is_possible(self, index, length) {
                return Err(CustomError::error(
                    "Modification incorrectly placed",
                    format!(
                        "Modification {modification} is not allowed on aminoacid {} index {index}",
                        self.aminoacid.char()
                    ),
                    Context::none(), // TODO: go and give the correct context here
                ));
            }
        }
        Ok(())
    }
}
