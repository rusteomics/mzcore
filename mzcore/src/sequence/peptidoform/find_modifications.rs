use std::collections::HashMap;

use itertools::Itertools;

use crate::{
    chemistry::{Chemical, MassMode, MolecularFormula},
    glycan::MonoSaccharide,
    ontology::{Ontologies, Ontology},
    quantities::{Tolerance, WithinTolerance},
    sequence::{
        AminoAcid, GnoComposition, Modification, Position, SimpleModification,
        SimpleModificationInner,
    },
    system::{Mass, Ratio, ratio::ppm},
};

use super::Peptidoform;

/// Search for modifications that fit the mass tolerance and optional position requirements. If the
/// `positions` is `None` it will not filter for possible positions. Otherwise only modifications
/// that are possible (see [`SimpleModificationInner::is_possible_aa`]) on any of the listed combinations
/// of amino acid and position. If the custom modifications are passed it will also search in them.
///
/// It returns the list of possible modifications.
pub fn modification_search_mass<'a>(
    mass: Mass,
    tolerance: Tolerance<Mass>,
    positions: Option<&'a [(Vec<AminoAcid>, Position)]>,
    mass_mode: MassMode,
    ontologies: &'a Ontologies,
) -> impl Iterator<Item = SimpleModification> + 'a {
    ontologies
        .data(&[])
        .filter(move |m| tolerance.within(&mass, &m.formula().mass(mass_mode)))
        .filter(move |m| {
            positions.is_none_or(|positions| {
                positions.iter().any(|(aas, p)| {
                    aas.iter()
                        .any(|aa| m.is_possible_aa(*aa, *p).any_possible())
                })
            })
        })
}

/// Search for modifications that have the exact same molecular formula as the given target.
///
/// It returns the list of possible modifications.
pub fn modification_search_formula<'a>(
    formula: &'a MolecularFormula,
    ontologies: &'a Ontologies,
) -> impl Iterator<Item = SimpleModification> + 'a {
    ontologies.data(&[]).filter(|m| *formula == m.formula())
}

/// Search for glycans in the GNOme database that have a similar composition. To detect similar
/// composition is converts all monosaccharides into molecular formulas then deduplicate this list.
/// This 'canonical composition' is then compared to the canonical composition for all GNOme
/// modification. Setting `search_topologies` to true allows any GNOme topology modification as
/// well as composition modification.
///
/// It returns the list of possible modifications.
pub fn modification_search_glycan(
    glycan: &[(MonoSaccharide, isize)],
    search_topologies: bool,
    ontologies: &Ontologies,
) -> impl Iterator<Item = SimpleModification> {
    let search = MonoSaccharide::search_composition(glycan);

    ontologies.gnome().data().filter(move |m| {
        if let SimpleModificationInner::Gno {
            composition: GnoComposition::Topology(structure),
            ..
        } = &**m
        {
            search_topologies
                && MonoSaccharide::search_composition(&structure.composition()) == *search
        } else if let SimpleModificationInner::Gno {
            composition: GnoComposition::Composition(composition),
            ..
        } = &**m
        {
            MonoSaccharide::search_composition(composition) == *search
        } else {
            false
        }
    })
}

/// Search for named modifications based on mass and/or chemical formula modifications in a peptide.
/// The struct is intended to be reused if multiple peptides need the same replacement strategy.
#[derive(Clone, Debug)]
pub struct PeptideModificationSearch<'ontologies> {
    /// If true forces the closest if there are multiple modifications within tolerance, if there are two as close it will still not provide any name
    force_closest: bool,
    /// If true also searches for named modifications for formula modifications
    replace_formulas: bool,
    /// Allow modifications specified as side chain on terminal amino acids to be redefined as terminal modifications
    allow_terminal_redefinition: bool,
    /// The mass mode of the mass modifications (ProForma defines mono isotopic)
    mass_mode: MassMode,
    /// The tolerance for mass based matching
    tolerance: Tolerance<Mass>,
    /// The list of manually given modifications
    modifications: Vec<SimpleModification>,
    /// The list of ontologies to search from
    selection: Vec<Ontology>,
    /// The custom modifications, if defined
    ontologies: Option<&'ontologies Ontologies>,
    /// The cache to speed up processing from mod + AA to the replacement mod
    cache: HashMap<(Position, Option<AminoAcid>, SimpleModification), Option<SimpleModification>>,
}

impl Default for PeptideModificationSearch<'_> {
    fn default() -> Self {
        Self {
            force_closest: false,
            replace_formulas: false,
            allow_terminal_redefinition: true,
            mass_mode: MassMode::Monoisotopic,
            tolerance: Tolerance::Relative(Ratio::new::<ppm>(10.0).into()),
            modifications: Vec::new(),
            selection: Vec::new(),
            ontologies: None,
            cache: HashMap::new(),
        }
    }
}

impl<'ontologies> PeptideModificationSearch<'ontologies> {
    /// Search in the specified list of modifications
    pub fn in_modifications(modifications: Vec<SimpleModification>) -> Self {
        Self {
            modifications,
            ..Self::default()
        }
    }

    /// Search in the given ontologies. Do not forget to add [`Ontology::Custom`] if you want to
    /// allow finding modification in the custom database.
    pub fn in_ontologies(selection: Vec<Ontology>, ontologies: &'ontologies Ontologies) -> Self {
        Self {
            selection,
            ontologies: Some(ontologies),
            ..Self::default()
        }
    }

    /// Set the tolerance of matches, default is 10 ppm
    #[must_use]
    pub fn tolerance(self, tolerance: Tolerance<Mass>) -> Self {
        Self {
            tolerance,
            cache: HashMap::new(),
            ..self
        }
    }

    /// Set the mass mode, all mass modifications will be interpreted as this mode, the default is [`MassMode::Monoisotopic`]
    #[must_use]
    pub fn mass_mode(self, mass_mode: MassMode) -> Self {
        Self {
            mass_mode,
            cache: HashMap::new(),
            ..self
        }
    }

    /// Set the strategy for handling one modification with multiple matches, on default (false) it
    /// will not provide any named modification if multiple are within the tolerance, on true the
    /// closest modification will be picked. If multiple modifications are just as close no
    /// modification will be picked.
    #[must_use]
    pub fn force_closest(self, force_closest: bool) -> Self {
        Self {
            force_closest,
            cache: HashMap::new(),
            ..self
        }
    }

    /// Also allow formula modifications to be replaced, defaults to false
    #[must_use]
    pub fn replace_formulas(self, replace_formulas: bool) -> Self {
        Self {
            replace_formulas,
            cache: HashMap::new(),
            ..self
        }
    }

    /// Allow modifications on the side chains of terminal amino acids to be redefined as terminal
    /// modifications, default true. It is very common to see such definitions as `Q[-17.02655]AA`
    /// which is a pyroglutamic acid on the Q at the N terminus, these are supposed to be defined
    /// as `[Gln->pyro-glu]-QAA` in ProForma.
    #[must_use]
    pub fn allow_terminal_redefinition(self, allow_terminal_redefinition: bool) -> Self {
        Self {
            allow_terminal_redefinition,
            cache: HashMap::new(),
            ..self
        }
    }

    /// Search for modifications that can be replaced by named modifications in this peptide.
    #[expect(clippy::similar_names)]
    pub fn search<Complexity>(
        &mut self,
        mut peptide: Peptidoform<Complexity>,
    ) -> Peptidoform<Complexity> {
        fn find_replacement_all_positions(
            settings: &mut PeptideModificationSearch,
            n_term: bool,
            c_term: bool,
            aminoacid: Option<AminoAcid>,
            in_place: &SimpleModification,
        ) -> Option<(SimpleModification, Position)> {
            if !settings.allow_terminal_redefinition || (!n_term && !c_term) {
                settings
                    .find_replacement(Position::Anywhere, aminoacid, in_place)
                    .map(|r| (r, Position::Anywhere))
            } else if n_term && c_term {
                settings
                    .find_replacement(Position::Anywhere, aminoacid, in_place)
                    .map(|r| (r, Position::Anywhere))
                    .or_else(|| {
                        settings
                            .find_replacement(Position::AnyNTerm, aminoacid, in_place)
                            .map(|r| (r, Position::AnyNTerm))
                    })
                    .or_else(|| {
                        settings
                            .find_replacement(Position::AnyCTerm, aminoacid, in_place)
                            .map(|r| (r, Position::AnyCTerm))
                    })
            } else if n_term {
                settings
                    .find_replacement(Position::Anywhere, aminoacid, in_place)
                    .map(|r| (r, Position::Anywhere))
                    .or_else(|| {
                        settings
                            .find_replacement(Position::AnyNTerm, aminoacid, in_place)
                            .map(|r| (r, Position::AnyNTerm))
                    })
            } else {
                // The case when c_term
                settings
                    .find_replacement(Position::Anywhere, aminoacid, in_place)
                    .map(|r| (r, Position::Anywhere))
                    .or_else(|| {
                        settings
                            .find_replacement(Position::AnyCTerm, aminoacid, in_place)
                            .map(|r| (r, Position::AnyCTerm))
                    })
            }
        }

        fn find_replacement_modification(
            settings: &mut PeptideModificationSearch,
            position: Position,
            aminoacid: Option<AminoAcid>,
            in_place: &Modification,
        ) -> Option<Modification> {
            match in_place {
                Modification::Simple(simple) => settings
                    .find_replacement(position, aminoacid, simple)
                    .map(Modification::Simple),
                Modification::CrossLink { .. } | Modification::Ambiguous { .. } => None, // TODO: potentially the cross-linker could be replaced?
            }
        }

        // Start with N and C terminal mods
        let mut n_term = peptide
            .get_n_term()
            .iter()
            .map(|m| {
                find_replacement_modification(
                    self,
                    Position::AnyNTerm,
                    peptide.sequence().first().map(|p| p.aminoacid.aminoacid()),
                    m,
                )
                .unwrap_or_else(|| m.clone())
            })
            .collect_vec();
        let mut c_term = peptide
            .get_c_term()
            .iter()
            .map(|m| {
                find_replacement_modification(
                    self,
                    Position::AnyCTerm,
                    peptide.sequence().last().map(|p| p.aminoacid.aminoacid()),
                    m,
                )
                .unwrap_or_else(|| m.clone())
            })
            .collect_vec();
        let len = peptide.len();

        // Go over all main stretch mods
        for (index, position) in peptide.sequence_mut().iter_mut().enumerate() {
            let is_n_term = index == 0;
            let is_c_term = index == len;
            let mut remove = None;
            for (i, m) in position.modifications.iter_mut().enumerate() {
                match m {
                    Modification::Simple(modification) => {
                        if let Some((replace, location)) = find_replacement_all_positions(
                            self,
                            is_n_term,
                            is_c_term,
                            Some(position.aminoacid.aminoacid()),
                            modification,
                        ) {
                            if location == Position::AnyNTerm {
                                n_term.push(Modification::Simple(replace));
                                remove = Some(i);
                            } else if location == Position::AnyCTerm {
                                c_term.push(Modification::Simple(replace));
                                remove = Some(i);
                            } else if location == Position::Anywhere {
                                *m = Modification::Simple(replace);
                            }
                            // If it can only be a terminal mod but there already is a terminal mod keep it in the original state
                        }
                    }
                    Modification::Ambiguous { .. } | Modification::CrossLink { .. } => (), //TODO: potentially the cross-linker could be replaced as well as the ambiguous mod, but that takes some more logic
                }
            }
            if let Some(remove) = remove.take() {
                position.modifications.remove(remove);
            }
        }
        for m in peptide.get_labile_mut_inner() {
            if let Some(replace) = self.find_replacement(Position::Anywhere, None, m) {
                *m = replace;
            }
        }
        peptide.n_term(n_term).c_term(c_term)
    }

    fn find_replacement(
        &mut self,
        position: Position,
        aminoacid: Option<AminoAcid>,
        in_place: &SimpleModification,
    ) -> Option<SimpleModification> {
        if matches!(&**in_place, SimpleModificationInner::Mass(_))
            || self.replace_formulas && matches!(&**in_place, SimpleModificationInner::Formula(_))
        {
            self.cache
                .entry((position, aminoacid, in_place.clone()))
                .or_insert_with(|| {
                    Self::find_replacement_uncached(
                        self.mass_mode,
                        self.tolerance,
                        self.replace_formulas,
                        self.force_closest,
                        &self.modifications,
                        &self.selection,
                        self.ontologies,
                        position,
                        aminoacid,
                        in_place,
                    )
                })
                .clone()
        } else {
            None
        }
    }

    #[expect(clippy::missing_panics_doc, clippy::too_many_arguments)]
    fn find_replacement_uncached(
        mass_mode: MassMode,
        tolerance: Tolerance<Mass>,
        replace_formulas: bool,
        force_closest: bool,
        modifications: &[SimpleModification],
        selection: &[Ontology],
        ontologies: Option<&Ontologies>,
        position: Position,
        aminoacid: Option<AminoAcid>,
        in_place: &SimpleModification,
    ) -> Option<SimpleModification> {
        let check_matches =
            |in_place: &SimpleModification, provided: &SimpleModification| match &**in_place {
                SimpleModificationInner::Mass(mass) => {
                    tolerance.within(&mass.into_inner(), &provided.formula().mass(mass_mode))
                }
                SimpleModificationInner::Formula(formula) if replace_formulas => {
                    *formula == provided.formula()
                }
                _ => false,
            };

        let options: Vec<_> = if modifications.is_empty() {
            ontologies
                .iter()
                .flat_map(|o| {
                    o.data(selection)
                        .filter(|modification| {
                            aminoacid.is_none_or(|aa| {
                                modification.is_possible_aa(aa, position).any_possible()
                            })
                        })
                        .filter(|modification| check_matches(in_place, modification))
                })
                .collect()
        } else {
            modifications
                .iter()
                .filter(|modification| {
                    aminoacid
                        .is_none_or(|aa| modification.is_possible_aa(aa, position).any_possible())
                })
                .filter(|modification| check_matches(in_place, modification))
                .cloned()
                .collect()
        };
        match options.len() {
            0 => None,
            1 => Some(options[0].clone()),
            _ if !force_closest => None,
            _ => {
                let distances: Vec<_> = options
                    .iter()
                    .map(|m| {
                        in_place
                            .formula()
                            .mass(mass_mode)
                            .ppm(m.formula().mass(mass_mode))
                            .value
                    })
                    .collect();
                let max: f64 = distances
                    .iter()
                    .copied()
                    .max_by(|a, b| (*a).total_cmp(b))
                    .unwrap(); // Guaranteed to always have a value
                let filtered: Vec<_> = distances
                    .iter()
                    .copied()
                    .enumerate()
                    .filter(|(_, d)| (*d - max).abs() < f64::EPSILON)
                    .collect();
                if filtered.len() == 1 {
                    Some(options[filtered[0].0].clone())
                } else {
                    None
                }
            }
        }
    }
}

#[test]
#[expect(clippy::missing_panics_doc)]
fn test_replacement() {
    let mut search = PeptideModificationSearch::in_ontologies(
        vec![Ontology::Unimod],
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .replace_formulas(true);
    let (peptide, _) = Peptidoform::pro_forma(
        "MSFNELT[79.9663]ESNKKSLM[+15.9949]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    let (expected, _) = Peptidoform::pro_forma(
        "MSFNELT[Phospho]ESNKKSLM[Oxidation]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    assert_eq!(search.search(peptide), expected);
    let (peptide, _) = Peptidoform::pro_forma(
        "Q[-17.02655]NKKSLM[+15.9949]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    let (expected, _) = Peptidoform::pro_forma(
        "[Gln->pyro-glu]-QNKKSLM[Oxidation]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    assert_eq!(search.search(peptide), expected);
    let (peptide, _) = Peptidoform::pro_forma(
        "M[Formula:O1]KSLM[+15.9949]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    let (expected, _) = Peptidoform::pro_forma(
        "M[Oxidation]KSLM[Oxidation]E",
        &crate::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap();
    assert_eq!(search.search(peptide), expected);
}
