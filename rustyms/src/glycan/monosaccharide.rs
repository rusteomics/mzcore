//! Handle monosaccharides

use itertools::Itertools;

use crate::{
    annotation::model::{FragmentationModel, GlycanModel},
    chemistry::{CachedCharge, Chemical, MolecularFormula},
    fragment::{DiagnosticPosition, Fragment, FragmentType},
    glycan::MonoSaccharide,
    quantities::Multi,
    sequence::{AminoAcid, SequencePosition},
    system::isize::Charge,
};

impl MonoSaccharide {
    /// Generate the composition used for searching on glycans
    pub(crate) fn search_composition(
        composition: &[(Self, isize)],
    ) -> Vec<(MolecularFormula, isize)> {
        // Sort on monosaccharide
        let mut composition = composition
            .iter()
            .filter(|(_, n)| *n != 0)
            .map(|(m, n)| (m.formula(), *n))
            .collect_vec();
        composition.sort_unstable_by(|a, b| a.0.cmp(&b.0));

        // Deduplicate
        let mut max = composition.len().saturating_sub(1);
        let mut index = 0;
        while index < max {
            let this = &composition[index];
            let next = &composition[index + 1];
            if this.0 == next.0 {
                composition[index].1 += next.1;
                composition.remove(index + 1);
                max = max.saturating_sub(1);
            } else {
                index += 1;
            }
        }
        composition.retain(|el| el.1 != 0);
        composition
    }

    /// Generate the theoretical fragments, if any monosaccharide is present a negative number of times no fragments are generated.
    pub(crate) fn theoretical_fragments(
        composition: &[(Self, isize)],
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
        let compositions = Self::composition_options(composition, model.glycan.compositional_range);

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
                        Self::composition_left_over(composition, &fragment_composition),
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
                sugar
                    .diagnostic_ions(
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

    fn composition_left_over(
        composition: &[(Self, isize)],
        remove: &[(Self, isize)],
    ) -> Vec<(Self, isize)> {
        let mut output = Vec::new();
        'outer: for part in composition {
            for other in remove {
                if part.0 == other.0 {
                    if part.1 - other.1 != 0 {
                        output.push((part.0.clone(), part.1 - other.1));
                    }
                    continue 'outer;
                }
            }
            if part.1 != 0 {
                output.push((part.0.clone(), part.1));
            }
        }
        output
    }

    /// Get all unique combinations of monosaccharides within the given range of number of monosaccharides used
    /// # Panics
    /// If any if the composition options has more then [`isize::MAX`] sugars.
    pub fn composition_options(
        composition: &[(Self, isize)],
        range: (Option<usize>, Option<usize>),
    ) -> Vec<Vec<(Self, isize)>> {
        if range.1 == Some(0) {
            return Vec::new();
        }
        let mut options: Vec<Vec<(Self, isize)>> = Vec::new();
        let mut result: Vec<Vec<(Self, isize)>> = Vec::new();
        for sugar in composition {
            let mut new_options = Vec::new();
            // Always start fresh
            for n in 1..=sugar.1 as usize {
                let new = vec![(sugar.0.clone(), isize::try_from(n).unwrap())];
                if range.0.is_none_or(|s| n >= s) && range.1.is_none_or(|s| n <= s) {
                    result.push(new.clone());
                }
                new_options.push(new);
            }
            // And build on previous combinations
            if !options.is_empty() {
                for n in 1..=sugar.1 as usize {
                    for o in &options {
                        let mut new = o.clone();
                        new.push((sugar.0.clone(), isize::try_from(n).unwrap()));
                        let size = new.iter().fold(0, |acc, (_, a)| acc + *a) as usize;
                        match range.1.map(|e| size.cmp(&e)) {
                            Some(std::cmp::Ordering::Greater) => (), // Ignore
                            None | Some(std::cmp::Ordering::Equal) => result.push(new), // Cannot get bigger
                            Some(std::cmp::Ordering::Less) => {
                                if range.0.is_none_or(|s| size >= s) {
                                    result.push(new.clone());
                                }
                                new_options.push(new); // Can get bigger
                            }
                        }
                    }
                }
            }
            options.extend_from_slice(&new_options);
        }
        result
    }

    /// Generate all uncharged diagnostic ions for this monosaccharide.
    pub(crate) fn diagnostic_ions(
        &self,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        position: DiagnosticPosition,
        add_base: bool,
        model: &GlycanModel,
    ) -> Vec<Fragment> {
        let base = Fragment::new(
            self.formula(),
            Charge::default(),
            peptidoform_ion_index,
            peptidoform_index,
            FragmentType::Diagnostic(position),
        );
        model
            .specific_neutral_losses
            .iter()
            .filter(|(ms, precise, _)| ms.equivalent(self, *precise))
            .flat_map(|(_, _, losses)| losses)
            .map(|loss| base.with_neutral_loss(loss))
            .chain(std::iter::repeat_n(base.clone(), usize::from(add_base)))
            .collect()
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use itertools::Itertools;

    use crate::{
        chemistry::Chemical,
        glycan::{
            glycan::{
                BaseSugar, Configuration, GlycanSubstituent, HexoseIsomer, MonoSaccharide,
                PentoseIsomer,
            },
            lists::GLYCAN_PARSE_LIST,
        },
        molecular_formula,
    };

    #[test]
    fn pro_forma_compliance() {
        let cases = &[
            ("hep", molecular_formula!(H 12 C 7 O 6)),
            ("phosphate", molecular_formula!(H 1 O 3 P 1)),
            ("a-hex", molecular_formula!(H 8 C 6 O 6)),
            ("sug", molecular_formula!(H 2 C 2 O 1)),
            ("hexn", molecular_formula!(H 11 C 6 N 1 O 4)),
            ("pen", molecular_formula!(H 8 C 5 O 4)),
            ("tet", molecular_formula!(H 6 C 4 O 3)),
            ("hexp", molecular_formula!(H 11 C 6 O 8 P 1)),
            ("neu5ac", molecular_formula!(H 17 C 11 N 1 O 8)),
            ("non", molecular_formula!(H 16 C 9 O 8)),
            ("hexnac(s)", molecular_formula!(H 13 C 8 N 1 O 8 S 1)),
            ("dec", molecular_formula!(H 18 C 10 O 9)),
            ("en,a-hex", molecular_formula!(H 6 C 6 O 5)),
            ("neu5gc", molecular_formula!(H 17 C 11 N 1 O 9)),
            ("neu", molecular_formula!(H 15 C 9 N 1 O 7)),
            ("hexnac", molecular_formula!(H 13 C 8 N 1 O 5)),
            ("fuc", molecular_formula!(H 10 C 6 O 4)),
            ("hexns", molecular_formula!(H 11 C 6 N 1 O 7 S 1)),
            ("tri", molecular_formula!(H 4 C 3 O 2)),
            ("oct", molecular_formula!(H 14 C 8 O 7)),
            ("sulfate", molecular_formula!(O 3 S 1)),
            ("d-hex", molecular_formula!(H 10 C 6 O 5)),
            ("hex", molecular_formula!(H 10 C 6 O 5)),
            ("hexs", molecular_formula!(H 10 C 6 O 8 S 1)),
        ];
        for (name, formula) in cases {
            assert_eq!(
                GLYCAN_PARSE_LIST
                    .iter()
                    .find(|p| p.0 == *name)
                    .unwrap_or_else(|| panic!("Assumed {name} would be defined"))
                    .1
                    .formula(),
                *formula,
                "{name}",
            );
        }
    }

    #[test]
    fn iupac_short_names() {
        let parse = |str: &str| {
            MonoSaccharide::from_short_iupac(str, 0, 0)
                .map(|(res, len)| {
                    assert_eq!(str.len(), len);
                    res
                })
                .unwrap()
        };
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Pyruvyl,
                ]
            )
        );
        assert_eq!(
            parse("GlcNAc"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[GlycanSubstituent::NAcetyl]
            )
        );
        assert_eq!(
            parse("Gal6S"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[GlycanSubstituent::Sulfate]
            )
        );
        assert_eq!(
            parse("GlcN2Gc"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[GlycanSubstituent::Amino, GlycanSubstituent::Glycolyl,]
            )
        );
        assert_eq!(
            parse("GalNAc3S"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[GlycanSubstituent::NAcetyl, GlycanSubstituent::Sulfate]
            )
        );
        assert_eq!(
            parse("GlcN2,6S2"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[
                    GlycanSubstituent::Amino,
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituent::Sulfate
                ]
            )
        );
        assert_eq!(
            parse("Tagf1,6P2"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)),
                &[GlycanSubstituent::Phosphate, GlycanSubstituent::Phosphate]
            )
            .furanose()
        );
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Pyruvyl,
                ]
            )
        );
        assert_eq!(
            parse("D-Araf"),
            MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)), &[])
                .furanose()
                .configuration(Configuration::D)
        );
        assert_eq!(
            parse("Xyl-onic"),
            MonoSaccharide::new(
                BaseSugar::Pentose(Some(PentoseIsomer::Xylose)),
                &[GlycanSubstituent::Acid]
            )
        );
        assert_eq!(
            parse("Glc2,3,4,6Ac4"),
            MonoSaccharide::new(
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                &[
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituent::Acetyl
                ]
            )
        );
    }

    #[test]
    fn composition_options() {
        let composition = &[
            (MonoSaccharide::new(BaseSugar::Hexose(None), &[]), 1),
            (MonoSaccharide::new(BaseSugar::Heptose(None), &[]), 2),
        ];
        let options_0 = MonoSaccharide::composition_options(composition, (Some(0), Some(0)));
        let options_1 = MonoSaccharide::composition_options(composition, (Some(1), Some(1)));
        let options_2 = MonoSaccharide::composition_options(composition, (Some(2), Some(2)));
        let options_3 = MonoSaccharide::composition_options(composition, (Some(3), Some(3)));
        let options_none_3 = MonoSaccharide::composition_options(composition, (None, Some(3)));
        let options_0_3 = MonoSaccharide::composition_options(composition, (Some(0), Some(3)));
        let human_readable = |options: &[Vec<(MonoSaccharide, isize)>]| {
            options
                .iter()
                .map(|option| {
                    option
                        .iter()
                        .map(|sug| format!("{}{}", sug.0, sug.1))
                        .join("&")
                })
                .join(",")
        };
        assert_eq!(
            human_readable(&options_0),
            "",
            "0 size composition should be empty"
        );
        assert_eq!(
            human_readable(&options_0_3.iter().cloned().sorted().collect_vec()),
            human_readable(
                &options_1
                    .iter()
                    .cloned()
                    .chain(options_2.iter().cloned())
                    .chain(options_3.iter().cloned())
                    .sorted()
                    .collect_vec()
            ),
            "0..=3 should be consistent with 0+1+2+3 separately"
        );
        assert_eq!(&options_0_3, &options_none_3,);
        assert_eq!(human_readable(&options_1), "Hex1,Hep1", "Options 1");
        assert_eq!(human_readable(&options_2), "Hep2,Hex1&Hep1", "Options 2");
        assert_eq!(human_readable(&options_3), "Hex1&Hep2", "Options 3");
    }
}
