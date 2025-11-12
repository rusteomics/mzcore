//! Handle monosaccharides

use itertools::Itertools;

use crate::{
    chemistry::{Chemical, MolecularFormula},
    glycan::MonoSaccharide,
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

    /// What is left over in a composition if the given set of monosaccharides is removed.
    pub fn composition_left_over(
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
    /// If any if the composition options has more than [`isize::MAX`] sugars.
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
            ("Hex", molecular_formula!(H 10 C 6 O 5)),
            ("HexNAc", molecular_formula!(H 13 C 8 N 1 O 5)),
            ("HexS", molecular_formula!(H 10 C 6 O 8 S 1)),
            ("HexP", molecular_formula!(H 11 C 6 O 8 P 1)),
            ("HexNAcS", molecular_formula!(H 13 C 8 N 1 O 8 S 1)),
            ("HexN", molecular_formula!(H 11 C 6 N 1 O 4)),
            ("HexNS", molecular_formula!(H 11 C 6 N 1 O 7 S 1)),
            ("dHex", molecular_formula!(H 10 C 6 O 4)),
            ("aHex", molecular_formula!(H 8 C 6 O 6)),
            ("en,aHex", molecular_formula!(H 6 C 6 O 5)),
            ("Neu", molecular_formula!(H 15 C 9 N 1 O 7)),
            ("NeuAc", molecular_formula!(H 17 C 11 N 1 O 8)),
            ("NeuGc", molecular_formula!(H 17 C 11 N 1 O 9)),
            ("Sug", molecular_formula!(H 2 C 2 O 1)),
            ("Tri", molecular_formula!(H 4 C 3 O 2)),
            ("Tet", molecular_formula!(H 6 C 4 O 3)),
            ("Pen", molecular_formula!(H 8 C 5 O 4)),
            ("Hep", molecular_formula!(H 12 C 7 O 6)),
            ("Oct", molecular_formula!(H 14 C 8 O 7)),
            ("Non", molecular_formula!(H 16 C 9 O 8)),
            ("Dec", molecular_formula!(H 18 C 10 O 9)),
            ("Fuc", molecular_formula!(H 10 C 6 O 4)),
            ("Sulfate", molecular_formula!(O 3 S 1)),
            ("Phosphate", molecular_formula!(H 1 O 3 P 1)),
        ];
        for (name, formula) in cases {
            let found = GLYCAN_PARSE_LIST
                .iter()
                .find(|p| p.0.iter().any(|n| n == name))
                .unwrap_or_else(|| panic!("Assumed {name} would be defined"));
            assert_eq!(found.1.formula(), *formula, "Formula incorrect: {name}",);
            assert_eq!(
                found.1.to_string().as_str(),
                *name,
                "ProForma name set wrong: {name}",
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

    #[test]
    fn out_of_spec() {
        use context_error::FullErrorContent;
        let (res, w) = MonoSaccharide::pro_forma_composition::<false>("Man2ManP").unwrap();
        assert_eq!(
            res,
            MonoSaccharide::pro_forma_composition::<true>("Hex3Phosphate1")
                .unwrap()
                .0
        );
        println!("{w:?}");
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].get_contexts().len(), 2);
        let (res, w) = MonoSaccharide::pro_forma_composition::<false>("Man2Man-1P").unwrap();
        assert_eq!(
            res,
            MonoSaccharide::pro_forma_composition::<true>("Hex1Phosphate1")
                .unwrap()
                .0
        );
        assert_eq!(w.len(), 1);
        // Maybe add warning that you are mixing single letter definitions with full names?
    }
}
