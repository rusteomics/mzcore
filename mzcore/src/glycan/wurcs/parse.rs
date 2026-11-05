use std::num::NonZeroU8;

use crate::{
    chemistry::Element,
    glycan::{
        BaseSugar, GlycanStructure, GlycanSubstituent, GlycanSubstituentLeavingGroup,
        GlycanSubstituentLocation, HexoseIsomer, MonoSaccharide,
        wurcs::structs::{
            BackBone, Carbon, GLIPOption, LIPOption, Linkage, Lip, MAPSymbol, Residue,
            TerminalCarbon, Wurcs,
        },
    },
};

/// The possible errors that can be returned from parsing a WURCS line
#[derive(Clone, Copy, Debug)]
pub enum WurcsParseError {
    /// The definition is empty
    Empty,
    /// It contains a repeating backbone, this is valid WURCS, but cannot be represented currently
    /// in a [`MonoSaccharide`].
    RepeatingBackbone,
    /// The backbone is too long for the current limitations of [`MonoSaccharide`].
    BackboneTooLong,
    /// The linkage is defined with a repeating pattern, this is not currently supported.
    RepeatedLinkage,
    /// The linkage is defined in an ambiguous manner, this is not currently supported
    AmbiguousLinkage,
    /// There is a cycle in the linkage
    Cycle,
    /// There are unconnected residues
    UnconnectedResidues,
    InvalidResidueCode,
}

impl Wurcs {
    /// Parse a tokenised WURCS line into [`MonoSaccharide`]s, note that this is very much work in
    /// progress.
    /// # Errors
    /// If these tokens cannot be parsed into a [`MonoSaccharide`].
    pub fn parse(self) -> Result<GlycanStructure, WurcsParseError> {
        let residues = self
            .residues
            .into_iter()
            .map(Residue::parse)
            .collect::<Result<Vec<_>, _>>()?;

        let mut connections = vec![Vec::new(); self.sequence.len()];

        for link in self.linkage {
            match link {
                Linkage::Known(link) => {
                    if link.lips.len() == 2 {
                        let a = match link.lips[0] {
                            GLIPOption::Known(glip) => glip.res_index,
                            _ => return Err(WurcsParseError::AmbiguousLinkage),
                        };
                        let b = match link.lips[1] {
                            GLIPOption::Known(glip) => glip.res_index,
                            _ => return Err(WurcsParseError::AmbiguousLinkage),
                        };
                        connections[a as usize].push(b);
                    }
                }
                Linkage::Repeated(..) => return Err(WurcsParseError::RepeatedLinkage),
            }
        }

        let mut added = vec![false; self.sequence.len()];

        Self::make_tree(0, &residues, &self.sequence, &connections, &mut added)
    }

    fn make_tree(
        residue: u8,
        residues: &[MonoSaccharide],
        sequence: &[u8],
        connections: &[Vec<u8>],
        added: &mut [bool],
    ) -> Result<GlycanStructure, WurcsParseError> {
        let residue = residue as usize;
        if residue >= added.len() {
            Err(WurcsParseError::InvalidResidueCode)
        } else if added[residue] {
            Err(WurcsParseError::Cycle)
        } else {
            Ok(GlycanStructure {
                sugar: residues[(sequence[residue] - 1) as usize].clone(),
                branches: connections[residue]
                    .iter()
                    .map(|c| Self::make_tree(*c, residues, sequence, connections, added))
                    .collect::<Result<thin_vec::ThinVec<_>, _>>()?,
            })
        }
    }
}

impl Residue {
    /// Parse a WURCS residue into a [`MonoSaccharide`].
    /// # Errors
    /// If this WURCS definition does not fit in the limitations of a [`MonoSaccharide`].
    pub(super) fn parse(self) -> Result<MonoSaccharide, WurcsParseError> {
        let (start, base, end, middle_length) = match self.backbone {
            BackBone::Defined(s, m, e) => (
                s,
                match m.len() {
                    0 => BaseSugar::Sugar,
                    1 => BaseSugar::Triose,
                    2 => BaseSugar::Tetrose(None),
                    3 => BaseSugar::Pentose(None),
                    4 => BaseSugar::Hexose(match m.as_slice() {
                        [
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                        ] => Some(HexoseIsomer::Glucose),
                        [
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                        ] => Some(HexoseIsomer::Mannose),
                        [
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                        ]
                        | [
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                        ] => Some(HexoseIsomer::Galactose),
                        [
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                        ] => Some(HexoseIsomer::Gulose),
                        [
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                        ] => Some(HexoseIsomer::Altrose),
                        [
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyRight,
                        ] => Some(HexoseIsomer::Allose),
                        [
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                        ] => Some(HexoseIsomer::Talose),
                        [
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                            Carbon::HydroxyRight,
                            Carbon::HydroxyLeft,
                        ] => Some(HexoseIsomer::Idose),
                        _ => None,
                    }),
                    5 => BaseSugar::Heptose(None),
                    6 => BaseSugar::Octose,
                    7 => BaseSugar::Nonose(None),
                    8 => BaseSugar::Decose,
                    _ => return Err(WurcsParseError::BackboneTooLong),
                },
                e,
                m.len() as u8,
            ),
            BackBone::Repeating(..) => return Err(WurcsParseError::RepeatingBackbone),
        };
        let mut res = MonoSaccharide::new(base, &[]);
        // Assumes a CHHO as standard
        match start {
            TerminalCarbon::CHHX
            | TerminalCarbon::AldehydeOrHemiacetal
            | TerminalCarbon::Hemiacetal
            | TerminalCarbon::CXH => (),
            TerminalCarbon::CHHH => res.substituents.push((
                GlycanSubstituent::Deoxy,
                GlycanSubstituentLocation::Single(
                    NonZeroU8::new(1),
                    GlycanSubstituentLeavingGroup::Hydroxy,
                ),
            )),
            o => println!("Not added yet: {o:?}"),
        }
        // Assumes a C=OH
        match end {
            TerminalCarbon::CXH => (),
            TerminalCarbon::CHHH => res.substituents.push((
                GlycanSubstituent::Deoxy,
                GlycanSubstituentLocation::Single(
                    NonZeroU8::new(2 + middle_length),
                    GlycanSubstituentLeavingGroup::Hydroxy,
                ),
            )),

            // TerminalCarbon::CHHX => res.substituents.push((
            //     GlycanSubstituent::Alcohol,
            //     GlycanSubstituentLocation::Single(
            //         NonZeroU8::new(2 + middle_length),
            //         GlycanSubstituentLeavingGroup::Hydroxy,
            //     ),
            // )),
            o => println!("Not added yet: {o:?}"),
        }
        let mut internal_cycle_count = 0;

        for m in self.mods {
            if m.lips.len() == 1 {
                match m.modification.as_slice() {
                    &[] => (),
                    // NCC/3=O
                    &[
                        MAPSymbol::Star(None),
                        MAPSymbol::Element(Element::N),
                        MAPSymbol::Element(Element::C),
                        MAPSymbol::Element(Element::C),
                        MAPSymbol::Branch(3),
                        MAPSymbol::DoubleBond,
                        MAPSymbol::Element(Element::O),
                    ] => res.substituents.push((GlycanSubstituent::NAcetyl, match m.lips[0] {
                        LIPOption::Known(Lip { position, .. }) => {
                            GlycanSubstituentLocation::Single(
                                position.and_then(NonZeroU8::new),
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            )
                        }
                        _ => GlycanSubstituentLocation::Single(
                            None,
                            GlycanSubstituentLeavingGroup::Hydroxy,
                        ),
                    })),
                    other => {
                        dbg!(other);
                        // todo!()
                    }
                }
            } else if m.lips.len() == 2 {
                if m.modification.is_empty() {
                    internal_cycle_count += 1;
                } else {
                    // todo!()
                }
            } else {
                // todo!()
            }
        }

        if internal_cycle_count == 0 {
            // res.substituents.push((GlycanSubstituent::Alcohol, None));
        }
        if internal_cycle_count > 1 {
            // todo!();
        }
        res.substituents.sort_unstable();
        #[cfg(debug_assertions)]
        {
            res.source_field = self.source_field;
        }
        Ok(res)
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc, clippy::missing_errors_doc)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::{
        chemistry::Molecule,
        glycan::{
            BaseSugar, GlycanStructure, HexoseIsomer, MonoSaccharide, tokenise_wurcs,
            wurcs::structs::Wurcs,
        },
        molecular_formula,
    };

    fn test_tokenise(value: &str) -> Result<Wurcs, BoxedError<'_, BasicKind>> {
        tokenise_wurcs(value, &Context::default().lines(0, value), ..)
    }

    #[ignore = "WURCS monosaccharides parsing is in its infancy"]
    #[test]
    fn monosaccharides() {
        // G58645HA,"Gal(b1-3)GalNAc-ol","WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/
        // a3-b1"

        for (wurcs, sugar, proforma, iupac, formula) in [
            // (
            //     "a2112h-1b_1-5",
            //     BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
            //     "?",
            //     "Gal(b1- ",
            //     molecular_formula!(C 6 H 12 O 6),
            // ),
            (
                "h2112h_2*NCC/3=O",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "?",
                "Gal2NAc-6-ol",
                molecular_formula!(C 8 H 17 N 1 O 6),
            ),
            (
                "h2122h_2*NCC/3=O",
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                "?",
                "Glc2NAc-ol",
                molecular_formula!(C 8 H 17 N 1 O 6),
            ),
            (
                "o3444h",
                BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
                "Hex",
                "?",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u1112h",
                BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
                "Hex",
                "Tal",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u1121h",
                BaseSugar::Hexose(Some(HexoseIsomer::Gulose)),
                "?",
                "L-Gul",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u1122h",
                BaseSugar::Hexose(Some(HexoseIsomer::Mannose)),
                "Hex",
                "Man",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u1221m",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "Fuc",
                "Fuc",
                molecular_formula!(C 6 H 12 O 5),
            ),
            (
                "u2111h",
                BaseSugar::Hexose(Some(HexoseIsomer::Altrose)),
                "Hex",
                "Alt",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u2112h",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "Hex",
                "Gal",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u2112m",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "?",
                "D-Fuc",
                molecular_formula!(C 6 H 12 O 5),
            ),
            (
                "u2121h",
                BaseSugar::Hexose(Some(HexoseIsomer::Idose)),
                "Hex",
                "Ido",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u2122h",
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                "?",
                "Glc",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u2122h_2*NCC/3=",
                BaseSugar::Hexose(Some(HexoseIsomer::Glucose)),
                "HexNAc",
                "GlcNAc",
                molecular_formula!(C 8 H 15 N 1 O 6),
            ),
            (
                "u2222h",
                BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
                "Hex",
                "All",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "u2222h_2-3",
                BaseSugar::Hexose(Some(HexoseIsomer::Allose)),
                "?",
                "2,3-Anhydro-All ",
                molecular_formula!(C 6 H 10 O 5),
            ),
            (
                "u3334h",
                BaseSugar::Hexose(Some(HexoseIsomer::Talose)),
                "Hex",
                "?",
                molecular_formula!(C 6 H 12 O 6),
            ),
        ] {
            let full_wurcs = format!("WURCS=2.0/1,1,0/[{wurcs}]/1/");
            let parsed = test_tokenise(&full_wurcs).unwrap().parse().unwrap();
            let monosaccharide = parsed.sugar.clone();
            assert_eq!(
                monosaccharide.base_sugar, sugar,
                "{wurcs} not detected as {sugar} but as {}",
                monosaccharide.base_sugar
            );
            if proforma != "?" {
                let v = monosaccharide.pro_forma_name();
                assert_eq!(v, proforma, "{wurcs} not detected as {proforma} but as {v}");
            }
            if iupac != "?" {
                let from_iupac = MonoSaccharide::from_short_iupac(iupac, 0, 0).unwrap().0;
                assert_eq!(
                    monosaccharide, from_iupac,
                    "{wurcs} not the same as {iupac} ({monosaccharide:?} != {from_iupac:?})"
                );
            }
            assert_eq!(
                monosaccharide.formula(),
                formula,
                "{wurcs} not detected as having formula {formula} instead it has {}",
                monosaccharide.formula()
            );
        }
    }

    #[test]
    #[cfg(feature = "glycan-render")]
    fn structure() {
        use std::io::Write;

        use crate::glycan::{GlycanModificationDetails, GlycanRenderSettings};

        let mut output = std::io::BufWriter::new(
            std::fs::File::create("/home/douwe/src/mzcore/glycan_test.html").unwrap(),
        );

        write!(
            &mut output,
            "<html><body style='color:white;background:black'>"
        )
        .unwrap();

        for (code, iupac, wurcs) in [
            (
                "G58645HA",
                "Gal(b1-3)GalNAc-ol",
                "WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1",
            ),
            (
                "G00037MO",
                "GlcNAc(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",
                "WURCS=2.0/2,3,2/[a2112h-1a_1-5_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O]/1-2-2/a3-b1_a6-c1",
            ),
            (
                "G00514XC",
                "Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-4)]Man(a1-3)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-4)Gal(b1-4)GlcNAc(b1-2)[Neu5Ac(a2-6)Gal(b1-4)GlcNAc(b1-6)]Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc(b1-",
                "WURCS=2.0/6,20,19/[a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-1-2-3-1-4-5-1-4-5-3-1-4-1-4-5-1-4-5-6/a4-b1_a6-t1_b4-c1_c3-d1_c6-k1_d2-e1_d4-h1_e4-f1_f6-g2_h4-i1_i6-j2_k2-l1_k6-q1_l4-m1_m4-n1_n4-o1_o6-p2_q4-r1_r6-s2",
            ),
            (
                "G46474UX",
                "Neu5Ac(a2-6)Gal(b1-?)GlcNAc(b1-2)Man(a1-3)[Gal(b1-?)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)[Fuc(a1-6)]GlcNAc-ol",
                "WURCS=2.0/7,11,10/[h2122h_2*NCC/3=O][a2122h-1b_1-5_2*NCC/3=O][a1122h-1b_1-5][a1122h-1a_1-5][a2112h-1b_1-5][Aad21122h-2a_2-6_5*NCC/3=O][a1221m-1a_1-5]/1-2-3-4-2-5-6-4-2-5-7/a4-b1_a6-k1_b4-c1_c3-d1_c6-h1_d2-e1_f6-g2_h2-i1_e?-f1_i?-j1",
            ),
            (
                "G71068GK",
                "QuiNAc(b1-2)[Oli3NAc(b1-3)][Oli3NAc(b1-3)Qui(b1-4)]Qui(?1-",
                "WURCS=2.0/4,5,4/[a2122m-1x_1-5][a2122m-1b_1-5_2*NCC/3=O][ad122m-1b_1-5_3*NCC/3=O][a2122m-1b_1-5]/1-2-3-4-3/a2-b1_a3-c1_a4-d1_d3-e1",
            ),
        ] {
            write!(&mut output, "<details open><summary>{code}</summary>").unwrap();
            let tokens = test_tokenise(wurcs).unwrap();

            let wurcs_tree = tokens.parse().unwrap();
            println!("{wurcs_tree}");
            let iupac_tree = GlycanStructure::from_short_iupac(iupac, 0..iupac.len(), 0).unwrap();
            println!("{iupac_tree}");

            let mut footnotes = Vec::new();
            let settings = GlycanRenderSettings {
                column_size: 30.0,
                sugar_size: 15.0,
                stroke_size: 1.5,
                foreground: [255, 255, 255],
                background: [0, 0, 0],
                modification_details: GlycanModificationDetails::Full,
            };
            let rendered_wurcs = wurcs_tree
                .render(
                    crate::glycan::GlycanRoot::Line,
                    crate::glycan::GlycanDirection::LeftToRight,
                    crate::glycan::GlycanSelection::FULL,
                    &settings,
                    &mut footnotes,
                )
                .unwrap();
            let rendered_iupac = iupac_tree
                .render(
                    crate::glycan::GlycanRoot::Line,
                    crate::glycan::GlycanDirection::LeftToRight,
                    crate::glycan::GlycanSelection::FULL,
                    &settings,
                    &mut footnotes,
                )
                .unwrap();

            write!(&mut output, "<p>WURCS</p>").unwrap();
            rendered_wurcs.to_svg(&mut output).unwrap();
            write!(&mut output, "<p>IUPAC</p>").unwrap();
            rendered_iupac.to_svg(&mut output).unwrap();
            if !footnotes.is_empty() {
                write!(
                    &mut output,
                    "<p>Footnotes</p><ol>{}</ol>",
                    footnotes.iter().map(|f| format!("<li>{f}</li>")).collect::<String>()
                )
                .unwrap();
            }

            let mut pairs = Vec::new();
            connect_trees(&wurcs_tree, &iupac_tree, &mut pairs);
            write!(
                &mut output,
                "<p>Pairs</p><ol>{}</ol>",
                pairs
                    .iter()
                    .map(|(a, b, c, d)| format!("<li>{a}: {b} vs {c}: {d}</li>"))
                    .collect::<String>()
            )
            .unwrap();
            write!(&mut output, "</details>").unwrap();
        }
        write!(&mut output, "</body></html>").unwrap();
    }

    fn connect_trees(
        a: &GlycanStructure,
        b: &GlycanStructure,
        pairs: &mut Vec<(String, String, String, String)>,
    ) {
        pairs.push((
            a.sugar.source_field.clone(),
            a.sugar.to_string(),
            b.sugar.source_field.clone(),
            b.sugar.to_string(),
        ));
        a.branches
            .iter()
            .zip(b.branches.iter())
            .for_each(|(a, b)| connect_trees(a, b, pairs));
    }
}
