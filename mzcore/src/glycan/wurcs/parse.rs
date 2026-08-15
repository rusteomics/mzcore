use crate::{
    chemistry::Element,
    glycan::{
        BaseSugar, GlycanSubstituent, HexoseIsomer, MonoSaccharide,
        wurcs::structs::{
            BackBone, Carbon, LIPOption, Lip, MAPSymbol, Residue, TerminalCarbon, Wurcs,
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
}

impl Wurcs {
    /// Parse a tokenised WURCS line into [`MonoSaccharide`]s, note that this is very much work in
    /// progress.
    /// # Errors
    /// If these tokens cannot be parsed into a [`MonoSaccharide`].
    pub fn parse(self) -> Result<Vec<MonoSaccharide>, WurcsParseError> {
        self.residues.into_iter().map(Residue::parse).collect()
    }
}

impl Residue {
    /// Parse a WURCS residue into a [`MonoSaccharide`].
    /// # Errors
    /// If this WURCS definition does not fit in the limitations of a [`MonoSaccharide`].
    pub(super) fn parse(self) -> Result<MonoSaccharide, WurcsParseError> {
        let (start, base, end) = match self.backbone {
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
            TerminalCarbon::CHHH => res.substituents.push((GlycanSubstituent::Deoxy, Some(1))),
            o => todo!("Not added yet: {:?}", o),
        }
        // Assumes a C=OH
        match end {
            TerminalCarbon::CXH => (),
            // TerminalCarbon::CHHH => res.substituents.push((GlycanSubstituent::Deoxy, None)), /*
            // TODO: get location */
            TerminalCarbon::CHHX => res.substituents.push((GlycanSubstituent::Alcohol, None)), /* TODO: get location */
            o => todo!("Not added yet: {:?}", o),
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
                        LIPOption::Known(Lip { position, .. }) => position,
                        _ => None,
                    })),
                    other => {
                        dbg!(other);
                        todo!()
                    }
                }
            } else if m.lips.len() == 2 {
                if m.modification.is_empty() {
                    internal_cycle_count += 1;
                } else {
                    todo!()
                }
            } else {
                todo!()
            }
        }

        if internal_cycle_count == 0 {
            // res.substituents.push((GlycanSubstituent::Alcohol, None));
        }
        if internal_cycle_count > 1 {
            todo!();
        }
        res.substituents.sort_unstable();
        Ok(res)
    }
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc, clippy::missing_errors_doc)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::{
        chemistry::Molecule,
        glycan::{BaseSugar, HexoseIsomer, MonoSaccharide, tokenise_wurcs, wurcs::structs::Wurcs},
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
            (
                "a2112h-1b_1-5",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "?",
                "Gal(b1- ",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "h2112h_2*NCC/3=O",
                BaseSugar::Hexose(Some(HexoseIsomer::Galactose)),
                "?",
                "Gal2NAc-ol",
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
            let monosaccharide = parsed[0].clone();
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
    fn structure() {
        // G58645HA,"Gal(b1-3)GalNAc-ol","WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/
        // a3-b1"

        let tokens =
            test_tokenise("WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1").unwrap();
        println!("{}", tokens.residues[0].to_structure().unwrap().to_dot());
        println!("{}", tokens.residues[1].to_structure().unwrap().to_dot());
        // todo!();
    }
}
