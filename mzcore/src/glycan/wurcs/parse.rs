use crate::{
    chemistry::Element,
    glycan::{
        BaseSugar, GlycanSubstituent, HexoseIsomer, MonoSaccharide,
        wurcs::structs::{BackBone, Carbon, LIP, LIPOption, MAPSymbol, Residue, Wurcs},
    },
};

#[derive(Debug)]
pub enum WurcsParseError {
    Empty,
    RepeatingBackbone,
    BackboneTooLong,
}

impl Wurcs {
    pub fn parse(self) -> Result<Vec<MonoSaccharide>, WurcsParseError> {
        self.residues.into_iter().map(Residue::parse).collect()
    }
}

impl Residue {
    pub fn parse(self) -> Result<MonoSaccharide, WurcsParseError> {
        let base = match self.backbone {
            BackBone::Defined(_s, m, _e) => match m.len() {
                0 => BaseSugar::Sugar,
                1 => BaseSugar::Triose,
                2 => BaseSugar::Tetrose(None),
                3 => BaseSugar::Pentose(None),
                4 => BaseSugar::Hexose(match m.as_slice() {
                    [
                        Carbon::HydroxyRight,
                        Carbon::HydroxyLeft,
                        Carbon::HydroxyLeft,
                        Carbon::HydroxyRight,
                    ] => Some(HexoseIsomer::Galactose), // TODO: also needs the cycle ether?
                    _ => None,
                }),
                5 => BaseSugar::Heptose(None),
                6 => BaseSugar::Octose,
                7 => BaseSugar::Nonose(None),
                8 => BaseSugar::Decose,
                _ => return Err(WurcsParseError::BackboneTooLong),
            },
            BackBone::Repeating(..) => return Err(WurcsParseError::RepeatingBackbone),
        };
        let mut res = MonoSaccharide::new(base, &[]);
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
                        LIPOption::Known(LIP { position, .. }) => position,
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
            res.substituents.push((GlycanSubstituent::Alcohol, None));
        }
        if internal_cycle_count > 1 {
            todo!();
        }
        res.substituents.sort_unstable();
        Ok(res)
    }
}

#[cfg(test)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::glycan::{
        BaseSugar, GlycanSubstituent, HexoseIsomer, tokenise_wurcs, wurcs::structs::Wurcs,
    };

    fn test_tokenise(value: &str) -> Result<Wurcs, BoxedError<'_, BasicKind>> {
        tokenise_wurcs(value, &Context::default().lines(0, value), ..)
    }

    #[test]
    fn composition() {
        // G58645HA,"Gal(b1-3)GalNAc-ol","WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/
        // a3-b1"

        let tokens =
            test_tokenise("WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1").unwrap();
        let parsed = dbg!(tokens.parse().unwrap());
        assert_eq!(parsed.len(), 2);
        assert_eq!(
            parsed[0].base_sugar,
            BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
        );
        assert_eq!(parsed[0].substituents[0].0, GlycanSubstituent::Alcohol);
        assert_eq!(parsed[0].substituents[1].0, GlycanSubstituent::NAcetyl);
        assert_eq!(parsed[0].pro_forma_name(), "{C8H15N1O5}");
        assert_eq!(
            parsed[1].base_sugar,
            BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
        );
        assert_eq!(parsed[1].pro_forma_name(), "Hex");
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
