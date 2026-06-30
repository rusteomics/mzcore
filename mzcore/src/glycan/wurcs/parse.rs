use crate::glycan::{
    BaseSugar, HexoseIsomer, MonoSaccharide,
    wurcs::structs::{BackBone, Carbon, Residue, Wurcs},
};

#[derive(Debug)]
pub enum WurcsParseError {
    Empty,
    RepeatingBackbone,
    BackboneTooLong,
}

impl Wurcs {
    pub fn parse(self) -> Result<Vec<MonoSaccharide>, WurcsParseError> {
        self.residues.into_iter().map(|r| r.parse()).collect()
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
                    ] => Some(HexoseIsomer::Galactose),
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

        Ok(MonoSaccharide::new(base, &[]))
    }
}

#[cfg(test)]
mod tests {
    use context_error::{BasicKind, BoxedError, Context};

    use crate::glycan::{BaseSugar, HexoseIsomer, tokenise_wurcs, wurcs::structs::Wurcs};

    fn test_tokenise(value: &str) -> Result<Wurcs, BoxedError<'_, BasicKind>> {
        tokenise_wurcs(value, &Context::default().lines(0, value), ..)
    }

    #[test]
    fn composition() {
        // G58645HA,"Gal(b1-3)GalNAc-ol","WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/
        // a3-b1" [h2112h_2*NCC/3=O]

        let tokens =
            test_tokenise("WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1").unwrap();
        let parsed = dbg!(tokens.parse().unwrap());
        assert_eq!(parsed.len(), 2);
        assert_eq!(
            parsed[0].base_sugar,
            BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
        );
        assert_eq!(
            parsed[1].base_sugar,
            BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
        );
    }

    #[test]
    fn structure() {
        // G58645HA,"Gal(b1-3)GalNAc-ol","WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1"

        let tokens =
            test_tokenise("WURCS=2.0/2,2,1/[h2112h_2*NCC/3=O][a2112h-1b_1-5]/1-2/a3-b1").unwrap();
        println!(
            "{}",
            tokens.residues[0].to_structure().unwrap().to_dot().unwrap()
        );
        println!(
            "{}",
            tokens.residues[1].to_structure().unwrap().to_dot().unwrap()
        );
        todo!();
    }
}
