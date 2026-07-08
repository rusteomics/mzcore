use crate::{
    chemistry::{Connection, Element, StructuralFormula},
    system::i8::Charge,
};

#[derive(Debug)]
pub struct Wurcs {
    pub residues: Vec<Residue>,
    pub sequence: Vec<u8>,
    pub linkage: Vec<Linkage>,
}

#[derive(Debug)]
pub struct Residue {
    pub backbone: BackBone,
    pub anomeric: Option<(Option<u8>, AnomericSymbol)>,
    pub mods: Vec<Mod>,
}

impl Residue {
    pub fn to_structure(&self) -> Result<StructuralFormula, String> {
        if let BackBone::Defined(s, m, e) = &self.backbone {
            let mut structure = StructuralFormula::default();
            let mut carbon_numbers = Vec::new();
            let mut last = s.add_to_structure(&mut structure)?;
            carbon_numbers.push(last.0);

            for part in m {
                last = part.add_to_structure(&mut structure, last)?;
                carbon_numbers.push(last.0);
            }

            let end = e.add_to_structure(&mut structure)?;
            carbon_numbers.push(end.0);

            if last.1 == end.1 {
                structure.connections.push((last.0, end.0, last.1));
                for m in &self.mods {
                    if m.modification.is_empty() && m.lips.len() == 2 {
                        let LIPOption::Known(lip1) = m.lips[0] else {
                            return Err("Complex LIP detected".to_string());
                        };
                        let LIPOption::Known(lip2) = m.lips[1] else {
                            return Err("Complex LIP detected".to_string());
                        };

                        if let Some(lip1) = lip1.position
                            && let Some(lip2) = lip2.position
                        {
                            // Assume ether
                            let i = structure.elements.len();
                            structure.elements.push((Some(Element::O), None, Charge::default()));
                            structure.connections.push((
                                carbon_numbers[lip1 as usize],
                                i,
                                Connection::SingleCovalent,
                            ));
                            structure.connections.push((
                                carbon_numbers[lip2 as usize],
                                i,
                                Connection::SingleCovalent,
                            ));
                        } else {
                            return Err("Undefined carbon numbers detected in LIP".to_string());
                        }
                    } else if m.lips.len() == 1 {
                        let LIPOption::Known(lip) = m.lips[0] else {
                            return Err("Complex LIP detected".to_string());
                        };

                        if let Some(lip) = lip.position {
                            let mut connection = None;
                            let mut base = Connection::SingleCovalent;
                            let mut aromatic_start = None;
                            let mut lut = Vec::new();
                            for sym in &m.modification {
                                match sym {
                                    MAPSymbol::Element(e) => {
                                        if let Some((i, c)) = connection.take() {
                                            structure.connections.push((
                                                i,
                                                structure.elements.len(),
                                                c,
                                            ));
                                        } else {
                                            structure.connections.push((
                                                structure.elements.len() - 1,
                                                structure.elements.len(),
                                                base,
                                            ));
                                        }
                                        lut.push(structure.elements.len());
                                        structure.elements.push((
                                            Some(*e),
                                            None,
                                            Charge::default(),
                                        ));
                                    }
                                    MAPSymbol::Star(i) => {
                                        connection = Some((
                                            carbon_numbers[i.unwrap_or(lip) as usize],
                                            connection
                                                .map(|c| c.1)
                                                .unwrap_or(Connection::SingleCovalent),
                                        ));
                                        lut.push(carbon_numbers[i.unwrap_or(lip) as usize]);
                                    }
                                    MAPSymbol::DoubleBond => {
                                        connection = Some((
                                            connection
                                                .map(|c| c.0)
                                                .unwrap_or(structure.elements.len() - 1),
                                            Connection::DoubleCovalent,
                                        ));
                                    }
                                    MAPSymbol::TripleBond => {
                                        connection = Some((
                                            connection
                                                .map(|c| c.0)
                                                .unwrap_or(structure.elements.len() - 1),
                                            Connection::TripleCovalent,
                                        ));
                                    }
                                    MAPSymbol::AromaticStart => {
                                        base = Connection::DoubleCovalent;
                                        aromatic_start = Some(structure.elements.len() - 1);
                                    }
                                    MAPSymbol::AromaticEnd => {
                                        base = Connection::SingleCovalent;
                                        if let Some(start) = aromatic_start {
                                            structure.connections.push((
                                                start,
                                                structure.elements.len() - 1,
                                                Connection::DoubleCovalent,
                                            ));
                                        }
                                    }
                                    MAPSymbol::Branch(i) | MAPSymbol::Cyclic(i) => {
                                        connection = Some((
                                            lut[*i as usize - 1],
                                            connection
                                                .map(|c| c.1)
                                                .unwrap_or(Connection::SingleCovalent),
                                        ));
                                    }
                                    MAPSymbol::Chirality(_) => (), // Ignore for now
                                }
                            }
                        } else {
                            return Err("Undefined carbon numbers detected in LIP".to_string());
                        }
                    } else {
                        dbg!(m);
                    }
                }

                structure.infer([
                    (vec![(0, Connection::SingleCovalent)], StructuralFormula {
                        elements: vec![
                            (Some(Element::O), None, Charge::default()),
                            (Some(Element::H), None, Charge::default()),
                        ],
                        connections: vec![(0, 1, Connection::SingleCovalent)],
                    }),
                    (vec![(0, Connection::DoubleCovalent)], StructuralFormula {
                        elements: vec![(Some(Element::O), None, Charge::default())],
                        connections: vec![],
                    }),
                    (
                        vec![
                            (0, Connection::SingleCovalent),
                            (2, Connection::SingleCovalent),
                            (3, Connection::SingleCovalent),
                        ],
                        StructuralFormula {
                            elements: vec![
                                (Some(Element::O), None, Charge::default()),
                                (Some(Element::H), None, Charge::default()),
                                (Some(Element::H), None, Charge::default()),
                                (Some(Element::H), None, Charge::default()),
                            ],
                            connections: vec![(0, 1, Connection::SingleCovalent)],
                        },
                    ),
                ]);

                Ok(structure)
            } else {
                Err("Invalid connection".to_string())
            }
        } else {
            Err("repeating backbone".to_string())
        }
    }
}

#[derive(Debug)]
pub enum BackBone {
    Defined(TerminalCarbon, Vec<Carbon>, TerminalCarbon),
    Repeating(Option<TerminalCarbon>, Vec<Carbon>, Option<TerminalCarbon>),
}

/// The carbon descriptors extended with the options in: https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.6b00650/suppl_file/ci6b00650_si_001.pdf section 2.8.
#[derive(Clone, Copy, Debug)]
pub enum Carbon {
    /// 'd' deoxy `H-C-H`
    Methylene,
    /// 'C' `X-C-X`
    Dual,
    /// '1' `X-C-H`
    HydroxyLeft,
    /// '2' `H-C-X`
    HydroxyRight,
    /// '3' `C(X)(H)`
    HydroxyOpposite,
    /// '4' `C(H)(X)`
    HydroxySame,
    /// 'x' one of 1 or 2
    HydroxyUnknown,
    /// '5' `X-C-Y`
    DualLeft,
    /// '6' `Y-C-X`
    DualRight,
    /// '7' `C(X)(Y)`
    DualOpposite,
    /// '8' `C(Y)(X)`
    DualSame,
    /// 'X' one of 5 or 6
    DualUnknown,
    /// 'O' `C=O`
    Ketone,
    /// `e`
    DoubleHydroxyEntgegen,
    /// `z`
    DoubleHydroxyZusammen,
    /// `n`
    DoubleHydroxyNoIsomer,
    /// `f`
    DoubleHydroxyUnknown,
    /// `E`
    DoubleEntgegen,
    /// `Z`
    DoubleZusammen,
    /// `N`
    DoubleNoIsomer,
    /// `F`
    DoubleUnknown,
    /// `K`
    DoubleBonded,
    /// `T`
    TripleBonded,
    /// 'a' anomeric
    Hemiketal,
    /// 'U'
    KetoneOrHemiketal,
    /// 'Q' can be any of the other carbon descriptors
    Unknown,
}

impl Carbon {
    fn add_to_structure(
        &self,
        structure: &mut StructuralFormula,
        last: (usize, Connection),
    ) -> Result<(usize, Connection), String> {
        let index = structure.elements.len();
        match self {
            Carbon::Methylene => {
                if last.1 != Connection::SingleCovalent {
                    return Err(format!("Expected single connection found {:?}", last.1));
                }
                structure.elements.extend_from_slice(&[
                    (Some(Element::C), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                ]);
                structure.connections.extend_from_slice(&[
                    (last.0, index, Connection::SingleCovalent),
                    (index, index + 1, Connection::SingleCovalent),
                    (index, index + 2, Connection::SingleCovalent),
                ]);
                Ok((index, Connection::SingleCovalent))
            }
            Carbon::HydroxyLeft
            | Carbon::HydroxyRight
            | Carbon::HydroxySame
            | Carbon::HydroxyOpposite
            | Carbon::HydroxyUnknown => {
                if last.1 != Connection::SingleCovalent {
                    return Err(format!("Expected single connection found {:?}", last.1));
                }
                structure.elements.extend_from_slice(&[
                    (Some(Element::C), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                ]);
                structure.connections.extend_from_slice(&[
                    (last.0, index, Connection::SingleCovalent),
                    (index, index + 1, Connection::SingleCovalent),
                ]);
                Ok((index, Connection::SingleCovalent))
            }
            Carbon::Ketone => {
                if last.1 != Connection::SingleCovalent {
                    return Err(format!("Expected single connection found {:?}", last.1));
                }
                structure.elements.extend_from_slice(&[
                    (Some(Element::C), None, Charge::default()),
                    (Some(Element::O), None, Charge::default()),
                ]); // TODO: maybe this needs to be empty and fall back to O if not defined
                structure.connections.extend_from_slice(&[
                    (last.0, index, Connection::SingleCovalent),
                    (index, index + 1, Connection::DoubleCovalent),
                ]);
                Ok((index, Connection::SingleCovalent))
            }
            c => Err(format!("Carbon symbol {c:?} not supported yet")),
        }
    }
}

#[derive(Clone, Copy, Debug)]
#[allow(non_camel_case_types)]
pub enum TerminalCarbon {
    /// 'm'
    CHHH,
    /// 'M'
    CXXX,
    /// 'h'
    CHHX,
    /// 'c'
    CXXH,
    /// 'C'
    CXXY,
    /// '1'
    CXYH,
    /// '2'
    CYXH,
    /// '3'
    CXYH_opposite,
    /// '4'
    CYXH_same,
    /// 'x'
    CXYH_unknown,
    /// '5'
    CXYZ,
    /// '6'
    CYXZ,
    /// '7'
    CXYZ_opposite,
    /// '8'
    CYXZ_same,
    /// 'X'
    CXYZ_unknown,
    /// 'o' -C=XH
    CXH,
    /// 'A' -C=XY
    CXY,
    /// 'n' =CHH
    CHH,
    /// 'N' =CXX
    CXX,
    /// 'e'
    Double_CXH_entgegen,
    /// 'z'
    Double_CXH_zusammen,
    /// 'f'
    Double_CXH_unknown,
    /// 'E'
    Double_CXY_entgegen,
    /// 'Z'
    Double_CXY_zusammen,
    /// 'F'
    Double_CXY_unknown,
    /// 'T' ≡C-X or -C≡X
    Triple_CX,
    /// 'K'
    Double_CX,
    /// 't'
    Triple_CH,
    /// 'a'
    Hemiacetal,
    /// 'u'
    AldehydeOrHemiacetal,
    /// 'Q'
    Unknown,
}

impl TerminalCarbon {
    fn add_to_structure(
        &self,
        structure: &mut StructuralFormula,
    ) -> Result<(usize, Connection), String> {
        let index = structure.elements.len();
        match self {
            TerminalCarbon::CHHH => {
                structure.elements.extend_from_slice(&[
                    (Some(Element::C), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                ]);
                structure.connections.extend_from_slice(&[
                    (index, index + 1, Connection::SingleCovalent),
                    (index, index + 2, Connection::SingleCovalent),
                    (index, index + 3, Connection::SingleCovalent),
                ]);
                Ok((index, Connection::SingleCovalent))
            }
            TerminalCarbon::CHHX => {
                structure.elements.extend_from_slice(&[
                    (Some(Element::C), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()),
                    (Some(Element::H), None, Charge::default()), // Needs to fall back to O
                ]);
                structure.connections.extend_from_slice(&[
                    (index, index + 1, Connection::SingleCovalent),
                    (index, index + 2, Connection::SingleCovalent),
                ]);
                Ok((index, Connection::SingleCovalent))
            }
            TerminalCarbon::Hemiacetal => {
                structure.elements.push((Some(Element::C), None, Charge::default()));
                Ok((index, Connection::SingleCovalent))
            }
            c => Err(format!("Terminal carbon symbol {c:?} not supported yet")),
        }
    }
}

#[derive(Clone, Copy, Debug)]
pub enum AnomericSymbol {
    /// 'a'
    Alpha,
    /// 'b'
    Beta,
    /// 'u'
    Up,
    /// 'd'
    Down,
    /// 'x'
    Unknown,
    /// 'o'
    None,
}

#[derive(Debug)]
pub struct Mod {
    pub lips: Vec<LIPOption>,
    pub modification: Vec<MAPSymbol>,
}

#[derive(Clone, Debug)]
pub enum LIPOption {
    Known(LIP),
    Statistic(bool, Probability, LIP),
    Alternative(Vec<LIP>),
}

#[derive(Clone, Copy, Debug)]
pub struct LIP {
    pub position: Option<u8>,
    pub direction: Direction,
    pub star_index: u8,
}

#[derive(Clone, Copy, Debug)]
pub enum Direction {
    /// 'u'
    Upside,
    /// 'd'
    Downside,
    /// 't'
    Tres,
    /// 'a'
    Same,
    /// 'b'
    Opposite,
    /// 'c'
    Third,
    /// 'x'
    Unknown,
    /// 'e'
    Entgegen,
    /// 'z'
    Zusammen,
    /// 'f'
    UnknownGeometricalIsomerism,
    /// 'n'
    Obvious,
}

#[derive(Clone, Copy, Debug)]
pub enum MAPSymbol {
    Element(Element),
    /// '*' or '*n'
    Star(Option<u8>),
    /// '/n'
    Branch(u8),
    /// '$n'
    Cyclic(u8),
    /// '='
    DoubleBond,
    /// '#'
    TripleBond,
    Chirality(Chirality),
    AromaticStart,
    AromaticEnd,
}

#[derive(Clone, Copy, Debug)]
pub enum Chirality {
    R,
    S,
    E,
    Z,
    Unknown,
}

#[derive(Debug)]
pub enum Linkage {
    Known(LIN),
    Repeated(Repeat, LIN),
}

#[derive(Debug)]
pub struct LIN {
    pub lips: Vec<GLIPOption>,
    pub modification: Vec<MAPSymbol>,
}

#[derive(Clone, Debug)]
pub enum GLIPOption {
    Known(GLIP),
    Statistic(bool, Probability, GLIP),
    Alternative(Vec<GLIP>),
    RESAlternative(bool, Vec<GLIP>),
}

#[derive(Clone, Copy, Debug)]
pub struct GLIP {
    /// Encoded using base 52
    pub res_index: u8,
    /// Position or unknown
    pub position: Option<u8>,
    pub direction: Direction,
    pub star_index: u8,
}

#[derive(Clone, Copy, Debug)]
pub enum Probability {
    Single(Option<f32>),
    Range(Option<f32>, Option<f32>),
}

#[derive(Clone, Copy, Debug)]
pub enum Repeat {
    Single(Option<u8>),
    /// Max, min
    Range(Option<u8>, Option<u8>),
}
