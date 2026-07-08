use std::{collections::HashMap, num::NonZeroU16, ops::Range};

use context_error::{BasicKind, BoxedError, Context, CreateError};

use crate::{
    chemistry::{Connection, ELEMENT_PARSE_LIST, Element, StructuralFormula},
    system::i8::Charge,
};

impl StructuralFormula {
    /// Parse a structural formula from an OpenSMILES v1.0 string.
    ///
    /// See <http://opensmiles.org/opensmiles.html> for the specification.
    ///
    /// # Errors
    /// If the string does not conform to the specification.
    pub fn from_smiles(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        // TODO: aromatic rings
        // TODO: cis/trans
        let tokens = tokenise_smiles(value)?;
        let mut structure = Self::default();
        let mut branches = Vec::new();
        let mut rings: HashMap<u8, (usize, Option<Connection>)> = HashMap::new();
        let mut last: Option<(usize, Option<Connection>)> = None;

        for (location, token) in tokens {
            match token {
                Token::Atom(_aromatic, isotope, element, _chiral, hcount, charge, _class) => {
                    let index = structure.elements.len();
                    structure.elements.push((element, isotope, charge));
                    if let Some((last_index, bond)) = last.take() {
                        structure.connections.push((last_index, index, bond.unwrap_or_default()));
                    }
                    last = Some((index, None));
                    for i in 0..hcount as usize {
                        structure.elements.push((Some(Element::H), None, Charge::default()));
                        structure.connections.push((
                            index,
                            index + i + 1,
                            Connection::SingleCovalent,
                        ));
                    }
                }
                Token::Bond(bond) => {
                    if let Some((last_index, last_bond)) = &mut last {
                        *last_bond = Some(match bond {
                            BondClass::Single => Connection::SingleCovalent,
                            BondClass::Double => Connection::DoubleCovalent,
                            BondClass::Triple => Connection::TripleCovalent,
                            BondClass::Quadruple => Connection::QuadrupleCovalent,
                            _ => todo!(),
                        });
                    }
                }
                Token::BranchOpen => {
                    if let Some(last) = last {
                        branches.push(last);
                        // Keep last in place for the next atom
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "Cannot open a branch if no atom was defined before",
                            Context::default().lines(0, value).add_highlight((0, location)),
                        ));
                    }
                }
                Token::BranchClose => {
                    if let Some(previous) = branches.pop() {
                        last = Some(previous);
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "Too many branches closed, this branch is missing the branch open bracket",
                            Context::default().lines(0, value).add_highlight((0, location)),
                        ));
                    }
                }
                Token::Dot => {
                    last = None;
                }
                Token::Reference(num) => {
                    if let Some(current) = last {
                        if let Some(previous) = rings.remove(&num) {
                            structure.connections.push((previous.0, current.0, match (previous.1, current.1) {
                                (None, a) | (a, None) => a.unwrap_or_default(),
                                (a, b) => {
                                    if a == b {
                                        a.unwrap_or_default()
                                    } else {
                                        return Err(BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid SMILES",
                                            "The bond symbol has to be the same on both locations of the bond number",
                                            Context::default().lines(0, value).add_highlight((0, location)),
                                        ));
                                    }
                                }
                            }));
                        } else {
                            rings.insert(num, current);
                        }
                        last = Some((current.0, None)); // Scrub the connection
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "A bond number can only be given after an atom is defined",
                            Context::default().lines(0, value).add_highlight((0, location)),
                        ));
                    }
                }
            }
        }

        if !branches.is_empty() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid SMILES",
                "No all branches are closed",
                Context::default().lines(0, value),
            ));
        }
        if !rings.is_empty() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid SMILES",
                "No all rings are closed",
                Context::default().lines(0, value),
            ));
        }
        // TODO: Validate that only one connection is made between each pair of atoms
        // TODO: Validate that no atom is bound to itself
        // TODO: Technically only the organic subset has to be inferred (not the square bracket
        // atoms)

        structure.infer_hydrogens();
        Ok(structure)
    }
}

fn tokenise_smiles(value: &str) -> Result<Vec<(Range<usize>, Token)>, BoxedError<'_, BasicKind>> {
    if !value.is_ascii() {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid SMILES",
            "SMILES can only contain ASCII characters",
            Context::default().lines(0, value),
        ));
    }
    let mut index = 0;
    let mut tokens = Vec::new();

    while let Some(c) = value.as_bytes().get(index) {
        match c {
            b'-' => tokens.push((index..index + 1, Token::Bond(BondClass::Single))),
            b'=' => tokens.push((index..index + 1, Token::Bond(BondClass::Double))),
            b'#' => tokens.push((index..index + 1, Token::Bond(BondClass::Triple))),
            b'$' => tokens.push((index..index + 1, Token::Bond(BondClass::Quadruple))),
            b'/' => tokens.push((index..index + 1, Token::Bond(BondClass::Up))),
            b'\\' => tokens.push((index..index + 1, Token::Bond(BondClass::Down))),
            b':' => tokens.push((index..index + 1, Token::Bond(BondClass::Aromatic))),
            b'(' => tokens.push((index..index + 1, Token::BranchOpen)),
            b')' => tokens.push((index..index + 1, Token::BranchClose)),
            b'.' => tokens.push((index..index + 1, Token::Dot)),
            b'%' => {
                let d1 = value.as_bytes().get(index + 1);
                let d2 = value.as_bytes().get(index + 2);
                if let Some(d1) = d1
                    && let Some(d2) = d2
                {
                    if d1.is_ascii_digit() && d2.is_ascii_digit() {
                        tokens.push((
                            index..index + 3,
                            Token::Reference((d1 - b'0') * 10 + (d2 - b'0')),
                        ));
                        index += 2;
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "Invalid digit after bond reference symbol '%'",
                            Context::default().lines(0, value).add_highlight((0, index, 3)),
                        ));
                    }
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid SMILES",
                        "Too few characters after bond reference symbol '%', at least two characters are required",
                        Context::default().lines(0, value).add_highlight((0, index, 3)),
                    ));
                }
            }
            b'[' => {
                let start = index;
                index += 1;
                // Isotope
                let mut isotope = 0;
                while let Some(c) = value.as_bytes().get(index)
                    && c.is_ascii_digit()
                {
                    isotope *= 10;
                    isotope += (c - b'0') as u16; // TODO: detect overflow?
                    index += 1;
                }

                // Element
                let mut aromatic = false;
                let mut element = None;
                for (a, option, el) in
                    ELEMENT_PARSE_LIST.iter().map(|(o, e)| (false, *o, *e)).chain([
                        (true, "b", Element::B),
                        (true, "c", Element::C),
                        (true, "n", Element::N),
                        (true, "o", Element::O),
                        (true, "p", Element::P),
                        (true, "se", Element::Se),
                        (true, "s", Element::S),
                        (true, "as", Element::As),
                    ])
                {
                    if value.as_bytes()[index..].starts_with(option.as_bytes()) {
                        element = Some(el);
                        index += option.len();
                        aromatic = a;
                        break;
                    }
                }
                if element.is_none() {
                    if value.as_bytes().get(index).copied() == Some(b'*') {
                        index += 1;
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "Invalid element, use either a known element, a lowercase aromatic element, or '*'",
                            Context::default().lines(0, value).add_highlight((0, index, 1)),
                        ));
                    }
                }

                // Chirality
                let mut chiral = None;
                if value.as_bytes().get(index).copied() == Some(b'@') {
                    index += 1;
                    chiral = match value.as_bytes().get(index) {
                        Some(b'@') => {
                            index += 1;
                            Some(ChiralClass::At)
                        }
                        Some(b'T' | b'A' | b'S' | b'O') => {
                            match value.as_bytes().get(index..=index + 1) {
                                Some(c @ (b"TH" | b"AL" | b"SP")) => {
                                    if let Some(d) = value.as_bytes().get(index + 2) {
                                        if d.is_ascii_digit() {
                                            let d = d - b'0';
                                            index += 3;
                                            Some(match c {
                                                b"TH" => ChiralClass::Th(d),
                                                b"AL" => ChiralClass::Al(d),
                                                b"SP" => ChiralClass::Sp(d),
                                                _ => unreachable!(),
                                            })
                                        } else {
                                            return Err(BoxedError::new(
                                                BasicKind::Error,
                                                "Invalid SMILES",
                                                "Invalid digit after chiral class",
                                                Context::default().lines(0, value).add_highlight((
                                                    0,
                                                    index + 2,
                                                    1,
                                                )),
                                            ));
                                        }
                                    } else {
                                        return Err(BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid SMILES",
                                            "Too few characters after chiral symbol '@', at least four characters are required with this chiral class",
                                            Context::default().lines(0, value).add_highlight((
                                                0,
                                                index - 1,
                                                4,
                                            )),
                                        ));
                                    }
                                }
                                Some(c @ (b"TB" | b"OH")) => {
                                    if let Some(d) = value.as_bytes().get(index + 2) {
                                        if d.is_ascii_digit() {
                                            let mut d = d - b'0';
                                            if let Some(d2 @ (b'0'..=b'9')) =
                                                value.as_bytes().get(index + 3)
                                            {
                                                d += (d2 - b'0') * 10;
                                                index += 1;
                                            }
                                            index += 3;
                                            Some(match c {
                                                b"TB" => ChiralClass::Tb(d),
                                                b"OH" => ChiralClass::Oh(d),
                                                _ => unreachable!(),
                                            })
                                        } else {
                                            return Err(BoxedError::new(
                                                BasicKind::Error,
                                                "Invalid SMILES",
                                                "Invalid digit after chiral class",
                                                Context::default().lines(0, value).add_highlight((
                                                    0,
                                                    index + 2,
                                                    1,
                                                )),
                                            ));
                                        }
                                    } else {
                                        return Err(BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid SMILES",
                                            "Too few characters after chiral symbol '@', at least four characters are required with this chiral class",
                                            Context::default().lines(0, value).add_highlight((
                                                0,
                                                index - 1,
                                                4,
                                            )),
                                        ));
                                    }
                                }
                                _ => {
                                    return Err(BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SMILES",
                                        "Invalid chiral class, has to be one of TH, AL, SP, TB, or OH",
                                        Context::default().lines(0, value).add_highlight((
                                            0,
                                            index - 1,
                                            3,
                                        )),
                                    ));
                                }
                            }
                        }
                        _ => Some(ChiralClass::None),
                    };
                }
                // Hcount
                let mut hcount = 0;
                if value.as_bytes().get(index).copied() == Some(b'H') {
                    index += 1;
                    hcount = 1;
                    if let Some(c @ (b'0'..=b'9')) = value.as_bytes().get(index) {
                        index += 1;
                        hcount = c - b'0';
                    }
                }
                // Charge
                let mut charge = 0;
                if value.as_bytes().get(index).copied() == Some(b'-')
                    || value.as_bytes().get(index).copied() == Some(b'+')
                {
                    // Note that the notation ++ and -- is a deprecated notation, so ignore until
                    // needed
                    let neg = value.as_bytes().get(index).copied() == Some(b'-');
                    charge = 1; // Handle cases where only the sign is given
                    index += 1;
                    if let Some(c @ (b'0'..=b'9')) = value.as_bytes().get(index) {
                        index += 1;
                        charge = (c - b'0') as i8;
                        if let Some(c @ (b'0'..=b'9')) = value.as_bytes().get(index) {
                            index += 1;
                            charge *= 10;
                            charge += (c - b'0') as i8;
                        }
                    }
                    if neg {
                        charge *= -1;
                    }
                }
                let charge = Charge::new::<crate::system::charge::e>(charge);
                // Class
                let mut class = None;
                if value.as_bytes().get(index).copied() == Some(b':') {
                    index += 1;
                    let mut class_num = 0;
                    while let Some(c) = value.as_bytes().get(index)
                        && c.is_ascii_digit()
                    {
                        class_num *= 10;
                        class_num += (c - b'0') as u16; // TODO: detect overflow?
                        index += 1;
                    }
                    class = Some(class_num);
                }

                // Closing bracket
                if value.as_bytes().get(index).copied() != Some(b']') {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid SMILES",
                        "Invalid atom, the closing square bracket ']' is missing",
                        Context::default().lines(0, value).add_highlight((0, index, 1)),
                    ));
                }

                tokens.push((
                    start..index + 1,
                    Token::Atom(
                        aromatic,
                        NonZeroU16::new(isotope),
                        element,
                        chiral,
                        hcount,
                        charge,
                        class,
                    ),
                ));
            }
            d @ b'0'..=b'9' => {
                tokens.push((index..index + 1, Token::Reference(d - b'0')));
            }
            b'b' => tokens.push((index..index + 1, Token::simple_atom(true, Element::B))),
            b'c' => tokens.push((index..index + 1, Token::simple_atom(true, Element::C))),
            b'n' => tokens.push((index..index + 1, Token::simple_atom(true, Element::N))),
            b'o' => tokens.push((index..index + 1, Token::simple_atom(true, Element::O))),
            b's' => tokens.push((index..index + 1, Token::simple_atom(true, Element::S))),
            b'p' => tokens.push((index..index + 1, Token::simple_atom(true, Element::P))),
            b'B' => {
                if value.as_bytes().get(index + 1).is_some_and(|b| *b == b'r') {
                    tokens.push((index..index + 2, Token::simple_atom(false, Element::Br)));
                    index += 1;
                } else {
                    tokens.push((index..index + 1, Token::simple_atom(false, Element::B)))
                }
            }
            b'C' => {
                if value.as_bytes().get(index + 1).is_some_and(|b| *b == b'l') {
                    tokens.push((index..index + 2, Token::simple_atom(false, Element::Cl)));
                    index += 1;
                } else {
                    tokens.push((index..index + 1, Token::simple_atom(false, Element::C)))
                }
            }
            b'N' => tokens.push((index..index + 1, Token::simple_atom(false, Element::N))),
            b'O' => tokens.push((index..index + 1, Token::simple_atom(false, Element::O))),
            b'S' => tokens.push((index..index + 1, Token::simple_atom(false, Element::S))),
            b'P' => tokens.push((index..index + 1, Token::simple_atom(false, Element::P))),
            b'F' => tokens.push((index..index + 1, Token::simple_atom(false, Element::F))),
            b'I' => tokens.push((index..index + 1, Token::simple_atom(false, Element::I))),
            b'*' => tokens.push((
                index..index + 1,
                Token::Atom(false, None, None, None, 0, Charge::default(), None),
            )),
            _ => {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SMILES",
                    "Invalid character encountered",
                    Context::default().lines(0, value).add_highlight((0, index, 1)),
                ));
            }
        }
        index += 1;
    }

    Ok(tokens)
}

#[derive(Debug)]
enum Token {
    Bond(BondClass),
    // Aromatic?, Isotope, Element, Chiral, Hses, charge, class
    Atom(
        bool,
        Option<NonZeroU16>,
        Option<Element>,
        Option<ChiralClass>,
        u8,
        Charge,
        Option<u16>,
    ),
    BranchOpen,
    BranchClose,
    Reference(u8),
    Dot,
}

impl Token {
    fn simple_atom(aromatic: bool, element: Element) -> Self {
        Self::Atom(
            aromatic,
            None,
            Some(element),
            None,
            0,
            Charge::default(),
            None,
        )
    }
}

#[derive(Debug)]
enum BondClass {
    Single,
    Double,
    Triple,
    Quadruple,
    /// `/`
    Up,
    /// `\`
    Down,
    /// `:`
    Aromatic,
}

#[derive(Debug)]
enum ChiralClass {
    None,
    At,
    Th(u8),
    Al(u8),
    Sp(u8),
    Tb(u8),
    Oh(u8),
}

#[cfg(test)]
mod tests {
    use crate::chemistry::{StructuralFormula, formula::smiles::tokenise_smiles};

    #[test]
    fn tokenise() {
        for o in [
            "Oc1c(*)cccc1",
            "C=C",
            "C#N",
            "CC#CC",
            "CCC=O",
            "[CH4:2]",
            "[Rh-](Cl)(Cl)(Cl)(Cl)$[Rh-](Cl)(Cl)(Cl)Cl",
            "N1CC2CCCCC2CC1",
            "[H]C([H])([H])[H]",
            "[238U]",
            "[12CH3+:234]",
            "Oc1cc(.NCCO)ccc1",
            "c1c2c3c4cc1.Br2.Cl3.Cl4",
            "N[C@](Br)(O)C",
            "N[C@@](Br)(O)C",
            "F/C=C/F",
            "NC(Br)=[C@]=C(O)C",
            "F[As@TB15](Cl)(S)(Br)N",
        ] {
            let tokens = tokenise_smiles(o).unwrap();
            let mut last = 0;
            for (location, _token) in tokens {
                assert_eq!(last, location.start);
                last = location.end;
            }
            assert_eq!(last, o.len());
        }
    }

    #[test]
    fn parse() {
        let structure = StructuralFormula::from_smiles("CCC=O").unwrap();
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(C 3 H 6 O 1)
        );
        // MOD:00064|N6-acetyl-L-lysine
        let structure = StructuralFormula::from_smiles("CC(=O)NCCCC[C@H](N-*)C(-*)=O").unwrap();
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(C 8 H 14 N 2 O 2)
        );
        let structure = StructuralFormula::from_smiles("N1CC2CCCCC2CC1").unwrap();
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(C 9 H 17 N 1)
        );
        let structure = StructuralFormula::from_smiles("C=1CCCCC=1").unwrap();
        let structure_a = StructuralFormula::from_smiles("C=1CCCCC=1").unwrap();
        let structure_b = StructuralFormula::from_smiles("C=1CCCCC=1").unwrap();
        assert_eq!(structure, structure_a);
        assert_eq!(structure, structure_b);
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(C 6 H 10)
        );
        let structure = StructuralFormula::from_smiles("C1CCCCC1C1CCCCC1").unwrap();
        let structure_a = StructuralFormula::from_smiles("C1CCCCC1C2CCCCC2").unwrap();
        let structure_b = StructuralFormula::from_smiles("C1CCCCC1C%42CCCCC%42").unwrap();
        assert_eq!(structure, structure_a);
        assert_eq!(structure, structure_b);
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(C 12 H 22)
        );
        let structure = StructuralFormula::from_smiles("[NH4+].[NH4+].[O-]S(=O)(=O)[S-]").unwrap();
        // println!("{}", structure.to_dot().unwrap());
        assert_eq!(
            structure.composition().unwrap(),
            molecular_formula!(H 8 N 2 O 3 S 2)
        );
    }
}
