use std::{num::NonZeroU16, ops::Range};

use context_error::{BasicKind, BoxedError, Context, CreateError};

use crate::{
    chemistry::{Connection, ELEMENT_PARSE_LIST, Element, StructuralFormula},
    system::isize::Charge,
};

impl StructuralFormula {
    /// Parse a structural formula from an OpenSMILES v1.0 string.
    ///
    /// See <http://opensmiles.org/opensmiles.html> for the specification.
    ///
    /// # Errors
    /// If the string does not conform to the specification.
    pub fn from_smiles(value: &str) -> Result<Self, BoxedError<'_, BasicKind>> {
        let tokens = tokenise_smiles(value)?;
        let mut structure = Self::default();

        for (location, token) in tokens {
            match token {
                Token::Atom(_, isotope, element, _, hcount, charge, _) => {
                    let index = structure.elements.len();
                    structure.elements.push((element, isotope, charge));
                    for i in 0..hcount as usize {
                        structure.elements.push((Some(Element::H), None, Charge::default()));
                        structure.connections.push((index, index + i, Connection::SingleCovalent));
                    }
                }
                _ => todo!(),
            }
        }

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
                    // Note that the notation ++ and -- is a deprecated notation, so ignore until needed
                    let neg = value.as_bytes().get(index).copied() == Some(b'-');
                    index += 1;
                    if let Some(c @ (b'0'..=b'9')) = value.as_bytes().get(index) {
                        index += 1;
                        charge = (c - b'0') as isize;
                        if let Some(c @ (b'0'..=b'9')) = value.as_bytes().get(index) {
                            index += 1;
                            charge *= 10;
                            charge += (c - b'0') as isize;
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
            b'b' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::B),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'c' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::C),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'n' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::N),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'o' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::O),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b's' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::S),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'p' => tokens.push((
                index..index + 1,
                Token::Atom(
                    true,
                    None,
                    Some(Element::P),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'B' => {
                if value.as_bytes().get(index + 1).is_some_and(|b| *b == b'r') {
                    tokens.push((
                        index..index + 2,
                        Token::Atom(
                            false,
                            None,
                            Some(Element::Br),
                            None,
                            0,
                            Charge::default(),
                            None,
                        ),
                    ));
                    index += 1;
                } else {
                    tokens.push((
                        index..index + 1,
                        Token::Atom(
                            false,
                            None,
                            Some(Element::B),
                            None,
                            0,
                            Charge::default(),
                            None,
                        ),
                    ))
                }
            }
            b'C' => {
                if value.as_bytes().get(index + 1).is_some_and(|b| *b == b'l') {
                    tokens.push((
                        index..index + 2,
                        Token::Atom(
                            false,
                            None,
                            Some(Element::Cl),
                            None,
                            0,
                            Charge::default(),
                            None,
                        ),
                    ));
                    index += 1;
                } else {
                    tokens.push((
                        index..index + 1,
                        Token::Atom(
                            false,
                            None,
                            Some(Element::C),
                            None,
                            0,
                            Charge::default(),
                            None,
                        ),
                    ))
                }
            }
            b'N' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::N),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'O' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::O),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'S' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::S),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'P' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::P),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'F' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::F),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
            b'I' => tokens.push((
                index..index + 1,
                Token::Atom(
                    false,
                    None,
                    Some(Element::I),
                    None,
                    0,
                    Charge::default(),
                    None,
                ),
            )),
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
    use crate::chemistry::formula::smiles::tokenise_smiles;

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
}
