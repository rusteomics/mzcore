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
        Self::from_smiles_inner(&Context::default().lines(0, value), value, 0..value.len())
    }

    /// Parse a structural formula from an OpenSMILES v1.0 string.
    ///
    /// See <http://opensmiles.org/opensmiles.html> for the specification.
    ///
    /// # Errors
    /// If the string does not conform to the specification.
    pub fn from_smiles_inner<'a>(
        base_context: &Context<'a>,
        line: &'a str,
        range: Range<usize>,
    ) -> Result<Self, BoxedError<'a, BasicKind>> {
        let tokens = tokenise_smiles(base_context, line, range)?;
        let mut structure = Self::default();
        let mut branches = Vec::new();
        let mut rings: HashMap<u8, (usize, Option<Connection>, bool)> = HashMap::new();
        let mut last: Option<(usize, Option<Connection>, bool)> = None;
        let mut infer_positions = Vec::new();
        let mut locations = Vec::new();

        for (location, token) in tokens {
            match token {
                Token::Atom(
                    inline,
                    aromatic,
                    isotope,
                    element,
                    _chiral,
                    hcount,
                    charge,
                    _class,
                ) => {
                    let index = structure.atoms.len();
                    locations.push(location.clone());
                    if inline {
                        infer_positions.push(index);
                    }
                    structure.atoms.push((element, isotope, charge));
                    if let Some((last_index, bond, last_aromatic)) = last.take() {
                        structure.connections.push((
                            last_index,
                            index,
                            bond.unwrap_or({
                                if aromatic && last_aromatic {
                                    Connection::Aromatic
                                } else {
                                    Connection::SingleCovalent
                                }
                            }),
                        ));
                    }
                    last = Some((index, None, aromatic));
                    for i in 0..hcount as usize {
                        locations.push(location.clone()); // Duplicate the location as they are inferred based on this position
                        structure.atoms.push((Some(Element::H), None, Charge::default()));
                        structure.connections.push((
                            index,
                            index + i + 1,
                            Connection::SingleCovalent,
                        ));
                    }
                }
                Token::Bond(bond) => {
                    if let Some((_, last_bond, _)) = &mut last {
                        *last_bond = Some(match bond {
                            BondClass::Single | BondClass::Down | BondClass::Up => {
                                Connection::SingleCovalent
                            }
                            BondClass::Double => Connection::DoubleCovalent,
                            BondClass::Triple => Connection::TripleCovalent,
                            BondClass::Quadruple => Connection::QuadrupleCovalent,
                            BondClass::Aromatic => Connection::Aromatic,
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
                            base_context.clone().add_highlight((0, location)),
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
                            base_context.clone().add_highlight((0, location)),
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
                                (None, a) | (a, None) => a,
                                (a, b) => {
                                    if a == b {
                                        a
                                    } else {
                                        return Err(BoxedError::new(
                                            BasicKind::Error,
                                            "Invalid SMILES",
                                            "The bond symbol has to be the same on both locations of the bond number",
                                            base_context.clone().add_highlight((0, location)),
                                        ));
                                    }
                                }
                            }.unwrap_or(if current.2 && previous.2 {
                                    Connection::Aromatic
                                } else {
                                    Connection::SingleCovalent
                                })));
                        } else {
                            rings.insert(num, current);
                        }
                        last = Some((current.0, None, current.2)); // Scrub the connection
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "A bond number can only be given after an atom is defined",
                            base_context.clone().add_highlight((0, location)),
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
                base_context.clone().add_highlights(
                    branches
                        .iter()
                        .filter_map(|(index, ..)| locations.get(*index).map(|l| (0, l.clone()))),
                ),
            ));
        }

        if !rings.is_empty() {
            return Err(BoxedError::new(
                BasicKind::Error,
                "Invalid SMILES",
                "No all rings are closed",
                base_context.clone().add_highlights(rings.iter().filter_map(
                    |(key, (index, ..))| {
                        locations.get(*index).map(|l| (0, l.clone(), key.to_string()))
                    },
                )),
            ));
        }

        // Normalising helps make the validation easier
        structure.normalise_connections();

        // Check for self bonds
        for connection in &structure.connections {
            if connection.0 == connection.1 {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SMILES",
                    "There is a bond from one atom to itself",
                    base_context.clone().add_highlights(
                        locations.get(connection.0).into_iter().map(|l| (0, l.clone())),
                    ),
                ));
            }
        }

        // Check for multiple bonds between the same atoms
        for window in structure.connections.windows(2) {
            let one = window[0];
            let two = window[1];
            if one.0 == two.0 && one.1 == two.1 {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SMILES",
                    "There are two bonds binding the same pair of atoms",
                    base_context.clone().add_highlights(
                        locations
                            .get(one.0)
                            .into_iter()
                            .chain(locations.get(one.1))
                            .map(|l| (0, l.clone())),
                    ),
                ));
            }
        }

        // Check for wrong numbers of aromatic bonds
        for index in 0..structure.atoms.len() {
            let mut num_aromatic = 0_usize;
            // This is a very nonlinear approach to this problem (Natoms * Nconnections) it should
            // be possible to get the complexity down with some hashmaps, but I did not bother to do
            // that yet.
            for connection in &structure.connections {
                if (connection.0 == index || connection.1 == index)
                    && connection.2 == Connection::Aromatic
                {
                    num_aromatic += 1;
                }
            }
            if num_aromatic == 1 {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SMILES",
                    "This atom only has one aromatic bond",
                    base_context
                        .clone()
                        .add_highlights(locations.get(index).into_iter().map(|l| (0, l.clone()))),
                ));
            }
        }

        // Infer hydrogens based
        structure.infer_hydrogens(&infer_positions);

        Ok(structure)
    }
}

/// Tokenise a SMILES string for ease of parsing later.
///
/// # Errors
/// If the string is not a valid SMILES string.
fn tokenise_smiles<'a>(
    base_context: &Context<'a>,
    line: &'a str,
    range: Range<usize>,
) -> Result<Vec<(Range<usize>, Token)>, BoxedError<'a, BasicKind>> {
    if !line[range.clone()].is_ascii() {
        return Err(BoxedError::new(
            BasicKind::Error,
            "Invalid SMILES",
            "SMILES can only contain ASCII characters",
            base_context.clone().add_highlight((0, range)),
        ));
    }
    let mut index = 0;
    let mut tokens = Vec::new();
    let bytes = line[range].as_bytes();

    while let Some(c) = bytes.get(index) {
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
                let d1 = bytes.get(index + 1);
                let d2 = bytes.get(index + 2);
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
                            base_context.clone().add_highlight((0, index, 3)),
                        ));
                    }
                } else {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid SMILES",
                        "Too few characters after bond reference symbol '%', at least two characters are required",
                        base_context.clone().add_highlight((0, index, 3)),
                    ));
                }
            }
            b'[' => {
                let start = index;
                index += 1;

                // Isotope
                let mut isotope: u16 = 0;
                while let Some(c) = bytes.get(index)
                    && c.is_ascii_digit()
                {
                    if let Some(i) =
                        isotope.checked_mul(10).and_then(|i| i.checked_add(u16::from(c - b'0')))
                    {
                        isotope = i;
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "The isotope number is too high",
                            base_context.clone().add_highlight((0, start..=index)),
                        ));
                    }
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
                    if bytes[index..].starts_with(option.as_bytes()) {
                        element = Some(el);
                        index += option.len();
                        aromatic = a;
                        break;
                    }
                }
                if element.is_none() {
                    if bytes.get(index).copied() == Some(b'*') {
                        index += 1;
                    } else {
                        return Err(BoxedError::new(
                            BasicKind::Error,
                            "Invalid SMILES",
                            "Invalid element, use either a known element, a lowercase aromatic element, or '*'",
                            base_context.clone().add_highlight((0, index, 1)),
                        ));
                    }
                }

                // Chirality
                let chiral = if bytes.get(index).copied() == Some(b'@') {
                    index += 1;
                    match bytes.get(index) {
                        Some(b'@') => {
                            index += 1;
                            Some(ChiralClass::Clockwise)
                        }
                        Some(b'T' | b'A' | b'S' | b'O') => match bytes.get(index..=index + 1) {
                            Some(c @ (b"TH" | b"AL" | b"SP")) => {
                                if let Some(d) = bytes.get(index + 2) {
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
                                            base_context.clone().add_highlight((0, index + 2, 1)),
                                        ));
                                    }
                                } else {
                                    return Err(BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SMILES",
                                        "Too few characters after chiral symbol '@', at least four characters are required with this chiral class",
                                        base_context.clone().add_highlight((0, index - 1, 4)),
                                    ));
                                }
                            }
                            Some(c @ (b"TB" | b"OH")) => {
                                if let Some(d) = bytes.get(index + 2) {
                                    if d.is_ascii_digit() {
                                        let mut d = d - b'0';
                                        if let Some(d2 @ (b'0'..=b'9')) = bytes.get(index + 3) {
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
                                            base_context.clone().add_highlight((0, index + 2, 1)),
                                        ));
                                    }
                                } else {
                                    return Err(BoxedError::new(
                                        BasicKind::Error,
                                        "Invalid SMILES",
                                        "Too few characters after chiral symbol '@', at least four characters are required with this chiral class",
                                        base_context.clone().add_highlight((0, index - 1, 4)),
                                    ));
                                }
                            }
                            _ => {
                                return Err(BoxedError::new(
                                    BasicKind::Error,
                                    "Invalid SMILES",
                                    "Invalid chiral class, has to be one of TH, AL, SP, TB, or OH",
                                    base_context.clone().add_highlight((0, index - 1, 3)),
                                ));
                            }
                        },
                        _ => Some(ChiralClass::Anticlockwise),
                    }
                } else {
                    None
                };

                // Hcount
                let mut hcount = 0;
                if bytes.get(index).copied() == Some(b'H') {
                    index += 1;
                    hcount = 1;
                    if let Some(c @ (b'0'..=b'9')) = bytes.get(index) {
                        index += 1;
                        hcount = c - b'0';
                    }
                }

                // Charge
                let mut charge = 0;
                if bytes.get(index).copied() == Some(b'-')
                    || bytes.get(index).copied() == Some(b'+')
                {
                    // Note that the notation ++ and -- is a deprecated notation, so ignore until
                    // needed
                    let neg = bytes.get(index).copied() == Some(b'-');
                    charge = 1; // Handle cases where only the sign is given
                    index += 1;
                    if let Some(c @ (b'0'..=b'9')) = bytes.get(index) {
                        index += 1;
                        charge = (c - b'0') as i8;
                        if let Some(c @ (b'0'..=b'9')) = bytes.get(index) {
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
                let class = if bytes.get(index).copied() == Some(b':') {
                    index += 1;
                    let mut class_num: u16 = 0;
                    while let Some(c) = bytes.get(index)
                        && c.is_ascii_digit()
                    {
                        if let Some(i) = class_num
                            .checked_mul(10)
                            .and_then(|i| i.checked_add(u16::from(c - b'0')))
                        {
                            class_num = i;
                        } else {
                            return Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid SMILES",
                                "The class number is too high",
                                base_context.clone().add_highlight((0, start..=index)),
                            ));
                        }
                        index += 1;
                    }
                    Some(class_num)
                } else {
                    None
                };

                // Closing bracket
                if bytes.get(index).copied() != Some(b']') {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid SMILES",
                        "Invalid atom, the closing square bracket ']' is missing",
                        base_context.clone().add_highlight((0, index, 1)),
                    ));
                }

                // Detect invalid hydrogens
                if element.is_some_and(|e| e == Element::H) && hcount != 0 {
                    return Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid SMILES",
                        "A hydrogen atom cannot have a hydrogen count, use explicit square bracket notation",
                        base_context.clone().add_highlight((0, start..index + 1)),
                    ));
                }

                tokens.push((
                    start..index + 1,
                    Token::Atom(
                        false,
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
                if bytes.get(index + 1).is_some_and(|b| *b == b'r') {
                    tokens.push((index..index + 2, Token::simple_atom(false, Element::Br)));
                    index += 1;
                } else {
                    tokens.push((index..index + 1, Token::simple_atom(false, Element::B)));
                }
            }
            b'C' => {
                if bytes.get(index + 1).is_some_and(|b| *b == b'l') {
                    tokens.push((index..index + 2, Token::simple_atom(false, Element::Cl)));
                    index += 1;
                } else {
                    tokens.push((index..index + 1, Token::simple_atom(false, Element::C)));
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
                Token::Atom(false, false, None, None, None, 0, Charge::default(), None),
                // Technically this is inline, but no hydrogen inference should be run so it can be
                // set to false here.
            )),
            _ => {
                return Err(BoxedError::new(
                    BasicKind::Error,
                    "Invalid SMILES",
                    "Invalid character encountered",
                    base_context.clone().add_highlight((0, index, 1)),
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
    // Inline, Aromatic?, Isotope, Element, Chiral, Hses, charge, class
    Atom(
        bool,
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
            true,
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

/// The SMILES bond class
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

/// The chirality class
#[allow(dead_code)] // Yes the class numbers are currently unused
#[derive(Debug)]
enum ChiralClass {
    /// Anticlockwise written neighbours `@`
    Anticlockwise,
    /// Clockwise written neighbours `@@`
    Clockwise,
    Th(u8),
    Al(u8),
    Sp(u8),
    Tb(u8),
    Oh(u8),
}

#[cfg(test)]
mod tests {
    use context_error::Context;

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
            let tokens = tokenise_smiles(&Context::default().lines(0, o), o, 0..o.len()).unwrap();
            let mut last = 0;
            for (location, _token) in tokens {
                assert_eq!(last, location.start);
                last = location.end;
            }
            assert_eq!(last, o.len());
        }
    }

    #[test]
    fn equivalence() {
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
    }

    #[test]
    fn parse() {
        let tests = [
            ("CCC=O", molecular_formula!(C 3 H 6 O 1)),
            (
                "CC(=O)NCCCC[C@H](N-*)C(-*)=O",
                molecular_formula!(C 8 H 14 N 2 O 2),
            ),
            ("N1CC2CCCCC2CC1", molecular_formula!(C 9 H 17 N 1)),
            (
                "[NH4+].[NH4+].[O-]S(=O)(=O)[S-]",
                molecular_formula!(H 8 N 2 O 3 S 2),
            ),
            // ("C1:C:C:C:C:C1", molecular_formula!(C 6 H 6)), // The wrong bond is inferred here
            // because neither are marked as aromatic, so even though this might be an example
            // in the spec I will ignore this for now.
            ("c1ccccc1", molecular_formula!(C 6 H 6)),
            ("n1ccccc1", molecular_formula!(C 5 H 5 N 1)),
            ("o1cccc1", molecular_formula!(C 4 H 4 O 1)),
            ("n1c[nH]cc1", molecular_formula!(C 3 H 4 N 2)),
            ("c1ccccc1-c2ccccc2", molecular_formula!(C 12 H 10)),
            ("c1ccccc1C", molecular_formula!(C 7 H 8)),
            ("N#N", molecular_formula!(N 2)),
            ("CN=C=O", molecular_formula!(C 2 H 3 N 1 O 1)),
            ("[Cu+2].[O-]S(=O)(=O)[O-]", molecular_formula!(Cu 1 S 1 O 4)),
            ("O=Cc1ccc(O)c(OC)c1", molecular_formula!(C 8 H 8 O 3)), // Vanillin
            ("COc1cc(C=O)ccc1O", molecular_formula!(C 8 H 8 O 3)),
            (
                "CC(=O)NCCC1=CNc2c1cc(OC)cc2",
                molecular_formula!(C 13 H 16 N 2 O 2),
            ),
            (
                "CC(=O)NCCc1c[nH]c2ccc(OC)cc12",
                molecular_formula!(C 13 H 16 N 2 O 2),
            ),
            (
                "CCc(c1)ccc2[n+]1ccc3c2[nH]c4c3cccc4",
                molecular_formula!(C 17 H 15 N 2 :z+1),
            ),
            (
                "CCc1c[n+]2ccc3c4ccccc4[nH]c3c2cc1",
                molecular_formula!(C 17 H 15 N 2 :z+1),
            ),
            ("CN1CCC[C@H]1c2cccnc2", molecular_formula!(C 10 H 14 N 2)),
            (
                r"CCC[C@@H](O)CC\C=C\C=C\C#CC#C\C=C\CO",
                molecular_formula!(C 17 H 22 O 2),
            ),
            (
                "CCC[C@@H](O)CC/C=C/C=C/C#CC#C/C=C/CO",
                molecular_formula!(C 17 H 22 O 2),
            ),
            (
                r"CC1=C(C(=O)C[C@@H]1OC(=O)[C@@H]2[C@H](C2(C)C)/C=C(\C)/C(=O)OC)C/C=C\C=C",
                molecular_formula!(C 22 H 28 O 5),
            ),
            (
                "O1C=C[C@H]([C@H]1O2)c3c2cc(OC)c4c3OC(=O)C5=C4CCC(=O)5",
                molecular_formula!(C 17 H 12 O 6),
            ),
            (
                "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)1",
                molecular_formula!(C 6 H 12 O 6),
            ),
            (
                "OC[C@@H](O1)[C@@H](O)[C@H](O)[C@@H]2[C@@H]1c3c(O)c(OC)c(O)cc3C(=O)O2",
                molecular_formula!(C 14 H 16 O 9),
            ),
            (
                "CC(C)[C@@]12C[C@@H]1[C@@H](C)C(=O)C2",
                molecular_formula!(C 10 H 16 O 1),
            ),
            (
                "OCCc1c(C)[n+](cs1)Cc2cnc(C)nc2N",
                molecular_formula!(C 12 H 17 N 4 O 1 S 1 :z+1),
            ),
            (
                "CC(C)(O1)C[C@@H](O)[C@@]1(O2)[C@@H](C)[C@@H]3CC=C4[C@]3(C2)C(=O)C[C@H]5[C@H]4CC[C@@H](C6)[C@]5(C)Cc(n7)c6nc(C[C@@]89(C))c7C[C@@H]8CC[C@@H]%10[C@@H]9C[C@@H](O)[C@@]%11(C)C%10=C[C@H](O%12)[C@]%11(O)[C@H](C)[C@]%12(O%13)[C@H](O)C[C@@]%13(C)CO",
                molecular_formula!(C 54 H 74 N 2 O 10 ),
            ),
        ];
        for (test, expected) in tests {
            let structure = StructuralFormula::from_smiles(test).unwrap();
            let composition = structure.composition().unwrap();
            if composition != expected {
                println!("{}", structure.to_dot());
                panic!(
                    "SMILES={test} was calculated to have composition {composition} but should have {expected}, difference: {}",
                    &expected - &composition
                );
            }
        }
    }

    #[test]
    fn invalid() {
        let tests = [
            "C-1CCCCC=1",
            "C12CCCCC12",
            "C12C2CCC1",
            "[HH1]",
            "C.1CCCCC.1",
            "C11",
            "CcccC",
        ];
        for test in tests {
            if let Ok(structure) = StructuralFormula::from_smiles(test) {
                println!("{}", structure.to_dot());
                panic!("SMILES={test} did not fail");
            }
        }
    }
}
