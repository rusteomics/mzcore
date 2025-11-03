use std::fmt::Write;

use crate::{
    chemistry::{AmbiguousLabel, MassMode, MolecularFormula},
    glycan::GlycanPosition,
    parse_json::{ParseJson, use_serde},
    system::{Mass, OrderedMass, Ratio, da, fraction},
};

use itertools::Itertools;

impl From<&MolecularFormula> for OrderedMass {
    /// Create an ordered mass from the monoisotopic mass (needed for [`Multi<MolecularFormula>`](crate::quantities::Multi))
    fn from(value: &MolecularFormula) -> Self {
        value.monoisotopic_mass().into()
    }
}

impl MolecularFormula {
    /// The mass of the molecular formula of this element, if all element species (isotopes) exists
    #[expect(clippy::missing_panics_doc)]
    pub fn monoisotopic_mass(&self) -> Mass {
        let mut mass = da(*self.additional_mass);
        for (e, i, n) in &self.elements {
            mass += e
                .mass(*i)
                .expect("An invalid molecular formula was created, please report this crash")
                * Ratio::new::<fraction>(f64::from(*n));
        }
        mass
    }

    /// The average weight of the molecular formula of this element, if all element species (isotopes) exists
    #[expect(clippy::missing_panics_doc)]
    pub fn average_weight(&self) -> Mass {
        let mut mass = da(*self.additional_mass); // Technically this is wrong, the additional mass is defined to be monoisotopic
        for (e, i, n) in &self.elements {
            mass += e
                .average_weight(*i)
                .expect("An invalid molecular formula was created, please report this crash")
                * Ratio::new::<fraction>(f64::from(*n));
        }
        mass
    }

    /// The most abundant mass, meaning the isotope that will have the highest intensity.
    /// It uses an averagine model for the isotopes so the mass will not reflect any isotopomer
    /// exact mass but will be in the form of `monoisotopic exact mass + n * δ`. Where n is the
    /// integer Dalton offset for that isomer and δ the difference between 13 C and 12 C
    /// (1.0033548353399997 Da).
    ///
    /// Only available with crate feature 'isotopes'.
    #[cfg(feature = "isotopes")]
    pub fn most_abundant_mass(&self) -> Mass {
        let isotopes = self.isotopic_distribution(0.01);
        let max = isotopes
            .iter()
            .enumerate()
            .max_by_key(|s| ordered_float::OrderedFloat(*s.1));
        self.monoisotopic_mass() + da(max.map_or(0.0, |f| f.0 as f64 * 1.0033548353399997)) // Offset of 12C to 13C
    }

    /// Get the mass in the given mode
    pub fn mass(&self, mode: MassMode) -> Mass {
        match mode {
            MassMode::Monoisotopic => self.monoisotopic_mass(),
            MassMode::Average => self.average_weight(),
            #[cfg(feature = "isotopes")]
            MassMode::MostAbundant => self.most_abundant_mass(),
        }
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from
    /// this collection of elements merged with the ProForma notation for specific isotopes. This
    /// only displays the elements part and not any potential charge.
    pub fn hill_notation_core(&self) -> String {
        self.hill_notation_generic(
            |element, buffer| {
                if let Some(isotope) = element.1 {
                    write!(buffer, "[{}{}{}]", isotope, element.0, element.2,).unwrap();
                } else {
                    write!(buffer, "{}{}", element.0, element.2,).unwrap();
                }
            },
            false,
            false,
        )
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements merged with the ProForma notation for specific isotopes
    pub fn hill_notation(&self) -> String {
        self.hill_notation_generic(
            |element, buffer| {
                if let Some(isotope) = element.1 {
                    write!(buffer, "[{}{}{}]", isotope, element.0, element.2,).unwrap();
                } else {
                    write!(buffer, "{}{}", element.0, element.2,).unwrap();
                }
            },
            true,
            true,
        )
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of
    /// elements merged with the ProForma notation for specific isotopes. Using fancy unicode characters for subscript
    /// and superscript numbers.
    pub fn hill_notation_fancy(&self) -> String {
        self.hill_notation_generic(
            |element, buffer| {
                if let Some(isotope) = element.1 {
                    write!(buffer, "{}", to_superscript_num(isotope.get())).unwrap();
                }
                write!(buffer, "{}", element.0,).unwrap();
                if element.2 != 1 {
                    write!(buffer, "{}", to_subscript_num(element.2 as isize)).unwrap();
                }
            },
            true,
            true,
        )
    }

    /// Create a [Hill notation](https://en.wikipedia.org/wiki/Chemical_formula#Hill_system) from this collections of elements encoded in HTML
    pub fn hill_notation_html(&self) -> String {
        self.hill_notation_generic(
            |element, buffer| {
                if let Some(isotope) = element.1 {
                    write!(buffer, "<sup>{isotope}</sup>",).unwrap();
                }
                write!(buffer, "{}", element.0,).unwrap();
                if element.2 != 1 {
                    write!(buffer, "<sub>{}</sub>", element.2).unwrap();
                }
            },
            true,
            true,
        )
    }

    /// Check if this formula contains a negative number of any element (ignores a negative number of electrons).
    pub fn contains_negative_amount(&self) -> bool {
        self.elements()
            .iter()
            .any(|e| e.0 != crate::chemistry::Element::Electron && e.2 < 0)
    }
}

impl ParseJson for MolecularFormula {
    fn from_json_value(
        value: serde_json::Value,
    ) -> Result<Self, context_error::BoxedError<'static, context_error::BasicKind>> {
        use_serde(value)
    }
}

impl std::fmt::Display for AmbiguousLabel {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::AminoAcid {
                option,
                sequence_index,
                peptidoform_index,
                peptidoform_ion_index,
            } => write!(
                f,
                "{option}@p{}.{}i{}",
                peptidoform_ion_index + 1,
                peptidoform_index + 1,
                sequence_index + 1
            ),
            Self::Modification {
                id,
                sequence_index,
                peptidoform_index,
                peptidoform_ion_index,
            } => write!(
                f,
                "\x23{id}@p{}.{}i{}",
                peptidoform_ion_index + 1,
                peptidoform_index + 1,
                sequence_index
            ),
            Self::ChargeCarrier(formula) => write!(f, "[{}]", formula.hill_notation()),
            Self::CrossLinkBound(name) => write!(f, "intact{name}"),
            Self::CrossLinkBroken(name, formula) => {
                write!(f, "broken{name}@{}", formula.hill_notation())
            }
            Self::GlycanFragment(bonds) => {
                write!(f, "Y{}", bonds.iter().map(GlycanPosition::label).join("Y"))
            }
            Self::GlycanFragmentComposition(composition) => write!(
                f,
                "Y{}",
                composition
                    .iter()
                    .map(|(sugar, amount)| format!("{sugar}{amount}"))
                    .join("")
            ),
        }
    }
}

#[expect(clippy::missing_panics_doc)] // Cannot panic
fn to_subscript_num(input: isize) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        if *c == b'-' {
            output.push('\u{208B}');
        } else {
            output.push(char::from_u32(u32::from(*c) + 0x2080 - 0x30).unwrap());
        }
    }
    output
}

#[expect(clippy::missing_panics_doc)] // Cannot panic
fn to_superscript_num(input: u16) -> String {
    let text = input.to_string();
    let mut output = String::new();
    for c in text.as_bytes() {
        // b'-' could be '\u{207B}' but that is useless when using u16
        if *c == b'1' {
            output.push('\u{00B9}');
        } else if *c == b'2' {
            output.push('\u{00B2}');
        } else if *c == b'3' {
            output.push('\u{00B3}');
        } else {
            output.push(char::from_u32(u32::from(*c) + 0x2070 - 0x30).unwrap());
        }
    }
    output
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use crate::chemistry::MolecularFormula;

    #[test]
    fn sorted() {
        assert_eq!(molecular_formula!(H 2 O 2), molecular_formula!(O 2 H 2));
        assert_eq!(
            molecular_formula!(H 6 C 2 O 1),
            molecular_formula!(O 1 C 2 H 6)
        );
        assert_eq!(
            molecular_formula!(H 6 C 2 O 1),
            molecular_formula!(O 1 H 6 C 2)
        );
    }

    #[test]
    fn simplified() {
        assert_eq!(molecular_formula!(H 2 O 1 O 1), molecular_formula!(O 2 H 2));
        assert_eq!(
            molecular_formula!(H 2 O 1 O 1 H 1 H -2 H 0 H -1 H 2),
            molecular_formula!(O 2 H 2)
        );
        assert_eq!(
            molecular_formula!(H 2 Sb 0 O 1 O 1 H 1 H -2 H 0 H -1 N 0 P 0 Na 0 H 2),
            molecular_formula!(O 2 H 2)
        );
    }

    #[test]
    fn add() {
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O 1) + molecular_formula!(H 1 O 1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O 3) + molecular_formula!(H 1 O -1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O -1) + molecular_formula!(H 1 O 3)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 1 O -1) + molecular_formula!(O 3 H 1)
        );
        assert_eq!(
            molecular_formula!(H 2 O 2),
            molecular_formula!(H 2 O -1) + molecular_formula!(O 3)
        );
    }

    #[test]
    fn xlmod() {
        assert_eq!(
            MolecularFormula::from_xlmod("C7 D10 H2 N4", ..).unwrap(),
            molecular_formula!(C 7 [2 H 10] H 2 N 4)
        );
        assert_eq!(
            MolecularFormula::from_xlmod("-C1 -H2 O1", ..).unwrap(),
            molecular_formula!(C -1 H -2 O 1)
        );
        assert_eq!(
            MolecularFormula::from_xlmod("13C6 H6 O2", ..).unwrap(),
            molecular_formula!([13 C 6] H 6 O 2)
        );
    }

    #[test]
    fn pro_forma_spaces() {
        assert_eq!(
            MolecularFormula::pro_forma::<false, false>("C1[13C1]H6"),
            MolecularFormula::pro_forma::<false, false>("C 1 [ 13 C 1 ] H 6")
        );
    }

    #[test]
    fn pro_forma_empty() {
        assert_eq!(
            MolecularFormula::pro_forma::<false, true>("(empty)"),
            MolecularFormula::pro_forma::<false, true>("H0")
        );
    }

    #[test]
    fn unimod() {
        assert_eq!(
            MolecularFormula::from_unimod("C(1) 13C(1) H(6)", ..),
            Ok(molecular_formula!(C 1 [13 C 1] H 6))
        );
        assert_eq!(
            MolecularFormula::from_unimod("H(25) C(8) 13C(7) N 15N(2) O(3)", ..),
            Ok(molecular_formula!(H 25 C 8 [13 C 7] N 1 [15 N 2] O 3))
        );
        assert_eq!(
            MolecularFormula::from_unimod("H(6) C(4) N(2) dHex", ..),
            Ok(molecular_formula!(C 10 H 16 N 2 O 4))
        );
    }
}
