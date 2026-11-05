//! Handle monosaccharides

use std::{borrow::Cow, fmt::Display, hash::Hash};

use context_error::{BasicKind, BoxedError, Context, CreateError};
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use thin_vec::ThinVec;

use crate::{
    chemistry::{MassOutputMode, MolecularFormula, Molecule},
    glycan::{
        BaseSugar, Configuration, GlycanSubstituent, GlycanSubstituentLeavingGroup,
        GlycanSubstituentLocation, HexoseIsomer,
    },
    parse_json::ParseJson,
    sequence::SequencePosition,
    space::{Space, UsedSpace},
};

/// A monosaccharide with all its complexity
#[derive(Clone, Debug, Deserialize, Ord, PartialOrd, Serialize)]
pub struct MonoSaccharide {
    pub(super) base_sugar: BaseSugar,
    pub(super) substituents: ThinVec<(GlycanSubstituent, GlycanSubstituentLocation)>,
    pub(super) furanose: bool,
    pub(super) configuration: Option<Configuration>,
    #[cfg(debug_assertions)]
    pub(super) source_field: String,
}

impl Molecule for MonoSaccharide {
    fn calculate_mass_inner<Mode: MassOutputMode>(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Mode::Output {
        self.base_sugar
            .calculate_mass_inner::<Mode>(sequence_index, peptidoform_index)
            + self
                .substituents
                .iter()
                .map(|s| {
                    s.0.calculate_mass_inner::<Mode>(sequence_index, peptidoform_index)
                        + s.1.calculate_mass_inner::<Mode>(sequence_index, peptidoform_index)
                })
                .sum::<Mode::Output>()
    }
}

impl Space for MonoSaccharide {
    fn space(&self) -> UsedSpace {
        (self.base_sugar.space()
            + self.substituents.space()
            + self.furanose.space()
            + self.configuration.space())
        .set_total::<Self>()
    }
}

impl MonoSaccharide {
    /// Check if this monosacharide is similar to another, this checks if the
    /// [`BaseSugar::equivalent`] is true (passing the precise flag there as well) and if the
    /// substituents are identical. If the precise flag is on it also checks if the furanose and
    /// configuration state are the same as well.
    pub fn equivalent(&self, other: &Self, precise: bool) -> bool {
        self.base_sugar.equivalent(&other.base_sugar, precise)
            && self.substituents == other.substituents
            && (!precise
                || (self.furanose == other.furanose && self.configuration == other.configuration))
    }
}

impl PartialEq for MonoSaccharide {
    fn eq(&self, other: &Self) -> bool {
        self.equivalent(other, true)
    }
}

impl Eq for MonoSaccharide {}

impl Hash for MonoSaccharide {
    fn hash<H: std::hash::Hasher>(&self, hasher: &mut H) {
        self.base_sugar.hash(hasher);
        self.substituents.hash(hasher);
        self.furanose.hash(hasher);
        self.configuration.hash(hasher);
    }
}

impl Display for MonoSaccharide {
    /// Display as valid ProForma glycan
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pro_forma_name())
    }
}

impl ParseJson for MonoSaccharide {
    fn from_json_value(value: Value) -> Result<Self, BoxedError<'static, BasicKind>> {
        if let Value::Object(mut map) = value {
            let context = |map: &serde_json::Map<String, Value>| {
                Context::default().lines(
                    0,
                    map.iter().map(|(k, v)| format!("\"{k}\": {v}")).join(","),
                )
            };
            Ok(Self {
                base_sugar: BaseSugar::from_json_value(map.remove("base_sugar").ok_or_else(
                    || {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid MonoSaccharide",
                            "The required property 'base_sugar' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                substituents: if let Value::Array(arr) =
                    map.remove("substituents").ok_or_else(|| {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid MonoSaccharide",
                            "The required property 'substituents' is missing",
                            context(&map),
                        )
                    })? {
                    arr.into_iter()
                        .map(|el| match el {
                            Value::Array(ref arr) if let [ref sub, ref num] = arr[..] => {
                                GlycanSubstituent::from_json_value(sub.clone()).and_then(|sub| {
                                    GlycanSubstituentLocation::from_json_value(num.clone()).map(|num| (sub, num))
                                })
                            }
                            sub @ Value::String(_) => {
                                GlycanSubstituent::from_json_value(sub).map(|sub| (sub, GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)))
                            }
                            _ => Err(BoxedError::new(
                                BasicKind::Error,
                                "Invalid MonoSaccharide",
                                "The substituents need to be an array with the glycansubsistuent and the location",
                                context(&map),
                            )),
                        })
                        .collect()
                } else {
                    Err(BoxedError::new(
                        BasicKind::Error,
                        "Invalid MonoSaccharide",
                        "The substituents need to be an array",
                        context(&map),
                    ))
                }?,
                configuration: Option::from_json_value(map.remove("configuration").ok_or_else(
                    || {
                        BoxedError::new(
                            BasicKind::Error,
                            "Invalid MonoSaccharide",
                            "The required property 'configuration' is missing",
                            context(&map),
                        )
                    },
                )?)?,
                furanose: bool::from_json_value(map.remove("furanose").ok_or_else(|| {
                    BoxedError::new(
                        BasicKind::Error,
                        "Invalid MonoSaccharide",
                        "The required property 'furanose' is missing",
                        context(&map),
                    )
                })?)?,
                #[cfg(debug_assertions)]
                source_field: String::new(),
            })
        } else {
            Err(BoxedError::new(
                BasicKind::Error,
                "Invalid MonoSaccharide",
                "The value has to be a map",
                Context::default().lines(0, value.to_string()),
            ))
        }

        // use_serde(value)
    }
}

impl MonoSaccharide {
    /// Get the name that should be used to represent this monosaccharide in ProForma glycan
    /// compositions
    pub fn pro_forma_name(&self) -> Cow<'static, str> {
        match self.base_sugar {
            BaseSugar::Sugar
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Sug")
            }
            BaseSugar::Triose
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Tri")
            }
            BaseSugar::Tetrose(_)
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Tet")
            }
            BaseSugar::Pentose(_)
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Pen")
            }
            BaseSugar::Hexose(isomer) if !self.furanose && self.configuration.is_none() => {
                match *self.substituents.as_slice() {
                    [] => Cow::Borrowed("Hex"),
                    [
                        (
                            GlycanSubstituent::Acid,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("aHex"),
                    [
                        (
                            GlycanSubstituent::Acid,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Deoxy,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Didehydro,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("en,aHex"),
                    [
                        (
                            GlycanSubstituent::Deoxy,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] if isomer == Some(HexoseIsomer::Galactose) => Cow::Borrowed("Fuc"),
                    [
                        (
                            GlycanSubstituent::Deoxy,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("dHex"),
                    [
                        (
                            GlycanSubstituent::NAcetyl,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Sulfate,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexNAcS"),
                    [
                        (
                            GlycanSubstituent::Amino,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Sulfate,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexNS"),
                    [
                        (
                            GlycanSubstituent::Amino,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexN"),
                    [
                        (
                            GlycanSubstituent::NAcetyl,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexNAc"),
                    [
                        (
                            GlycanSubstituent::Sulfate,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexS"),
                    [
                        (
                            GlycanSubstituent::Phosphate,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("HexP"),
                    _ => Cow::Owned(format!("{{{}}}", self.formula())),
                }
            }
            BaseSugar::Heptose(_)
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Hep")
            }
            BaseSugar::Octose
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Oct")
            }
            BaseSugar::Nonose(None) if !self.furanose && self.configuration.is_none() => {
                match *self.substituents.as_slice() {
                    [] => Cow::Borrowed("Non"),
                    [
                        (
                            GlycanSubstituent::Acetyl,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Acid,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Amino,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("NeuAc"),
                    [
                        (
                            GlycanSubstituent::Acid,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Amino,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Glycolyl,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("NeuGc"),
                    [
                        (
                            GlycanSubstituent::Acid,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Amino,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                        (
                            GlycanSubstituent::Deoxy,
                            GlycanSubstituentLocation::Single(
                                _,
                                GlycanSubstituentLeavingGroup::Hydroxy,
                            ),
                        ),
                    ] => Cow::Borrowed("Neu"),
                    _ => Cow::Owned(format!("{{{}}}", self.formula())),
                }
            }
            BaseSugar::Decose
                if self.substituents.is_empty()
                    && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Dec")
            }
            BaseSugar::Custom(_)
                if matches!(self.substituents.as_slice(), [(
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(_, GlycanSubstituentLeavingGroup::Hydroxy,),
                )]) && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Phosphate")
            }
            BaseSugar::Custom(_)
                if matches!(self.substituents.as_slice(), [(
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(_, GlycanSubstituentLeavingGroup::Hydroxy,),
                )]) && !self.furanose
                    && self.configuration.is_none() =>
            {
                Cow::Borrowed("Sulfate")
            }
            _ => Cow::Owned(format!("{{{}}}", self.formula())),
        }
    }

    /// Create a IUPAC condensed name for this monosaccharide
    pub fn iupac_condensed_name(&self) -> String {
        // TODO: recognise Fucose
        use std::fmt::Write;
        let mut name = String::new();

        match self.configuration {
            Some(Configuration::L) => name.push_str("L-"),
            Some(Configuration::D) => name.push_str("D-"),
            _ => (),
        }

        let mut substitutions = self
            .substituents
            .iter()
            .copied()
            .map(|(s, l)| {
                (
                    (s, match l {
                        GlycanSubstituentLocation::Single(_, l)
                        | GlycanSubstituentLocation::Double(_, l, _) => l,
                    }),
                    l,
                )
            })
            .into_group_map();
        let write_locs = |locations: &[GlycanSubstituentLocation],
                          symbol: &str,
                          leaving: GlycanSubstituentLeavingGroup,
                          name: &mut String| {
            let mut following = false;
            let mut double = false;
            for l in locations {
                if following {
                    name.push(',');
                } else {
                    following = true;
                }
                match l {
                    GlycanSubstituentLocation::Single(None, _) => {
                        name.push('?');
                    }
                    GlycanSubstituentLocation::Single(Some(l), _) => {
                        let _ = write!(name, "{l}");
                    }
                    GlycanSubstituentLocation::Double(None, ..) => {
                        double = true;
                        name.push_str("?-?");
                    }
                    GlycanSubstituentLocation::Double(Some((a, b)), ..) => {
                        double = true;
                        let _ = write!(name, "{a}-{b}");
                    }
                }
            }
            if leaving == GlycanSubstituentLeavingGroup::Hydrogen {
                name.push('C');
            }
            name.push_str(symbol);
            if !double && locations.len() > 1 {
                let _ = write!(name, "{}", locations.len());
            }
        };

        for locations in substitutions.values_mut() {
            locations.sort_unstable();
        }

        if let Some(locations) = substitutions.remove(&(
            GlycanSubstituent::Deoxy,
            GlycanSubstituentLeavingGroup::Hydroxy,
        )) {
            write_locs(
                &locations,
                "d",
                GlycanSubstituentLeavingGroup::Hydroxy,
                &mut name,
            );
        }

        if let Some(sugar) = self.base_sugar.iupac_name() {
            name.push_str(sugar);
        } else {
            name.push_str("???");
        }
        if self.furanose {
            name.push('f');
        }

        for ((substitution, leaving), locations) in
            substitutions.into_iter().sorted_by_key(|(_, l)| match l.first() {
                Some(
                    GlycanSubstituentLocation::Single(Some(n), _)
                    | GlycanSubstituentLocation::Double(Some((n, _)), ..),
                ) => Some(*n),
                _ => None,
            })
        {
            write_locs(&locations, substitution.notation(), leaving, &mut name);
        }

        name
    }

    /// Get the base sugar of this monosaccharide
    pub const fn base_sugar(&self) -> &BaseSugar {
        &self.base_sugar
    }

    /// Check if this is a fucose
    pub fn is_fucose(&self) -> bool {
        self.base_sugar == BaseSugar::Hexose(Some(HexoseIsomer::Galactose))
            && self.substituents.iter().any(|(s, ..)| *s == GlycanSubstituent::Deoxy)
    }

    /// Create a new monosaccharide
    pub fn new(
        sugar: BaseSugar,
        substituents: &[(GlycanSubstituent, GlycanSubstituentLocation)],
    ) -> Self {
        Self {
            base_sugar: sugar,
            substituents: substituents.iter().copied().sorted().collect(),
            furanose: false,
            configuration: None,
            #[cfg(debug_assertions)]
            source_field: String::new(),
        }
    }

    /// Set this saccharide up as to be a furanose
    #[must_use]
    pub fn furanose(self) -> Self {
        Self {
            furanose: true,
            ..self
        }
    }

    /// Set this saccharide up to be a certain configuration
    #[must_use]
    pub fn configuration(self, configuration: Configuration) -> Self {
        Self {
            configuration: Some(configuration),
            ..self
        }
    }

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

    /// Get all unique combinations of monosaccharides within the given range of number of
    /// monosaccharides used
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
                            None | Some(std::cmp::Ordering::Equal) => result.push(new), /* Cannot get bigger */
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
    use std::num::NonZeroU8;

    use itertools::Itertools;

    use crate::{
        chemistry::Molecule,
        glycan::{
            glycan::{
                BaseSugar, Configuration, GlycanSubstituent, GlycanSubstituentLeavingGroup,
                GlycanSubstituentLocation, HexoseIsomer, PentoseIsomer,
            },
            lists::GLYCAN_PARSE_LIST,
            monosaccharide::MonoSaccharide,
        },
        molecular_formula,
        parse_json::ParseJson,
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
            assert_eq!(found.1.formula(), *formula, "Formula incorrect: {name}");
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
                    assert_eq!(str.len(), len, "Did not fully consume: '{str}'");
                    res
                })
                .unwrap()
        };
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Pyruvyl,
                    GlycanSubstituentLocation::Double(
                        None,
                        GlycanSubstituentLeavingGroup::Hydroxy,
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Pyruvyl,
                    GlycanSubstituentLocation::Double(
                        None,
                        GlycanSubstituentLeavingGroup::Hydroxy,
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
            ])
        );
        assert_eq!(
            parse("Gal2,3Ac2S"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
                ),
            ])
        );
        assert_eq!(
            parse("GlcNAc"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[(
                GlycanSubstituent::NAcetyl,
                GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
            )])
        );
        assert_eq!(
            parse("Gal6S"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[(
                GlycanSubstituent::Sulfate,
                GlycanSubstituentLocation::Single(
                    NonZeroU8::new(6),
                    GlycanSubstituentLeavingGroup::Hydroxy
                )
            )])
        );
        assert_eq!(
            parse("GlcN2Gc"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
                (
                    GlycanSubstituent::Amino,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
                ),
                (
                    GlycanSubstituent::Glycolyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
            ])
        );
        assert_eq!(
            parse("GalNAc3S"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
                (
                    GlycanSubstituent::NAcetyl,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
                ),
                (
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                )
            ])
        );
        assert_eq!(
            parse("GlcN2,6S2"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
                (
                    GlycanSubstituent::Amino,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
                ),
                (
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(6),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                )
            ])
        );
        assert_eq!(
            parse("Tagf1,6P2"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Tagatose)), &[
                (
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(1),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Phosphate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(6),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                )
            ])
            .furanose()
        );
        assert_eq!(
            parse("Gal2,3Ac24-1,6-1Py"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Pyruvyl,
                    GlycanSubstituentLocation::Double(
                        None,
                        GlycanSubstituentLeavingGroup::Hydroxy,
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Pyruvyl,
                    GlycanSubstituentLocation::Double(
                        None,
                        GlycanSubstituentLeavingGroup::Hydroxy,
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
            ])
        );
        assert_eq!(
            parse("D-Araf"),
            MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Arabinose)), &[])
                .furanose()
                .configuration(Configuration::D)
        );
        assert_eq!(
            parse("Xyl-onic"),
            MonoSaccharide::new(BaseSugar::Pentose(Some(PentoseIsomer::Xylose)), &[(
                GlycanSubstituent::Acid,
                GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
            )])
        );
        assert_eq!(
            parse("Glc2,3,4,6Ac4"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Glucose)), &[
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(2),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(4),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                ),
                (
                    GlycanSubstituent::Acetyl,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(6),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                )
            ])
        );
        // This is an 'alternative LIP' (in WURCS terms) so the sulfate is supposed to be at 3 OR 4,
        // but just force it as 3 for now.
        assert_eq!(
            parse("GalNAc3/4S"),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[
                (
                    GlycanSubstituent::NAcetyl,
                    GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
                ),
                (
                    GlycanSubstituent::Sulfate,
                    GlycanSubstituentLocation::Single(
                        NonZeroU8::new(3),
                        GlycanSubstituentLeavingGroup::Hydroxy
                    )
                )
            ])
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
                .map(|option| option.iter().map(|sug| format!("{}{}", sug.0, sug.1)).join("&"))
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
            MonoSaccharide::pro_forma_composition::<true>("Hex3Phosphate1").unwrap().0
        );
        println!("{w:?}");
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].get_contexts().len(), 2);
        let (res, w) = MonoSaccharide::pro_forma_composition::<false>("Man2Man-1P").unwrap();
        assert_eq!(
            res,
            MonoSaccharide::pro_forma_composition::<true>("Hex1Phosphate1").unwrap().0
        );
        assert_eq!(w.len(), 1);
        // Maybe add warning that you are mixing single letter definitions with full names?
    }

    #[test]
    fn parse_json() {
        let hex = r#"{"base_sugar": {"Hexose": null},"substituents": [],"furanose": false,"configuration": null,"proforma_name": "Hex"}"#;
        let hexnac = r#"{"base_sugar": {"Hexose": null},"substituents": ["NAcetyl"],"furanose": false,"configuration": null,"proforma_name": "HexNAc"}"#;
        let fuc = r#"{"base_sugar": {"Hexose": "Galactose"},"substituents": ["Deoxy"],"furanose": false,"configuration": null,"proforma_name": "Fuc"}"#;

        assert_eq!(
            MonoSaccharide::from_json(hex).unwrap(),
            MonoSaccharide::new(BaseSugar::Hexose(None), &[])
        );
        assert_eq!(
            MonoSaccharide::from_json(hexnac).unwrap(),
            MonoSaccharide::new(BaseSugar::Hexose(None), &[(
                GlycanSubstituent::NAcetyl,
                GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
            )])
        );
        assert_eq!(
            MonoSaccharide::from_json(fuc).unwrap(),
            MonoSaccharide::new(BaseSugar::Hexose(Some(HexoseIsomer::Galactose)), &[(
                GlycanSubstituent::Deoxy,
                GlycanSubstituentLocation::Single(None, GlycanSubstituentLeavingGroup::Hydroxy)
            )])
        );
    }
}
