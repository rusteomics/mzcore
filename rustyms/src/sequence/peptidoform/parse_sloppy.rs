use std::sync::{Arc, LazyLock};

use custom_error::*;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    glycan::GLYCAN_PARSE_LIST,
    helper_functions::{ResultExtensions, end_of_enclosure, parse_named_counter},
    ontology::{CustomDatabase, Ontology},
    sequence::{
        AminoAcid, CheckedAminoAcid, Modification, PeptideModificationSearch, Peptidoform,
        SemiAmbiguous, SequenceElement, SequencePosition, SimpleModification,
        SimpleModificationInner, peptidoform::parse_modification,
    },
    system::Mass,
};

/// Parameters to control the parsing of 'sloppy' ProForma sequences.
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[serde(bound(deserialize = "'de: 'static"))]
pub struct SloppyParsingParameters {
    /// Ignore a prefix lowercase n as in `n[211]GC[779]RQSSEEK` as this indicates an N terminal modification in MSFragger
    pub ignore_prefix_lowercase_n: bool,
    /// Allow AA+12AA instead of AA[+12]AA as used by Casanovo
    pub allow_unwrapped_modifications: bool,
    /// Allow Xmod as modification indication as used by DeepNovo
    pub mod_indications: (Option<&'static str>, Vec<(AminoAcid, SimpleModification)>),
    /// Support for custom encodings, e.g. `AAAmAAA` instead of `AAAM[oxidation]AAA` as used by NovoB
    pub custom_alphabet: Vec<(u8, SequenceElement<SemiAmbiguous>)>,
    /// Replacing mass mods with known predefined mods, e.g. `AAA(+79.97)AAA` instead of `AAA[phospho]AAA` as used by InstaNovo
    pub replace_mass_modifications: Option<Vec<SimpleModification>>,
}

impl Peptidoform<SemiAmbiguous> {
    /// Read a sequence as a ProForma sequence, or if that fails try to read it as using
    /// [`Self::sloppy_pro_forma`].
    ///
    /// # Errors
    /// If both parsers fail. It returns an error that combines the feedback from both parsers.
    pub fn pro_forma_or_sloppy(
        line: &str,
        location: std::ops::Range<usize>,
        custom_database: Option<&CustomDatabase>,
        parameters: &SloppyParsingParameters,
    ) -> Result<Self, CustomError> {
        Peptidoform::pro_forma(&line[location.clone()], custom_database).and_then(|p| p.into_semi_ambiguous().ok_or_else(|| 
            CustomError::error(
                "Peptidoform too complex",
                "A peptidoform as used here should not contain any complex parts of the ProForma specification, only amino acids and simple placed modifications are allowed",
                Context::line_range(None, line, location.clone()),
            ))).or_else(|pro_forma_error| 
                Self::sloppy_pro_forma(line, location.clone(), custom_database, parameters)
                .map_err(|sloppy_error| 
                    CustomError::error(
                        "Invalid peptidoform", 
                        "The sequence could not be parsed as a ProForma nor as a more loosly defined peptidoform, see the underlying errors for details", 
                        Context::line_range(None, line, location.clone())).with_underlying_errors(vec![pro_forma_error, sloppy_error])))
    }

    /// Read sloppy ProForma like sequences. Defined by the use of square or round braces to indicate
    /// modifications and missing any particular method of defining the N or C terminal modifications.
    /// Additionally, any underscores will be ignored both on the ends and inside the sequence.
    ///
    /// All modifications follow the same definitions as the strict ProForma syntax, if it cannot be
    /// parsed as a strict ProForma modification it falls back to [`Modification::sloppy_modification`].
    ///
    /// # Errors
    /// If it does not fit the above description.
    #[expect(clippy::missing_panics_doc)] // Cannot panic
    pub fn sloppy_pro_forma<'a>(
        line: &'a str,
        location: std::ops::Range<usize>,
        custom_database: Option<&CustomDatabase>,
        parameters: &SloppyParsingParameters,
    ) -> Result<Self, BoxedError<'a>> {
        if line[location.clone()].trim().is_empty() {
            return Err(BoxedError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, location.start, 1),
            ));
        }

        let mut peptide = Self::default();
        let chars: &[u8] = line[location.clone()].as_bytes();
        peptide
            .sequence_mut()
            .reserve(chars.iter().map(u8::is_ascii_uppercase).count()); // Reserve approximately the right length for the vector, this will overestimate in some cases but not by a lot
        let mut index = 0;

        while index < chars.len() {
            match chars[index] {
                b'n' if parameters.ignore_prefix_lowercase_n && index == 0 => index += 1, //ignore
                b',' | b'_' => index += 1,                                                //ignore
                b'[' | b'(' => {
                    let (open, close) = if chars[index] == b'[' {
                        (b'[', b']')
                    } else {
                        (b'(', b')')
                    };
                    let end_index =
                        end_of_enclosure(&line[location.clone()], index + 1, open, close)
                            .ok_or_else(|| {
                                BoxedError::error(
                                    "Invalid modification",
                                    "No valid closing delimiter",
                                    Context::line(None, line, location.start + index, 1),
                                )
                            })?;
                    let modification = Modification::sloppy_modification(
                        line,
                        location.start + index + 1..location.start + end_index,
                        peptide.sequence().last(),
                        custom_database,
                    )
                    .map(Modification::Simple)?;
                    index = end_index + 1;

                    let pep_len = peptide.len();
                    let n_term_empty = peptide.get_n_term().is_empty();
                    match peptide.sequence_mut().last_mut() {
                        Some(aa) => {
                            if pep_len == 1
                                && !modification
                                    .is_possible(aa, SequencePosition::Index(0))
                                    .any_possible()
                                && modification
                                    .is_possible(aa, SequencePosition::NTerm)
                                    .any_possible()
                                && n_term_empty
                            {
                                peptide.add_simple_n_term(
                                    modification
                                        .simple()
                                        .expect(
                                            "Can only put a simple modification on an N terminus.",
                                        )
                                        .clone(),
                                );
                            } else {
                                aa.modifications.push(modification);
                            }
                        }
                        None => {
                            peptide.add_simple_n_term(
                                modification
                                    .simple()
                                    .expect("Can only put a simple modification on an N terminus.")
                                    .clone(),
                            );
                        }
                    }
                }
                _ if parameters.mod_indications.0.is_some_and(|pattern| {
                    line[location.start + index..location.end].starts_with(pattern)
                }) =>
                {
                    index += parameters
                        .mod_indications
                        .0
                        .map(str::len)
                        .unwrap_or_default();

                    match peptide.sequence_mut().last_mut() {
                        Some(seq) => parameters
                            .mod_indications
                            .1
                            .iter()
                            .find(|(aa, _)| *aa == seq.aminoacid.aminoacid())
                            .map(|(_, m)| seq.modifications.push(Modification::Simple(m.clone())))
                            .ok_or_else(|| {
                                BoxedError::error(
                                    "Invalid mod indication",
                                    "There is no given mod for this amino acid.",
                                    Context::line(None, line, location.start + index - 4, 4),
                                )
                            })?,
                        None => {
                            return Err(BoxedError::error(
                                "Invalid mod indication",
                                "A mod indication should always follow an amino acid.",
                                Context::line(None, line, location.start + index - 3, 3),
                            ));
                        }
                    }
                }

                ch if parameters.allow_unwrapped_modifications
                    && (ch == b'-' || ch == b'+' || ch.is_ascii_digit()) =>
                {
                    let length = 1 + chars[index + 1..]
                        .iter()
                        .take_while(|c| c.is_ascii_digit() || **c == b'.')
                        .count();
                    let modification = SimpleModificationInner::Mass(Mass::new::<crate::system::dalton>(
                        line[location.start + index..location.start + index + length]
                        .parse::<f64>()
                        .map_err(|err|
                            BoxedError::error(
                                "Invalid mass shift modification", 
                                format!("Mass shift modification must be a valid number but this number is invalid: {err}"), 
                                Context::line(None, line, location.start + index, length))
                            )?).into()).into();
                    match peptide.sequence_mut().last_mut() {
                        Some(aa) => aa.modifications.push(Modification::Simple(modification)),
                        None => {
                            peptide.add_simple_n_term(modification);
                        }
                    }
                    index += length;
                }
                ch => {
                    if let Some(seq) = parameters
                        .custom_alphabet
                        .iter()
                        .find_map(|(c, seq)| (*c == ch).then_some(seq))
                    {
                        peptide.sequence_mut().push(seq.clone());
                    } else {
                        peptide.sequence_mut().push(SequenceElement::new(
                            ch.try_into().map_err(|()| {
                                BoxedError::error(
                                    "Invalid amino acid",
                                    "This character is not a valid amino acid",
                                    Context::line(None, line, location.start + index, 1),
                                )
                            })?,
                            None,
                        ));
                    }
                    index += 1;
                }
            }
        }
        if peptide.is_empty() {
            return Err(BoxedError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, location.start, location.len()),
            ));
        }
        peptide.enforce_modification_rules()?;
        Ok(
            if let Some(modifications) = parameters.replace_mass_modifications.clone() {
                PeptideModificationSearch::in_modifications(modifications)
                    .tolerance(crate::quantities::Tolerance::Absolute(crate::system::da(
                        0.05,
                    )))
                    .search(peptide)
            } else {
                peptide
            },
        )
    }
}

static SLOPPY_MOD_OPAIR_REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(?:[^:]+:)?(.*) (?:(?:on)|(?:from)) ([A-Z])").unwrap());
static SLOPPY_MOD_ON_REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(.*)\s*\([- @a-zA-Z]+\)").unwrap());
static SLOPPY_MOD_NUMERIC_END_REGEX: LazyLock<Regex> =
    LazyLock::new(|| Regex::new(r"(.*)\d+").unwrap());

impl Modification {
    /// Parse a modification defined by sloppy names
    /// # Errors
    /// If the name is not in Unimod, PSI-MOD, the custom database, or the predefined list of common trivial names.
    /// Or if this is the case when the modification follows a known structure (eg `mod (AAs)`).
    pub fn sloppy_modification<'a>(
        line: &'a str,
        location: std::ops::Range<usize>,
        position: Option<&SequenceElement<SemiAmbiguous>>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<SimpleModification, BoxedError<'a>> {
        let full_context = Context::line(None, line, location.start, location.len());
        let name = &line[location];

        Self::find_name(name, position, custom_database)
            .or_else( || {
                match name.trim().to_lowercase().split_once(':') {
                    Some(("u", tail)) => Ontology::Unimod.find_name(tail, None),
                    Some(("unimod", tail)) => Ontology::Unimod.find_id(tail.parse::<usize>().ok()?, None),
                    Some(("m", tail)) => Ontology::Psimod.find_name(tail, None),
                    Some(("c", tail)) => Ontology::Custom.find_name(tail, custom_database),
                    _ => None
                }
            })
            .or_else( || {
                name.trim().split_ascii_whitespace().next().and_then(|head| Self::find_name::<SemiAmbiguous>(head, position, custom_database))
            })
            .or_else(|| {SLOPPY_MOD_OPAIR_REGEX
                .captures(name)
                .and_then(|capture| {
                    let pos = capture[2].chars().next().and_then(|a| AminoAcid::try_from(a).ok().map(|a| SequenceElement::new(CheckedAminoAcid::new(a), None)));
                    Self::find_name::<SemiAmbiguous>(&capture[1], position.or(pos.as_ref()), custom_database)
                        .ok_or_else(|| {
                            parse_named_counter(
                                &capture[1].to_ascii_lowercase(),
                                &GLYCAN_PARSE_LIST,
                                false,
                            )
                            .map(|g| Arc::new(SimpleModificationInner::Glycan(g)))
                        })
                        .flat_err()
                        .ok()
                })
                .or_else(|| {
                    // Common sloppy naming: `modification (AAs)` also accepts `modification (Protein N-term)`
                    SLOPPY_MOD_ON_REGEX
                        .captures(name)
                        .and_then(|capture| {
                            Self::find_name(&capture[1], position, custom_database)
                        })
                })
                .or_else(|| {
                    // Common sloppy naming: `modification1`
                    SLOPPY_MOD_NUMERIC_END_REGEX
                        .captures(name)
                        .and_then(|capture| {
                            Self::find_name(&capture[1], position, custom_database)
                        })
                })
            }).ok_or_else(|| {
                BoxedError::error(
                    "Could not interpret modification",
                    "Modifications have to be defined as a number, Unimod, or PSI-MOD name, if this is a custom modification make sure to add it to the database",
                    full_context,
                ).suggestions(
                    Ontology::find_closest_many(
                        &[Ontology::Unimod, Ontology::Psimod],
                        &name.trim().to_lowercase(),
                        custom_database).get_suggestions().iter().map(ToString::to_string))
            })
    }

    fn find_name<T>(
        name: &str,
        position: Option<&SequenceElement<T>>,
        custom_database: Option<&CustomDatabase>,
    ) -> Option<SimpleModification> {
        let name = name.trim().to_lowercase();
        match name.as_str() {
            "o" | "ox" => Ontology::Unimod.find_id(35, None), // oxidation
            "cam" | "carbamidomethylation" => Ontology::Unimod.find_id(4, None), // carbamidomethyl
            "nem" => Ontology::Unimod.find_id(108, None),     // Nethylmaleimide
            "deamidation" => Ontology::Unimod.find_id(7, None), // deamidated
            "formylation" => Ontology::Unimod.find_id(122, None), // formyl
            "pyro-glu" => Ontology::Unimod.find_id(
                if position.is_some_and(|p| p.aminoacid.aminoacid() == AminoAcid::GlutamicAcid) {
                    27
                } else {
                    28
                },
                None,
            ), // pyro Glu with the logic to pick the correct modification based on the amino acid it is placed on
            "sub a" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Arginine => Some(1189),
                    AminoAcid::Asparagine => Some(1155),
                    AminoAcid::AsparticAcid => Some(553),
                    AminoAcid::Cysteine => Some(1055),
                    AminoAcid::Glutamine => Some(1177),
                    AminoAcid::GlutamicAcid => Some(560),
                    AminoAcid::Glycine => Some(571),
                    AminoAcid::Histidine => Some(1113),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1125)
                    }
                    AminoAcid::Lysine => Some(1131),
                    AminoAcid::Methionine => Some(1142),
                    AminoAcid::Phenylalanine => Some(1090),
                    AminoAcid::Proline => Some(624),
                    AminoAcid::Serine => Some(648),
                    AminoAcid::Threonine => Some(659),
                    AminoAcid::Tryptophan => Some(1224),
                    AminoAcid::Tyrosine => Some(1237),
                    AminoAcid::Valine => Some(667),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub c" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(1044),
                    AminoAcid::Arginine => Some(644),
                    AminoAcid::Asparagine => Some(1156),
                    AminoAcid::AsparticAcid => Some(1067),
                    AminoAcid::Glutamine => Some(1178),
                    AminoAcid::GlutamicAcid => Some(1078),
                    AminoAcid::Glycine => Some(577),
                    AminoAcid::Histidine => Some(1114),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1126)
                    }
                    AminoAcid::Lysine => Some(1132),
                    AminoAcid::Methionine => Some(1143),
                    AminoAcid::Phenylalanine => Some(567),
                    AminoAcid::Proline => Some(1166),
                    AminoAcid::Serine => Some(654),
                    AminoAcid::Threonine => Some(1203),
                    AminoAcid::Tryptophan => Some(674),
                    AminoAcid::Tyrosine => Some(683),
                    AminoAcid::Valine => Some(1213),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub d" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(542),
                    AminoAcid::Arginine => Some(1190),
                    AminoAcid::Asparagine => Some(621),
                    AminoAcid::Cysteine => Some(1056),
                    AminoAcid::Glutamine => Some(1179),
                    AminoAcid::GlutamicAcid => Some(562),
                    AminoAcid::Glycine => Some(576),
                    AminoAcid::Histidine => Some(349),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1127)
                    }
                    AminoAcid::Lysine => Some(1133),
                    AminoAcid::Methionine => Some(1144),
                    AminoAcid::Phenylalanine => Some(1091),
                    AminoAcid::Proline => Some(1167),
                    AminoAcid::Serine => Some(1196),
                    AminoAcid::Threonine => Some(1204),
                    AminoAcid::Tryptophan => Some(1225),
                    AminoAcid::Tyrosine => Some(682),
                    AminoAcid::Valine => Some(670),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub e" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(545),
                    AminoAcid::Arginine => Some(1191),
                    AminoAcid::Asparagine => Some(1157),
                    AminoAcid::AsparticAcid => Some(558),
                    AminoAcid::Cysteine => Some(1057),
                    AminoAcid::Glutamine => Some(632),
                    AminoAcid::Glycine => Some(574),
                    AminoAcid::Histidine => Some(1115),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1128)
                    }
                    AminoAcid::Lysine => Some(596),
                    AminoAcid::Methionine => Some(1145),
                    AminoAcid::Phenylalanine => Some(1092),
                    AminoAcid::Proline => Some(1168),
                    AminoAcid::Serine => Some(1197),
                    AminoAcid::Threonine => Some(1205),
                    AminoAcid::Tryptophan => Some(1226),
                    AminoAcid::Tyrosine => Some(1238),
                    AminoAcid::Valine => Some(668),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub f" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(1045),
                    AminoAcid::Arginine => Some(1195),
                    AminoAcid::Asparagine => Some(1158),
                    AminoAcid::AsparticAcid => Some(1068),
                    AminoAcid::Cysteine => Some(547),
                    AminoAcid::Glutamine => Some(1180),
                    AminoAcid::GlutamicAcid => Some(1079),
                    AminoAcid::Glycine => Some(1103),
                    AminoAcid::Histidine => Some(1116),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(602)
                    }
                    AminoAcid::Lysine => Some(1134),
                    AminoAcid::Methionine => Some(1146),
                    AminoAcid::Proline => Some(1169),
                    AminoAcid::Serine => Some(647),
                    AminoAcid::Threonine => Some(1206),
                    AminoAcid::Tryptophan => Some(1227),
                    AminoAcid::Tyrosine => Some(678),
                    AminoAcid::Valine => Some(666),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub g" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(544),
                    AminoAcid::Arginine => Some(646),
                    AminoAcid::Asparagine => Some(1159),
                    AminoAcid::AsparticAcid => Some(556),
                    AminoAcid::Cysteine => Some(552),
                    AminoAcid::Glutamine => Some(1181),
                    AminoAcid::GlutamicAcid => Some(564),
                    AminoAcid::Histidine => Some(1117),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1129)
                    }
                    AminoAcid::Lysine => Some(1135),
                    AminoAcid::Methionine => Some(1147),
                    AminoAcid::Phenylalanine => Some(1093),
                    AminoAcid::Proline => Some(1170),
                    AminoAcid::Serine => Some(657),
                    AminoAcid::Threonine => Some(1207),
                    AminoAcid::Tryptophan => Some(676),
                    AminoAcid::Tyrosine => Some(1239),
                    AminoAcid::Valine => Some(672),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub h" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(1046),
                    AminoAcid::Arginine => Some(641),
                    AminoAcid::Asparagine => Some(620),
                    AminoAcid::AsparticAcid => Some(554),
                    AminoAcid::Cysteine => Some(1058),
                    AminoAcid::Glutamine => Some(1181),
                    AminoAcid::GlutamicAcid => Some(633),
                    AminoAcid::Glycine => Some(1104),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(606)
                    }
                    AminoAcid::Lysine => Some(1136),
                    AminoAcid::Methionine => Some(1148),
                    AminoAcid::Phenylalanine => Some(1094),
                    AminoAcid::Proline => Some(625),
                    AminoAcid::Serine => Some(1198),
                    AminoAcid::Threonine => Some(1208),
                    AminoAcid::Tryptophan => Some(1228),
                    AminoAcid::Tyrosine => Some(681),
                    AminoAcid::Valine => Some(1214),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub i" | "sub l" | "sub j" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Alanine => Some(1047),
                    AminoAcid::Arginine => Some(645),
                    AminoAcid::Asparagine => Some(622),
                    AminoAcid::AsparticAcid => Some(1069),
                    AminoAcid::Cysteine => Some(1059),
                    AminoAcid::Glutamine => Some(635),
                    AminoAcid::GlutamicAcid => Some(1081),
                    AminoAcid::Histidine => Some(585),
                    AminoAcid::Glycine => Some(1105),
                    AminoAcid::Lysine => Some(600),
                    AminoAcid::Methionine => Some(614),
                    AminoAcid::Phenylalanine => Some(568),
                    AminoAcid::Proline => Some(629),
                    AminoAcid::Serine => Some(656),
                    AminoAcid::Threonine => Some(664),
                    AminoAcid::Tryptophan => Some(677),
                    AminoAcid::Tyrosine => Some(1248),
                    AminoAcid::Valine => Some(671),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub k" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1229),
                    AminoAcid::Tyrosine => Some(1240),
                    AminoAcid::Arginine => Some(640),
                    AminoAcid::Phenylalanine => Some(1095),
                    AminoAcid::Histidine => Some(1119),
                    AminoAcid::Methionine => Some(613),
                    AminoAcid::GlutamicAcid => Some(563),
                    AminoAcid::Glutamine => Some(631),
                    AminoAcid::AsparticAcid => Some(1070),
                    AminoAcid::Asparagine => Some(618),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(590)
                    }
                    AminoAcid::Cysteine => Some(1060),
                    AminoAcid::Threonine => Some(661),
                    AminoAcid::Valine => Some(1215),
                    AminoAcid::Proline => Some(1171),
                    AminoAcid::Serine => Some(1199),
                    AminoAcid::Alanine => Some(1048),
                    AminoAcid::Glycine => Some(1106),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub m" => position // Some mod options are available
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1230),
                    AminoAcid::Tyrosine => Some(1241),
                    AminoAcid::Arginine => Some(643),
                    AminoAcid::Phenylalanine => Some(1096),
                    AminoAcid::Histidine => Some(1120),
                    AminoAcid::GlutamicAcid => Some(1082),
                    AminoAcid::Lysine => Some(598),
                    AminoAcid::Glutamine => Some(1182),
                    AminoAcid::AsparticAcid => Some(1071),
                    AminoAcid::Asparagine => Some(1160),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(608)
                    }
                    AminoAcid::Cysteine => Some(1061),
                    AminoAcid::Threonine => Some(663),
                    AminoAcid::Valine => Some(669),
                    AminoAcid::Proline => Some(1172),
                    AminoAcid::Serine => Some(1200),
                    AminoAcid::Alanine => Some(1049),
                    AminoAcid::Glycine => Some(1107),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub n" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1231),
                    AminoAcid::Tyrosine => Some(680),
                    AminoAcid::Arginine => Some(1192),
                    AminoAcid::Phenylalanine => Some(1097),
                    AminoAcid::Histidine => Some(348),
                    AminoAcid::Methionine => Some(1149),
                    AminoAcid::GlutamicAcid => Some(1083),
                    AminoAcid::Lysine => Some(595),
                    AminoAcid::Glutamine => Some(1183),
                    AminoAcid::AsparticAcid => Some(555),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(589)
                    }
                    AminoAcid::Cysteine => Some(1062),
                    AminoAcid::Threonine => Some(660),
                    AminoAcid::Valine => Some(1216),
                    AminoAcid::Proline => Some(1173),
                    AminoAcid::Serine => Some(651),
                    AminoAcid::Alanine => Some(1050),
                    AminoAcid::Glycine => Some(1108),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub p" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1232),
                    AminoAcid::Tyrosine => Some(1242),
                    AminoAcid::Arginine => Some(639),
                    AminoAcid::Phenylalanine => Some(1098),
                    AminoAcid::Histidine => Some(580),
                    AminoAcid::Methionine => Some(1150),
                    AminoAcid::GlutamicAcid => Some(1084),
                    AminoAcid::Lysine => Some(1137),
                    AminoAcid::Glutamine => Some(630),
                    AminoAcid::AsparticAcid => Some(1072),
                    AminoAcid::Asparagine => Some(1161),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(604)
                    }
                    AminoAcid::Cysteine => Some(1063),
                    AminoAcid::Threonine => Some(662),
                    AminoAcid::Valine => Some(1217),
                    AminoAcid::Serine => Some(652),
                    AminoAcid::Alanine => Some(543),
                    AminoAcid::Glycine => Some(1109),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub q" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1233),
                    AminoAcid::Tyrosine => Some(1243),
                    AminoAcid::Arginine => Some(642),
                    AminoAcid::Phenylalanine => Some(1099),
                    AminoAcid::Histidine => Some(582),
                    AminoAcid::Methionine => Some(1151),
                    AminoAcid::GlutamicAcid => Some(561),
                    AminoAcid::Lysine => Some(597),
                    AminoAcid::AsparticAcid => Some(1073),
                    AminoAcid::Asparagine => Some(1162),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(607)
                    }
                    AminoAcid::Cysteine => Some(1064),
                    AminoAcid::Threonine => Some(1209),
                    AminoAcid::Valine => Some(1218),
                    AminoAcid::Proline => Some(626),
                    AminoAcid::Serine => Some(1201),
                    AminoAcid::Alanine => Some(1051),
                    AminoAcid::Glycine => Some(1110),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub r" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(675),
                    AminoAcid::Tyrosine => Some(1244),
                    AminoAcid::Phenylalanine => Some(1100),
                    AminoAcid::Histidine => Some(584),
                    AminoAcid::Methionine => Some(611),
                    AminoAcid::GlutamicAcid => Some(1085),
                    AminoAcid::Lysine => Some(599),
                    AminoAcid::Glutamine => Some(634),
                    AminoAcid::AsparticAcid => Some(1074),
                    AminoAcid::Asparagine => Some(1163),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(609)
                    }
                    AminoAcid::Cysteine => Some(551),
                    AminoAcid::Threonine => Some(665),
                    AminoAcid::Valine => Some(1219),
                    AminoAcid::Proline => Some(628),
                    AminoAcid::Serine => Some(655),
                    AminoAcid::Alanine => Some(1052),
                    AminoAcid::Glycine => Some(578),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub s" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(673),
                    AminoAcid::Tyrosine => Some(679),
                    AminoAcid::Arginine => Some(636),
                    AminoAcid::Phenylalanine => Some(566),
                    AminoAcid::Histidine => Some(1121),
                    AminoAcid::Methionine => Some(1152),
                    AminoAcid::GlutamicAcid => Some(1086),
                    AminoAcid::Lysine => Some(1138),
                    AminoAcid::Glutamine => Some(1184),
                    AminoAcid::AsparticAcid => Some(1075),
                    AminoAcid::Asparagine => Some(616),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(601)
                    }
                    AminoAcid::Cysteine => Some(548),
                    AminoAcid::Threonine => Some(658),
                    AminoAcid::Valine => Some(1220),
                    AminoAcid::Proline => Some(623),
                    AminoAcid::Alanine => Some(540),
                    AminoAcid::Glycine => Some(572),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub t" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1234),
                    AminoAcid::Tyrosine => Some(1245),
                    AminoAcid::Arginine => Some(638),
                    AminoAcid::Phenylalanine => Some(1101),
                    AminoAcid::Histidine => Some(1122),
                    AminoAcid::Methionine => Some(610),
                    AminoAcid::GlutamicAcid => Some(1087),
                    AminoAcid::Lysine => Some(594),
                    AminoAcid::Glutamine => Some(1185),
                    AminoAcid::AsparticAcid => Some(1076),
                    AminoAcid::Asparagine => Some(617),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(588)
                    }
                    AminoAcid::Cysteine => Some(1065),
                    AminoAcid::Valine => Some(1221),
                    AminoAcid::Proline => Some(627),
                    AminoAcid::Serine => Some(650),
                    AminoAcid::Alanine => Some(541),
                    AminoAcid::Glycine => Some(1111),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub v" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1235),
                    AminoAcid::Tyrosine => Some(1246),
                    AminoAcid::Arginine => Some(1193),
                    AminoAcid::Phenylalanine => Some(570),
                    AminoAcid::Histidine => Some(1123),
                    AminoAcid::Methionine => Some(615),
                    AminoAcid::GlutamicAcid => Some(565),
                    AminoAcid::Lysine => Some(1139),
                    AminoAcid::Glutamine => Some(1186),
                    AminoAcid::AsparticAcid => Some(559),
                    AminoAcid::Asparagine => Some(1164),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(605)
                    }
                    AminoAcid::Cysteine => Some(1066),
                    AminoAcid::Threonine => Some(1210),
                    AminoAcid::Proline => Some(1174),
                    AminoAcid::Serine => Some(1202),
                    AminoAcid::Alanine => Some(546),
                    AminoAcid::Glycine => Some(575),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub w" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tyrosine => Some(1247),
                    AminoAcid::Arginine => Some(637),
                    AminoAcid::Phenylalanine => Some(1102),
                    AminoAcid::Histidine => Some(1124),
                    AminoAcid::Methionine => Some(1153),
                    AminoAcid::GlutamicAcid => Some(1088),
                    AminoAcid::Lysine => Some(1140),
                    AminoAcid::Glutamine => Some(1187),
                    AminoAcid::AsparticAcid => Some(1077),
                    AminoAcid::Asparagine => Some(1165),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(603)
                    }
                    AminoAcid::Cysteine => Some(549),
                    AminoAcid::Valine => Some(1211),
                    AminoAcid::Threonine => Some(1222),
                    AminoAcid::Proline => Some(1175),
                    AminoAcid::Serine => Some(649),
                    AminoAcid::Alanine => Some(1053),
                    AminoAcid::Glycine => Some(573),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "sub y" => position
                .filter(|p| p.modifications.is_empty())
                .map(|p| p.aminoacid.aminoacid())
                .and_then(|aa| match aa {
                    AminoAcid::Tryptophan => Some(1236),
                    AminoAcid::Arginine => Some(1194),
                    AminoAcid::Phenylalanine => Some(569),
                    AminoAcid::Histidine => Some(581),
                    AminoAcid::Methionine => Some(1154),
                    AminoAcid::GlutamicAcid => Some(1089),
                    AminoAcid::Lysine => Some(1141),
                    AminoAcid::Glutamine => Some(1188),
                    AminoAcid::AsparticAcid => Some(557),
                    AminoAcid::Asparagine => Some(619),
                    AminoAcid::Isoleucine | AminoAcid::AmbiguousLeucine | AminoAcid::Leucine => {
                        Some(1130)
                    }
                    AminoAcid::Cysteine => Some(550),
                    AminoAcid::Valine => Some(1212),
                    AminoAcid::Threonine => Some(1223),
                    AminoAcid::Proline => Some(1176),
                    AminoAcid::Serine => Some(653),
                    AminoAcid::Alanine => Some(1054),
                    AminoAcid::Glycine => Some(1112),
                    _ => None,
                })
                .and_then(|i| Ontology::Unimod.find_id(i, None)),
            "ala->ile" | "ala->leu" | "ala->xle" => Ontology::Unimod.find_id(1125, None),
            "cys->ile" | "cys->leu" | "cys->xle" => Ontology::Unimod.find_id(1126, None),
            "asp->ile" | "asp->leu" | "asp->xle" => Ontology::Unimod.find_id(1127, None),
            "glu->ile" | "glu->leu" | "glu->xle" => Ontology::Unimod.find_id(1128, None),
            "phe->ile" | "phe->leu" | "phe->xle" => Ontology::Unimod.find_id(602, None),
            "gly->ile" | "gly->leu" | "gly->xle" => Ontology::Unimod.find_id(1129, None),
            "his->ile" | "his->leu" | "his->xle" => Ontology::Unimod.find_id(606, None),
            "lys->ile" | "lys->leu" | "lys->xle" => Ontology::Unimod.find_id(590, None),
            "met->ile" | "met->leu" | "met->xle" => Ontology::Unimod.find_id(608, None),
            "asn->ile" | "asn->leu" | "asn->xle" => Ontology::Unimod.find_id(589, None),
            "pro->ile" | "pro->leu" | "pro->xle" => Ontology::Unimod.find_id(604, None),
            "gln->ile" | "gln->leu" | "gln->xle" => Ontology::Unimod.find_id(607, None),
            "arg->ile" | "arg->leu" | "arg->xle" => Ontology::Unimod.find_id(609, None),
            "ser->ile" | "ser->leu" | "ser->xle" => Ontology::Unimod.find_id(601, None),
            "thr->ile" | "thr->leu" | "thr->xle" => Ontology::Unimod.find_id(588, None),
            "val->ile" | "val->leu" | "val->xle" => Ontology::Unimod.find_id(605, None),
            "trp->ile" | "trp->leu" | "trp->xle" => Ontology::Unimod.find_id(603, None),
            "tyr->ile" | "tyr->leu" | "tyr->xle" => Ontology::Unimod.find_id(1130, None),
            "ile->ala" | "leu->ala" | "xle->ala" => Ontology::Unimod.find_id(1047, None),
            "ile->arg" | "leu->arg" | "xle->arg" => Ontology::Unimod.find_id(645, None),
            "ile->asn" | "leu->asn" | "xle->asn" => Ontology::Unimod.find_id(622, None),
            "ile->asp" | "leu->asp" | "xle->asp" => Ontology::Unimod.find_id(1069, None),
            "ile->cys" | "leu->cys" | "xle->cys" => Ontology::Unimod.find_id(1059, None),
            "ile->gln" | "leu->gln" | "xle->gln" => Ontology::Unimod.find_id(635, None),
            "ile->glu" | "leu->glu" | "xle->glu" => Ontology::Unimod.find_id(1081, None),
            "ile->his" | "leu->his" | "xle->his" => Ontology::Unimod.find_id(585, None),
            "ile->gly" | "leu->gly" | "xle->gly" => Ontology::Unimod.find_id(1105, None),
            "ile->lys" | "leu->lys" | "xle->lys" => Ontology::Unimod.find_id(600, None),
            "ile->met" | "leu->met" | "xle->met" => Ontology::Unimod.find_id(614, None),
            "ile->phe" | "leu->phe" | "xle->phe" => Ontology::Unimod.find_id(568, None),
            "ile->pro" | "leu->pro" | "xle->pro" => Ontology::Unimod.find_id(629, None),
            "ile->ser" | "leu->ser" | "xle->ser" => Ontology::Unimod.find_id(656, None),
            "ile->thr" | "leu->thr" | "xle->thr" => Ontology::Unimod.find_id(664, None),
            "ile->trp" | "leu->trp" | "xle->trp" => Ontology::Unimod.find_id(677, None),
            "ile->tyr" | "leu->tyr" | "xle->tyr" => Ontology::Unimod.find_id(1248, None),
            "ile->val" | "leu->val" | "xle->val" => Ontology::Unimod.find_id(671, None),
            _ => parse_modification::numerical_mod(&name)
                .ok()
                .or_else(|| Ontology::Unimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Psimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Custom.find_name(&name, custom_database)),
        }
    }
}
