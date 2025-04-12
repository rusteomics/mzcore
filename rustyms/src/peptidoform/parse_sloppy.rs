use std::sync::{Arc, LazyLock};

use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::{
    checked_aminoacid::CheckedAminoAcid,
    error::{Context, CustomError},
    glycan::GLYCAN_PARSE_LIST,
    helper_functions::{end_of_enclosure, parse_named_counter, ResultExtensions},
    modification::{Modification, Ontology, SimpleModification, SimpleModificationInner},
    ontologies::CustomDatabase,
    peptidoform::*,
    system::Mass,
    AminoAcid, SequenceElement,
};

/// Parameters to control the parsing of 'sloppy' ProForma sequences.
#[derive(Debug, Clone, Default, Eq, PartialEq, Hash, Serialize, Deserialize)]
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
    /// Read sloppy ProForma like sequences. Defined by the use of square or round braces to indicate
    /// modifications and missing any particular method of defining the N or C terminal modifications.
    /// Additionally any underscores will be ignored both on the ends and inside the sequence.
    ///
    /// All modifications follow the same definitions as the strict ProForma syntax, if it cannot be
    /// parsed as a strict ProForma modification it falls back to [`Modification::sloppy_modification`].
    ///
    /// # Errors
    /// If it does not fit the above description.
    #[expect(clippy::missing_panics_doc)] // Cannot panic
    pub fn sloppy_pro_forma(
        line: &str,
        location: std::ops::Range<usize>,
        custom_database: Option<&CustomDatabase>,
        parameters: &SloppyParsingParameters,
    ) -> Result<Self, CustomError> {
        if line[location.clone()].trim().is_empty() {
            return Err(CustomError::error(
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
                                CustomError::error(
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
                                    .is_possible(aa, crate::SequencePosition::Index(0))
                                    .any_possible()
                                && modification
                                    .is_possible(aa, crate::SequencePosition::NTerm)
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
                                CustomError::error(
                                    "Invalid mod indication",
                                    "There is no given mod for this amino acid.",
                                    Context::line(None, line, location.start + index - 4, 4),
                                )
                            })?,
                        None => {
                            return Err(CustomError::error(
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
                            CustomError::error(
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
                                CustomError::error(
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
            return Err(CustomError::error(
                "Peptide sequence is empty",
                "A peptide sequence cannot be empty",
                Context::line(None, line, location.start, location.len()),
            ));
        }
        peptide.enforce_modification_rules()?;
        Ok(
            if let Some(modifications) = parameters.replace_mass_modifications.clone() {
                PeptideModificationSearch::in_modifications(modifications)
                    .tolerance(crate::Tolerance::Absolute(crate::system::da(0.05)))
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
    #[expect(clippy::missing_panics_doc)]
    pub fn sloppy_modification(
        line: &str,
        location: std::ops::Range<usize>,
        position: Option<&SequenceElement<SemiAmbiguous>>,
        custom_database: Option<&CustomDatabase>,
    ) -> Result<SimpleModification, CustomError> {
        let full_context = Context::line(None, line, location.start, location.len());
        let name = &line[location];

        Self::find_name(name, position, custom_database)
            .or_else( || {
                match name.trim().to_lowercase().split_once(':') {
                    Some(("u", tail)) => Ontology::Unimod.find_name(tail, None),
                    Some(("m", tail)) => Ontology::Psimod.find_name(tail, None),
                    Some(("c", tail)) => Ontology::Custom.find_name(tail, custom_database),
                    _ => None
                }
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
                CustomError::error(
                    "Could not interpret modification",
                    "Modifications have to be defined as a number, Unimod, or PSI-MOD name, if this is a custom modification make sure to add it to the database",
                    full_context,
                ).with_suggestions(
                    Ontology::find_closest_many(
                        &[Ontology::Unimod, Ontology::Psimod],
                        &name.trim().to_lowercase(),
                        custom_database).suggestions())
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
            _ => crate::peptidoform::parse_modification::numerical_mod(&name)
                .ok()
                .or_else(|| Ontology::Unimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Psimod.find_name(&name, custom_database))
                .or_else(|| Ontology::Custom.find_name(&name, custom_database)),
        }
    }
}
