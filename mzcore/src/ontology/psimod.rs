//! Code to handle the PSI-MOD ontology
use std::{borrow::Cow, str::FromStr};

use context_error::{
    BoxedError, Context, CreateError, FullErrorContent, StaticErrorContent, combine_error,
};
use itertools::Itertools;

use mzcv::{
    CVData, CVError, CVFile, CVSource, CVVersion, HashBufReader, OboIdentifier, OboOntology,
    OboStanzaType, SynonymScope,
};

use crate::{
    chemistry::MolecularFormula,
    helper_functions::UnwrapInfallible,
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        LinkerSpecificity, ModificationId, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};

/// PSI-MOD modifications
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct PsiMod {}

impl CVData for SimpleModificationInner {
    type Index = u32;
    fn index(&self) -> Option<u32> {
        self.description().and_then(ModificationId::id)
    }
    fn name(&self) -> Option<Cow<'_, str>> {
        self.description().map(|d| Cow::Borrowed(d.name.as_ref()))
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.description().into_iter().flat_map(|d| {
            d.synonyms
                .iter()
                .filter(|(s, _)| *s == SynonymScope::Exact)
                .map(|(_, s)| s.as_ref())
        })
    }
}

impl CVSource for PsiMod {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "PSI-MOD"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "PSI-MOD",
            extension: "obo",
            url: Some(
                "https://raw.githubusercontent.com/HUPO-PSI/psi-mod-CV/refs/heads/master/PSI-MOD.obo",
            ),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/psimod.dat"),
                Configuration::default(),
            )
            .unwrap()
            .0;
            Some(cache)
        }
        #[cfg(feature = "internal-no-data")]
        None
    }

    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<
        (
            CVVersion,
            Self::Structure,
            Vec<BoxedError<'static, CVError>>,
        ),
        Vec<BoxedError<'static, CVError>>,
    > {
        let reader = reader.next().ok_or_else(|| {
            vec![BoxedError::new(
                CVError::MissingReader,
                "Missing reader",
                "One file reader should be given for parsing PSI-MOD files",
                Context::default(),
            )]
        })?;
        let obo = OboOntology::from_raw(reader).map_err(|e| {
            vec![
                BoxedError::small(
                    CVError::FileCouldNotBeParsed,
                    e.get_short_description(),
                    e.get_long_description(),
                )
                .add_contexts(e.get_contexts().iter().cloned()),
            ]
        })?;
        let mut mods: Vec<SimpleModification> = Vec::new();
        let mut errors = Vec::new();

        'stanza: for obj in &obo.objects {
            if obj.stanza_type != OboStanzaType::Term
                || obj.id == (Some("MOD".into()), "00000".into())
                || obj.id == (Some("MOD".into()), "00004".into())
                || obj.id == (Some("MOD".into()), "00008".into())
            {
                continue;
            }
            let Ok(id) = obj.id.1.parse() else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        CVError::ItemError,
                        "Invalid ID",
                        "A PSI-MOD ID should be numerical",
                        Context::show(obj.id.1.to_string()),
                    ),
                );
                continue;
            };
            let mut modification = OntologyModification {
                id,
                name: obj.lines["name"][0].0.clone(),
                ontology: Ontology::Psimod,
                obsolete: obj.obsolete,
                ..OntologyModification::default()
            };
            if let Some((description, cross_ids, _, _)) = &obj.definition {
                modification.description = description.clone();
                modification.cross_ids = cross_ids.clone().into();
            }
            for synonym in &obj.synonyms {
                modification
                    .synonyms
                    .push((synonym.scope, synonym.synonym.clone()));
            }

            let mut origins = Vec::new();
            let mut term = None;
            for (id, _values, _comment) in &obj.xref {
                match (id.0.as_deref(), &id.1) {
                    (Some("DiffFormula"), s) if s.as_ref() != "\"none\"" => {
                        match MolecularFormula::psi_mod(
                            s.trim_start_matches('\"').trim_end_matches('\"'),
                        ) {
                            Ok(formula) => modification.formula = formula,
                            Err(err) => combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Invalid formula",
                                    "The formula was invalid",
                                    Context::show(obj.id.to_string()),
                                )
                                .add_underlying_error(
                                    err.to_owned()
                                        .convert::<CVError, BoxedError<'static, CVError>>(|_| {
                                            CVError::ItemError
                                        }),
                                ),
                            ),
                        }
                    }
                    (Some("Origin"), value) => {
                        let commas = value.chars().filter(|c| *c == ',').count();
                        let cross_link_numbers = obj.lines.get("comment").map_or(Vec::new(), |v| {
                            v[0].0
                                .split(';')
                                .flat_map(|s| s.split('.'))
                                .filter_map(|tag| {
                                    tag.trim()
                                        .to_ascii_lowercase()
                                        .strip_prefix("cross-link")
                                        .map(|v| {
                                            v.trim()
                                                .trim_end_matches('.')
                                                .parse::<usize>()
                                                .map_err(|e| {
                                                    format!("Could not parse '{tag}': {e}")
                                                })
                                                .unwrap()
                                        })
                                })
                                .collect()
                        });
                        if cross_link_numbers.len() > 1 {
                            combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Multiple cross-link numbers",
                                    format!(
                                        "This modification was defined to be a cross-linker with the following number of sites: {}",
                                        cross_link_numbers.iter().join(", ")
                                    ),
                                    Context::show(obj.id.to_string()),
                                ),
                            );
                            continue 'stanza;
                        }

                        if let Some(amount) = cross_link_numbers.first() {
                            if commas + 1 != *amount {
                                combine_error(
                                    &mut errors,
                                    BoxedError::new(
                                        CVError::ItemError,
                                        "Different cross-link numbers",
                                        format!(
                                            "This modification was defined be a cross-linker with {amount} sites, but has {} orgins",
                                            commas + 1
                                        ),
                                        Context::show(obj.id.to_string()),
                                    ),
                                );
                                continue 'stanza;
                            }
                        } else if value.contains(',') {
                            combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Missing cross-link number",
                                    format!(
                                        "This modification has {} orgins, but is not defined to be a cross-linker in the comment",
                                        commas + 1
                                    ),
                                    Context::show(obj.id.to_string()),
                                ),
                            );
                            continue 'stanza;
                        }
                        if value.as_ref() != "\"none\"" {
                            origins = value
                                .trim_start_matches('\"')
                                .trim_end_matches('\"')
                                .split(',')
                                .map(|s| s.trim().to_string())
                                .collect();
                        }
                    }
                    (Some("TermSpec"), value) => {
                        let v = value.trim_start_matches('\"').trim_end_matches('\"');
                        if v == "N-term" {
                            term = Some(Position::AnyNTerm);
                        } else if v == "C-term" {
                            term = Some(Position::AnyCTerm);
                        } else if v == "none" {
                            term = Some(Position::Anywhere);
                        } else {
                            combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Invalid TermSpec",
                                    "The termSpec should be 'N-term', 'C-term', or 'none'",
                                    Context::show(format!("{}: '{v}'", obj.id)),
                                ),
                            );
                            continue 'stanza;
                        }
                    }
                    _ => (), // ignore
                }
            }
            if origins.len() <= 1 {
                let mut rules = Vec::new();
                if let Some(origin) = origins.first() {
                    match parse_rule(origin, term) {
                        Ok(rule) => rules.push((vec![rule], Vec::new(), Vec::new())),
                        Err(err) => {
                            combine_error(&mut errors, err);
                            continue;
                        }
                    }
                } else if let Some(term) = term {
                    rules.push((vec![PlacementRule::Terminal(term)], Vec::new(), Vec::new()));
                }
                modification.data = ModData::Mod {
                    specificities: rules,
                };
            } else if origins.len() == 2 {
                let left = match parse_rule(&origins[0], term) {
                    Ok(rule) => rule,
                    Err(err) => {
                        combine_error(&mut errors, err);
                        continue;
                    }
                };
                let right = match parse_rule(&origins[1], term) {
                    Ok(rule) => rule,
                    Err(err) => {
                        combine_error(&mut errors, err);
                        continue;
                    }
                };
                modification.data = ModData::Linker {
                    length: crate::sequence::LinkerLength::Unknown,
                    specificities: vec![LinkerSpecificity::Asymmetric {
                        rules: (vec![left], vec![right]),
                        stubs: Vec::new(),
                        neutral_losses: Vec::new(),
                        diagnostic: Vec::new(),
                    }],
                };
            } else {
                combine_error(
                    &mut errors,
                    BoxedError::new(
                        CVError::ItemWarning,
                        "Higher order cross-linker",
                        format!(
                            "This cross-linker links {} sites, but only cross-linkers with two sites are supported for now, this modification will be ignored",
                            origins.len()
                        ),
                        Context::show(obj.id.to_string()),
                    ),
                );
            }
            mods.push(modification.into());
        }

        Ok((obo.version(), mods, errors))
    }
}

/// Parse a Origin rule, either a single amino acid, or a PSI-MOD ID
/// # Errors
/// If this is not a valid rule
fn parse_rule(
    rule: &str,
    term: Option<Position>,
) -> Result<PlacementRule, BoxedError<'static, CVError>> {
    if rule.len() == 1 {
        if rule == "X" {
            Ok(term
                .filter(|t| *t != Position::Anywhere)
                .map_or_else(|| PlacementRule::Anywhere, PlacementRule::Terminal))
        } else if let Ok(aa) = rule.try_into() {
            Ok(PlacementRule::AminoAcid(
                vec![aa].into(),
                term.unwrap_or(Position::Anywhere),
            ))
        } else {
            Err(BoxedError::new(
                CVError::ItemError,
                "Invalid amino acid",
                "A single letter origin is assumed to be an amino acid but this is not",
                Context::show(rule.to_string()),
            ))
        }
    } else {
        let id = OboIdentifier::from_str(rule).unwrap_infallible();

        if id.0.is_none_or(|ns| ns.as_ref() != "MOD") {
            Err(BoxedError::new(
                CVError::ItemError,
                "Invalid ID",
                "A modification as an origin should come from within PSI-MOD",
                Context::show(rule.to_string()),
            ))
        } else if let Ok(id) = id.1.parse() {
            Ok(PlacementRule::PsiModification(
                id,
                term.unwrap_or(Position::Anywhere),
            ))
        } else {
            Err(BoxedError::new(
                CVError::ItemError,
                "Invalid ID",
                "The ID should be numeric",
                Context::show(rule.to_string()),
            ))
        }
    }
}
