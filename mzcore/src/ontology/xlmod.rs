//! Code to handle the XL-MOD ontology
use context_error::{
    BoxedError, Context, CreateError, FullErrorContent, StaticErrorContent, combine_error,
    combine_errors,
};
use itertools::Itertools;
use std::{collections::HashMap, sync::Arc};
use thin_vec::ThinVec;

use mzcv::{
    CVError, CVFile, CVSource, CVVersion, HashBufReader, OboIdentifier, OboOntology, OboStanzaType,
    OboValue,
};

use crate::{
    chemistry::{DiagnosticIon, MolecularFormula, NeutralLoss},
    ontology::{
        Ontology,
        ontology_modification::{ModData, OntologyModification},
    },
    sequence::{
        AminoAcid, LinkerLength, LinkerSpecificity, PlacementRule, Position, SimpleModification,
        SimpleModificationInner,
    },
};

/// XL-MOD modifications
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct XlMod {}

impl CVSource for XlMod {
    type Data = SimpleModificationInner;
    type Structure = Vec<SimpleModification>;
    fn cv_name() -> &'static str {
        "XLMOD"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "XLMOD",
            extension: "obo",
            url: Some(
                "https://raw.githubusercontent.com/HUPO-PSI/xlmod-CV/refs/heads/main/XLMOD.obo",
            ),
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("../databases/xlmod.dat"),
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
        let reader = reader.next().unwrap();
        let mut errors = Vec::new();
        OboOntology::from_raw(reader)
            .map_err(|e| {
                vec![
                    BoxedError::small(
                        CVError::FileCouldNotBeParsed,
                        e.get_short_description(),
                        e.get_long_description(),
                    )
                    .add_contexts(e.get_contexts().iter().cloned()),
                ]
            })
            .map(|obo| {
                (
                    obo.version(),
                    {
                        let mut mods: Vec<SimpleModification> = Vec::new();

                        for obj in &obo.objects {
                            if obj.stanza_type != OboStanzaType::Term {
                                continue;
                            }
                            let id: u32 = obj
                                .id
                                .1
                                .parse()
                                .expect("Incorrect XLMOD id, should be numerical");
                            let name = obj.lines["name"][0].0.clone();

                            let description = obj
                                .definition
                                .as_ref()
                                .map_or_else(Box::default, |d| d.0.clone());
                            let cross_ids = obj
                                .definition
                                .as_ref()
                                .map_or_else(ThinVec::new, |d| d.1.clone().into());
                            let synonyms = obj
                                .synonyms
                                .iter()
                                .map(|s| (s.scope, s.synonym.clone()))
                                .collect();

                            // Get all properties from any ancestors in the tree then get all properties from this definition
                            let mut properties = Properties::default();
                            let mut stack = Vec::new();
                            stack.extend(obj.is_a.clone());

                            while let Some(id) = stack.pop() {
                                if let Some(obj) = obo.objects.iter().find(|o| o.id == id) {
                                    combine_errors(
                                        &mut errors,
                                        parse_property_values(
                                            &obj.property_values,
                                            &mut properties,
                                            id,
                                        ),
                                    );
                                    stack.extend(obj.is_a.clone());
                                }
                            }

                            combine_errors(
                                &mut errors,
                                parse_property_values(
                                    &obj.property_values,
                                    &mut properties,
                                    obj.id.clone(),
                                ),
                            );

                            if properties.origins.0.is_empty() {
                                properties.origins.0 = vec![PlacementRule::Anywhere];
                            }

                            if properties.sites == Some(2) {
                                mods.push(
                                    OntologyModification {
                                        formula: properties.formula.unwrap_or_default(),
                                        name,
                                        description,
                                        cross_ids,
                                        synonyms,
                                        id,
                                        ontology: Ontology::Xlmod,
                                        obsolete: obj.obsolete,
                                        data: ModData::Linker {
                                            length: properties.length,
                                            specificities: vec![
                                                if properties.origins.1.is_empty() {
                                                    LinkerSpecificity::Symmetric {
                                                        rules: properties.origins.0,
                                                        stubs: properties.stubs,
                                                        neutral_losses: properties.neutral_losses,
                                                        diagnostic: properties.diagnostic_ions,
                                                    }
                                                } else {
                                                    LinkerSpecificity::Asymmetric {
                                                        rules: (
                                                            properties.origins.0,
                                                            properties.origins.1,
                                                        ),
                                                        stubs: properties.stubs,
                                                        neutral_losses: properties.neutral_losses,
                                                        diagnostic: properties.diagnostic_ions,
                                                    }
                                                },
                                            ],
                                        },
                                    }
                                    .into(),
                                );
                            } else if properties.sites == Some(1) {
                                mods.push(Arc::new(
                                    OntologyModification {
                                        formula: properties.formula.unwrap_or_default(),
                                        name,
                                        description,
                                        cross_ids,
                                        synonyms,
                                        ontology: Ontology::Xlmod,
                                        obsolete: obj.obsolete,
                                        id,
                                        data: ModData::Mod {
                                            specificities: vec![(
                                                properties.origins.0,
                                                properties.neutral_losses,
                                                properties.diagnostic_ions,
                                            )],
                                        },
                                    }
                                    .into(),
                                ));
                            } else if !properties.dna_linker && !obj.property_values.is_empty() {
                                if let Some(sites) = properties.sites {
                                    combine_error(&mut errors, BoxedError::new(
                                        CVError::ItemWarning,
                                        "Higher order cross-linker", 
                                        format!("This cross-linker links {sites} sites, but only cross-linkers with two sites are supported for now, this modification will be ignored"), 
                                        Context::default().lines(0, format!("XLMOD:{id:05}"))));
                                } else {
                                    combine_error(&mut errors, BoxedError::new(
                                        CVError::ItemWarning,
                                        "Undefined modification", 
                                        "This modification has no definition of the number of sites it links, this modification will be ignored", 
                                        Context::default().lines(0, format!("XLMOD:{id:05}"))));
                                }
                            }
                        }

                        mods
                    },
                    errors,
                )
            })
    }
}

#[derive(Default)]
struct Properties {
    dna_linker: bool,
    sites: Option<u8>,
    length: LinkerLength,
    formula: Option<MolecularFormula>,
    origins: (Vec<PlacementRule>, Vec<PlacementRule>),
    diagnostic_ions: Vec<DiagnosticIon>,
    neutral_losses: Vec<NeutralLoss>,
    stubs: Vec<(MolecularFormula, MolecularFormula)>,
}

fn parse_property_values(
    property_values: &HashMap<
        Box<str>,
        Vec<(OboValue, Vec<(Box<str>, Box<str>)>, Option<Box<str>>)>,
    >,
    properties: &mut Properties,
    id: OboIdentifier,
) -> Vec<BoxedError<'static, CVError>> {
    let mut mass = None;
    let mut errors = Vec::new();
    for (property, value) in property_values {
        match property.as_ref() {
            "reactionSites" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple reaction sites",
                            "More than 1 'reactionSites` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                properties.sites = if let OboValue::Integer(n) = value[0].0 {
                    u8::try_from(n).map_or_else(
                        |_| {
                            combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Out of range number of sites",
                                    "The number of sites can only be in range 0—255",
                                    Context::default()
                                        .lines(0, value[0].0.to_string())
                                        .to_owned(),
                                ),
                            );
                            None
                        },
                        Some,
                    )
                } else {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid type",
                            "Invalid item type",
                            Context::none().lines(
                                0,
                                format!(
                                    "{id}: reactionSites: type: {}, expected: integer",
                                    value[0].0.datatype()
                                ),
                            ),
                        ),
                    );
                    None
                }
            }
            "spacerLength" => {
                for (def, _, _) in value {
                    let length = if let OboValue::Float(n) = def {
                        *n
                    } else {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                CVError::ItemError,
                                "Invalid type",
                                "Invalid item type",
                                Context::none().lines(
                                    0,
                                    format!(
                                        "{id}: spacerLength: type: {}, expected: float",
                                        def.datatype()
                                    ),
                                ),
                            ),
                        );
                        0.0
                    };
                    match &mut properties.length {
                        LinkerLength::Discreet(options) => {
                            options.push(length.into());
                        }
                        l => *l = LinkerLength::Discreet(vec![length.into()]),
                    }
                }
            }
            "minSpacerLength" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple minSpacerLength",
                            "More than 1 'minSpacerLength` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                let length = if let OboValue::Float(n) = value[0].0 {
                    n
                } else {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid type",
                            "Invalid item type",
                            Context::none().lines(
                                0,
                                format!(
                                    "{id}: minSpacerLength: type: {}, expected: float",
                                    value[0].0.datatype()
                                ),
                            ),
                        ),
                    );
                    0.0
                };
                match &mut properties.length {
                    LinkerLength::InclusiveRange(start, _) => {
                        *start = length.into();
                    }
                    l => {
                        *l = LinkerLength::InclusiveRange(length.into(), length.into());
                    }
                }
            }
            "maxSpacerLength" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple maxSpacerLength",
                            "More than 1 'maxSpacerLength` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                let length = if let OboValue::Float(n) = value[0].0 {
                    n
                } else {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid type",
                            "Invalid item type",
                            Context::none().lines(
                                0,
                                format!(
                                    "{id}: maxSpacerLength: type: {}, expected: float",
                                    value[0].0.datatype()
                                ),
                            ),
                        ),
                    );
                    0.0
                };
                match &mut properties.length {
                    LinkerLength::InclusiveRange(_, end) => {
                        *end = length.into();
                    }
                    l => {
                        *l = LinkerLength::InclusiveRange(length.into(), length.into());
                    }
                }
            }
            "monoIsotopicMass" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple monoIsotopicMass",
                            "More than 1 'monoIsotopicMass` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                mass = if let OboValue::Float(n) = value[0].0 {
                    Some(ordered_float::OrderedFloat(n))
                } else {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemError,
                            "Invalid type",
                            "Invalid item type",
                            Context::none().lines(
                                0,
                                format!(
                                    "{id}: monoIsotopicMass: type: {}, expected: float",
                                    value[0].0.datatype()
                                ),
                            ),
                        ),
                    );
                    None
                }
            }
            "deadEndFormula" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple deadEndFormula",
                            "More than 1 'deadEndFormula` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                properties.sites = Some(1);
                properties.formula = MolecularFormula::xlmod(&value[0].0.to_string())
                    .map_err(|e| {
                        combine_error(&mut errors, e.to_owned().convert(|_| CVError::ItemError));
                    })
                    .ok();
            }
            "neutralLossFormula" => {
                for (def, _, _) in value {
                    properties
                        .neutral_losses
                        .push(MolecularFormula::xlmod(&def.to_string()).unwrap().into());
                }
            }
            "bridgeFormula" => {
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple bridgeFormula",
                            "More than 1 'bridgeFormula` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                properties.sites = properties.sites.or(Some(2));
                properties.formula =
                    Some(MolecularFormula::xlmod(&value[0].0.to_string()).unwrap());
            }
            "specificities" => {
                // specificities: "(C,U)" xsd:string
                // specificities: "(K,N,Q,R,Protein N-term)&(E,D,Protein C-term)" xsd:string
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple specificities",
                            "More than 1 'specificities` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                if let Some((l, r)) = value[0].0.to_string().split_once('&') {
                    properties.sites = properties.sites.or(Some(2));
                    properties
                        .origins
                        .0
                        .extend(l.trim_matches(['(', ')']).split(',').filter_map(|s| {
                            read_placement_rule(s.trim())
                                .map_err(|err| combine_error(&mut errors, err))
                                .ok()
                        }));
                    properties
                        .origins
                        .1
                        .extend(r.trim_matches(['(', ')']).split(',').filter_map(|s| {
                            read_placement_rule(s.trim())
                                .map_err(|err| combine_error(&mut errors, err))
                                .ok()
                        }));
                } else {
                    properties.origins.0.extend(
                        value[0]
                            .0
                            .to_string()
                            .trim_matches(['(', ')'])
                            .split(',')
                            .filter_map(|s| {
                                read_placement_rule(s.trim())
                                    .map_err(|err| combine_error(&mut errors, err))
                                    .ok()
                            }),
                    );
                }
            }
            "secondarySpecificities" => {
                // TODO: keep track that these are 'secondary' somewhere, similarly to the Unimod hidden state
                // secondarySpecificities: "(S,T,Y)" xsd:string
                if value.len() > 1 {
                    combine_error(
                        &mut errors,
                        BoxedError::new(
                            CVError::ItemWarning,
                            "Multiple secondarySpecificities",
                            "More than 1 'secondarySpecificities` definitions for this entry",
                            Context::none().lines(0, id.to_string()),
                        ),
                    );
                }
                properties.origins.0.extend(
                    value[0]
                        .0
                        .to_string()
                        .trim_matches(['(', ')'])
                        .split(',')
                        .filter_map(|s| {
                            read_placement_rule(s.trim())
                                .map_err(|err| combine_error(&mut errors, err))
                                .ok()
                        }),
                );
            }
            "baseSpecificities" | "secondarybaseSpecificities" => {
                properties.dna_linker = true;
            }
            "reporterMass" | "CID_Fragment" => {
                // reporterMass: "555.2481" xsd:double
                // CID_Fragment: "828.5" xsd:double
                for (def, _, _) in value {
                    properties.diagnostic_ions.push(DiagnosticIon(
                        MolecularFormula::with_additional_mass(if let OboValue::Float(n) = def {
                            *n
                        } else {
                            combine_error(
                                &mut errors,
                                BoxedError::new(
                                    CVError::ItemError,
                                    "Invalid type",
                                    "Invalid item type",
                                    Context::none().lines(
                                        0,
                                        format!(
                                            "{id}: {property}: type: {}, expected: float",
                                            value[0].0.datatype()
                                        ),
                                    ),
                                ),
                            );
                            0.0
                        }),
                    ));
                }
            }
            "reporterFormula" => {
                for (def, _, _) in value {
                    properties.diagnostic_ions.push(DiagnosticIon(
                        MolecularFormula::xlmod(&def.to_string()).unwrap(),
                    ));
                }
            }
            "stubDefinition" | "stubFormula" => {
                // CID = H4 C3 O2 S1, -H2 -O1 : H2 C3 O1
                // ETD = -H1 :
                // TODO: Extend the stub logic to handle the techniques and neutral losses
                for (def, _, _) in value {
                    if let OboValue::String(definition) = &def {
                        let (techniques, definition) = definition.split_once('=').unwrap();
                        let (first, second) = definition.split_once(':').unwrap();
                        let mut split_first = first.split(',');
                        let mut split_second = second.split(',');
                        let formula_first =
                            split_first.next().map_or(MolecularFormula::default(), |v| {
                                MolecularFormula::xlmod(v).unwrap()
                            });
                        let losses_first: Vec<NeutralLoss> = split_first
                            .map(|v| MolecularFormula::xlmod(v).unwrap().into())
                            .collect();
                        let formula_second = split_second
                            .next()
                            .map_or(MolecularFormula::default(), |v| {
                                MolecularFormula::xlmod(v).unwrap()
                            });
                        let losses_second: Vec<NeutralLoss> = split_second
                            .map(|v| MolecularFormula::xlmod(v).unwrap().into())
                            .collect();
                        properties.stubs.push((formula_first, formula_second));
                    } else {
                        combine_error(
                            &mut errors,
                            BoxedError::new(
                                CVError::ItemError,
                                "Invalid type",
                                "Invalid item type",
                                Context::none().lines(
                                    0,
                                    format!(
                                        "{id}: {property}: type: {}, expected: string",
                                        value[0].0.datatype()
                                    ),
                                ),
                            ),
                        );
                    }
                }
            }
            _ => {} // TODO: handle: hydrophilicPEGchain?, maxAbsorption?, waveLengthRange?
        }
    }

    // Ignore the mass if a formula is set
    properties.formula = properties
        .formula
        .clone()
        .or_else(|| mass.map(|v| MolecularFormula::with_additional_mass(v.0)));

    errors
}

/// Read a placement rule
/// # Errors
/// If not a valid placement rule
fn read_placement_rule(brick: &str) -> Result<PlacementRule, BoxedError<'static, CVError>> {
    if brick.len() == 1
        && let Ok(aa) = AminoAcid::try_from(brick)
    {
        Ok(PlacementRule::AminoAcid(
            vec![aa].into(),
            Position::Anywhere,
        ))
    } else if brick.eq_ignore_ascii_case("Protein N-term") {
        Ok(PlacementRule::Terminal(Position::ProteinNTerm))
    } else if brick.eq_ignore_ascii_case("Protein C-term") {
        Ok(PlacementRule::Terminal(Position::ProteinCTerm))
    } else {
        Err(BoxedError::new(
            CVError::ItemError,
            "Invalid placement rule",
            "Placement rule has to be an amino acid, protein N-term, or protein C-term.",
            Context::default().lines(0, brick).to_owned(),
        ))
    }
}
