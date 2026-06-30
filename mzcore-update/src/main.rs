//! Update the ontologies used in mzcore

use context_error::{self as _, BasicKind, BoxedError, CreateError};
use mzcore::{
    chemistry::Chemical,
    ontology::*,
    sequence::{CrossId, PlacementRule, Position, SimpleModificationInner},
};
use mzcv::CVIndex;

fn main() {
    let mut local = false;
    let mut link = false;
    let args: Vec<String> = std::env::args()
        .skip(1)
        .map(|v| v.to_ascii_lowercase())
        .filter(|v| {
            if v == "local" {
                local = true;
                false
            } else if v == "link" {
                link = true;
                false
            } else {
                true
            }
        })
        .collect();
    let psimod = if args.is_empty()
        || args.contains(&"psimod".to_string())
        || args.contains(&"psi-mod".to_string())
    {
        // PSI-MOD
        let mut index = CVIndex::<PsiMod>::empty();
        let errs = if local {
            index.update_from_path(None, false).unwrap()
        } else {
            index.update_from_url(&[]).unwrap()
        };
        for err in errs {
            println!("{err}");
        }
        println!(
            "PSI-MOD version: {}, last updated: {}, modifications: {}",
            index.version().version.as_deref().unwrap_or("-"),
            index.version().last_updated().as_deref().unwrap_or("-"),
            index.data().len()
        );
        index
            .save_to_cache_at(std::path::Path::new("mzcore/src/databases/psimod.dat"))
            .unwrap();
        Some(index)
    } else {
        None
    };
    let resid = if args.is_empty() || args.contains(&"resid".to_string()) {
        // RESID
        let path = std::path::Path::new("mzcore-update/data/RESID.xml");
        if path.exists() {
            let mut index = CVIndex::<Resid>::empty();
            let errs = index.update_from_path([Some(path)], false).unwrap();
            for err in errs {
                println!("{err}");
            }
            println!(
                "RESID version: {}, last updated: {}, modifications: {}",
                index.version().version.as_deref().unwrap_or("-"),
                index.version().last_updated().as_deref().unwrap_or("-"),
                index.data().len()
            );
            index
                .save_to_cache_at(std::path::Path::new("mzcore/src/databases/resid.dat"))
                .unwrap();
            Some(index)
        } else {
            println!(
                "RESID is ignored, to update download the file to `mzcore-update/data/RESID.xml` (because FTP is not supported)"
            );
            None
        }
    } else {
        None
    };
    let xlmod = if args.is_empty()
        || args.contains(&"xlmod".to_string())
        || args.contains(&"xl-mod".to_string())
    {
        // XLMOD
        let mut index = CVIndex::<XlMod>::empty();
        let errs = if local {
            index.update_from_path(None, false).unwrap()
        } else {
            index.update_from_url(&[]).unwrap()
        };
        for err in errs {
            println!("{err}");
        }
        println!(
            "XLMOD version: {}, last updated: {}, modifications: {}",
            index.version().version.as_deref().unwrap_or("-"),
            index.version().last_updated().as_deref().unwrap_or("-"),
            index.data().len()
        );
        index
            .save_to_cache_at(std::path::Path::new("mzcore/src/databases/xlmod.dat"))
            .unwrap();
        Some(index)
    } else {
        None
    };
    let mut unimod = if args.is_empty() || args.contains(&"unimod".to_string()) {
        // Unimod
        let mut index = CVIndex::<Unimod>::empty();
        let errs = if local {
            index.update_from_path(None, false).unwrap()
        } else {
            index.update_from_url(&[]).unwrap()
        };
        for err in errs {
            println!("{err}");
        }
        println!(
            "Unimod version: {}, last updated: {}, modifications: {}",
            index.version().version.as_deref().unwrap_or("-"),
            index.version().last_updated().as_deref().unwrap_or("-"),
            index.data().len()
        );
        Some(index)
    } else {
        None
    };
    if args.is_empty() || args.contains(&"gno".to_string()) || args.contains(&"gnome".to_string()) {
        // GNOme
        let mut index = CVIndex::<Gnome>::empty();
        let errs = if local {
            index.update_from_path(None, false).unwrap()
        } else {
            index.update_from_url(&[]).unwrap()
        };
        for err in errs {
            println!("{err}");
        }
        println!(
            "GNOme version: {}, last updated: {}, modifications: {}",
            index.version().version.as_deref().unwrap_or("-"),
            index.version().last_updated().as_deref().unwrap_or("-"),
            index.data().len()
        );
        index
            .save_to_cache_at(std::path::Path::new("mzcore/src/databases/gnome.dat"))
            .unwrap();
    }
    if link
        && let Some(psimod) = psimod
        && let Some(unimod) = &mut unimod
        && let Some(xlmod) = xlmod
        && let Some(resid) = resid
    {
        validate_linking(unimod, &psimod, &xlmod, &resid);
    }
    if let Some(index) = unimod {
        index
            .save_to_cache_at(std::path::Path::new("mzcore/src/databases/unimod.dat"))
            .unwrap();
    }
}

fn validate_linking(
    unimod: &mut CVIndex<Unimod>,
    psimod: &CVIndex<PsiMod>,
    xlmod: &CVIndex<XlMod>,
    resid: &CVIndex<Resid>,
) {
    let mut new_unimod = unimod.data().iter().map(|m| (**m).clone()).collect::<Vec<_>>();
    let mut errors = Vec::new();
    for m in psimod.data() {
        if let Some(psimod_id) = m.description() {
            // let mut unimod_link = false;
            for cross_id in &psimod_id.cross_ids {
                if let CrossId::Mod(Ontology::Unimod, code, rule) = cross_id {
                    // unimod_link = true;
                    if let Some(unimod_m) = new_unimod
                        .iter_mut()
                        .find(|m| m.description().is_some_and(|c| c.id() == *code))
                    {
                        if valid_link(m, unimod_m, rule, &mut errors) {
                            println!("Valid link to {unimod_m}");
                            if let SimpleModificationInner::Database { id, .. }
                            | SimpleModificationInner::Linker { id, .. }
                            | SimpleModificationInner::Gno { id, .. } = unimod_m
                                && !id.cross_ids.contains(&CrossId::Mod(
                                    Ontology::Psimod,
                                    id.id(),
                                    PlacementRule::Position(Position::Anywhere),
                                ))
                            {
                                id.cross_ids.push(CrossId::Mod(
                                    Ontology::Psimod,
                                    id.id(),
                                    PlacementRule::Position(Position::Anywhere),
                                ));
                            }
                        }
                    } else {
                        context_error::combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid Unimod code",
                                "Unimod modification does not exist",
                                context_error::Context::default()
                                    .lines(0, format!("MOD:{} {code}", psimod_id.id())),
                            ),
                        );
                    }
                } else if let CrossId::Mod(Ontology::Xlmod, code, rule) = cross_id {
                    if let Some(other) = xlmod.get_by_index(code) {
                        valid_link(m, &other, rule, &mut errors);
                    } else {
                        context_error::combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid XLMOD code",
                                "XLMOD modification does not exist",
                                context_error::Context::default()
                                    .lines(0, format!("MOD:{} {code}", psimod_id.id())),
                            ),
                        );
                    }
                } else if let CrossId::Mod(Ontology::Resid, code, rule) = cross_id {
                    if let Some(other) = resid.get_by_index(code) {
                        valid_link(m, &other, rule, &mut errors);
                    } else {
                        context_error::combine_error(
                            &mut errors,
                            BoxedError::new(
                                BasicKind::Error,
                                "Invalid RESID code",
                                "RESID modification does not exist",
                                context_error::Context::default()
                                    .lines(0, format!("MOD:{} {code}", psimod_id.id())),
                            ),
                        );
                    }
                }
            }

            // TODO: check at some point but be sure to first inherit any unimod links from parents
            // and maybe also propose a link? if !unimod_link {
            //     errors
            //         .entry("Not linked to Unimod")
            //         .or_default()
            //         .push((id.id(), String::new()))
            // }
        }
    }

    unimod
        .update(
            unimod.version().clone(),
            new_unimod.into_iter().map(std::sync::Arc::new),
        )
        .unwrap();

    for error in errors {
        println!("{error}");
    }
}

fn valid_link(
    one: &SimpleModificationInner,
    two: &SimpleModificationInner,
    rule: &PlacementRule,
    warnings: &mut Vec<BoxedError<'static, BasicKind>>,
) -> bool {
    let mut suspicious = false;
    let name_one = one.description().map_or_else(
        || one.to_string(),
        |o| format!("{}:{}", o.ontology.name(), o.id()),
    );
    let name_two = two.description().map_or_else(
        || two.to_string(),
        |o| format!("{}:{}", o.ontology.name(), o.id()),
    );
    if one.formula() != two.formula() {
        context_error::combine_error(
            warnings,
            BoxedError::new(
                BasicKind::Warning,
                "Suspicious link",
                format!(
                    "{} and {} do not have the same formula",
                    one.description().map_or(Ontology::Custom, |d| d.ontology),
                    two.description().map_or(Ontology::Custom, |d| d.ontology)
                ),
                context_error::Context::default().lines(
                    0,
                    format!("{name_one}={} {name_two}={}", one.formula(), two.formula()),
                ),
            ),
        );
    }
    if *rule != PlacementRule::Position(Position::Anywhere)
        && let SimpleModificationInner::Database { specificities, .. } = two
        && !specificities.iter().any(|s| s.0.iter().any(|r| rule.is_subset(r)))
    {
        context_error::combine_error(
            warnings,
            BoxedError::new(
                BasicKind::Warning,
                "Suspicious link",
                format!(
                    "{} does not allow this modification on this location",
                    one.description().map_or(Ontology::Custom, |d| d.ontology),
                ),
                context_error::Context::default()
                    .lines(0, format!("{name_one} {name_two} {rule:?}")),
            ),
        );
        suspicious = true;
    }
    !suspicious
}
