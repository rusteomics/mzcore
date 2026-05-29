//! Update the ontologies used in mzcore

use context_error::{self as _, BasicKind, BoxedError, CreateError};
use mzcore::{
    chemistry::Chemical,
    ontology::*,
    sequence::{CrossId, PlacementRule, SimpleModificationInner},
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
    if args.is_empty() || args.contains(&"resid".to_string()) {
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
        } else {
            println!(
                "RESID is ignored, to update download the file to `mzcore-update/data/RESID.xml` (because FTP is not supported)"
            );
        }
    }
    if args.is_empty()
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
    }
    let unimod = if args.is_empty() || args.contains(&"unimod".to_string()) {
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
        index
            .save_to_cache_at(std::path::Path::new("mzcore/src/databases/unimod.dat"))
            .unwrap();
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
        && let Some(unimod) = unimod
    {
        validate_linking(&unimod, &psimod);
    }
}

fn validate_linking(unimod: &CVIndex<Unimod>, psimod: &CVIndex<PsiMod>) {
    let mut errors = Vec::new();
    for m in psimod.data() {
        if let Some(id) = m.description() {
            // let mut unimod_link = false;
            for cross_id in &id.cross_ids {
                if let CrossId::Mod(Ontology::Unimod, code, rule) = cross_id {
                    // unimod_link = true;
                    if let Some(unimod_m) = unimod.get_by_index(&code) {
                        if m.formula() != unimod_m.formula() {
                            context_error::combine_error(
                                &mut errors,
                                BoxedError::new(
                                    BasicKind::Warning,
                                    "Suspicious link",
                                    "Unimod and PSI-MOD do not have the same formula",
                                    context_error::Context::default().lines(
                                        0,
                                        format!(
                                            "MOD:{}={} UNIMOD:{code}={}",
                                            id.id(),
                                            m.formula(),
                                            unimod_m.formula()
                                        ),
                                    ),
                                ),
                            );
                        }
                        if *rule != PlacementRule::Anywhere {
                            match &*unimod_m {
                                SimpleModificationInner::Database { specificities, .. } => {
                                    if !specificities
                                        .iter()
                                        .any(|s| s.0.iter().any(|r| rule.is_subset(r)))
                                    {
                                        context_error::combine_error(
                                            &mut errors,
                                            BoxedError::new(
                                                BasicKind::Warning,
                                                "Suspicious link",
                                                "Unimod does not allow this modification on this location",
                                                context_error::Context::default().lines(
                                                    0,
                                                    format!(
                                                        "MOD:{} UNIMOD:{code} {rule:?}",
                                                        id.id()
                                                    ),
                                                ),
                                            ),
                                        );
                                    }
                                }
                                _ => (),
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
                                    .lines(0, format!("MOD:{} {code}", id.id())),
                            ),
                        );
                    }
                }
            }

            // TODO: check at some point but be sure to first inherit any unimod links from parents and maybe also propose a link?
            // if !unimod_link {
            //     errors
            //         .entry("Not linked to Unimod")
            //         .or_default()
            //         .push((id.id(), String::new()))
            // }
        }
    }

    for error in errors {
        println!("{error}");
    }
}
