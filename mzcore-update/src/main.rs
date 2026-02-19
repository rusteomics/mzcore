//! Update the ontologies used in mzcore
use mzcore::ontology::*;
use mzcv::CVIndex;

fn main() {
    let args: Vec<String> = std::env::args()
        .skip(1)
        .map(|v| v.to_ascii_lowercase())
        .collect();
    if args.is_empty()
        || args.contains(&"psimod".to_string())
        || args.contains(&"psi-mod".to_string())
    {
        // PSI-MOD
        let mut index = CVIndex::<PsiMod>::empty();
        let errs = index.update_from_url(&[]).unwrap();
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
    }
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
        let errs = index.update_from_url(&[]).unwrap();
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
    if args.is_empty() || args.contains(&"unimod".to_string()) {
        // Unimod
        let mut index = CVIndex::<Unimod>::empty();
        let errs = index.update_from_url(&[]).unwrap();
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
    }
    if args.is_empty() || args.contains(&"gno".to_string()) || args.contains(&"gnome".to_string()) {
        // GNOme
        let mut index = CVIndex::<Gnome>::empty();
        let errs = index.update_from_url(&[]).unwrap();
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
}
