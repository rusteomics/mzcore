//! Update the ontologies used in mzcore
use mzcore::ontology::*;
use mzcv::CVIndex;

fn main() {
    // PSI-MOD
    // let mut index = CVIndex::<PsiMod>::empty();
    // index.update_from_url(&[]).unwrap();
    // println!(
    //     "PSI-MOD version: {}, last updated: {}, modifications: {}",
    //     index.version().version.as_deref().unwrap_or("-"),
    //     index.version().last_updated().as_deref().unwrap_or("-"),
    //     index.data().len()
    // );
    // index
    //     .save_to_cache_at(std::path::Path::new("mzcore/src/databases/psimod.dat"))
    //     .unwrap();
    // // RESID
    // let path = std::path::Path::new("mzcore-update-ontologies/data/RESID.xml");
    // if path.exists() {
    //     let mut index = CVIndex::<Resid>::empty();
    //     index.update_from_path([Some(path)], false).unwrap();
    //     println!(
    //         "RESID version: {}, last updated: {}, modifications: {}",
    //         index.version().version.as_deref().unwrap_or("-"),
    //         index.version().last_updated().as_deref().unwrap_or("-"),
    //         index.data().len()
    //     );
    //     index
    //         .save_to_cache_at(std::path::Path::new("mzcore/src/databases/resid.dat"))
    //         .unwrap();
    // } else {
    //     println!(
    //         "RESID is ignored, to update download the file to `mzcore-update-ontologies/data/RESID.xml` (because FTP is not supported)"
    //     );
    // }
    // XLMOD
    let mut index = CVIndex::<XlMod>::empty();
    index
        .update_from_path([Some(std::path::Path::new("../xlmod-CV/XLMOD.obo"))], false)
        .unwrap();
    println!(
        "XLMOD version: {}, last updated: {}, modifications: {}",
        index.version().version.as_deref().unwrap_or("-"),
        index.version().last_updated().as_deref().unwrap_or("-"),
        index.data().len()
    );
    // index
    //     .save_to_cache_at(std::path::Path::new("mzcore/src/databases/xlmod.dat"))
    //     .unwrap();
    // Unimod
    // let mut index = CVIndex::<Unimod>::empty();
    // index.update_from_url(&[]).unwrap();
    // println!(
    //     "Unimod version: {}, last updated: {}, modifications: {}",
    //     index.version().version.as_deref().unwrap_or("-"),
    //     index.version().last_updated().as_deref().unwrap_or("-"),
    //     index.data().len()
    // );
    // index
    //     .save_to_cache_at(std::path::Path::new("mzcore/src/databases/unimod.dat"))
    //     .unwrap();
    // // GNOme
    // let mut index = CVIndex::<Gnome>::empty();
    // index.update_from_url(&[]).unwrap();
    // println!(
    //     "GNOme version: {}, last updated: {}, modifications: {}",
    //     index.version().version.as_deref().unwrap_or("-"),
    //     index.version().last_updated().as_deref().unwrap_or("-"),
    //     index.data().len()
    // );
    // index
    //     .save_to_cache_at(std::path::Path::new("mzcore/src/databases/gnome.dat"))
    //     .unwrap();
}
