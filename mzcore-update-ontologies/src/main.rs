//! Update the ontologies used in mzcore
use mzcore::ontology::PsiMod;
use mzcv::CVIndex;

fn main() {
    let mut index = CVIndex::<PsiMod>::empty();
    index.update_from_url(&[None]).unwrap();
    println!(
        "PSI-MOD version: {}, last updated: {}, modifications: {}",
        index.version().version.as_deref().unwrap_or("-"),
        index.version().last_updated().as_deref().unwrap_or("-"),
        index.data().len()
    );
    index
        .save_to_cache_at(std::path::Path::new("mzcore/src/databases/psimod_new.dat"))
        .unwrap();
}
