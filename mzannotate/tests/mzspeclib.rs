//! Test all present mzspeclib files
use std::collections::{HashMap, HashSet};

use mzannotate::mzspeclib::MzSpecLibParser;

#[test]
fn read_all_files() {
    let mut files = 0;
    let mut errors: HashMap<std::path::PathBuf, Vec<_>> = HashMap::new();
    let mut spectrum_unused_attributes = HashSet::new();
    let mut interpretation_unused_attributes = HashSet::new();

    for entry in std::fs::read_dir("../data").unwrap().flatten() {
        if entry
            .path()
            .to_string_lossy()
            .to_ascii_lowercase()
            .ends_with(".mzspeclib.txt")
        {
            files += 1;
            let spectra = MzSpecLibParser::new(std::io::BufReader::new(
                std::fs::File::open(entry.path()).unwrap(),
            ))
            .unwrap();
            for spectrum in spectra {
                match spectrum {
                    Ok(spectrum) => {
                        spectrum_unused_attributes.extend(
                            spectrum
                                .attributes
                                .iter()
                                .filter(|a| a.name != mzannotate::term!(MS:1003254|peak attribute))
                                .map(|a| a.name.clone()),
                        );
                        interpretation_unused_attributes.extend(
                            spectrum
                                .interpretations
                                .iter()
                                .flat_map(|a| a.attributes.iter())
                                .map(|a| a.name.clone()),
                        );
                    }
                    Err(e) => {
                        errors.entry(entry.path()).or_default().push(e);
                    }
                }
            }
        }
    }

    for (path, errors) in &errors {
        println!("[{}] errors: {}", path.display(), errors.len());
        for error in errors {
            println!("{error:?}");
        }
        println!();
    }

    println!("Unused spectrum attributes:");
    for attribute in &spectrum_unused_attributes {
        println!("{attribute}");
    }
    println!();

    println!("Unused interpretation attributes:");
    for attribute in &interpretation_unused_attributes {
        println!("{attribute}");
    }
    println!();

    let total_errors = errors.values().map(Vec::len).sum::<usize>();

    println!(
        "files: {files}, files with errors: {}, total errors: {total_errors}, ignored spectrum attributes: {}, ignored interpretation attributes: {}",
        errors.len(),
        spectrum_unused_attributes.len(),
        interpretation_unused_attributes.len(),
    );

    assert!(total_errors == 0, "some errors encountered");
    assert!(
        spectrum_unused_attributes.is_empty(),
        "Some spectrum attributes unused"
    );
    assert!(
        interpretation_unused_attributes.is_empty(),
        "Some interpretation attributes unused"
    );
}
