//! Test all present mzSpecLib files
#![allow(clippy::missing_panics_doc)]
use std::collections::{HashMap, HashSet};

use context_error as _;
use indexmap as _;
use itertools as _;
use mzannotate::mzspeclib::{MzSpecLibTextParser, MzSpecLibTextWriter};
use ordered_float as _;
use serde as _;
use serde_json as _;
use thin_vec as _;
use uom as _;

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
            let mut parsed_spectra = Vec::new();
            files += 1;
            let spectra = MzSpecLibTextParser::open(
                std::io::BufReader::new(std::fs::File::open(entry.path()).unwrap()),
                Some(entry.path()),
                &mzcore::ontology::STATIC_ONTOLOGIES,
            )
            .unwrap();
            let header = spectra.header().clone();
            for spectrum in spectra {
                match spectrum {
                    Ok(spectrum) => {
                        spectrum_unused_attributes.extend(
                            spectrum
                                .attributes
                                .iter()
                                .flatten()
                                .filter(|a| a.name != mzcv::term!(MS:1003254|peak attribute))
                                .map(|a| a.name.clone()),
                        );
                        interpretation_unused_attributes.extend(
                            spectrum
                                .interpretations
                                .iter()
                                .flat_map(|a| a.attributes.iter())
                                .flatten()
                                .filter(|a| {
                                    a.name != mzcv::term!(MS:1002252|Comet:xcorr)
                                        && a.name != mzcv::term!(MS:1002354|PSM-level q-value)
                                })
                                .map(|a| a.name.clone()),
                        );
                        parsed_spectra.push(spectrum);
                    }
                    Err(e) => {
                        errors.entry(entry.path()).or_default().push(e);
                    }
                }
            }

            let rewrite_path = entry.path().with_extension("txt.out");
            let mut writer = MzSpecLibTextWriter::new(std::io::BufWriter::new(
                std::fs::File::create(&rewrite_path).unwrap(),
            ));
            *writer.header_mut() = header;
            let mut writer = writer.write_header().unwrap();
            writer.write_spectra(&parsed_spectra).unwrap();
            drop(writer);

            let reparsed_spectra = MzSpecLibTextParser::open(
                std::io::BufReader::new(std::fs::File::open(&rewrite_path).unwrap()),
                Some(rewrite_path),
                &mzcore::ontology::STATIC_ONTOLOGIES,
            )
            .unwrap()
            .collect::<Result<Vec<_>, _>>()
            .unwrap();

            assert_eq!(parsed_spectra.len(), reparsed_spectra.len());

            for (p, r) in parsed_spectra.iter().zip(reparsed_spectra.iter()) {
                assert_eq!(p.peaks, r.peaks);
                // The other information fields have too many issues with sort order and equivalence
                // that equality is not a good test
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
    // assert!(
    //     spectrum_unused_attributes.is_empty(),
    //     "Some spectrum attributes unused"
    // );
    // assert!(
    //     interpretation_unused_attributes.is_empty(),
    //     "Some interpretation attributes unused"
    // );
}

#[test]
fn a_second_and_a_minute_declaration_of_one_instant_agree() {
    /// The same instant, declared twice: 1719.50946 s is 28.658491 min.
    const LIBRARY: &str = "<mzSpecLib>
MS:1003186|library format version=1.0
<Spectrum=1>
[1]MS:1000894|retention time=1719.509460
[1]UO:0000000|unit=UO:0000010|second
<Peaks>

<Spectrum=2>
[1]MS:1000894|retention time=28.658491
[1]UO:0000000|unit=UO:0000031|minute
<Peaks>
";

    let spectra: Vec<_> = MzSpecLibTextParser::open(
        LIBRARY.as_bytes(),
        None,
        &mzcore::ontology::STATIC_ONTOLOGIES,
    )
    .unwrap()
    .map(|s| s.unwrap())
    .collect();

    let times: Vec<f64> = spectra
        .iter()
        .map(|x| x.description.acquisition.scans[0].start_time)
        .collect();

    assert_eq!(times[0], times[1], "1719.50946 s is 28.658491 min");
    // mzdata documents start_time as minutes (mzdata::/src/spectrum/scan_properties.rs:138). The
    // tolerance is for the to_f32 the parser applies on the way in.
    assert!((times[0] - 28.658_491).abs() < 1e-4, "got {}", times[0]);
}
