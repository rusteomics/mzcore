#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter},
    sync::{atomic::AtomicUsize, Arc},
};

use clap::Parser;
use directories::ProjectDirs;
use fragment::FragmentType;
use itertools::Itertools;
use mzdata::io::{MZFileReader, SpectrumSource};
use rayon::prelude::*;
use rustyms::{
    identification::IdentifiedPeptideSource,
    model::MatchingParameters,
    spectrum::{Score, Scores},
    *,
};
use spectrum::{AnnotatedPeak, PeakSpectrum};

#[derive(Parser)]
struct Cli {
    /// The input csv file, should have the following columns: 'raw_file' (full path), 'scan_index', 'z', 'sequence', and can have 'mode' (etd/td_etd/ethcd/etcad/eacid/ead/hcd/cid/all/none, defaults to the global model)
    #[arg(short, long)]
    in_path: String,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: String,
    /// Global mode, will be overruled by line specific modes (etd/td_etd/ethcd/etcad/eacid/ead/hcd/cid/all/none)
    #[arg(long, default_value_t = String::from("all"))]
    mode: String,
    /// Turns on reporting of glycan Y-ions in a charge independent manner
    #[arg(long)]
    charge_independent_Y: bool,
    /// Turns on reporting of intensity statistics, the fraction of total TIC that could be annotated as well as the TIC for each spectrum
    #[arg(long)]
    report_intensity: bool,
    /// Turns on reporting of I/L coverage by satellite ions, returns a list with a 0 (not covered) or 1 (covered) for each I or L in the peptide
    #[arg(long)]
    report_IL_satellite_coverage: bool,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
}

fn select_model(text: &str, default: &'static FragmentationModel) -> &'static FragmentationModel {
    match text.to_ascii_lowercase().as_str() {
        "etd" => FragmentationModel::etd(),
        "td_etd" => FragmentationModel::td_etd(),
        "ethcd" | "etcad" => FragmentationModel::ethcd(),
        "eacid" => FragmentationModel::eacid(),
        "ead" => FragmentationModel::ead(),
        "hcd" | "cid" => FragmentationModel::cid_hcd(),
        "all" => FragmentationModel::all(),
        "none" => FragmentationModel::none(),
        _ => default,
    }
}

fn main() {
    let args = Cli::parse();
    let model = select_model(&args.mode, FragmentationModel::all());
    let parameters = MatchingParameters::default();
    let path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .unwrap()
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(serde_json::from_reader(BufReader::new(File::open(path).unwrap())).unwrap())
    };
    let files = rustyms::identification::BasicCSVData::parse_file(
        args.in_path,
        custom_database.as_ref(),
        true,
    )
    .expect("Invalid input file")
    .filter_map(|a| a.ok())
    .into_group_map_by(|l| l.raw_file.clone());
    let out_file = BufWriter::new(File::create(args.out_path).unwrap());
    let total_peptides = files.values().map(|f| f.len()).sum::<usize>();
    let peptides_counter = AtomicUsize::default();
    let raw_file_counter = AtomicUsize::default();
    println!("Raw files: 0/{}, Peptides: 0/{total_peptides}", files.len());

    let out_data: Vec<_> =  files.par_iter().flat_map(|(file_name, lines)| {
        let mut file = mzdata::io::MZReaderType::open_path(file_name).unwrap();

        let rows = lines
            .iter()
            .filter_map(|line| {
                let selected_model = line.mode.as_ref()
                    .map_or(model, |text| select_model(text, model));
                if let Some(spectrum) = file.get_spectrum_by_index(line.scan_index)
                {
                    let fragments = line.sequence.generate_theoretical_fragments(line.z, selected_model);
                    let annotated = spectrum.annotate(
                        line.sequence.clone(),
                        &fragments,
                        &parameters,
                        MassMode::Monoisotopic,
                    );
                    let scores: &Scores = &annotated
                        .scores(&fragments, &parameters, MassMode::Monoisotopic)
                        .1[0][0];

                    let mut row: BTreeMap<Arc<String>, String> = line.full_csv_line().unwrap_or(&[]).iter().cloned().collect();

                    if args.charge_independent_Y {
                        let unique_Y = fragments
                            .iter()
                            .filter_map(|fragment| {
                                if let FragmentType::Y(pos) = &fragment.ion {
                                    Some(pos)
                                } else {
                                    None
                                }
                            })
                            .unique()
                            .count();
                        let unique_Y_found = annotated
                            .spectrum()
                            .flat_map(|peak| &peak.annotation)
                            .filter_map(|fragment| {
                                if let FragmentType::Y(pos) = &fragment.ion {
                                    Some(pos)
                                } else {
                                    None
                                }
                            })
                            .unique()
                            .count();
                        row.insert(
                            Arc::new("ion_Y_charge_independent".to_string()),
                            format!("{}", (unique_Y_found as f64 / unique_Y as f64),),
                        );
                    }
                    if args.report_intensity {
                        row.insert(
                            Arc::new("intensity_combined".to_string()),
                            match scores.score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.fraction().to_string()
                                }
                            },
                        );
                        row.insert(
                            Arc::new("total_ion_current".to_string()),
                            match scores.score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.total.to_string()
                                }
                            },
                        );
                        for (ion, score) in &scores.ions {
                            row.insert(Arc::new(format!("intensity_{ion}")), match score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.fraction().to_string()
                                }
                            });
                        }
                    }
                    if args.report_IL_satellite_coverage {
                        row.insert(
                            Arc::new("IL_satellite_coverage".to_string()),
                            annotated.peptide.clone().singular_peptide().map_or(
                                String::new(),
                                |p| {
                                    p.sequence()
                                        .iter()
                                        .enumerate()
                                        .filter(|(_, s)| {
                                            s.aminoacid.aminoacid() == AminoAcid::Isoleucine
                                                || s.aminoacid.aminoacid() == AminoAcid::Leucine
                                        })
                                        .map(|(i, _)| {
                                            if annotated.spectrum().any(|p: &AnnotatedPeak| {
                                                p.annotation.iter().any(|a: &Fragment| {
                                                    matches!(a.ion, FragmentType::w(s, _, 0, _, _) | FragmentType::d(s, _, 0, _, _) if s.sequence_index
                                                    == SequencePosition::Index(i))
                                                })
                                            }) {
                                                '1'
                                            } else {
                                                '0'
                                            }
                                        })
                                        .collect()
                                },
                            ),
                        );
                    }
                    row.insert(
                        Arc::new("ion_combined".to_string()),
                        match scores.score {
                            Score::Position {
                                theoretical_positions,
                                ..
                            } => format!("{}", theoretical_positions.fraction(),),
                            Score::UniqueFormulas {
                                unique_formulas, ..
                            } => format!("{}", unique_formulas.fraction(),),
                        },
                    );

                    for (ion, score) in &scores.ions {
                        let recovered = match score {
                            Score::Position {
                                theoretical_positions,
                                ..
                            } => theoretical_positions,
                            Score::UniqueFormulas {
                                unique_formulas, ..
                            } => unique_formulas,
                        };
                        row.insert(Arc::new(format!("ion_{ion}")), format!("{}", recovered.fraction(),));
                    }
                    Some(row)
                } else {
                    eprintln!("Could not find scan index {} for file {}", line.scan_index, file_name.to_string_lossy());
                    None
                }
            })
            .collect::<Vec<_>>();
        println!(
            "Raw files: {}/{}, Peptides: {}/{total_peptides}",
            raw_file_counter.fetch_add(1, std::sync::atomic::Ordering::Relaxed) + 1,
            files.len(),
            peptides_counter.fetch_add(lines.len(), std::sync::atomic::Ordering::Relaxed) + lines.len(),
        );
        rows.into_iter().map(|row| row.into_iter().collect_vec()).collect_vec()
    }).collect();

    rustyms::csv::write_csv(
        out_file,
        out_data.into_iter().map(|r| {
            r.into_iter()
                .map(|(k, v)| (std::sync::Arc::<std::string::String>::unwrap_or_clone(k), v))
        }),
        ',',
    )
    .unwrap();
}
