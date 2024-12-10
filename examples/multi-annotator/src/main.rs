#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter},
};

use clap::Parser;
use directories::ProjectDirs;
use fragment::FragmentType;
use itertools::Itertools;
use rayon::prelude::*;
use rustyms::{
    spectrum::{Score, Scores},
    system::{e, usize::Charge, Mass},
    *,
};
use spectrum::{AnnotatedPeak, PeakSpectrum};

#[derive(Parser)]
struct Cli {
    /// The input csv file, should have the following columns: 'mgf_path', 'scan_number', 'z', 'sequence', and can have 'fragmentation' (etd/td_etd/ethcd/etcad/hot eacid/eacid/ead/hcd/cid/all/none, defaults to the global model)
    #[arg(short, long)]
    in_path: String,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: String,
    /// The tolerance for matching fragments, use `<x>ppm` or `<x>da` to control the unit, e.g. `10.0ppm` or `2.3da`
    #[arg(short, long, default_value_t = Tolerance::new_ppm(20.0), value_parser=mass_tolerance_parse)]
    pub tolerance: Tolerance<Mass>,
    /// Global model, will be overruled by line specific models (etd/td_etd/ethcd/etcad/hot eacid/eacid/ead/hcd/cid/all/none)
    #[arg(long, default_value_t = String::from("all"))]
    model: String,
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

fn mass_tolerance_parse(input: &str) -> Result<Tolerance<Mass>, &'static str> {
    input.parse().map_err(|()| "Invalid tolerance parameter")
}

fn select_model(text: &str, default: &Model) -> Model {
    match text.to_ascii_lowercase().as_str() {
        "etd" => Model::etd(),
        "td_etd" => Model::td_etd(),
        "ethcd" | "etcad" => Model::ethcd(),
        "hot eacid" | "eacid" => Model::hot_eacid(),
        "ead" => Model::ead(),
        "hcd" | "cid" => Model::cid_hcd(),
        "all" => Model::all(),
        "none" => Model::none(),
        _ => default.clone(),
    }
}

fn main() {
    let args = Cli::parse();
    let model = select_model(&args.model, &Model::all());
    let path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .unwrap()
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(serde_json::from_reader(BufReader::new(File::open(path).unwrap())).unwrap())
    };
    let files = rustyms::csv::parse_csv(args.in_path, b',', None)
        .unwrap()
        .filter_map(|a| a.ok())
        .into_group_map_by(|l| l.index_column("mgf_path").unwrap().0.to_string());
    let out_file = BufWriter::new(File::create(args.out_path).unwrap());
    let mut out_data = Vec::new();
    let total_peptides = files.values().map(|f| f.len()).sum::<usize>();
    let mut done_peptides = 0;
    println!("Raw files: 0/{}, Peptides: 0/{total_peptides}", files.len());

    for (index, (file_name, lines)) in files.iter().enumerate() {
        let file = rustyms::rawfile::mgf::open(file_name).unwrap();

        let rows = lines
            .par_iter()
            .filter_map(|line| {
                let scan_number = line
                    .index_column("scan_number")
                    .unwrap()
                    .0
                    .parse::<usize>()
                    .unwrap();
                let z = line.index_column("z").unwrap().0.parse::<usize>().unwrap();
                let peptide = CompoundPeptidoform::pro_forma(
                    line.index_column("sequence").unwrap().0,
                    custom_database.as_ref(),
                )
                .unwrap();
                let selected_model = line
                    .index_column("fragmentation")
                    .map_or_else(|_| model.clone(), |(text, _)| select_model(text, &model));
                if let Some(spectrum) = file.iter().find(|s| s.raw_scan_number == Some(scan_number))
                {
                    let fragments =
                        peptide.generate_theoretical_fragments(Charge::new::<e>(z), &model);
                    let annotated = spectrum.annotate(
                        peptide,
                        &fragments,
                        &selected_model,
                        MassMode::Monoisotopic,
                    );
                    let scores: &Scores = &annotated
                        .scores(&fragments, &selected_model, MassMode::Monoisotopic)
                        .1[0][0];

                    let mut row: BTreeMap<_, _> = line.into();

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
                            "ion_Y_charge_independent".to_string(),
                            format!("{}", (unique_Y_found as f64 / unique_Y as f64),),
                        );
                    }
                    if args.report_intensity {
                        row.insert(
                            "annotated_intensity".to_string(),
                            match scores.score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.fraction().to_string()
                                }
                            },
                        );
                        row.insert(
                            "total_ion_current".to_string(),
                            match scores.score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.total.to_string()
                                }
                            },
                        );
                    }
                    if args.report_IL_satellite_coverage {
                        row.insert(
                            "IL_satellite_coverage".to_string(),
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
                                                    matches!(a.ion, FragmentType::w(s) | FragmentType::d(s) if s.sequence_index
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
                        "ion_combined".to_string(),
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
                        row.insert(format!("ion_{ion}"), format!("{}", recovered.fraction(),));
                    }
                    Some(row)
                } else {
                    eprintln!("Could not find scan number {scan_number} for file {file_name}");
                    None
                }
            })
            .collect::<Vec<_>>();
        out_data.extend_from_slice(&rows);
        done_peptides += lines.len();
        println!(
            "Raw files: {}/{}, Peptides: {done_peptides}/{total_peptides}",
            index + 1,
            files.len()
        );
    }

    rustyms::csv::write_csv(out_file, out_data).unwrap();
}
