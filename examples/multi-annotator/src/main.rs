//! Annotate many peptides at once
#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation

use mzsignal as _; // To aid vesion selection

use std::{
    collections::BTreeMap,
    fs::File,
    io::BufWriter,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use clap::Parser;
use directories::ProjectDirs;
use fragment::FragmentType;
use itertools::Itertools;
use mzdata::io::{MZFileReader, SpectrumSource};
use rayon::prelude::*;
use rustyms::{
    annotation::{
        AnnotatableSpectrum, AnnotatedPeak, Score, Scores,
        model::{FragmentationModel, MatchingParameters},
    },
    chemistry::MassMode,
    fragment::{DiagnosticPosition, Fragment},
    glycan::MonoSaccharide,
    identification::{BasicCSVData, IdentifiedPeptidoformSource, csv::write_csv},
    sequence::{
        AminoAcid, GnoComposition, SequencePosition, SimpleModificationInner,
        parse_custom_modifications,
    },
    spectrum::PeakSpectrum,
    *,
};

/// The command line interface arguments
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
struct Cli {
    /// The input csv file, should have the following columns: 'raw_file' (full path), 'scan_index' (0-based index in the raw file), 'z', 'sequence', and can have 'mode' (etd/td_etd/ethcd/etcad/eacid/ead/hcd/cid/all/none, defaults to the global model)
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
    /// Turns on reporting of ambiguous amino acids, returns a column for each ambiguous amino acid
    /// (J/B/Z). This column it contains `{found_a}/{total_a}|{found_b}/{total_b}` if there are
    /// multiple ambiguous amino acids the values are separated by semicolons and the options are
    /// ordered by location in the peptidoform(s). For J the order is always I/L, for B D/N, for Z
    /// E/Q.
    ///
    /// If `report_intensity` is also turned on an additional column is added per amino acid which
    /// contains `{fraction_tic_annotated_a}|{fraction_tic_annotated_b}`. This follows the same
    /// structure as the fragments column.
    #[arg(long)]
    report_ambiguous_amino_acids: bool,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
    /// Turns on specific counting for glycan fragments. This collects statistics for any fragment
    /// that contains the given monosaccharide. It returns two columns for every given monosaccharide
    /// one containing B numbers `{found}/{total}` and one containing Y numbers `{found}/{total}`.
    ///
    /// This parameters should be specified with monosaccharide names separated by commas.
    #[arg(long, value_parser=monosaccharide_parser, value_delimiter = ',')]
    glycan_buckets: Vec<MonoSaccharide>,
}

fn monosaccharide_parser(value: &str) -> Result<MonoSaccharide, String> {
    MonoSaccharide::from_short_iupac(value, 0, 0)
        .map(|s| s.0)
        .map_err(|e| e.to_string())
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
        .expect("Could not generate Annotator configurationpath (needed to check if custom modifications are defined)")
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(parse_custom_modifications(&path).expect("Could not parse custom modifications file, if you do not need these you can skip parsing them using the appropriate flag"))
    };
    let files = BasicCSVData::parse_file(args.in_path, custom_database.as_ref(), true, None)
        .expect("Invalid input file")
        .filter_map(Result::ok)
        .into_group_map_by(|l| l.raw_file.clone());
    let out_file =
        BufWriter::new(File::create(args.out_path).expect("Could not create out CSV file"));
    let total_peptides = files.values().map(Vec::len).sum::<usize>();
    let peptides_counter = AtomicUsize::default();
    let raw_file_counter = AtomicUsize::default();
    println!("Raw files: 0/{}, Peptides: 0/{total_peptides}", files.len());

    let precursor_mass = Arc::new("precursor_mass".to_string());
    let column_y_independent = Arc::new("ion_Y_charge_independent".to_string());
    let intensity_combined = Arc::new("intensity_combined".to_string());
    let tic = Arc::new("total_ion_current".to_string());

    let out_data: Vec<_> =  files.par_iter().flat_map(|(file_name, lines)| {
        let mut file = mzdata::io::MZReaderType::open_path(file_name).unwrap_or_else(|err| {eprintln!("Could not open raw file: {}\nError: {err}", file_name.to_string_lossy()); std::process::exit(2)});

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

                    row.insert(precursor_mass.clone(), line.sequence.formulas().unique().iter().map(|f| f.monoisotopic_mass().value.to_string()).join(";"));

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
                            column_y_independent.clone(),
                            format!("{}", (unique_Y_found as f64 / unique_Y as f64),),
                        );
                    }
                    if args.report_intensity {
                        row.insert(
                            intensity_combined.clone(),
                            match scores.score {
                                Score::Position { intensity, .. }
                                | Score::UniqueFormulas { intensity, .. } => {
                                    intensity.fraction().to_string()
                                }
                            },
                        );
                        row.insert(
                            tic.clone(),
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
                            annotated.peptide.clone().singular_peptidoform().map_or(
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
                    if args.report_ambiguous_amino_acids {
                        let stats = annotated.ambigous_statistics(&fragments, &parameters, MassMode::Monoisotopic);
                        row.insert(
                            Arc::new("ambiguous_J_found".to_string()),
                            stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::Isoleucine).map(|aa| format!("{}/{}|{}/{}", aa.optiona_a.1.found, aa.optiona_a.1.total, aa.optiona_b.1.found, aa.optiona_b.1.total)).join(";")
                        );
                        if args.report_intensity {
                            row.insert(
                                Arc::new("ambiguous_J_intensity".to_string()),
                                stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::Isoleucine).map(|aa| format!("{}|{}", aa.optiona_a.2.fraction(), aa.optiona_b.2.fraction())).join(";")
                            );
                        }
                        row.insert(
                            Arc::new("ambiguous_B_found".to_string()),
                            stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::AsparticAcid).map(|aa| format!("{}/{}|{}/{}", aa.optiona_a.1.found, aa.optiona_a.1.total, aa.optiona_b.1.found, aa.optiona_b.1.total)).join(";")
                        );
                        if args.report_intensity {
                            row.insert(
                                Arc::new("ambiguous_B_intensity".to_string()),
                                stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::AsparticAcid).map(|aa| format!("{}|{}", aa.optiona_a.2.fraction(), aa.optiona_b.2.fraction())).join(";")
                            );
                        }
                        row.insert(
                            Arc::new("ambiguous_Z_found".to_string()),
                            stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::GlutamicAcid).map(|aa| format!("{}/{}|{}/{}", aa.optiona_a.1.found, aa.optiona_a.1.total, aa.optiona_b.1.found, aa.optiona_b.1.total)).join(";")
                        );
                        if args.report_intensity {
                            row.insert(
                                Arc::new("ambiguous_Z_intensity".to_string()),
                                stats.aminoacids.iter().filter(|aa| aa.optiona_a.0 == AminoAcid::GlutamicAcid).map(|aa| format!("{}|{}", aa.optiona_a.2.fraction(), aa.optiona_b.2.fraction())).join(";")
                            );
                        }
                    }
                    if !args.glycan_buckets.is_empty() {
                        #[derive(Debug)]
                        struct Match<'a> {
                            target: &'a MonoSaccharide,
                            found_B: usize,
                            total_B: usize,
                            found_Y: usize,
                            total_Y: usize,
                        }

                        let mut buckets = args.glycan_buckets.iter().map(|m| Match {target: m,found_B: 0, total_B: 0, found_Y: 0, total_Y: 0}).collect_vec();
                        for (theoretical, f) in fragments.iter().map(|f| (true, f)).chain(annotated.spectrum().flat_map(|p| p.annotation.iter().map(|a| (false, a)))) {
                            match &f.ion {
                                FragmentType::Y(pos) => {
                                    if let Some((_, seq)) = pos.first().and_then(|p| p.attachment) {
                                        let element = &annotated.peptide.peptidoform_ions()[f.peptidoform_ion_index.unwrap_or_default()].peptidoforms()[f.peptidoform_index.unwrap_or_default()][seq];
                                        if let Some(glycan) = element.modifications.iter().find_map(|m| match (*m).clone().into_simple().as_deref() {
                                            Some(SimpleModificationInner::Gno { composition: GnoComposition::Topology(structure), .. } | SimpleModificationInner::GlycanStructure(structure)) => Some(structure.clone()), _ => None}) {
                                                for bucket in &mut buckets {
                                                    if glycan.contains(bucket.target, false, None, pos).unwrap_or_default() {
                                                        if theoretical {
                                                            bucket.total_Y += 1; }
                                                        else {
                                                            bucket.found_Y += 1 ;
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                },
                                FragmentType::YComposition(composition, _) => {
                                    for bucket in &mut buckets {
                                        if composition.iter().any(|c| c.0.equivalent(bucket.target, false)) {
                                            if theoretical {
                                                bucket.total_Y += 1; }
                                            else {
                                                bucket.found_Y += 1 ;
                                            }
                                        }
                                    }
                                }
                                FragmentType::B{b, y,end: _} => {
                                    if let Some((_, seq)) = b.attachment.as_ref() {
                                        let element = &annotated.peptide.peptidoform_ions()[f.peptidoform_ion_index.unwrap_or_default()].peptidoforms()[f.peptidoform_index.unwrap_or_default()][*seq];
                                        if let Some(glycan) = element.modifications.iter().find_map(|m| match (*m).clone().into_simple().as_deref() {
                                            Some(SimpleModificationInner::Gno { composition: GnoComposition::Topology(structure), .. } | SimpleModificationInner::GlycanStructure(structure)) => Some(structure.clone()), _ => None}) {
                                                for bucket in &mut buckets {
                                                    if glycan.contains(bucket.target, false, Some(b), y).unwrap_or_default() {
                                                        if theoretical {
                                                            bucket.total_B += 1; }
                                                        else {
                                                            bucket.found_B += 1 ;
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                },
                                FragmentType::BComposition(composition, _) => {
                                    for bucket in &mut buckets {
                                        if composition.iter().any(|c| c.0.equivalent(bucket.target, false)) {
                                            if theoretical {
                                                bucket.total_B += 1; }
                                            else {
                                                bucket.found_B += 1 ;
                                            }
                                        }
                                    }
                                },
                                FragmentType::Diagnostic(DiagnosticPosition::Glycan(_, sug) | DiagnosticPosition::GlycanCompositional(sug, _)) => {
                                    for bucket in &mut buckets {
                                        if sug.equivalent(bucket.target, false) {
                                            if theoretical {
                                                bucket.total_B += 1 ;}
                                            else {
                                                bucket.found_B += 1 ;
                                            }
                                        }
                                    }
                                }
                                _ => ()
                            }
                        }

                        for bucket in buckets {
                            row.insert(Arc::new(format!("glycan_bucket_{}_B", bucket.target)), format!("{}/{}", bucket.found_B, bucket.total_B));
                            row.insert(Arc::new(format!("glycan_bucket_{}_Y", bucket.target)), format!("{}/{}", bucket.found_Y, bucket.total_Y));
                        }
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
            }).map(|row| row.into_iter().collect_vec()).collect_vec();
        println!(
            "Raw files: {}/{}, Peptides: {}/{total_peptides}",
            raw_file_counter.fetch_add(1, Ordering::Relaxed) + 1,
            files.len(),
            peptides_counter.fetch_add(lines.len(), Ordering::Relaxed) + lines.len(),
        );
        rows
    }).collect();

    write_csv(
        out_file,
        out_data
            .into_iter()
            .map(|r| r.into_iter().map(|(k, v)| (Arc::unwrap_or_clone(k), v))),
        ',',
    )
    .expect("Could not write results in out CSV file");
}
