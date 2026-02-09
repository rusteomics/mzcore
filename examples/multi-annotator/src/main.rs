//! Annotate many peptides at once
#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation

use core::f32;
use std::{
    collections::BTreeMap,
    fs::File,
    io::BufWriter,
    sync::{
        Arc,
        atomic::{AtomicUsize, Ordering},
    },
};

use context_error::CombineErrorsExtender;

use clap::Parser;
use directories::ProjectDirs;
use itertools::Itertools;
use mzannotate::{
    annotation::{Score, Scores, model::parse_custom_models},
    fragment::{DiagnosticPosition, FragmentType},
    mzspeclib::AnalyteTarget,
    prelude::*,
};
use mzcore::{
    chemistry::MassMode,
    csv::write_csv,
    glycan::MonoSaccharide,
    ontology::Ontologies,
    quantities::Tolerance,
    sequence::{AminoAcid, GnoComposition, SequencePosition, SimpleModificationInner},
    system::MassOverCharge,
};
use mzdata::{
    io::{MZFileReader, SpectrumSource},
    mzpeaks::PeakCollection,
    mzsignal::PeakPicker,
    spectrum::{SignalContinuity, SpectrumLike},
};
use mzident::{BasicCSVPSM, PSMSource};
use rayon::prelude::*;

/// The command line interface arguments
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
struct Cli {
    /// The input csv file, should have the following columns: `raw_file` (full path), `scan_index` (0-based index in the raw file), `z`, `sequence`, and can have `mode` (`etd`/`td_etd`/`ethcd`/`etcad`/`eacid`/`ead`/`hcd`/`cid`/`all`/`none`, defaults to the global model)
    #[arg(short, long)]
    in_path: String,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: String,
    /// Global mode, will be overruled by line specific modes (`etd`/`td_etd`/`ethcd`/`etcad`/`eacid`/`ead`/`hcd`/`cid`/`all`/`none`)
    #[arg(long, default_value_t = String::from("all"))]
    mode: String,
    /// Turns on reporting of glycan Y-ions in a charge independent manner
    #[arg(long)]
    charge_independent_Y: bool,
    /// Turns on reporting of intensity statistics, the fraction of total TIC that could be annotated as well as the TIC for each spectrum
    #[arg(long)]
    report_intensity: bool,
    /// Turns on reporting of false match rate statistics
    #[arg(long)]
    report_false_match_rate: bool,
    /// Turns on reporting of I/L coverage by satellite ions, returns a list with a 0 (not covered) or 1 (covered) for each I or L in the peptide
    #[arg(long)]
    report_IL_satellite_coverage: bool,
    /// Turns on reporting of fragments found and peaks annotated statistics
    #[arg(long)]
    report_basic_statistics: bool,
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
    /// To turn off loading the custom models from the Annotator (if installed)
    #[arg(long)]
    no_custom_models: bool,
    /// Turns on specific counting for glycan fragments. This collects statistics for any fragment
    /// that contains the given monosaccharide. It returns two columns for every given monosaccharide
    /// one containing B numbers `{found}/{total}` and one containing Y numbers `{found}/{total}`.
    ///
    /// This parameters should be specified with monosaccharide names separated by commas.
    #[arg(long, value_parser=monosaccharide_parser, value_delimiter = ',')]
    glycan_buckets: Vec<MonoSaccharide>,
    /// Filter the MS2 spectra to only contain peaks above this fraction of the intensity value of the TIC.
    /// First the TIC filter is applied, then the base peak, then the absolute.
    #[arg(long)]
    tic_noise_threshold: Option<f32>,
    /// Filter the MS2 spectra to only contain peaks above this fraction of the intensity value of the base peaks.
    /// First the TIC filter is applied, then the base peak, then the absolute.
    #[arg(long)]
    basepeak_noise_threshold: Option<f32>,
    /// Filter the MS2 spectra to only contain peaks above this intensity value.
    /// First the TIC filter is applied, then the base peak, then the absolute.
    #[arg(long)]
    absolute_noise_threshold: Option<f32>,
    /// MS2 fragment match tolerance. Use a number followed by `ppm`, `Th`, `m/z`, or `mz` (capitalisation is ignored).
    #[arg(short, long, default_value_t = Tolerance::new_ppm(20.0.into()), value_parser=mass_tolerance_parse)]
    tolerance: Tolerance<MassOverCharge>,
}

/// Parse a monosaccharide from the string.
/// # Errors
/// If the string does not contain a well formed monosacchride.
fn monosaccharide_parser(value: &str) -> Result<MonoSaccharide, String> {
    MonoSaccharide::from_short_iupac(value, 0, 0)
        .map(|s| s.0)
        .map_err(|e| e.to_string())
}
/// # Errors
/// If the string is not `xxppm` or `xxth`.
fn mass_tolerance_parse(input: &str) -> Result<Tolerance<MassOverCharge>, &'static str> {
    input.parse().map_err(|()| "Invalid tolerance parameter")
}

fn select_model<'a>(
    text: &str,
    default: &'a FragmentationModel,
    custom: Option<&'a [(String, FragmentationModel)]>,
) -> &'a FragmentationModel {
    match text.to_ascii_lowercase().as_str() {
        "etd" => FragmentationModel::etd(),
        "td_etd" => FragmentationModel::td_etd(),
        "ethcd" | "etcad" => FragmentationModel::etcid(),
        "eacid" => FragmentationModel::eacid(),
        "ead" => FragmentationModel::ead(),
        "hcd" | "cid" => FragmentationModel::cid(),
        "all" => FragmentationModel::all(),
        "none" => FragmentationModel::none(),
        _ => custom.map_or(default, |m| {
            m.iter()
                .find_map(|(name, model)| name.eq_ignore_ascii_case(text).then_some(model))
                .unwrap_or(default)
        }),
    }
}

/// Run the multi annotator program
/// # Panics
/// * If the projectdirs could not be made.
/// * If the custom modifications file could not be parsed.
/// * If a selected spectrum could not be peak picked.
/// * If the input file does not exist, or could not be recognised as a supported format.
/// * If the output csv file could not be created, or written to.
#[allow(clippy::cognitive_complexity)]
fn main() {
    let args = Cli::parse();

    let parameters = MatchingParameters::default().tolerance(args.tolerance);
    let custom_models_path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .expect("Could not generate Annotator configuration path (needed to check if custom models are defined)")
        .config_dir()
        .join("../custom_models.json");
    let custom_models = if args.no_custom_models || !custom_models_path.exists() {
        None
    } else {
        Some(parse_custom_models(&custom_models_path).expect("Could not parse custom models file, if you do not need these you can skip parsing them using the appropriate flag"))
    };
    let model = select_model(
        &args.mode,
        FragmentationModel::all(),
        custom_models.as_deref(),
    );
    let ontologies = Ontologies::init().0;
    let mut peptidoforms = BasicCSVPSM::parse_file(args.in_path, &ontologies, true, None)
        .expect("Invalid input file")
        .combine_errors();
    let files = peptidoforms.into_group_map_by(|l| l.raw_file.clone());
    if !peptidoforms.errors().is_empty() {
        for e in peptidoforms.errors() {
            eprintln!("{e}");
        }
        eprintln!(
            "Errors were found while parsing the peptidoform CSV file, the program will continue but will ignore all failed lines"
        );
    }
    let out_file =
        BufWriter::new(File::create(args.out_path).expect("Could not create out CSV file"));
    let total_peptides = files.values().map(Vec::len).sum::<usize>();
    let peptides_counter = AtomicUsize::default();
    let raw_file_counter = AtomicUsize::default();
    println!("Raw files: 0/{}, Peptides: 0/{total_peptides}", files.len());

    let precursor_mass = Arc::new("precursor_mass".to_string());
    let column_y_independent = Arc::new("ion_Y_charge_independent".to_string());
    let intensity_combined = Arc::new("intensity_combined".to_string());
    let fragments_found = Arc::new("fragments_found".to_string());
    let peaks_annotated = Arc::new("peaks_annotated".to_string());
    let peaks_fdr_column = Arc::new("peaks_fdr".to_string());
    let intensity_fdr_column = Arc::new("intensity_fdr".to_string());
    let tic = Arc::new("total_ion_current".to_string());
    let base_peak = Arc::new("base_peak_intensity".to_string());

    let out_data: Vec<_> =  files.par_iter().flat_map(|(file_name, lines)| {
        let mut file = mzdata::io::MZReaderType::open_path(file_name).unwrap_or_else(|err| {eprintln!("Could not open raw file: {}\nError: {err}", file_name.to_string_lossy()); std::process::exit(2)});

        let rows = lines
            .iter()
            .filter_map(|line| {
                let selected_model = line.mode.as_ref()
                    .map_or(model, |text| select_model(text, model, custom_models.as_deref(),));
                if let Some(mut spectrum) = file.get_spectrum_by_index(line.scan_index)
                {
                    if spectrum.signal_continuity() == SignalContinuity::Profile {
                        spectrum
                            .pick_peaks_with(&PeakPicker::default())
                            .unwrap_or_else(|err| {
                                panic!(
                                    "Spectrum could not be peak picked: {err}",
                                )
                            });
                    } else if spectrum.arrays.is_some() && spectrum.peaks.is_none() {
                        // USI spectra are mostly loaded as the binary array maps instead of peaks regardless of the signal continuity level
                        spectrum.peaks = spectrum.arrays.as_ref().map(Into::into);
                    }
                    if let Some(threshold) = args.tic_noise_threshold && let Some(peaks) = spectrum.peaks.as_mut() {
                        let threshold = peaks.total_ion_current() * threshold;
                        peaks.peaks.retain(|p: &mzdata::mzpeaks::CentroidPeak| p.intensity > threshold);
                    }
                    if let Some(threshold) = args.basepeak_noise_threshold && let Some(peaks) = spectrum.peaks.as_mut() {
                        let threshold = peaks.base_peak().map_or(0.0, |v| v.intensity  * threshold);
                        peaks.peaks.retain(|p: &mzdata::mzpeaks::CentroidPeak| p.intensity > threshold);
                    }
                    if let Some(threshold) = args.absolute_noise_threshold && let Some(peaks) = spectrum.peaks.as_mut() {
                            peaks.peaks.retain(|p: &mzdata::mzpeaks::CentroidPeak| p.intensity > threshold);
                    }
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
                            .peaks.iter()
                            .flat_map(|peak| &peak.annotations)
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
                        row.insert(
                            base_peak.clone(),
                            annotated.peaks.iter().map(|p| p.intensity).max_by(f32::total_cmp).unwrap_or(0.0).to_string(),
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
                    if args.report_basic_statistics {
                        row.insert(
                            fragments_found.clone(),
                            match scores.score {
                                Score::Position { fragments, .. }
                                | Score::UniqueFormulas { fragments, .. } => {
                                    fragments.fraction().to_string()
                                }
                            },
                        );
                        row.insert(
                            peaks_annotated.clone(),
                            match scores.score {
                                Score::Position { peaks, .. }
                                | Score::UniqueFormulas { peaks, .. } => {
                                    peaks.fraction().to_string()
                                }
                            },
                        );
                    }
                    if args.report_false_match_rate {
                        let (fdr, _) = annotated.fdr(&fragments, &parameters, MassMode::Monoisotopic);
                        row.insert(
                            peaks_fdr_column.clone(),
                            fdr.peaks_fdr().to_string(),
                        );
                        row.insert(
                            intensity_fdr_column.clone(),
                            fdr.intensity_fdr().to_string(),
                        );
                    }
                    if args.report_IL_satellite_coverage {
                        row.insert(
                            Arc::new("IL_satellite_coverage".to_string()),
                            annotated.analytes.first().and_then(|a| match &a.target {
                                    AnalyteTarget::PeptidoformIon(pep) => pep.singular_ref(),
                                    _ => None,
                                }).map_or(
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
                                            if annotated.peaks.iter().any(|p| {
                                                p.annotations.iter().any(|a: &Fragment| {
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
                            ));
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
                        for (theoretical, f) in fragments.iter().map(|f| (true, f)).chain(annotated.peaks.iter().flat_map(|p| p.annotations.iter().map(|a| (false, a)))) {
                            match &f.ion {
                                FragmentType::Y(pos) => {
                                    if let Some((_, seq)) = pos.first().and_then(|p| p.attachment) {
                                        let element = &annotated.analytes.iter().find(|a| a.id as usize == f.peptidoform_ion_index.unwrap_or_default()).and_then(|a| match &a.target {
                                    AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                                    _ => None,
                                }).unwrap().peptidoforms()[f.peptidoform_index.unwrap_or_default()][seq];
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
                                        let element = &annotated.analytes.iter().find(|a| a.id as usize == f.peptidoform_ion_index.unwrap_or_default()).and_then(|a| match &a.target {
                                    AnalyteTarget::PeptidoformIon(pep) => Some(pep),
                                    _ => None,
                                }).unwrap().peptidoforms()[f.peptidoform_index.unwrap_or_default()][*seq];
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
