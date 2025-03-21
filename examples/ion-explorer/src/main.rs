#![allow(non_snake_case)] // charge_independent_Y needs the capital as it means the glycan fragmentation
use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter},
    path::Path,
};

use clap::Parser;
use directories::ProjectDirs;
use itertools::Itertools;
use mzdata::{
    io::{MZFileReader, SpectrumSource},
    mzpeaks::{CentroidPeak, PeakCollection},
    prelude::SpectrumLike,
    spectrum::{MultiLayerSpectrum, RefPeakDataLevel},
};
use rustyms::{
    fragment::{FragmentKind, FragmentType},
    identification::{SpectrumId, SpectrumIds},
    model::{ChargeRange, PrimaryIonSeries, SatelliteIonSeries, SatelliteLocation},
    modification::{Ontology, SimpleModification},
    *,
};

#[derive(Parser)]
struct Cli {
    /// The input identified peptides file
    #[arg(short, long)]
    in_path: String,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: String,
    /// The raw file to use for any file without a raw file
    #[arg(long)]
    raw_file: Option<String>,
    /// The directory where to find raw files that are named in the peptides file
    #[arg(long)]
    raw_file_directory: String,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
    /// The bin width of the mz bins
    #[arg(long, default_value = "0.25")]
    resolution: f64,
    /// The high end of the range for looking at diagnostic/immonium ions, selects all peaks in 0..=`max_start`
    #[arg(long, default_value = "200.0")]
    max_start: f64,
    /// The number of Thompson to select before a fragment ion
    #[arg(long, default_value = "100.0")]
    before_fragment: f64,
    /// The number of Thompson to select after a fragment ion
    #[arg(long, default_value = "100.0")]
    after_fragment: f64,
}

fn main() {
    let args = Cli::parse();

    let model = Model::none()
        .b(PrimaryIonSeries::default())
        .d(SatelliteIonSeries::default().location(SatelliteLocation {
            rules: Vec::new(),
            base: Some(6),
        }))
        .v(SatelliteIonSeries::default().location(SatelliteLocation {
            rules: Vec::new(),
            base: Some(6),
        }))
        .w(SatelliteIonSeries::default().location(SatelliteLocation {
            rules: Vec::new(),
            base: Some(6),
        }))
        .y(PrimaryIonSeries::default())
        .precursor(Vec::new(), Vec::new(), (0, None), ChargeRange::PRECURSOR);

    let path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .unwrap()
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(serde_json::from_reader(BufReader::new(File::open(path).unwrap())).unwrap())
    };
    let peptides = rustyms::identification::open_identified_peptides_file(
        &args.in_path,
        custom_database.as_ref(),
    )
    .expect("Could not open identified peptides file")
    .filter_map(|a| a.ok())
    .into_group_map_by(|l| match l.scans() {
        SpectrumIds::FileKnown(spectra) => spectra.first().map(|s| s.0.clone()),
        _ => None,
    });

    let mut stack = Stack::default();

    for (file, peptides) in peptides {
        let mut file = mzdata::io::MZReaderType::open_path(file.map(|f|Path::new(&args.raw_file_directory).join(f)).or(args.raw_file.as_ref().map(|p| p.into())).expect("The raw file parameter has to be defined if there are peptides without a defined raw file")).unwrap();

        for peptide in peptides {
            if peptide.charge().is_none() {
                continue;
            }
            if let Some(cpi) = peptide.peptide() {
                let id = match peptide.scans() {
                    SpectrumIds::FileKnown(spectra) => {
                        spectra.first().and_then(|s| s.1.first().cloned())
                    }
                    SpectrumIds::FileNotKnown(ids) => ids.first().cloned(),
                    _ => None,
                };
                if let Some(spectrum) = match id {
                    Some(SpectrumId::Index(i)) => file.get_spectrum_by_index(i),
                    Some(SpectrumId::Native(n)) => file.get_spectrum_by_id(&n),
                    _ => continue,
                } {
                    let cpi = cpi.compound_peptidoform();
                    let fragments =
                        cpi.generate_theoretical_fragments(peptide.charge().unwrap(), &model);
                    extract_and_merge(
                        &mut stack,
                        &spectrum,
                        &fragments,
                        &[(
                            (AminoAcid::Asparagine, Vec::new()),
                            (
                                AminoAcid::Asparagine,
                                vec![Ontology::Unimod.find_id(7, None).unwrap()],
                            ),
                        )],
                        &cpi,
                        &args,
                    );
                }
            }
        }
    }
    stack.store(Path::new(&args.out_path));
}

type ComparisonKey = (AminoAcid, Vec<SimpleModification>);

fn extract_and_merge(
    stack: &mut Stack,
    spectrum: &MultiLayerSpectrum,
    fragments: &[Fragment],
    comparisons: &[(ComparisonKey, ComparisonKey)],
    peptidoform: &CompoundPeptidoformIon,
    args: &Cli,
) {
    let spectrum = match spectrum.peaks() {
        RefPeakDataLevel::Centroid(c) => c,
        _ => return,
    };
    for fragment in fragments {
        if let Some(mz) = fragment.mz(MassMode::Monoisotopic) {
            let key = match fragment.ion {
                FragmentType::d(_, _, d, _)
                | FragmentType::v(_, _, d, _)
                | FragmentType::w(_, _, d, _) => (fragment.ion.kind(), Some(d)),
                _ => (fragment.ion.kind(), None),
            };
            let low = mz.value - args.before_fragment;
            let high = mz.value + args.after_fragment;
            let sub_spectrum = spectrum.between(low, high, mzdata::prelude::Tolerance::Da(0.0));
            merge_stack(
                stack.fragments.entry(key).or_default(),
                sub_spectrum,
                mz.value,
                args.resolution,
            );
            // Comparison
            let (kind, pos) = match fragment.ion {
                FragmentType::a(pos, _)
                | FragmentType::b(pos, _)
                | FragmentType::c(pos, _)
                | FragmentType::x(pos, _)
                | FragmentType::y(pos, _)
                | FragmentType::z(pos, _) => (fragment.ion.kind(), pos),
                _ => continue,
            };
            let seq = &peptidoform.peptidoform_ions()
                [fragment.peptidoform_ion_index.unwrap_or_default()]
            .peptidoforms()[fragment.peptidoform_index.unwrap_or_default()][pos.sequence_index];
            let key: ComparisonKey = (
                seq.aminoacid.aminoacid(),
                seq.modifications
                    .iter()
                    .filter_map(|m| m.simple())
                    .cloned()
                    .collect(),
            );
            for (a, b) in comparisons {
                if key == *a {
                    merge_comparison_stack(
                        stack
                            .comparison
                            .entry((kind, a.clone(), b.clone()))
                            .or_default(),
                        sub_spectrum,
                        mz.value,
                        args.resolution,
                        true,
                    );
                } else if key == *b {
                    merge_comparison_stack(
                        stack
                            .comparison
                            .entry((kind, a.clone(), b.clone()))
                            .or_default(),
                        sub_spectrum,
                        mz.value,
                        args.resolution,
                        false,
                    );
                }
            }
        }
    }
    // Get start
    let sub_spectrum = spectrum.between(0.0, args.max_start, mzdata::prelude::Tolerance::Da(0.0));
    merge_stack(&mut stack.start, sub_spectrum, 0.0, args.resolution);
}

fn merge_stack(points: &mut Vec<Point>, slice: &[CentroidPeak], center: f64, resolution: f64) {
    for found_peak in slice {
        let normalised_mz = ((found_peak.mz - center) / resolution).round() * resolution;
        match points.binary_search_by(|p| p.mz.total_cmp(&normalised_mz)) {
            Ok(index) => {
                points[index].count += 1;
                points[index].total_intensity += found_peak.intensity as f64;
            }
            Err(index) => {
                points.insert(
                    index,
                    Point {
                        mz: normalised_mz,
                        count: 1,
                        total_intensity: found_peak.intensity as f64,
                    },
                );
            }
        }
    }
}

fn merge_comparison_stack(
    points: &mut Vec<ComparisonPoint>,
    slice: &[CentroidPeak],
    center: f64,
    resolution: f64,
    is_a: bool,
) {
    for found_peak in slice {
        let normalised_mz = ((found_peak.mz - center) / resolution).round() * resolution;
        match points.binary_search_by(|p| p.mz.total_cmp(&normalised_mz)) {
            Ok(index) => {
                if is_a {
                    points[index].count_a += 1;
                    points[index].total_intensity_a += found_peak.intensity as f64;
                } else {
                    points[index].count_b += 1;
                    points[index].total_intensity_b += found_peak.intensity as f64;
                }
            }
            Err(index) => {
                points.insert(
                    index,
                    if is_a {
                        ComparisonPoint {
                            mz: normalised_mz,
                            count_a: 1,
                            count_b: 0,
                            total_intensity_a: found_peak.intensity as f64,
                            total_intensity_b: 0.0,
                        }
                    } else {
                        ComparisonPoint {
                            mz: normalised_mz,
                            count_a: 0,
                            count_b: 1,
                            total_intensity_a: 0.0,
                            total_intensity_b: found_peak.intensity as f64,
                        }
                    },
                );
            }
        }
    }
}

#[derive(Default, Debug)]
struct Stack {
    start: Vec<Point>,
    fragments: BTreeMap<(FragmentKind, Option<u8>), Vec<Point>>,
    comparison: BTreeMap<(FragmentKind, ComparisonKey, ComparisonKey), Vec<ComparisonPoint>>,
}

impl Stack {
    fn store(&self, base_path: &Path) {
        write_stack(&base_path.join("start.csv"), &self.start);
        for (key, stack) in self.fragments.iter() {
            let path = match key {
                (i, None) => base_path.join(format!("fragment_{i}.csv")),
                (i, Some(d)) => base_path.join(format!("fragment_{i}_{d}.csv")),
            };
            write_stack(&path, stack);
        }
        for ((fragment, a, b), stack) in self.comparison.iter() {
            let a = format!("{}{}", a.0, a.1.iter().map(|m| format!("[{m}]")).join(""));
            let b = format!("{}{}", b.0, b.1.iter().map(|m| format!("[{m}]")).join(""));
            let path = base_path.join(format!("comparison_{fragment}_{a}_{b}.csv"));
            write_comparison_stack(&path, stack);
        }
    }
}

fn write_stack(path: &Path, points: &[Point]) {
    let out_file = BufWriter::new(File::create(path).unwrap());
    rustyms::csv::write_csv(
        out_file,
        points.iter().map(|p| {
            [
                ("mz".to_string(), p.mz.to_string()),
                ("count".to_string(), p.count.to_string()),
                (
                    "avg_intensity".to_string(),
                    (p.total_intensity / p.count as f64).to_string(),
                ),
            ]
        }),
    )
    .unwrap();
}

fn write_comparison_stack(path: &Path, points: &[ComparisonPoint]) {
    let out_file = BufWriter::new(File::create(path).unwrap());
    rustyms::csv::write_csv(
        out_file,
        points.iter().map(|p| {
            [
                ("mz".to_string(), p.mz.to_string()),
                ("count_a".to_string(), p.count_a.to_string()),
                ("count_b".to_string(), p.count_b.to_string()),
                (
                    "avg_intensity_a".to_string(),
                    if p.count_a == 0 {0.0} else {p.total_intensity_a / p.count_a as f64}.to_string(),
                ),
                (
                    "avg_intensity_b".to_string(),
                    if p.count_b == 0 {0.0} else {p.total_intensity_b / p.count_b as f64}.to_string(),
                ),
            ]
        }),
    )
    .unwrap();
}

#[derive(Debug)]
struct Point {
    mz: f64,
    count: usize,
    total_intensity: f64,
}

#[derive(Debug)]
struct ComparisonPoint {
    mz: f64,
    count_a: usize,
    count_b: usize,
    total_intensity_a: f64,
    total_intensity_b: f64,
}
