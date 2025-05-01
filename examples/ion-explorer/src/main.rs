//! Explore common masses around fragments to detect which fragments actually occur

use std::{
    collections::BTreeMap,
    fs::File,
    io::{BufReader, BufWriter},
    ops::RangeInclusive,
    path::{Path, PathBuf},
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
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rustyms::{
    annotation::model::{
        ChargeRange, FragmentationModel, PrimaryIonSeries, SatelliteIonSeries, SatelliteLocation,
    },
    chemistry::MassMode,
    fragment::{Fragment, FragmentKind, FragmentType},
    identification::{SpectrumId, SpectrumIds},
    sequence::{AminoAcid, CompoundPeptidoformIon, SimpleModification},
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
    raw_file: Option<PathBuf>,
    /// The directory where to find non absolute raw files that are named in the peptides file
    #[arg(long)]
    raw_file_directory: Option<PathBuf>,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
    /// The bin width of the mass bins
    #[arg(long, default_value = "0.25")]
    resolution: f64,
    /// The high end of the range for looking at diagnostic/immonium ions, selects all peaks in 0..=`max_start`
    #[arg(long, default_value = "200.0")]
    max_start: f64,
    /// The number of Dalton to select before a fragment ion
    #[arg(long, default_value = "100.0")]
    before_fragment: f64,
    /// The number of Dalton to select after a fragment ion
    #[arg(long, default_value = "100.0")]
    after_fragment: f64,
}

fn main() {
    let args = Cli::parse();

    let model = FragmentationModel::none()
        .clone()
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
        false,
    )
    .expect("Could not open identified peptides file")
    .filter_map(Result::ok)
    .into_group_map_by(|l| match l.scans() {
        SpectrumIds::FileKnown(spectra) => spectra.first().map(|s| s.0.clone()),
        _ => None,
    });

    let stack: Stack = peptides.par_iter().map(|(file, peptides)| {
        let mut stack = Stack::default();
        let path = match file {
            Some(file) => {
                if file.is_absolute() {
                    file.clone()
                } else {
                    args.raw_file_directory.as_ref().expect("The raw file directory parameter has to be defined if there are peptides with a non absolute path").join(file)
                }
            }
            None => args.raw_file.clone().expect("The raw file parameter has to be defined if there are peptides without a defined raw file"),
        };
        let mut file = mzdata::io::MZReaderType::open_path(path).unwrap();

        for peptide in peptides {
            if peptide.charge().is_none() {
                continue;
            }
            if let Some(cpi) = peptide.peptidoform() {
                let id = match peptide.scans() {
                    SpectrumIds::FileKnown(spectra) => {
                        spectra.first().and_then(|s| s.1.first().cloned())
                    }
                    SpectrumIds::FileNotKnown(ids) => ids.first().cloned(),
                    SpectrumIds::None => None,
                };
                if let Some(spectrum) = match id {
                    Some(SpectrumId::Index(i)) => file.get_spectrum_by_index(i),
                    Some(SpectrumId::Native(n)) => file.get_spectrum_by_id(&n),
                    _ => continue,
                } {
                    let cpi = cpi.compound_peptidoform_ion();
                    let fragments =
                        cpi.generate_theoretical_fragments(peptide.charge().unwrap(), &model);
                    extract_and_merge(
                        &mut stack,
                        &spectrum,
                        &fragments,
                        &cpi,
                        &args,
                        peptide.mode().map(ToString::to_string),
                    );
                }
            }
        }
        stack
    }).collect();
    stack.store(Path::new(&args.out_path));
}

fn extract_and_merge(
    stack: &mut Stack,
    spectrum: &MultiLayerSpectrum,
    fragments: &[Fragment],
    peptidoform: &CompoundPeptidoformIon,
    args: &Cli,
    mode: Option<String>,
) {
    let spectrum = match spectrum.peaks() {
        RefPeakDataLevel::Centroid(c) => c,
        _ => return,
    };
    for fragment in fragments {
        if let Some(mz) = fragment.mz(MassMode::Monoisotopic) {
            let seq = fragment.ion.position().map(|pos| {
                let seq = &peptidoform.peptidoform_ions()
                    [fragment.peptidoform_ion_index.unwrap_or_default()]
                .peptidoforms()[fragment.peptidoform_index.unwrap_or_default()][pos.sequence_index];
                (
                    seq.aminoacid.aminoacid(),
                    seq.modifications
                        .iter()
                        .filter_map(|m| m.clone().into_simple())
                        .collect(),
                )
            });

            let key = match fragment.ion {
                FragmentType::d(_, _, d, _, _)
                | FragmentType::v(_, _, d, _)
                | FragmentType::w(_, _, d, _, _) => {
                    (fragment.ion.kind(), Some(d), seq, mode.clone())
                }
                _ => (fragment.ion.kind(), None, seq, mode.clone()),
            };
            let low = mz.value - args.before_fragment / fragment.charge.value as f64;
            let high = mz.value + args.after_fragment / fragment.charge.value as f64;
            let sub_spectrum = spectrum.between(low, high, mzdata::prelude::Tolerance::Da(0.0));

            merge_stack(
                stack.fragments.entry(key).or_default(),
                sub_spectrum,
                fragment.charge.value,
                mz.value,
                args.resolution,
                low..=high,
            );
        }
    }
    // Get start
    let sub_spectrum = spectrum.between(0.0, args.max_start, mzdata::prelude::Tolerance::Da(0.0));
    merge_stack(
        &mut stack.start,
        sub_spectrum,
        1,
        0.0,
        args.resolution,
        0.0..=args.max_start,
    );
}

fn merge_stack(
    points: &mut Vec<Point>,
    slice: &[CentroidPeak],
    charge: usize,
    center: f64,
    resolution: f64,
    range: RangeInclusive<f64>,
) {
    for found_peak in slice {
        if !range.contains(&found_peak.mz) {
            continue;
        }
        let normalised_mass =
            ((found_peak.mz - center) / resolution).round() * resolution * charge as f64;
        match points.binary_search_by(|p| p.mass.total_cmp(&normalised_mass)) {
            Ok(index) => {
                points[index].count += 1;
                points[index].total_intensity += f64::from(found_peak.intensity);
            }
            Err(index) => {
                points.insert(
                    index,
                    Point {
                        mass: normalised_mass,
                        count: 1,
                        total_intensity: found_peak.intensity as f64,
                    },
                );
            }
        }
    }
}

type ItemKey = (
    FragmentKind,
    Option<u8>,
    Option<(AminoAcid, Vec<SimpleModification>)>,
    Option<String>,
);

#[derive(Debug, Default)]
struct Stack {
    start: Vec<Point>,
    fragments: BTreeMap<ItemKey, Vec<Point>>,
}

impl Stack {
    fn store(&self, base_path: &Path) {
        write_stack(&base_path.join("start.csv"), &self.start);
        for (key, stack) in &self.fragments {
            let path = match key {
                (i, None, el, mode) => base_path.join(format!(
                    "fragment_{i}_{}_{}.csv",
                    el.as_ref().map_or("-".to_string(), |(aa, mods)| format!(
                        "{aa}{}",
                        mods.iter().map(|m| format!("[{m}]")).join("")
                    )),
                    mode.as_ref().map_or("-".to_string(), ToString::to_string)
                )),
                (i, Some(d), el, mode) => base_path.join(format!(
                    "fragment_{d}{i}_{}_{}.csv",
                    el.as_ref().map_or("-".to_string(), |(aa, mods)| format!(
                        "{aa}{}",
                        mods.iter().map(|m| format!("[{m}]")).join("")
                    )),
                    mode.as_ref().map_or("-".to_string(), ToString::to_string)
                )),
            };
            write_stack(&path, stack);
        }
    }

    fn merge(mut self, other: &Self) -> Self {
        merge_points(&mut self.start, &other.start);
        for (key, points) in &other.fragments {
            let entry = self.fragments.entry(key.clone()).or_default();
            merge_points(entry, points);
        }
        self
    }
}

fn merge_points(points: &mut Vec<Point>, other: &[Point]) {
    for point in other {
        match points.binary_search_by(|p| p.mass.total_cmp(&point.mass)) {
            Ok(index) => {
                points[index].count += point.count;
                points[index].total_intensity += point.total_intensity;
            }
            Err(index) => {
                points.insert(index, point.clone());
            }
        }
    }
}

impl rayon::iter::FromParallelIterator<Self> for Stack {
    fn from_par_iter<I>(par_iter: I) -> Self
    where
        I: rayon::prelude::IntoParallelIterator<Item = Self>,
    {
        let par_iter = par_iter.into_par_iter();
        par_iter.reduce(Self::default, |acc, s| acc.merge(&s))
    }
}

fn write_stack(path: &Path, points: &[Point]) {
    let out_file =
        BufWriter::new(File::create(path).unwrap_or_else(|e| {
            panic!("Could not create file: '{}'\n{e}", path.to_string_lossy())
        }));
    rustyms::identification::csv::write_csv(
        out_file,
        points.iter().map(|p| {
            [
                ("mass".to_string(), p.mass.to_string()),
                ("count".to_string(), p.count.to_string()),
                ("total_intensity".to_string(), p.total_intensity.to_string()),
            ]
        }),
        ',',
    )
    .unwrap();
}

#[derive(Clone, Debug)]
struct Point {
    mass: f64,
    count: usize,
    total_intensity: f64,
}
