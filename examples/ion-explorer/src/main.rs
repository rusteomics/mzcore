//! Explore common masses around fragments to detect which fragments actually occur

use std::{
    collections::BTreeMap,
    fs::File,
    io::BufWriter,
    ops::RangeInclusive,
    path::{Path, PathBuf},
};

use clap::Parser;
use itertools::Itertools;
use mzannotate::{
    fragment::{FragmentKind, FragmentType},
    prelude::*,
};
use mzcore::{
    chemistry::MassMode,
    ontology::Ontologies,
    sequence::{Modification, PeptidoformIonSet},
};
use mzdata::{
    io::{MZFileReader, SpectrumSource},
    mzpeaks::PeakCollection,
    prelude::SpectrumLike,
    spectrum::RefPeakDataLevel,
};
use mzident::{SpectrumId, SpectrumIds, prelude::*};
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

#[derive(Parser)]
struct Cli {
    /// The input PSM file
    #[arg(short, long)]
    in_path: PathBuf,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: PathBuf,
    /// The raw file to use for any file without a raw file
    #[arg(long)]
    raw_file: Option<PathBuf>,
    /// The directory where to find non absolute raw files that are named in the peptides file
    #[arg(long)]
    raw_file_directory: Option<PathBuf>,
    /// The bin width of the mass bins
    #[arg(long, default_value = "0.25")]
    resolution: f64,
    /// The high end of the range for looking at diagnostic/immonium ions, selects all peaks in
    /// 0..=`max_start`
    #[arg(long, default_value = "200.0")]
    max_start: f64,
    /// The number of Dalton to select before a fragment ion
    #[arg(long, default_value = "100.0")]
    before_fragment: f64,
    /// The number of Dalton to select after a fragment ion
    #[arg(long, default_value = "100.0")]
    after_fragment: f64,
}

#[allow(clippy::missing_panics_doc)]
fn main() {
    let args = Cli::parse();

    let model = FragmentationModel::cid();
    let mut parameters = MatchingParameters::default();
    parameters.match_isotopes = true;
    let parameters = parameters;

    let peptides: std::collections::HashMap<PathBuf, Vec<PSM<mzcore::sequence::Linked, mzident::MaybePeptidoform>>> =
        open_psm_file(&args.in_path, &Ontologies::init().0, false)
            .expect("Could not open PSM file")
            .filter_map(Result::ok)
            .into_group_map_by(|l| match l.scans() {
                SpectrumIds::FileKnown(spectra) => spectra
                    .first()
                    .map(|s| {
                        if s.0.is_absolute() {
                            s.0.clone()
                        } else {
                            args.raw_file_directory
                                .as_ref()
                                .expect("The raw file directory parameter has to be defined if there are peptides with a non absolute path")
                                .join(&s.0)
                                .with_extension("mzML")
                        }
                    })
                    .expect("A file known spectra ref should have at least one file"),
                SpectrumIds::FileNotKnown(_) | SpectrumIds::None => args
                    .raw_file
                    .clone()
                    .expect("The raw file parameter has to be defined if there are peptides without a defined raw file"),
            });

    let stack: Stack = peptides
        .par_iter()
        .map(|(path, peptides)| {
            let mut stack = Stack::default();
            let mut file = mzdata::io::MZReaderType::open_path(path)
                .unwrap_or_else(|err| panic!("Could not open '{}' {err}", path.display()));

            println!("{}", peptides.len());
            for peptide in peptides {
                if peptide.charge().is_none() {
                    continue;
                }
                if let Some(cpi) = peptide.peptidoform_ion_set().map(std::borrow::Cow::into_owned) {
                    let id = match peptide.scans() {
                        SpectrumIds::FileKnown(spectra) => {
                            spectra.first().and_then(|s| s.1.first().cloned())
                        }
                        SpectrumIds::FileNotKnown(ids) => ids.first().cloned(),
                        SpectrumIds::None => None,
                    };
                    if let Some(spectrum) = match id {
                        Some(SpectrumId::Index(i)) => file.get_spectrum_by_index(i),
                        Some(SpectrumId::Number(i)) => file.get_spectrum_by_index(i - 1),
                        Some(SpectrumId::Native(n)) => file.get_spectrum_by_id(&n),
                        _ => continue,
                    } {
                        let fragments =
                            cpi.generate_theoretical_fragments(peptide.charge().unwrap(), model);
                        let mut annotated = spectrum.annotate(
                            cpi.clone(),
                            &fragments,
                            &parameters,
                            MassMode::Monoisotopic,
                        );
                        annotated.peaks.peaks.retain(|p| p.annotations.is_empty());
                        extract_and_merge(&mut stack, &annotated, &fragments, &cpi, &args);
                    }
                }
            }
            stack
        })
        .collect();
    stack.store(Path::new(&args.out_path));
}

fn extract_and_merge(
    stack: &mut Stack,
    spectrum: &AnnotatedSpectrum,
    fragments: &[Fragment],
    peptidoform: &PeptidoformIonSet,
    args: &Cli,
) {
    let RefPeakDataLevel::Centroid(spectrum) = spectrum.peaks() else {
        return;
    };

    let mods_present = peptidoform
        .peptidoform_ions()
        .iter()
        .flat_map(mzcore::sequence::PeptidoformIon::peptidoforms)
        .flat_map(mzcore::sequence::Peptidoform::sequence)
        .flat_map(|s| s.modifications.iter())
        .unique()
        .cloned()
        .map(Some)
        .chain([None])
        .collect::<Vec<_>>();
    for fragment in fragments {
        // Filter to only get the fragments of interest
        if matches!(
            fragment.ion,
            FragmentType::b(_, 0) | FragmentType::y(_, 0) | FragmentType::Precursor
        ) && fragment.neutral_loss.is_empty()
            && fragment.isotope.is_empty()
            && let Some(mz) = fragment.mz(MassMode::Monoisotopic)
        {
            let low = mz.value - args.before_fragment / fragment.charge.value as f64;
            let high = mz.value + args.after_fragment / fragment.charge.value as f64;
            let sub_spectrum = spectrum.between(low, high, mzdata::prelude::Tolerance::Da(0.0));

            for m in &mods_present {
                let key = (fragment.ion.kind(), m.clone());
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
    }
    // Get start
    let sub_spectrum = spectrum.between(0.0, args.max_start, mzdata::prelude::Tolerance::Da(0.0));

    for m in &mods_present {
        let key = (FragmentKind::diagnostic, m.clone());
        merge_stack(
            stack.fragments.entry(key).or_default(),
            sub_spectrum,
            1,
            0.0,
            args.resolution,
            0.0..=args.max_start,
        );
    }
}

fn merge_stack(
    points: &mut Vec<Point>,
    slice: &[AnnotatedPeak<Fragment>],
    charge: isize,
    center: f64,
    resolution: f64,
    range: RangeInclusive<f64>,
) {
    for found_peak in slice {
        if !range.contains(&found_peak.mz.value) {
            continue;
        }
        let normalised_mass =
            ((found_peak.mz.value - center) / resolution).round() * resolution * charge as f64;
        match points.binary_search_by(|p| p.mass.total_cmp(&normalised_mass)) {
            Ok(index) => {
                points[index].count += 1;
                points[index].total_intensity += f64::from(found_peak.intensity);
            }
            Err(index) => {
                points.insert(index, Point {
                    mass: normalised_mass,
                    count: 1,
                    total_intensity: f64::from(found_peak.intensity),
                });
            }
        }
    }
}

type ItemKey = (FragmentKind, Option<Modification>);

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
                (i, Some(m)) => base_path.join(format!("fragment_{i}_{m}.csv")),
                (i, None) => base_path.join(format!("fragment_{i}.csv")),
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

/// # Panics
/// If the out file could not be created or could not be written to.
fn write_stack(path: &Path, points: &[Point]) {
    let out_file =
        BufWriter::new(File::create(path).unwrap_or_else(|e| {
            panic!("Could not create file: '{}'\n{e}", path.to_string_lossy())
        }));
    mzcore::csv::write_csv(
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
