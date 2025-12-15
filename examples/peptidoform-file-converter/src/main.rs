//! Convert a PSM file

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use clap::Parser;
use context_error::{BasicKind, BoxedError, Context, CreateError, combine_error};
use mzcore::ontology::Ontologies;
use mzident::{FastaIdentifier, PSMMetaData, SpectrumId, SpectrumIds, open_psm_file};

/// The command line interface arguments
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
struct Cli {
    /// The input file
    #[arg(short, long)]
    in_path: PathBuf,
    /// The output path to output the resulting csv or mztab file
    #[arg(short, long)]
    out_path: PathBuf,
    /// If needed to provide the raw file that was used in the creation of this file
    #[arg(long)]
    raw_file: Option<PathBuf>,
}

fn main() {
    let args = Cli::parse();
    let out_file =
        BufWriter::new(File::create(&args.out_path).expect("Could not create output file"));
    let extension = args
        .out_path
        .extension()
        .map(|e| e.to_string_lossy().to_ascii_lowercase());
    let mut errors = Vec::new();
    let psms = open_psm_file(&args.in_path, &Ontologies::init().0, false)
        .expect("Invalid input file")
        .filter_map(|f| match f {
            Ok(v) => Some(v),
            Err(error) => {
                combine_error(&mut errors, error, ());
                None
            }
        })
        .collect::<Vec<_>>();

    if errors.is_empty() {
        println!("No errors, enjoy the new file!");
    } else {
        for e in &errors {
            println!("{e}");
        }
        println!(
            "Errors were found while parsing the peptidoform file. Output is still generated but all above lines are ignored."
        );
    }

    match extension.as_deref() {
        Some("csv") => save_csv(out_file, &psms, &args.in_path, args.raw_file.as_deref()),
        Some("mztab") => {
            mzident::mztab_writer::MzTabWriter::write::<_, mzident::MzTabProtein>(
                out_file,
                &[],
                &[],
                &psms,
                mzident::mztab_writer::MSRun {
                    location: args.raw_file.unwrap_or_default(),
                    ..Default::default()
                },
            )
            .unwrap();
        }
        Some(_) | None => println!("The output file has to be .csv or .mztab"),
    }
}

fn save_csv<PSM: PSMMetaData>(
    mut out_file: BufWriter<File>,
    psms: &[PSM],
    in_path: &std::path::Path,
    raw_file: Option<&std::path::Path>,
) {
    writeln!(&mut out_file, "sequence,scan_index,raw_file,z,score,class").unwrap();
    let mut warnings = Vec::new();
    for (peptidoform_index, psm) in psms.iter().enumerate() {
        let indices = match psm.scans() {
                SpectrumIds::None => Ok(Vec::new()),
                SpectrumIds::FileNotKnown(ids) => ids
                    .iter()
                    .map(|id| match id {
                        SpectrumId::Index(id) => {
                            Ok((*id, raw_file.map(PathBuf::from).ok_or_else(|| BoxedError::new(BasicKind::Error,
                            "Missing raw file",
                            "This format does not store the raw file so this should be given via the command line arguments",
                            Context::default()
                                .source(in_path.to_string_lossy())
                                .line_index(peptidoform_index as u32),
                        ))?))
                        }
                        SpectrumId::Number(id) => {
                            Ok((id - 1, raw_file.map(PathBuf::from).ok_or_else(|| BoxedError::new(BasicKind::Error,
                            "Missing raw file",
                            "This format does not store the raw file so this should be given via the command line arguments",
                            Context::default()
                                .source(in_path.to_string_lossy())
                                .line_index(peptidoform_index as u32),
                        ))?))
                        }
                        _ => {
                            combine_error(&mut warnings, BoxedError::new(BasicKind::Warning,
                            "Invalid spectrum id",
                            "Only spectrum indexes and spectrum numbers can be used, an empty scan number is used instead",
                            Context::default()
                                .source(in_path.to_string_lossy())
                                .line_index(peptidoform_index as u32),
                        ), ());
                        Ok((0, PathBuf::new()))
                    }
                    })
                    .collect(),
                SpectrumIds::FileKnown(files) => files
                    .iter()
                    .flat_map(|(file, ids)| {
                        ids.iter().map(|id| match id {
                            SpectrumId::Index(id) => Ok((*id, file.clone())),
                            SpectrumId::Number(id) => Ok((id - 1, file.clone())),
                            _ => Err(BoxedError::new(BasicKind::Error,
                                "Invalid spectrum id",
                                "Only spectrum indexes and spectrum numbers can be used",
                                Context::default()
                                    .source(in_path.to_string_lossy())
                                    .line_index(peptidoform_index as u32),
                            )),
                        })
                    })
                    .collect(),
            }.unwrap();
        let sequence = psm
            .compound_peptidoform_ion()
            .ok_or_else(|| {
                BoxedError::new(
                    BasicKind::Error,
                    "Missing sequence",
                    "This peptidoform misses a sequence",
                    Context::default()
                        .source(in_path.to_string_lossy())
                        .line_index(peptidoform_index as u32),
                )
            })
            .unwrap();
        let charge = psm.charge().map(|v| v.value).unwrap_or_default();
        let score = psm
            .original_confidence()
            .map(|(v, _)| v)
            .unwrap_or_default();
        let class = psm.protein_names().map_or("Unknown", |ids| {
            if ids.iter().any(FastaIdentifier::decoy) {
                "Decoy"
            } else {
                "Target"
            }
        });
        for (index, raw_file) in indices {
            writeln!(
                &mut out_file,
                "\"{sequence}\",{index},\"{}\",{charge},{score},{class}",
                raw_file.display()
            )
            .unwrap();
        }
    }

    for e in &warnings {
        println!("{e}");
    }
}
