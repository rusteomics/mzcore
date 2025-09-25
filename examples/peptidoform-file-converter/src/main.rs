//! Convert an identified peptidoform file

use std::{
    fs::File,
    io::{BufWriter, Write},
    path::PathBuf,
};

use clap::Parser;
use context_error::{BasicKind, BoxedError, Context, CreateError};
use directories::ProjectDirs;
use rustyms::{
    identification::{
        FastaIdentifier, MetaData, SpectrumId, SpectrumIds, open_identified_peptidoforms_file,
    },
    sequence::parse_custom_modifications,
};

/// The command line interface arguments
#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Parser)]
struct Cli {
    /// The input file
    #[arg(short, long)]
    in_path: String,
    /// The output path to output the resulting csv file
    #[arg(short, long)]
    out_path: String,
    /// If needed to provide the raw file that was used in the creation of this file
    #[arg(long)]
    raw_file: Option<String>,
    /// To turn off loading the custom modifications database from the Annotator (if installed)
    #[arg(long)]
    no_custom_mods: bool,
}

fn main() {
    let args = Cli::parse();
    let path = ProjectDirs::from("com", "com.snijderlab.annotator", "")
        .expect("Could not generate Annotator configurationpath (needed to check if custom modifications are defined)")
        .config_dir()
        .join("../custom_modifications.json");
    let custom_database = if args.no_custom_mods || !path.exists() {
        None
    } else {
        Some(parse_custom_modifications(&path).expect("Could not parse custom modifications file, if you do not need these you can skip parsing them using the appropriate flag"))
    };
    let mut out_file =
        BufWriter::new(File::create(args.out_path).expect("Could not create out CSV file"));
    writeln!(&mut out_file, "sequence,scan_index,raw_file,z,score,class").unwrap();
    let mut errors = Vec::new();
    let mut warnings = Vec::new();
    for (peptidoform_index, peptidoform) in
        open_identified_peptidoforms_file(&args.in_path, custom_database.as_ref(), false)
            .expect("Invalid input file")
            .enumerate()
    {
        if let Err(error) = peptidoform.and_then(|p| {
            let indices = match p.scans() {
                SpectrumIds::None => Ok(Vec::new()),
                SpectrumIds::FileNotKnown(ids) => ids
                    .iter()
                    .map(|id| match id {
                        SpectrumId::Index(id) => {
                            Ok((*id, args.raw_file.as_ref().map(PathBuf::from).ok_or_else(|| BoxedError::new(BasicKind::Error,
                            "Missing raw file",
                            "This format does not store the raw file so this should be given via the command line arguments",
                            Context::default()
                                .source(args.in_path.clone())
                                .line_index(peptidoform_index as u32),
                        ))?))
                        }
                        SpectrumId::Number(id) => {
                            Ok((id - 1, args.raw_file.as_ref().map(PathBuf::from).ok_or_else(|| BoxedError::new(BasicKind::Error,
                            "Missing raw file",
                            "This format does not store the raw file so this should be given via the command line arguments",
                            Context::default()
                                .source(args.in_path.clone())
                                .line_index(peptidoform_index as u32),
                        ))?))
                        }
                        _ => {
                            warnings.push(BoxedError::new(BasicKind::Warning,
                            "Invalid spectrum id",
                            "Only spectrum indexes and spectrum numbers can be used, an empty scan number is used instead",
                            Context::default()
                                .source(args.in_path.clone())
                                .line_index(peptidoform_index as u32),
                        ));
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
                                    .source(args.in_path.clone())
                                    .line_index(peptidoform_index as u32),
                            )),
                        })
                    })
                    .collect(),
            }?;
            let sequence = p
                .compound_peptidoform_ion().ok_or_else(|| BoxedError::new(BasicKind::Error,
                    "Missing sequence",
                    "This peptidoform misses a sequence",
                    Context::default()
                        .source(args.in_path.clone())
                        .line_index(peptidoform_index as u32),
                ))?;
            let charge = p.charge().map(|v| v.value).unwrap_or_default();
            let score = p.original_confidence().unwrap_or_default();
            let class = p.protein_names().map_or("Unknown",|ids| if ids.iter().any(FastaIdentifier::decoy) {"Decoy"} else {"Target"});
            for (index, raw_file) in indices {
                writeln!(&mut out_file, "\"{sequence}\",{index},\"{}\",{charge},{score},{class}",raw_file.display()).unwrap();
            }
            Ok((p
                .compound_peptidoform_ion()
                .map(std::borrow::Cow::into_owned),))
        }) {
            context_error::combine_error(&mut errors, error, ());
        }
    }
    if errors.is_empty() {
        for e in &warnings {
            println!("{e}");
        }
        println!("No errors, enjoy the new file!");
    } else {
        for e in &errors {
            println!("{e}");
        }
        for e in &warnings {
            println!("{e}");
        }
        println!(
            "Errors were found while parsing the peptidoform file. Output is still generated but all above lines are ignored."
        );
    }
}
