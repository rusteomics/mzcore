use std::path::Path;

use context_error::*;
use mzcore::{
    ontology::CustomDatabase,
    sequence::{Linked, SemiAmbiguous, SimpleLinear},
};

use crate::identification::*;

// TODO:
// * Merge multiple annotations for the same spectrum (e.g. all candidates peaks export, take care not to lose info on chimeric spectra)
// * Merge identical (or similar?) peptide sequences (for faster processing)

/// Open the selected path and automatically determine the filetype. It will decompress gzipped
/// files automatically.
///
/// # Errors
/// It errors if the filetype could not be determined or if opening the file errors.
pub fn open_identified_peptidoforms_file<'a>(
    path: impl AsRef<Path>,
    custom_database: Option<&'a CustomDatabase>,
    keep_all_columns: bool,
) -> Result<GeneralIdentifiedPeptidoforms<'a>, BoxedError<'static, BasicKind>> {
    let path = path.as_ref();
    let actual_extension = path
        .extension()
        .map(|ex| {
            (ex == "gz")
                .then_some(path)
                .and_then(|p| p.file_stem())
                .and_then(|p| Path::new(p).extension())
                .unwrap_or(ex)
        })
        .map(|ex| ex.to_string_lossy().to_lowercase());
    match actual_extension.as_deref() {
        Some("csv") => PeaksData::parse_file(path, custom_database, keep_all_columns, None)
            .map(IdentifiedPeptidoformIter::into_box)
            .or_else(|pe| {
                NovorData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|ne| (pe, ne))
            })
            .or_else(|(pe, ne)| {
                InstaNovoData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|ie| (pe, ne, ie))
            })
            .or_else(|(pe, ne, ie)| {
                PLinkData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|le| (pe, ne, ie, le))
            }).or_else(|(pe, ne, ie, le)| {
                PowerNovoData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|pne| (pe, ne, ie, le, pne))
            }).or_else(|(pe, ne, ie, le, pne)| {
                PLGSData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|ple| (pe, ne, ie, le, pne, ple))
            }).or_else(|(pe, ne, ie, le, pne, ple)| {
                BasicCSVData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|be| (pe, ne, ie, le, pne, ple, be))
            }).or_else(|(pe, ne, ie, le, pne, ple, be)| {
                PUniFindData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|pue| (pe, ne, ie, le, pne, ple, be, pue))
            }).map_err(|(pe, ne, ie, le, pne, ple, be, pue)| {
                BoxedError::new(BasicKind::Error,
                    "Unknown file format",
                    "Could not be recognised as either a Peaks, Novor, InstaNovo, pLink, PowerNovo, PLGS, pUniFind, or basic file",
                    Context::default().source(path.to_string_lossy()).to_owned(),
                )
                .add_underlying_errors(vec![pe, ne, ie, le, pne, ple, be, pue])
            }),
        Some("tsv") => SageData::parse_file(path, custom_database, keep_all_columns, None)
            .map(IdentifiedPeptidoformIter::into_box)
            .or_else(|se| {
                PepNetData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|pe| (se, pe))
            })
            .or_else(|(se, pe)| {
                MSFraggerData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|mfe| (se, pe, mfe))
            })
            .or_else(|(se, pe, mfe)| {
                ProteoscapeData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|pse| (se, pe, mfe, pse))
            })
            .or_else(|(se, pe, mfe, pse)| {
                NovoBData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|ne| (se, pe, mfe, pse, ne))
            })
            .or_else(|(se, pe, mfe, pse,  ne)| {
                PiPrimeNovoData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|pne| (se, pe, mfe, pse, ne, pne))
            })
            .map_err(|(se, pe, mfe, pse, ne, pne)| {
                BoxedError::new(BasicKind::Error,
                    "Unknown file format",
                    "Could not be recognised a Sage, PepNet, MSFragger, NovoB, π-PrimeNovo or Proteoscape file",
                    Context::default().source(path.to_string_lossy()).to_owned(),
                )
                .add_underlying_errors(vec![se, pe, mfe, pse, ne, pne])
            }),
        Some("psmtsv") => {
            OpairData::parse_file(path, custom_database, keep_all_columns, None).map(IdentifiedPeptidoformIter::into_box)
        }
        Some("fasta" | "fas" | "fa" | "faa" | "mpfa") => FastaData::parse_file(path).map(|peptides| {
            let a: Box<dyn Iterator<Item = Result<IdentifiedPeptidoform<Linked, MaybePeptidoform>, BoxedError<'static, BasicKind>>> + 'a>
                = Box::new(peptides.into_iter().map(|p| Ok(IdentifiedPeptidoform::<SemiAmbiguous, PeptidoformPresent>::from(p).cast())));
            a
        }),
        Some("txt") => {
            MaxQuantData::parse_file(path, custom_database, keep_all_columns, None)
            .map(IdentifiedPeptidoformIter::into_box)
            .or_else(|me| {
                PiHelixNovoData::parse_file(path, custom_database, keep_all_columns, None)
                    .map(IdentifiedPeptidoformIter::into_box)
                    .map_err(|hne| (me, hne))
            })
            .map_err(|(me, he)| {
                BoxedError::new(BasicKind::Error,
                    "Unknown file format",
                    "Could not be recognised as either a MaxQuant or π-HelixNovo file",
                    Context::default().source(path.to_string_lossy()).to_owned(),
                )
                .add_underlying_errors(vec![me, he])
            })
        }
        Some("mztab") => MZTabData::parse_file(path, custom_database).map(|peptides| {
            let a: Box<dyn Iterator<Item = Result<IdentifiedPeptidoform<Linked, MaybePeptidoform>, BoxedError<'static, BasicKind>>> + 'a>
                = Box::new(peptides.into_iter().map(|p| p.map(|p| {
                        IdentifiedPeptidoform::<SimpleLinear, MaybePeptidoform>::from(p).cast()
                    })));
            a
        }),
        Some("deepnovo_denovo") => {
            DeepNovoFamilyData::parse_file(path, custom_database, keep_all_columns, None).map(IdentifiedPeptidoformIter::into_box)
        },
        Some("ssl") => {
            SpectrumSequenceListData::parse_file(path, custom_database, keep_all_columns, None).map(IdentifiedPeptidoformIter::into_box)
        }
        _ => Err(BoxedError::new(BasicKind::Error,
            "Unknown extension",
            "Use CSV, SSL, TSV, TXT, PSMTSV, deepnovo_denovo, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::default().source(path.to_string_lossy()).to_owned(),
        )),
    }
}

#[expect(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use super::*;
    use std::fs::File;
    use std::io::BufReader;

    #[test]
    fn open_sage() {
        match test_format::<SageData>(
            BufReader::new(File::open("src/identification/test_files/sage_v0_14.tsv").unwrap()),
            None,
            true,
            false,
            Some(SageVersion::V0_14),
        ) {
            Ok(n) => assert_eq!(n, 19),
            Err(e) => {
                println!("{e}");
                panic!("Failed identified peptides test");
            }
        }
    }

    #[test]
    fn open_msfragger() {
        match test_format::<MSFraggerData>(
            BufReader::new(File::open("src/identification/test_files/msfragger_v21.tsv").unwrap()),
            None,
            true,
            false,
            Some(MSFraggerVersion::FragPipeV20Or21),
        ) {
            Ok(n) => assert_eq!(n, 19),
            Err(e) => {
                println!("{e}");
                panic!("Failed identified peptides test");
            }
        }
    }
}
