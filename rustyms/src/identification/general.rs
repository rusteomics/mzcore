use std::path::Path;

use super::{
    error::{Context, CustomError},
    ontologies::CustomDatabase,
    DeepNovoFamilyData, FastaData, IdentifiedPeptide, IdentifiedPeptideIter,
    IdentifiedPeptideSource, InstaNovoData, MSFraggerData, MZTabData, MaxQuantData, NovoBData,
    NovorData, OpairData, PLGSData, PLinkData, PeaksData, PepNetData, PowerNovoData, SageData,
    SpectrumSequenceListData,
};

// TODO:
// * Merge multiple annotations for the same spectrum (e.g. all candidates peaks export, take care not to lose info on chimeric spectra)
// * Merge identical (or similar?) peptide sequences (for faster processing)

/// Open the selected path and automatically determine the file type. It will uncompress gzipped
/// files automatically.
///
/// # Errors
/// It errors if the file type could not be determined or if opening the file errors.
pub fn open_identified_peptides_file<'a>(
    path: impl AsRef<Path>,
    custom_database: Option<&'a CustomDatabase>,
) -> Result<Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>, CustomError> {
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
        Some("csv") => PeaksData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|pe| {
                NovorData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|ne| (pe, ne))
            })
            .or_else(|(pe, ne)| {
                InstaNovoData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|ie| (pe, ne, ie))
            })
            .or_else(|(pe, ne, ie)| {
                PLinkData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|le| (pe, ne, ie, le))
            }).or_else(|(pe, ne, ie, le)| {
                PowerNovoData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|pne| (pe, ne, ie, le, pne))
            }).or_else(|(pe, ne, ie, le, pne)| {
                PLGSData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|ple| (pe, ne, ie, le, pne, ple))
            }).map_err(|(pe, ne, ie, le, pne, ple)| {
                CustomError::error(
                    "Unknown file format",
                    "Could not be recognised as either a Peaks, Novor, InstaNovo, pLink, PowerNovo, or PLGS file",
                    Context::show(path.to_string_lossy()),
                )
                .with_underlying_errors(vec![pe, ne, ie, le, pne, ple])
            }),
        Some("tsv") => MSFraggerData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|me| {
                SageData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|se| (me, se))
            })
            .or_else(|(me, se)| {
                PepNetData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|pe| (me, se, pe))
            })
            .map_err(|(me, se, pe)| {
                CustomError::error(
                    "Unknown file format",
                    "Could not be recognised a MSFragger, PepNet or Sage file",
                    Context::show(path.to_string_lossy()),
                )
                .with_underlying_errors(vec![me, se, pe])
            }),
        Some("psmtsv") => {
            OpairData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        }
        Some("fasta") => FastaData::parse_file(path).map(|peptides| {
            Box::new(peptides.into_iter().map(|p| Ok(p.into())))
                as Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>
        }),
        Some("txt") => {
            MaxQuantData::parse_file(path, custom_database)
            .map(IdentifiedPeptideIter::into_box)
            .or_else(|me| {
                NovoBData::parse_file(path, custom_database)
                    .map(IdentifiedPeptideIter::into_box)
                    .map_err(|ne| (me, ne))
            })
            .map_err(|(me, ne)| {
                CustomError::error(
                    "Unknown file format",
                    "Could not be recognised as either a MaxQuant or NovoB file",
                    Context::show(path.to_string_lossy()),
                )
                .with_underlying_errors(vec![me, ne])
            })
        }
        Some("mztab") => MZTabData::parse_file(path, custom_database).map(|peptides| {
            Box::new(peptides.into_iter().map(|p| p.map(Into::into)))
                as Box<dyn Iterator<Item = Result<IdentifiedPeptide, CustomError>> + 'a>
        }),
        Some("deepnovo_denovo") => {
            DeepNovoFamilyData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        },
        Some("ssl") => {
            SpectrumSequenceListData::parse_file(path, custom_database).map(IdentifiedPeptideIter::into_box)
        }
        _ => Err(CustomError::error(
            "Unknown extension",
            "Use CSV, SSL, TSV, TXT, PSMTSV, deepnovo_denovo, or Fasta, or any of these as a gzipped file (eg csv.gz).",
            Context::show(path.to_string_lossy()),
        )),
    }
}

#[allow(clippy::missing_panics_doc)]
#[cfg(test)]
mod tests {
    use super::*;
    use crate::identification::{test_format, MSFraggerVersion, SageVersion};
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
            Some(MSFraggerVersion::V21),
        ) {
            Ok(n) => assert_eq!(n, 19),
            Err(e) => {
                println!("{e}");
                panic!("Failed identified peptides test");
            }
        }
    }
}
