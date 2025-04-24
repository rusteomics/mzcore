//! Handle MGF reader reading
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use ordered_float::OrderedFloat;
use regex::Regex;
use uom::num_traits::Zero;

use crate::{
    error::{Context, CustomError},
    helper_functions::check_extension,
    spectrum::{PeakSpectrum, RawPeak, RawSpectrum},
    system::{
        charge::e,
        f64::{Mass, MassOverCharge, Time},
        mass::dalton,
        mass_over_charge::mz,
        time::s,
        usize::Charge,
    },
};
use flate2::bufread::GzDecoder;

/// Open a MGF file and return the contained spectra.
///
/// # Errors
/// It returns an error when:
/// * The file could not be opened
/// * Any line in the file could not be read
/// * When any expected number in the file is not a number
/// * When there is only one column (separated by space or tab) on a data row
pub fn open(path: impl AsRef<Path>) -> Result<Vec<RawSpectrum>, CustomError> {
    let path = path.as_ref();
    let file = File::open(path).map_err(|err| {
        CustomError::error(
            "Could not open file",
            format!("Additional info: {err}"),
            Context::show(path.display()),
        )
    })?;
    if check_extension(path, "gz") {
        open_raw(GzDecoder::new(BufReader::new(file)))
    } else {
        open_raw(file)
    }
}

/// Open a MGF file and return the contained spectra. Open it from a raw buffered reader.
///
/// # Errors
/// It returns an error when:
/// * The file could not be opened
/// * Any line in the file could not be read
/// * When any expected number in the file is not a number
/// * When there is only one column (separated by space or tab) on a data row
#[expect(clippy::missing_panics_doc)]
pub fn open_raw<T: std::io::Read>(reader: T) -> Result<Vec<RawSpectrum>, CustomError> {
    let reader = BufReader::new(reader);
    let mut current = RawSpectrum::default();
    let mut output = Vec::new();
    for (line_index, line) in reader.lines().enumerate() {
        let line = line.map_err(|err| {
            CustomError::error(
                "Could not read mgf file",
                format!("Error while reading line: {err}"),
                Context::show(format!("Line number {}", line_index + 1)),
            )
        })?;
        let base_error = CustomError::error(
            "Could not read mgf file",
            "..",
            Context::full_line(line_index, line.clone()),
        );
        match line.as_str() {
            "BEGIN IONS" | "" => (),
            "END IONS" => {
                output.push(current);
                current = RawSpectrum::default();
            }
            t if t.contains('=') => {
                // THe previous line made sure it will always contain an equals sign
                let (key, value) = t.split_once('=').unwrap();
                match key {
                    "PEPMASS" => match value.split_once(' ') {
                        None => {
                            current.mass =
                                Some(Mass::new::<dalton>(value.parse().map_err(|_| {
                                    base_error.with_long_description(format!(
                                        "Not a number {key} for PEPMASS"
                                    ))
                                })?));
                        }
                        Some((mass, intensity)) => {
                            current.mass =
                                Some(Mass::new::<dalton>(mass.parse().map_err(|_| {
                                    base_error.with_long_description(format!(
                                        "Not a number {key} for PEPMASS"
                                    ))
                                })?));
                            current.intensity = Some(intensity.parse().map_err(|_| {
                                base_error.with_long_description(format!(
                                    "Not a number {key} for PEPMASS"
                                ))
                            })?);
                        }
                    },
                    "CHARGE" => {
                        current.charge = Some(parse_charge(value).map_err(|()| {
                            base_error
                                .with_long_description(format!("Not a number {key} for CHARGE"))
                        })?);
                    }
                    "RT" => {
                        current.rt = Some(Time::new::<s>(value.parse().map_err(|_| {
                            base_error.with_long_description(format!("Not a number {key} for RT"))
                        })?));
                    }
                    "RTINSECONDS" => {
                        current.rt = Some(Time::new::<s>(value.parse().map_err(|_| {
                            base_error.with_long_description(format!("Not a number {key} for RT"))
                        })?));
                    }
                    "TITLE" => parse_title(value, &mut current),
                    "SEQUENCE" => current.sequence = Some(value.to_owned()),
                    "NUM_SCANS" => {
                        current.num_scans = value.parse().map_err(|_| {
                            base_error
                                .with_long_description(format!("Not a number {key} for NUM_SCANS"))
                        })?;
                    }
                    _ => (),
                }
            }
            t if t.contains(' ') || t.contains('\t') => {
                let split = if t.contains(' ') {
                    t.split(' ').collect::<Vec<_>>()
                } else {
                    t.split('\t').collect::<Vec<_>>()
                };
                let mut peak = RawPeak {
                    mz: MassOverCharge::zero(),
                    intensity: OrderedFloat(0.0),
                };
                if split.len() < 2 {
                    return Err(base_error.with_long_description("Not enough columns"));
                }
                peak.mz = MassOverCharge::new::<mz>(split[0].parse().map_err(|_| {
                    base_error.with_long_description(format!("Not a number {} for MZ", split[0]))
                })?);
                peak.intensity = split[1].parse().map_err(|_| {
                    base_error
                        .with_long_description(format!("Not a number {} for INTENSITY", split[1]))
                })?;
                if split.len() >= 3 {
                    _ = parse_charge(split[2]).map_err(|()| {
                        base_error
                            .with_long_description(format!("Not a number {} for CHARGE", split[2]))
                    })?;
                }
                current.add_peak(peak);
            }
            _ => {}
        }
    }
    Ok(output)
}

/// # Errors
/// When the charge could not be properly parsed. For example if it has a negative charge.
fn parse_charge(input: &str) -> Result<Charge, ()> {
    if input.ends_with('+') {
        Ok(Charge::new::<e>(
            input.trim_end_matches('+').parse().map_err(|_| ())?,
        ))
    } else if input.ends_with('-') {
        Err(())
    } else {
        Ok(Charge::new::<e>(input.parse().map_err(|_| ())?))
    }
}

#[expect(clippy::missing_panics_doc)]
fn parse_title(title: &str, spectrum: &mut RawSpectrum) {
    // basic structure: <name>.<scan>.<scan>.<experiment?>? File:"<name>", NativeID:"(<header>) +"
    let ms_convert_format: Regex =
        Regex::new(r#"(.+)\.(\d+)\.\d+\.\d* File:".*", NativeID:"(.+)""#).unwrap();
    // other structure: <name>.ScanId;v=<num>;d1=<scan>.<scan>.<experiment?>_INDEX<index>
    let other_format: Regex =
        Regex::new(r"(.+)\.ScanId;v=\d+;d1=(\d+)\.\d+\.\d*_INDEX(\d+)").unwrap();

    spectrum.title = title.to_string();
    if let Some(ms_convert) = ms_convert_format.captures(title) {
        spectrum.raw_file = Some(ms_convert[1].to_string());
        spectrum.raw_scan_number = ms_convert[2].parse().ok(); // By definition will always work thanks to the regex
        for header in ms_convert[3].split(' ') {
            match header.split_once('=') {
                Some(("sample", n)) => spectrum.sample = n.parse().ok(),
                Some(("period", n)) => spectrum.period = n.parse().ok(),
                Some(("cycle", n)) => spectrum.cycle = n.parse().ok(),
                Some(("experiment", n)) => spectrum.experiment = n.parse().ok(),
                Some(("controllerType", n)) => spectrum.controller_type = n.parse().ok(),
                Some(("controllerNumber", n)) => spectrum.controller_number = n.parse().ok(),
                None | Some(_) => (),
            }
        }
    } else if let Some(other) = other_format.captures(title) {
        spectrum.raw_file = Some(other[1].to_string());
        spectrum.raw_scan_number = other[2].parse().ok(); // By definition will always work thanks to the regex
        spectrum.raw_index = other[3].parse().ok(); // By definition will always work thanks to the regex
    }
    // Else just ignore
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::*;
    #[test]
    fn test_open() {
        let spectra =
            open(std::env::var("CARGO_MANIFEST_DIR").unwrap() + "/data/example.mgf").unwrap();
        assert_eq!(spectra.len(), 1);
        assert_eq!(spectra[0].spectrum().len(), 5);
        assert!(spectra[0][0].mz < spectra[0][1].mz);
    }

    #[test]
    fn test_titles() {
        assert_eq!(
            test_title_helper(
                r#"230629 1-Trypsin CID 06.93060.93060. File:"", NativeID:"sample=1 period=1 cycle=2024 experiment=2""#
            ),
            (
                Some("230629 1-Trypsin CID 06".to_string()),
                Some(93060),
                Some(1),
                Some(1),
                Some(2024),
                Some(2),
                None,
                None,
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"IgAmix_10KE_CE25_02-IgA mix 10KE 7000nA CE 25.130.130. File:"IgAmix_10KE_CE25_02.wiff", NativeID:"sample=1 period=1 cycle=1966 experiment=2""#
            ),
            (
                Some("IgAmix_10KE_CE25_02-IgA mix 10KE 7000nA CE 25".to_string()),
                Some(130),
                Some(1),
                Some(1),
                Some(1966),
                Some(2),
                None,
                None,
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"230629 1-Trypsin EAD 02.52990.52990. File:"", NativeID:"sample=1 period=1 cycle=2039 experiment=2""#
            ),
            (
                Some("230629 1-Trypsin EAD 02".to_string()),
                Some(52990),
                Some(1),
                Some(1),
                Some(2039),
                Some(2),
                None,
                None,
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"26072023_MS2_TZB_1µm_HiRes_Isol1284mzNarrow_ETD_RT10ms_1e5_IT200 Scan[1].1.1.3 File:"26072023_MS2_TZB_1µm_HiRes_Isol1284mzNarrow_ETD_RT10ms_1e5_IT200 Scan[1].raw", NativeID:"controllerType=0 controllerNumber=1 scan=1""#
            ),
            (
                Some(
                    "26072023_MS2_TZB_1µm_HiRes_Isol1284mzNarrow_ETD_RT10ms_1e5_IT200 Scan[1]"
                        .to_string()
                ),
                Some(1),
                None,
                None,
                None,
                None,
                Some(0),
                Some(1),
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"20230505_L1_UM3_5858593_SA_EXT00_totalIgAplasma_GT_ETHCD_20230513002852.990.990.2 File:"20230505_L1_UM3_5858593_SA_EXT00_totalIgAplasma_GT_ETHCD_20230513002852.raw", NativeID:"controllerType=0 controllerNumber=1 scan=990""#
            ),
            (
                Some(
                    "20230505_L1_UM3_5858593_SA_EXT00_totalIgAplasma_GT_ETHCD_20230513002852"
                        .to_string()
                ),
                Some(990),
                None,
                None,
                None,
                None,
                Some(0),
                Some(1),
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"20210728_L1_Vq1_Tamar002_MGUS_Fab_TCEP_ETHCD_22min_rep1_01.36.36.15 File:"20210728_L1_Vq1_Tamar002_MGUS_Fab_TCEP_ETHCD_22min_rep1_01.raw", NativeID:"controllerType=0 controllerNumber=1 scan=36""#
            ),
            (
                Some("20210728_L1_Vq1_Tamar002_MGUS_Fab_TCEP_ETHCD_22min_rep1_01".to_string()),
                Some(36),
                None,
                None,
                None,
                None,
                Some(0),
                Some(1),
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"20191211_F1_Ag5_peng0013_SA_her_Arg_C.2824.2824.3 File:"20191211_F1_Ag5_peng0013_SA_her_Arg_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=2824""#
            ),
            (
                Some("20191211_F1_Ag5_peng0013_SA_her_Arg_C".to_string()),
                Some(2824),
                None,
                None,
                None,
                None,
                Some(0),
                Some(1),
                None
            )
        );
        assert_eq!(
            test_title_helper(
                r#"20190502_F1_Ag5_3117030_SA_Herceptin_tryp_HCD_01.188.188.3 File:"20190502_F1_Ag5_3117030_SA_Herceptin_tryp_HCD_01.raw", NativeID:"controllerType=0 controllerNumber=1 scan=188""#
            ),
            (
                Some("20190502_F1_Ag5_3117030_SA_Herceptin_tryp_HCD_01".to_string()),
                Some(188),
                None,
                None,
                None,
                None,
                Some(0),
                Some(1),
                None
            )
        );
        assert_eq!(
            test_title_helper(
                "20230512_l1_um3_rodri078_sa_ext00_pro1_1_50_4h_tr.ScanId;v=1;d1=782.782.2_INDEX0"
            ),
            (
                Some("20230512_l1_um3_rodri078_sa_ext00_pro1_1_50_4h_tr".to_string()),
                Some(782),
                None,
                None,
                None,
                None,
                None,
                None,
                Some(0)
            )
        );
        assert_eq!(
            test_title_helper(
                "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a.ScanId;v=1;d1=3120.3120.3_INDEX0"
            ),
            (
                Some(
                    "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a"
                        .to_string()
                ),
                Some(3120),
                None,
                None,
                None,
                None,
                None,
                None,
                Some(0)
            )
        );
        assert_eq!(
            test_title_helper(
                "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a.ScanId;v=1;d1=10287.10287.4_INDEX6510"
            ),
            (
                Some(
                    "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a"
                        .to_string()
                ),
                Some(10287),
                None,
                None,
                None,
                None,
                None,
                None,
                Some(6510)
            )
        );
        assert_eq!(
            test_title_helper(
                "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a.ScanId;v=1;d1=14660.14660.3_INDEX10481"
            ),
            (
                Some(
                    "20230512_l1_um3_rodri078_sa_ext00_pro1_1_200_1h_a"
                        .to_string()
                ),
                Some(14660),
                None,
                None,
                None,
                None,
                None,
                None,
                Some(10481)
            )
        );
    }

    #[expect(clippy::type_complexity)]
    fn test_title_helper(
        title: &str,
    ) -> (
        Option<String>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
        Option<usize>,
    ) {
        let mut spectrum = RawSpectrum::default();
        parse_title(title, &mut spectrum);
        (
            spectrum.raw_file,
            spectrum.raw_scan_number,
            spectrum.sample,
            spectrum.period,
            spectrum.cycle,
            spectrum.experiment,
            spectrum.controller_type,
            spectrum.controller_number,
            spectrum.raw_index,
        )
    }
}
