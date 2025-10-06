#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::{PiHelixNovoData, PiHelixNovoVersion, test_format};

#[test]
fn pihelixnovo() {
    match test_format::<PiHelixNovoData>(
        BufReader::new(PIHELIXNOVO_V1_1.as_bytes()),
        None,
        true,
        false,
        Some(PiHelixNovoVersion::V1_1),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const PIHELIXNOVO_V1_1: &str = r#"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.5.5.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=5"	LLLLLLLLLLLLLLLHLLLLK	0.43
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.29.29.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=29"	AEELLLLLLLLLLLLLLLLLLLHG+42.011N+0.984N+0.984	0.41
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.55.55.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=55"	A+42.011+43.006-17.027+43.006-17.027EE+43.006-17.027NEE+43.006-17.027-17.027C+57.021+43.006-17.027EC+57.021+43.006-17.027EK	0.49
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.73.73.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=73"	TDAM+15.995LLN+0.984P+43.006-17.027DN+0.984	0.43
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.102.102.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=102"	HLHDDHPEEEK	0.3
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.134.134.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=134"	VAL+43.006-17.027-17.027EEE-17.027+42.011N+0.984	0.45
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.155.155.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=155"	HAEELLLLLLEE-17.027+43.006-17.027LKLLELLEEK	0.24
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.188.188.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=188"	THLHLLLLLLLK	0.3
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.213.213.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=213"	+43.006-17.027-17.027GHGM+15.995GM+15.995DN+0.984	0.29
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.246.246.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=246"	TQ+0.984LLLLLLLLLK	0.39
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.253.253.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=253"	+43.006-17.027-17.027-17.027-17.027-17.027-17.027-17.027-17.027-17.027-17.027N+0.984-17.027N+0.984	0.35
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.273.273.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=273"	N+0.984L+43.006-17.027-17.027EEE-17.027+42.011N+0.984	0.51
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.293.293.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=293"	LH+43.006-17.027M+15.995M+15.995HN+0.984	0.34
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.308.308.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=308"	DHHHHHH+42.011E+42.011-17.027-17.027EEEK	0.61
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.328.328.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=328"	L+42.011H+43.006-17.027C+57.021N+0.984EN+0.984	0.37
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.336.336.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=336"	HM-17.027-17.027-17.027EEE+42.011-17.027N+0.984	0.55
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.354.354.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=354"	Q+0.984+42.011HHLLLLLLK	0.35
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.384.384.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=384"	HEEN+0.984NEE-17.027+42.011N+0.984	0.52
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.405.405.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=405"	VL+43.006-17.027+43.006-17.027AKKK	0.65
20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.424.424.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_L.raw", NativeID:"controllerType=0 controllerNumber=1 scan=424"	I+42.011+43.006-17.027+43.006-17.027-17.027N+0.984-17.027N+0.984-17.027N+0.984	0.34"#;
