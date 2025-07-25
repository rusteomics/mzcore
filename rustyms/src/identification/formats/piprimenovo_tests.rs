#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::identification::{PiPrimeNovoData, PiPrimeNovoVersion, test_format};

#[test]
fn piprimenovo() {
    match test_format::<PiPrimeNovoData>(
        BufReader::new(PIPRIMENOVO_V0_1.as_bytes()),
        None,
        true,
        false,
        Some(PiPrimeNovoVersion::V0_1),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const PIPRIMENOVO_V0_1: &str = r#"label	prediction	charge	score
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.39.39.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=39"	VHLLLK	2	3.0026685635675676e-06
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.41.41.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=41"	Q[+0.984]DPDDLDDDDDN[+0.984]	2	1.855278456634312e-14
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.43.43.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=43"	ELEN[+0.984]N[+0.984]N[+0.984]	2	0.00035882205702364445
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.57.57.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=57"	[+43.006-17.027]-+[43.006-17.027]-[17.027]N[+0.984]DN[+0.984]	2	0.00010693447256926447
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.65.65.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=65"	LLLLLLK	2	2.1655953332810896e-06
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.71.71.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=71"	QDQDQN[+0.984]QN[+0.984]	2	4.976878287266118e-08
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.76.76.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=76"	E+[43.006-17.027]LLLLLLLL-[17.027]LLEEEEEK	3	2.374853571177096e-17
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.82.82.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=82"	LLHDDN[+0.984]	2	9.080250310944393e-05
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.93.93.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=93"	DPLLLLHLLELLEPLPLLHLHLHELK	3	1.4704309213157036e-22
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.107.107.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=107"	HDQN[+0.984]PQN[+0.984]QPDQLPQPQN[+0.984]	2	1.4770137017746862e-15
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.121.121.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=121"	LLHM[+15.995]DM[+15.995]K	2	1.1705211363732815e-05
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.124.124.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=124"	TN[+0.984]LELEN[+0.984]-[17.027]-[17.027]N[+0.984]	3	1.0969644392844202e-07
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.133.133.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=133"	M[+15.995]M[+15.995]LLLLDLLELK	3	9.197373418423638e-11
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.136.136.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=136"	VLLLLIDLDDDN[+0.984]	2	1.2713721134899325e-12
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.141.141.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=141"	VEEEEEEEEEEEK	2	8.02034083591252e-09
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.160.160.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=160"	LLDDK	2	1.7250285964109935e-05
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.167.167.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=167"	[-17.027]-QGLM[+15.995]LDN[+0.984]	2	6.090159558880259e-07
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.180.180.2 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=180"	EGGA+[43.006-17.027]LLLN[+0.984]	2	7.738479723684577e-08
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.182.182.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=182"	A+[43.006-17.027]LM[+15.995]LH+[43.006-17.027]LLDLK	3	5.53576628981034e-10
20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.189.189.3 File:"20250515_EX1_UM1_plitt001_SA_EXT00_d5_C.raw", NativeID:"controllerType=0 controllerNumber=1 scan=189"	HAEEEEEN[+0.984]-[17.027]N[+0.984]	3	7.733673328402801e-07"#;
