#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::{PsmUtilsPSM, PsmUtilsVersion, test_format};

#[test]
fn punifind() {
    match test_format::<PsmUtilsPSM>(
        BufReader::new(DATA.as_bytes()),
        &mzcore::ontology::STATIC_ONTOLOGIES,
        crate::TestSettings::default(),
        Some(PsmUtilsVersion::V1_5),
    ) {
        Ok(n) => assert_eq!(n, 3),
        Err(e) => {
            println!("{e}");
            panic!("Failed PSMs test");
        }
    }
}

const DATA: &str = r"peptidoform	spectrum_id	run	collection	is_decoy	score	qvalue	pep	precursor_mz	retention_time	ion_mobility	protein_list	rank	source	provenance:alphadia_filename	rescoring:rt_observed	rescoring:mobility_observed	rescoring:mz_observed	rescoring:charge	rescoring:delta_rt
SSYGSSSNDDSYGSSNNDDSYGSSNK/3	71876	LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01		False	136.160126	0.0	0.0	894.33783	2800.518555	1e-06	['P18899']	1	AlphaDIA	C:\Users\ralfg\git\psm_utils\example_files\alphadia.precursors.tsv	2800.518555	1e-06	894.33783	3.0	452.909424
SSQGSSSSTQSAPSETASASK/2	41978	LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01		False	122.27832	0.0	0.0	986.440491	1647.208252	1e-06	['Q9ULU4']	1	AlphaDIA	C:\Users\ralfg\git\psm_utils\example_files\alphadia.precursors.tsv	1647.208252	1e-06	986.440491	2.0	-23.25415
SSQTSGTNEQSSAIVSAR/2	68554	LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01		False	152.012512	0.0	0.0	905.432312	2678.317139	1e-06	['O60763']	1	AlphaDIA	C:\Users\ralfg\git\psm_utils\example_files\alphadia.precursors.tsv	2678.317139	1e-06	905.432312	2.0	31.525879
";
