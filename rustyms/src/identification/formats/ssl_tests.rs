#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::identification::{test_format, SpectrumSequenceListData, SpectrumSequenceListVersion};

#[test]
fn cascadia_v0_0_5() {
    match test_format::<SpectrumSequenceListData>(
        BufReader::new(CASCADIA_V0_0_5.as_bytes()),
        None,
        false,
        false,
        Some(SpectrumSequenceListVersion::SSL),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn small_molecule_example() {
    match test_format::<SpectrumSequenceListData>(
        BufReader::new(SMALL_MOLECULES_EXAMPLE.as_bytes()),
        None,
        false,
        false,
        Some(SpectrumSequenceListVersion::SSL),
    ) {
        Ok(n) => assert_eq!(n, 6),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn peptide_example() {
    match test_format::<SpectrumSequenceListData>(
        BufReader::new(PEPTIDE_EXAMPLE.as_bytes()),
        None,
        true,
        false,
        Some(SpectrumSequenceListVersion::SSL),
    ) {
        Ok(n) => assert_eq!(n, 12),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const CASCADIA_V0_0_5: &str = "file	scan	charge	sequence	score-type	score	retention-time	start-time	end-time
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	140	4.0	VNHKPSNTKVDKK	Cascadia Score	0.8716757	13.830939	13.700380273173717	13.961498312641712
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	143	4.0	LANVNHKPSNTKVDK	Cascadia Score	0.9688815	13.73074	13.600180573771862	13.861298613239857
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	284	4.0	LANVNHKPSNTKVDK	Cascadia Score	0.91546655	14.112868	13.982309289286999	14.243427328754994
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	350	4.0	TYLANVNHKPSNTK	Cascadia Score	0.9532166	14.182865	14.052306123088268	14.313424162556263
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	374	4.0	KENAGEDPGLARQAPKPR	Cascadia Score	0.9884338	14.223199	14.092639870952038	14.353757910420033
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	428	4.0	LANVNHKPSNTKVDK	Cascadia Score	0.8682205	14.452167	14.321607537578014	14.58272557704601
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	559	3.0	LTPPSREEMTK	Cascadia Score	0.94662637	14.79058	14.660020776103405	14.9211388155714
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	624	2.0	KDVDQYMTK	Cascadia Score	0.9709009	14.849122	14.718563027690319	14.979681067158314
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	832	3.0	LTPPSREEMTK	Cascadia Score	0.98843235	15.338381	15.207821793864635	15.46893983333263
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	939	2.0	YVDGVEVHNAK	Cascadia Score	0.8673074	15.534778	15.40421862156239	15.665336661030384
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1003	3.0	KAVGGLGKLGKDA	Cascadia Score	0.91004914	15.7533655	15.6228064969286	15.883924536396595
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1064	4.0	HKVYASEVTHQGLSSPVTK	Cascadia Score	0.98483497	15.964101	15.833541817973522	16.094659857441517
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1189	3.0	QAPGKGLESVAR	Cascadia Score	0.91872495	16.191717	16.06115812809315	16.322276167561146
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1195	3.0	LTPPSREEMTK	Cascadia Score	0.9785851	16.206034	16.07547468693104	16.336592726399036
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1212	2.0	KVLPVPQK	Cascadia Score	0.9976683	16.322432	16.191872544597057	16.452990584065052
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1219	3.0	EVTHQGLSSPVTK	Cascadia Score	0.99756336	16.327246	16.196686692546276	16.45780473201427
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1320	2.0	LTVDGVSR	Cascadia Score	0.97489554	16.50676	16.37620062382069	16.637318663288685
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1339	3.0	LTPPSREEMTK	Cascadia Score	0.9637745	16.652891	16.52233213932362	16.783450178791615
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	26981	4.0	THTC[Carbamidomethyl]PPC[Carbamidomethyl]PAPELLGGPSVFLFPPKPK	Cascadia Score	0.85148484	75.06698	74.93641943485585	75.19753747432384
../test_data/test/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp.mzML	1397	4.0	HKVYAGEVTHQGLSSPVTK	Cascadia Score	0.99907994	16.68329	16.552731461833385	16.81384950130138";

const SMALL_MOLECULES_EXAMPLE: &str = r#"file	scan	charge	adduct	inchikey	chemicalformula	moleculename	otherkeys
dexcaf_051017.mzML	01369	-1	[M-H]	ZXPLRDFHBYIQOX-BTBVOZEKSA-N	C24H44O21N0	Glc04Reduced	 
dexcaf_051017.mzML	01639	-1	[M-H]	NBVGBCYERZIRIP-JAMOUWTMSA-N	C30H54O26N0	Glc05Reduced	
dexcaf_051017.mzML	01855	-1	[M-H]	PNHJKLJIDNHXFR-ZGJYWSOBSA-N	C36H64O31N0	Glc06Reduced	
dexcaf_051017.mzML	02029	-1	[M-H]	NVKJDLBVRSXYRE-BMFDHOHESA-N	C42H74O36N0	Glc07Reduced	
dexcaf_051017.mzML	02179	-1	[M-H]	YMRGEPQWJZHXFF-MGQBKJSVSA-N	C48H84O41N0	Glc08Reduced	
dexcaf_051017.mzML	01079	-1	[M-H]	RYYVLZVUVIJVGH-UHFFFAOYSA-N	C8H10N4O2	Caffeine	"InChI:1S/C8H10N4O2/c1-10-4-9-6-5(10)7(13)12(3)8(14)11(6)2/h4H,1-3H3	HMDB:01847	CAS:58-08-2	SMILES:Cn1cnc2n(C)c(=O)n(C)c(=O)c12""#;

const PEPTIDE_EXAMPLE: &str = "file	scan	charge	sequence
demo.ms2	8	3	VGAGAPVYLAAVLEYLAAEVLELAGNAAR
demo.ms2	1806	2	LAESITIEQGK
demo.ms2	2572	2	ELAEDGC[+57.0]SGVEVR
demo.ms2	3088	2	TTAGAVEATSEITEGK
demo.ms2	3266	2	DC[+57.0]EEVGADSNEGGEEEGEEC[+57.0]
demo.ms2	9734	3	IWELEFPEEAADFQQQPVNAQ[-17.0]PQN
demo.ms2	20919	3	VHINIVVIGHVDSGK
../elsewhere/spec.mzXML	00497	2	LKEPAQNTADNAK
../elsewhere/spec.mzXML	00680	2	ALEGPGPGEDAAHSENNPPR
../elsewhere/spec.mzXML	00965	2	FFSHEAEQK
../elsewhere/spec.mzXML	01114	2	C[+57.0]GPSQPLK
../elsewhere/spec.mzXML	01382	2	AVHVQVTDAEAGK";
