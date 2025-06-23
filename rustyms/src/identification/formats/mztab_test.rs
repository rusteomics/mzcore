#![allow(clippy::missing_panics_doc)]
use std::io::{BufRead, BufReader};

use crate::{
    error::CustomError,
    identification::{IdentifiedPeptidoform, MZTabData, test_identified_peptidoform},
    sequence::SemiAmbiguous,
};

#[test]
fn pride_exp_excerpt_ac_1643() {
    assert_eq!(
        open_file(BufReader::new(PRIDE_EXP_EXCERPT_AC_1643.as_bytes())).unwrap(),
        25
    );
}

#[test]
fn pride_exp_excerpt_ac_16649() {
    assert_eq!(
        open_file(BufReader::new(PRIDE_EXP_EXCERPT_AC_16649.as_bytes())).unwrap(),
        15
    );
}

#[test]
fn silac_cqi() {
    assert_eq!(open_file(BufReader::new(SILAC_CQI.as_bytes())).unwrap(), 30);
}

#[test]
fn itraq_cqi() {
    assert_eq!(open_file(BufReader::new(ITRAQ_CQI.as_bytes())).unwrap(), 36);
}

#[test]
fn itraq_sqi() {
    assert_eq!(open_file(BufReader::new(ITRAQ_SQI.as_bytes())).unwrap(), 28);
}

#[test]
fn labelfree_cqi() {
    assert_eq!(
        open_file(BufReader::new(LABELFREE_CQI.as_bytes())).unwrap(),
        58
    );
}

#[test]
fn labelfree_sqi() {
    assert_eq!(
        open_file(BufReader::new(LABELFREE_SQI.as_bytes())).unwrap(),
        58
    );
}

#[test]
fn casanovo_v3_2_0() {
    assert_eq!(
        open_file(BufReader::new(CASANOVO_V3_2_0_A.as_bytes())).unwrap(),
        41
    );
    assert_eq!(
        open_file(BufReader::new(CASANOVO_V3_2_0_B.as_bytes())).unwrap(),
        41
    );
}

#[test]
fn casanovo_v4_2_1() {
    assert_eq!(
        open_file(BufReader::new(CASANOVO_V4_2_1.as_bytes())).unwrap(),
        39
    );
}

#[test]
fn contranovo_v1_0_0() {
    assert_eq!(
        open_file(BufReader::new(CONTRANOVO_V1_0_0.as_bytes())).unwrap(),
        20
    );
}

/// Open a MZTab file from the given reader.
/// # Errors
/// If any part of the process errors.
fn open_file(reader: impl BufRead) -> Result<usize, CustomError> {
    let mut peptides = 0;
    for read in MZTabData::parse_reader(reader, None) {
        let peptide: IdentifiedPeptidoform<SemiAmbiguous> = read?.into();
        peptides += 1;

        test_identified_peptidoform(&peptide, false, false).unwrap();
    }
    Ok(peptides)
}

const PRIDE_EXP_EXCERPT_AC_1643: &str = r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	title	COFRADIC N-terminal proteome of unstimulated human blood platelets, identified and unidentified spectra
MTD	description	date of export: Mon Jun 16 11:52:25 BST 2014
MTD	instrument[1]-name	[PRIDE, PRIDE:0000131, Instrument model, Micromass Q-TOF I]
MTD	instrument[1]-source	[PSI, PSI:1000008, Ionization Type, ESI]
MTD	instrument[1]-analyzer[1]	[PSI, PSI:1000010, Analyzer Type, Quadrupole-TOF]
MTD	instrument[1]-detector	[PSI, PSI:1000026, Detector Type, MultiChannelPlate]
MTD	software[1]	[MS, MS:1001456, analysis software, MassLynx v3.5]
MTD	protein_search_engine_score[1]	[MS, MS:1001171, Mascot:score, ]
MTD	psm_search_engine_score[1]	[MS, MS:1001153, search engine specific score, ]
MTD	publication[1]	pubmed:16038019
MTD	publication[2]	pubmed:12665801
MTD	publication[3]	pubmed:16518876
MTD	contact[1]-name	Kristian Flikka
MTD	contact[1]-affiliation	Computational Biology Unit, Bergen Center for Computational Science, University of Bergen
MTD	contact[1]-email	flikka@ii.uib.no
MTD	contact[2]-name	Lennart Martens
MTD	contact[2]-affiliation	Department of Medical Protein Research (GE07, VIB09), Faculty of Medicine and Health Sciences, Ghent University and Flanders Interuniversitary Institute for Biotechnology (VIB)
MTD	contact[2]-email	lennart.martens@UGent.be
MTD	uri[1]	http://www.ebi.ac.uk/pride/archive/assays/1643
MTD	fixed_mod[1]	[MS, MS:1002453, No fixed modifications searched, ]
MTD	variable_mod[1]	[MS, MS:1002454, No variable modifications searched, ]
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]
MTD	sample[1]-species[1]	[NEWT, 9606, Homo sapiens (Human), ]
MTD	sample[1]-cell_type[1]	[CL, CL:0000233, platelet, ]
MTD	sample[1]-description	Unstimulated human blood platelets
MTD	sample[1]-custom[1]	[MeSH, D001792, blood_platelets, ]
MTD	assay[1]-sample_ref	sample[1]
MTD	assay[1]-ms_run_ref	ms_run[1]

COM	Only variable modifications can be reported when the original source is a PRIDE XML file

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	search_engine_score[1]_ms_run[1]	num_psms_ms_run[1]	num_peptides_distinct_ms_run[1]	num_peptides_unique_ms_run[1]	ambiguity_members	modifications	protein_coverage
PRT	IPI00025512	null	9606	Homo sapiens (Human)	IPI human	2.31	[MS, MS:1001207, Mascot, ]	60.67	60.67	3	1	null	null	null	0.0
PRT	IPI00298497	null	9606	Homo sapiens (Human)	IPI human	2.31	[MS, MS:1001207, Mascot, ]	69.37	69.37	49	10	null	null	null	0.0

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end
PSM	LFDQAFGLPR	17699	IPI00025512	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=17699	null	null	28	37
PSM	LFDQAFGLPR	17815	IPI00025512	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=17815	null	null	28	37
PSM	LFDQAFGLPR	17587	IPI00025512	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=17587	null	null	28	37
PSM	VNDNEEGFFSAR	5927	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=5927	null	null	33	44
PSM	LRPAPPPISGGGYR	14345	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14345	null	null	59	72
PSM	QGVNDNEEGFFSAR	14083	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14083	null	null	31	44
PSM	KIQKLESDVSAQMEYCR	2689	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=2689	null	null	208	224
PSM	QGVNDNEEGFFSAR	14423	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14423	null	null	31	44
PSM	QGVNDNEEGFFSAR	7451	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=7451	null	null	31	44
PSM	QGVNDNEEGFFSAR	14219	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14219	null	null	31	44
PSM	QGVNDNEEGFFSAR	7615	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=7615	null	null	31	44
PSM	QGVNDNEEGFFSAR	14461	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14461	null	null	31	44
PSM	AAATQKKVER	12129	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=12129	null	null	78	87
PSM	LRPAPPPISGGGYR	7069	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=7069	null	null	59	72
PSM	LRPAPPPISGGGYR	14193	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14193	null	null	59	72
PSM	QGVNDNEEGFFSAR	6145	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=6145	null	null	31	44
PSM	QGVNDNEEGFFSAR	6045	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=6045	null	null	31	44
PSM	LRPAPPPISGGGYR	7223	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=7223	null	null	59	72
PSM	QGVNDNEEGFFSAR	7677	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=7677	null	null	31	44
PSM	LRPAPPPISGGGYR	14495	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14495	null	null	59	72
PSM	GHRPLDKKR	12699	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=12699	null	null	45	53
PSM	QGVNDNEEGFFSAR	14269	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14269	null	null	31	44
PSM	QGVNDNEEGFFSAR	15077	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=15077	null	null	31	44
PSM	QGVNDNEEGFFSAR	14371	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14371	null	null	31	44
PSM	EAPSLRPAPPPISGGGYR	14059	IPI00298497	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[1]:spectrum=14059	null	null	55	72";

const PRIDE_EXP_EXCERPT_AC_16649: &str = r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Summary
MTD	mzTab-type	Quantification
MTD	mzTab-ID	16649
MTD	title	The Synaptic Proteome during Development and Plasticity of the Mouse Visual Cortex
MTD	description	date of export: Mon Jun 16 10:57:56 BST 2014
MTD	sample_processing[1]	[PRIDE, PRIDE:0000025, Reduction, DTT]
MTD	sample_processing[2]	[MOD, MOD:00110, L-cysteine methyl disulfide, ]
MTD	sample_processing[3]	[PRIDE, PRIDE:0000160, Enzyme, Trypsin]
MTD	instrument[1]-name	[PRIDE, PRIDE:0000131, Instrument model, ABI 4800]
MTD	instrument[1]-source	[MS, MS:1000075, matrix-assisted laser desorption ionization, ]
MTD	instrument[1]-analyzer[1]	[MS, MS:1000084, time-of-flight, ]
MTD	instrument[1]-detector	[MS, MS:1000116, photomultiplier, ]
MTD	software[1]	[MS, MS:1001456, analysis software, Matrix Science Mascot v2.3.01]
MTD	protein_search_engine_score[1]	[MS, MS:1001171, Mascot:score, ]
MTD	psm_search_engine_score[1]	[PRIDE, PRIDE:0000069, Mascot score, ]
MTD	publication[1]	pubmed:21398567
MTD	contact[1]-name	August B. Smit
MTD	contact[1]-affiliation	Department of Molecular and Cellular Neurobiology, Center for Neurogenomics and Cognitive Research, Neuroscience Campus Amsterdam, VU University, De Boelelaan 1085, 1081 HV Amsterdam, The Netherlands
MTD	contact[1]-email	guus.smit@cncr.vu.nl
MTD	contact[2]-name	Christiaan N. Levelt
MTD	contact[2]-affiliation	Molecular Visual Plasticity group, Netherlands Institute for Neuroscience, an institute of the Royal Netherlands Academy of Arts and Sciences, Meibergdreef 47, 1105 BA Amsterdam, The Netherlands
MTD	contact[2]-email	c.levelt@nin.knaw.nl
MTD	contact[3]-name	Ka Wan Li
MTD	contact[3]-affiliation	Department of Molecular and Cellular Neurobiology, Center for Neurogenomics and Cognitive Research, Neuroscience Campus Amsterdam, VU University, De Boelelaan 1085, 1081 HV Amsterdam, The Netherlands
MTD	contact[3]-email	ka.wan.li@cncr.vu.nl
MTD	contact[4]-name	Pim van Nierop
MTD	contact[4]-affiliation	Department of Molecular and Cellular Neurobiology, Center for Neurogenomics and Cognitive Research, Neuroscience Campus Amsterdam, VU University, De Boelelaan 1085, 1081 HV Amsterdam, The Netherlands
MTD	contact[4]-email	pim.van.nierop@cncr.vu.nl
MTD	uri[1]	http://www.ebi.ac.uk/pride/archive/assays/16649
MTD	fixed_mod[1]	[MS, MS:1002453, No fixed modifications searched, ]
MTD	variable_mod[1]	[MOD, MOD:01499, iTRAQ4plex-116 reporter+balance reagent acylated residue, 144.102062]
MTD	variable_mod[2]	[MOD, MOD:00425, monohydroxylated residue, 15.994915]
MTD	quantification_method	[PRIDE, PRIDE:0000313, iTRAQ, ]
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000395, Ratio, ]
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2011/03/PRD000397/PRIDE_Exp_Complete_Ac_16649.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]
MTD	sample[1]-species[1]	[NEWT, 10090, Mus musculus (Mouse), subsample1]
MTD	sample[1]-description	Biological replicate 1, P30
MTD	sample[1]-custom[1]	[GO, GO:0044456, synapse part, subsample1]
MTD	sample[1]-custom[2]	[EFO, EFO:0000916, visual cortex, subsample1]
MTD	sample[2]-species[1]	[NEWT, 10090, Mus musculus (Mouse), subsample2]
MTD	sample[2]-description	Biological replicate 1, P30 right eye occlusion
MTD	sample[2]-custom[1]	[GO, GO:0044456, synapse part, subsample2]
MTD	sample[2]-custom[2]	[EFO, EFO:0000916, visual cortex, subsample2]
MTD	sample[3]-species[1]	[NEWT, 10090, Mus musculus (Mouse), subsample3]
MTD	sample[3]-description	Biological replicate 1, P46
MTD	sample[3]-custom[1]	[GO, GO:0044456, synapse part, subsample3]
MTD	sample[3]-custom[2]	[EFO, EFO:0000916, visual cortex, subsample3]
MTD	sample[4]-species[1]	[NEWT, 10090, Mus musculus (Mouse), subsample4]
MTD	sample[4]-description	Biological replicate 1, P46 dark rearing
MTD	sample[4]-custom[1]	[GO, GO:0044456, synapse part, subsample4]
MTD	sample[4]-custom[2]	[EFO, EFO:0000916, visual cortex, subsample4]
MTD	assay[1]-quantification_reagent	[PRIDE, PRIDE:0000114, iTRAQ reagent 114, subsample1]
MTD	assay[1]-sample_ref	sample[1]
MTD	assay[1]-ms_run_ref	ms_run[1]
MTD	assay[2]-quantification_reagent	[PRIDE, PRIDE:0000115, iTRAQ reagent 115, subsample2]
MTD	assay[2]-sample_ref	sample[2]
MTD	assay[2]-ms_run_ref	ms_run[1]
MTD	assay[3]-quantification_reagent	[PRIDE, PRIDE:0000116, iTRAQ reagent 116, subsample3]
MTD	assay[3]-sample_ref	sample[3]
MTD	assay[3]-ms_run_ref	ms_run[1]
MTD	assay[4]-quantification_reagent	[PRIDE, PRIDE:0000117, iTRAQ reagent 117, subsample4]
MTD	assay[4]-sample_ref	sample[4]
MTD	assay[4]-ms_run_ref	ms_run[1]

COM	Only variable modifications can be reported when the original source is a PRIDE XML file

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	search_engine_score[1]_ms_run[1]	num_psms_ms_run[1]	num_peptides_distinct_ms_run[1]	num_peptides_unique_ms_run[1]	ambiguity_members	modifications	protein_coverage	protein_abundance_assay[1]	protein_abundance_assay[2]	protein_abundance_assay[3]	protein_abundance_assay[4]
PRT	223462890	Spna2 protein [Mus musculus]	10090	Mus musculus (Mouse)	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	6539.67	6539.67	157	92	null	null	null	0	1	0.853	0.864	0.791

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end
PSM	QQVLDR	1661	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	37.76	0-MOD:01499	null	1	902.482117	902.518133	ms_run[1]:spectrum=1661	R	Y	20	25
PSM	LVQYLR	2280	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	44.64	0-MOD:01499	null	1	935.577454	935.580003	ms_run[1]:spectrum=2280	K	E	151	156
PSM	LVQYLR	2281	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	44.76	0-MOD:01499	null	1	935.583252	935.580003	ms_run[1]:spectrum=2281	K	E	151	156
PSM	LQQLFR	2537	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	45.41	0-MOD:01499	null	1	948.595581	948.575253	ms_run[1]:spectrum=2537	R	D	786	791
PSM	EAGSVSLR	2809	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	55.05	0-MOD:01499	null	1	962.509827	962.539253	ms_run[1]:spectrum=2809	K	M	1058	1065
PSM	LSILSEER	5465	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	39.82	0-MOD:01499	null	1	1090.623169	1090.622983	ms_run[1]:spectrum=5465	K	T	442	449
PSM	IDGITIQAR	6333	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	65.64	0-MOD:01499	null	1	1130.69458	1130.665513	ms_run[1]:spectrum=6333	R	Q	734	742
PSM	LFGAAEVQR	6418	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	57.63	0-MOD:01499	null	1	1134.634888	1134.639293	ms_run[1]:spectrum=6418	K	F	251	259
PSM	LFGAAEVQR	6420	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	48.86	0-MOD:01499	null	1	1134.64624	1134.639293	ms_run[1]:spectrum=6420	K	F	251	259
PSM	DLTGVQNLR	6915	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	49.74	0-MOD:01499	null	1	1159.643921	1159.655683	ms_run[1]:spectrum=6915	R	K	1800	1808
PSM	LGDSHDLQR	7536	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	36.23	0-MOD:01499	null	1	1184.487671	1184.614543	ms_run[1]:spectrum=7536	K	F	1335	1343
PSM	KQEALVAR	7964	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	47.21	0-MOD:01499,1-MOD:01499	null	1	1202.668945	1202.746446	ms_run[1]:spectrum=7964	K	Y	758	765
PSM	QGQIDNQTR	7979	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	53.35	0-MOD:01499	null	1	1203.583984	1203.620373	ms_run[1]:spectrum=7979	R	I	1046	1054
PSM	EFSMMFK	8068	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	42.5	0-MOD:01499,7-MOD:01499	null	1	1207.596069	1207.609896	ms_run[1]:spectrum=8068	K	H	2332	2338
PSM	SLQQLAEER	8285	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	73.16	0-MOD:01499	null	1	1217.661011	1217.661163	ms_run[1]:spectrum=8285	R	S	1217	1225";

const SILAC_CQI: &str = r#"COM	This	line	serves	as	a	size	and	separator	hint	for	spreadsheet	applications.	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
COM	Report of a "Complete Quantification report" SILAC experiment, quantification on 2 study variables (control/treatment), 3+3 assays (replicates) reported, identifications reported.
MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Complete
MTD	mzTab-type	Quantification
MTD	description	mzTab example file for reporting a summary report of quantification data quantified on the protein level
MTD	ms_run[1]-location	file://C:/path/to/my/file1.mzML
MTD	ms_run[2]-location	file://C:/path/to/my/file2.mzML
MTD	ms_run[3]-location	file://C:/path/to/my/file3.mzML
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000393, Relative quantification unit,]
MTD	protein_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	psm_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	software[1]	[MS, MS:1001583, MaxQuant,]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
COM	Specifiying site and position of search modifications is not mandatory but good practice
MTD	fixed_mod[1]-site	C
MTD	fixed_mod[1]-position	Anywhere
MTD	variable_mod[1]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[1]-site	M
MTD	variable_mod[1]-position	Anywhere
MTD	quantification_method	[MS, MS:1001835, SILAC, ]
MTD	assay[1]-quantification_reagent	[PRIDE, PRIDE:0000326, SILAC light,]
MTD	assay[2]-quantification_reagent	[PRIDE, PRIDE:0000325, SILAC heavy,]
MTD	assay[3]-quantification_reagent	[PRIDE, PRIDE:0000326, SILAC light,]
MTD	assay[4]-quantification_reagent	[PRIDE, PRIDE:0000325, SILAC heavy,]
MTD	assay[5]-quantification_reagent	[PRIDE, PRIDE:0000326, SILAC light,]
MTD	assay[6]-quantification_reagent	[PRIDE, PRIDE:0000325, SILAC heavy,]
MTD	assay[2]-quantification_mod[1]	[UNIMOD, UNIMOD:267, Label:13C(6)15N(4),]
COM	Specifiying details of quantification modifications is not mandatory but good practice
MTD	assay[2]-quantification_mod[1]-site	R
MTD	assay[2]-quantification_mod[1]-position	Anywhere
MTD	assay[2]-quantification_mod[2]	[UNIMOD, UNIMOD:259, Label:13C(6)15N(2),]
MTD	assay[2]-quantification_mod[2]-site	L
MTD	assay[2]-quantification_mod[2]-position	Anywhere
MTD	assay[4]-quantification_mod[1]	[UNIMOD, UNIMOD:267, Label:13C(6)15N(4),]
MTD	assay[4]-quantification_mod[1]-site	R
MTD	assay[4]-quantification_mod[1]-position	Anywhere
MTD	assay[4]-quantification_mod[2]	[UNIMOD, UNIMOD:259, Label:13C(6)15N(2),]
MTD	assay[4]-quantification_mod[2]-site	L
MTD	assay[4]-quantification_mod[2]-position	Anywhere
MTD	assay[6]-quantification_mod[1]	[UNIMOD, UNIMOD:267, Label:13C(6)15N(4),]
MTD	assay[6]-quantification_mod[1]-site	R
MTD	assay[6]-quantification_mod[1]-position	Anywhere
MTD	assay[6]-quantification_mod[2]	[UNIMOD, UNIMOD:259, Label:13C(6)15N(2),]
MTD	assay[6]-quantification_mod[2]-site	L
MTD	assay[6]-quantification_mod[2]-position	Anywhere
MTD	assay[1]-ms_run_ref	ms_run[1]
MTD	assay[2]-ms_run_ref	ms_run[1]
MTD	assay[3]-ms_run_ref	ms_run[2]
MTD	assay[4]-ms_run_ref	ms_run[2]
MTD	assay[5]-ms_run_ref	ms_run[3]
MTD	assay[6]-ms_run_ref	ms_run[3]
MTD	study_variable[1]-assay_refs	assay[1],assay[3],assay[5]
MTD	study_variable[2]-assay_refs	assay[2],assay[4],assay[6]
MTD	study_variable[1]-description	heat shock response of control
MTD	study_variable[2]-description	heat shock response of treatment

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	search_engine_score[1]_ms_run[1]	search_engine_score[1]_ms_run[2]	search_engine_score[1]_ms_run[3]	num_psms_ms_run[1]	num_psms_ms_run[2]	num_psms_ms_run[3]	num_peptides_distinct_ms_run[1]	num_peptides_distinct_ms_run[2]	num_peptides_distinct_ms_run[3]	num_peptides_unique_ms_run[1]	num_peptides_unique_ms_run[2]	num_peptides_unique_ms_run[3]	ambiguity_members	modifications	protein_coverage	protein_abundance_assay[1]	protein_abundance_assay[2]	protein_abundance_assay[3]	protein_abundance_assay[4]	protein_abundance_assay[5]	protein_abundance_assay[6]	protein_abundance_study_variable[1]	protein_abundance_stdev_study_variable[1]	protein_abundance_std_error_study_variable[1]	protein_abundance_study_variable[2]	protein_abundance_stdev_study_variable[2]	protein_abundance_std_error_study_variable[2]
COM	Accession	Description	Taxonomie ID	Species	Database	Version	Search Engine	best Mascot score	Mascot score (HSPRep1)	Mascot score (HSPRep2)	Mascot score (HSPRep3)	PSMs (HSPRep1)	PSMs (HSPRep2)	PSMs (HSPRep3)	Distinct Peptides (HSPRep1)	Distinct Peptides (HSPRep2)	Distinct Peptides (HSPRep3)	Unique Peptides (HSPRep1)	Unique Peptides (HSPRep2)	Unique Peptides (HSPRep3)	Ambiguity Members	Modifications	Protein Coverage (fraction)	Abundance (HSPControlRep1)	Abundance (HSPTreatmentRep1)	Abundance (HSPControlRep2)	Abundance (HSPTreatmentRep2)	Abundance (HSPControlRep3)	Abundance (HSPTreatmentRep3)	Abundance (HSPControl)	Standard Deviation (HSPControl)	Standard Error (HSPControl)	Abundance (HSPTreatment)	Standard Deviation (HSPTreatment)	Stdandard Error (HSPTreatment)
PRT	P63017	Heat shock cognate 71 kDa protein	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	46	26	40	1	1	1	1	1	1	1	1	1	null	0	0.34	44.60127454	1930.618571	37.99134855	1659.350761	26.32448593	1408.464016	36.30570301	9.25425854	5.342948659	1666.14445	261.1435628	150.7713063
PRT	P14602	Heat shock protein beta-1	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	100	100	12	3	3	1	3	3	1	2	2	1	Q340U4,Q5K0U2,P8L901	0	0.12	400.250675	1201.136763	329.8484513	992.5417895	286.5974881	81341.37681	338.8988715	57.36457362	33.11945202	27845.01845	46329.32275	26748.24696
PRT	Q8K0U4	Heat shock 70 kDa protein 12A	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	120	null	22	1	0	1	1	0	1	1	0	1	null	0	0.14	3703.367578	1270364.298	4524.070045	1099450.723	2616.197927	112524.8008	3614.545183	957.0324281	552.54293	827446.6074	625010.2036	360849.8093
PRT	Q61699	Heat shock protein 105 kDa	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	30	8	14	1	1	1	1	1	1	1	1	1	null	0	0.08	43.91120283	13.33293551	44.37102327	9.046028182	29.17909251	11.67165754	39.15377287	8.64138559	4.989106297	11.35020708	2.161455852	1.247917118
PRT	P07901	Heat shock protein HSP 90-alpha	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	45	-5	9	4	3	4	4	3	4	4	3	4	null	12-UNIMOD:35, 98-UNIMOD:35,727-UNIMOD:35	0.21	3.203884813	0.3813336733	3.917537412	9.307189513	2.009982987	0.8119629157	3.043801737	0.9638002458	0.5564503313	3.500162034	5.033640481	2.906173687

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	spectra_ref	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	pre	post	start	end 
COM	Sequence	PSM identifier	accession	Unqiue	Database	Database Version	Search Engine	Mascot score	Modifications	Spectra Reference	Retention Time	Charge	Experimental m/z	Calculated m/z	Pre	Post	Start	end 
PSM	QTQTFTTYSDNQPGVL	1	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	null	ms_run[1]:scan=1296	1336.62	3	600.6474638	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	2	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	null	ms_run[1]:scan=1300	1327.08	2	956.9883817	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	3	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	9-UNIMOD:4	ms_run[1]:scan=845	885.62	3	527.6578728	527.6362	R	K	262	274
PSM	DWYPAHSR	4	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	4	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	4	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	5	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	null	ms_run[1]:scan=1155	1195.62	3	550.9477381	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	6	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=1064	1104.62	3	542.9762487	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	7	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=2849	2876.08	2	1974.396543	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	8	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	23-UNIMOD:35	ms_run[1]:scan=2584	2611.08	2	1788.285098	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	9	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	76	11-UNIMOD:35	ms_run[1]:scan=1092	1132.62	3	450.6037522	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	10	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	87	0-UNIMOD:35	ms_run[1]:scan=3157	3184.08	2	2405.606618	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	11	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	null	ms_run[2]:scan=1530	1336.62	3	600.6109415	600.6197	K	I	424	439
PSM	ALLRLHQECEKLK	12	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	9-UNIMOD:4	ms_run[2]:scan=1079	885.62	3	527.6036448	527.6362	R	K	262	274
PSM	DWYPAHSR	13	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	13	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	13	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	14	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	13	null	ms_run[2]:scan=1389	1195.62	3	550.9312114	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	15	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	33	null	ms_run[2]:scan=1298	1104.62	3	542.9670575	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	16	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-5	null	ms_run[2]:scan=3083	2876.08	2	1974.436352	1974.3984	R	M	692	728
PSM	TLTIVDTGIGMTK	17	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	28	null	ms_run[2]:scan=1326	1132.62	3	450.583395	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	18	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	33	null	ms_run[2]:scan=3391	3184.08	2	2405.631884	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	19	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	null	ms_run[3]:scan=1062	1336.62	3	600.629071	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	20	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	22	null	ms_run[3]:scan=1066	1327.08	2	956.9932461	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	21	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	9-UNIMOD:4	ms_run[3]:scan=611	885.62	3	527.6168147	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	22	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	12	null	ms_run[3]:scan=921	1195.62	3	550.9193765	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	23	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	9	null	ms_run[3]:scan=2615	2876.08	2	1974.397537	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	24	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	9	23-UNIMOD:35	ms_run[3]:scan=2350	2611.08	2	1788.276486	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	25	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	null	ms_run[3]:scan=858	1132.62	3	450.5823742	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	26	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	6	12-UNIMOD:35	ms_run[3]:scan=2923	3184.08	2	2405.596441	2405.6084	-	E	1	41"#;

const ITRAQ_CQI: &str = r#"COM	This	line	serves	as	a	size	and	separator	hint	for	spreadsheet	applications.	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
COM	Report of a minimal "Complete Quantification report" iTRAQ experiment, quantification on 4 study variables (t=0, t=1, t=2, t=3), 4*4 assays (4 replicate experiments) reported, identifications reported.
MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Complete
MTD	mzTab-type	Quantification
MTD	description	mzTab example file for reporting a summary report of quantification data quantified on the protein level
MTD	protein_search_engine_score[1]	[MS, MS:1001171,Mascot:score, ]
MTD	psm_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	ms_run[1]-location	file://C:/path/to/my/file1.mzML
MTD	ms_run[2]-location	file://C:/path/to/my/file2.mzML
MTD	ms_run[3]-location	file://C:/path/to/my/file3.mzML
MTD	ms_run[4]-location	file://C:/path/to/my/file4.mzML
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000393, Relative quantification unit, ]
MTD	software[1]	[MS, MS:1000752, TOPP software, ]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
COM	Specifiying site and position of search modifications is not mandatory but good practice
MTD	fixed_mod[1]-site	C
MTD	fixed_mod[1]-position	Anywhere
MTD	fixed_mod[2]	[UNIMOD, UNIMOD:214, iTRAQ4plex, ]
MTD	fixed_mod[2]-site	K
MTD	fixed_mod[2]-position	Anywhere
MTD	fixed_mod[3]	[UNIMOD, UNIMOD:214, iTRAQ4plex, ]
MTD	fixed_mod[3]-site	N-term
MTD	fixed_mod[3]-position	Any N-term
MTD	variable_mod[1]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[1]-site	M
MTD	variable_mod[1]-position	Anywhere
MTD	quantification_method	[PRIDE,PRIDE:0000313,iTRAQ, ]
MTD	assay[1]-quantification_reagent	[PRIDE, PRIDE:0000114, iTRAQ reagent 114, ]
MTD	assay[2]-quantification_reagent	[PRIDE, PRIDE:0000115, iTRAQ reagent 115, ]
MTD	assay[3]-quantification_reagent	[PRIDE, PRIDE:0000116, iTRAQ reagent 116, ]
MTD	assay[4]-quantification_reagent	[PRIDE, PRIDE:0000117, iTRAQ reagent 117, ]
MTD	assay[5]-quantification_reagent	[PRIDE, PRIDE:0000114, iTRAQ reagent 114, ]
MTD	assay[6]-quantification_reagent	[PRIDE, PRIDE:0000115, iTRAQ reagent 115, ]
MTD	assay[7]-quantification_reagent	[PRIDE, PRIDE:0000116, iTRAQ reagent 116, ]
MTD	assay[8]-quantification_reagent	[PRIDE, PRIDE:0000117, iTRAQ reagent 117, ]
MTD	assay[9]-quantification_reagent	[PRIDE, PRIDE:0000114, iTRAQ reagent 114, ]
MTD	assay[10]-quantification_reagent	[PRIDE, PRIDE:0000115, iTRAQ reagent 115, ]
MTD	assay[11]-quantification_reagent	[PRIDE, PRIDE:0000116, iTRAQ reagent 116, ]
MTD	assay[12]-quantification_reagent	[PRIDE, PRIDE:0000117, iTRAQ reagent 117, ]
MTD	assay[13]-quantification_reagent	[PRIDE, PRIDE:0000114, iTRAQ reagent 114, ]
MTD	assay[14]-quantification_reagent	[PRIDE, PRIDE:0000115, iTRAQ reagent 115, ]
MTD	assay[15]-quantification_reagent	[PRIDE, PRIDE:0000116, iTRAQ reagent 116, ]
MTD	assay[16]-quantification_reagent	[PRIDE, PRIDE:0000117, iTRAQ reagent 117, ]
MTD	assay[1]-ms_run_ref	ms_run[1]
MTD	assay[2]-ms_run_ref	ms_run[1]
MTD	assay[3]-ms_run_ref	ms_run[1]
MTD	assay[4]-ms_run_ref	ms_run[1]
MTD	assay[5]-ms_run_ref	ms_run[2]
MTD	assay[6]-ms_run_ref	ms_run[2]
MTD	assay[7]-ms_run_ref	ms_run[2]
MTD	assay[8]-ms_run_ref	ms_run[2]
MTD	assay[9]-ms_run_ref	ms_run[3]
MTD	assay[10]-ms_run_ref	ms_run[3]
MTD	assay[11]-ms_run_ref	ms_run[3]
MTD	assay[12]-ms_run_ref	ms_run[3]
MTD	assay[13]-ms_run_ref	ms_run[4]
MTD	assay[14]-ms_run_ref	ms_run[4]
MTD	assay[15]-ms_run_ref	ms_run[4]
MTD	assay[16]-ms_run_ref	ms_run[4]
MTD	study_variable[1]-assay_refs	assay[1],assay[6],assay[11],assay[16]
MTD	study_variable[2]-assay_refs	assay[2],assay[5],assay[12],assay[15]
MTD	study_variable[3]-assay_refs	assay[3],assay[7],assay[9],assay[14]
MTD	study_variable[4]-assay_refs	assay[4],assay[8],assay[10],assay[13]
MTD	study_variable[1]-description	Incubation at t=0
MTD	study_variable[2]-description	Incubation at t=1
MTD	study_variable[3]-description	Incubation at t=2
MTD	study_variable[4]-description	Incubation at t=3

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	search_engine_score[1]_ms_run[1]	search_engine_score[1]_ms_run[2]	search_engine_score[1]_ms_run[3]	search_engine_score[1]_ms_run[4]	num_psms_ms_run[1]	num_psms_ms_run[2]	num_psms_ms_run[3]	num_psms_ms_run[4]	num_peptides_distinct_ms_run[1]	num_peptides_distinct_ms_run[2]	num_peptides_distinct_ms_run[3]	num_peptides_distinct_ms_run[4]	num_peptides_unique_ms_run[1]	num_peptides_unique_ms_run[2]	num_peptides_unique_ms_run[3]	num_peptides_unique_ms_run[4]	ambiguity_members	modifications	protein_coverage	protein_abundance_assay[1]	protein_abundance_assay[2]	protein_abundance_assay[3]	protein_abundance_assay[4]	protein_abundance_assay[5]	protein_abundance_assay[6]	protein_abundance_assay[7]	protein_abundance_assay[8]	protein_abundance_assay[9]	protein_abundance_assay[10]	protein_abundance_assay[11]	protein_abundance_assay[12]	protein_abundance_assay[13]	protein_abundance_assay[14]	protein_abundance_assay[15]	protein_abundance_assay[16]	protein_abundance_study_variable[1]	protein_abundance_stdev_study_variable[1]	protein_abundance_std_error_study_variable[1]	protein_abundance_study_variable[2]	protein_abundance_stdev_study_variable[2]	protein_abundance_std_error_study_variable[2]	protein_abundance_study_variable[3]	protein_abundance_stdev_study_variable[3]	protein_abundance_std_error_study_variable[3]	protein_abundance_study_variable[4]	protein_abundance_stdev_study_variable[4]	protein_abundance_std_error_study_variable[4]
COM	Accession	Description	Taxonomie ID	Species	Database	Version	Search Engine	best Mascot score	Mascot score (Run1)	Mascot score (Run2)	Mascot score (Run3)	Mascot score (Run4)	PSMs (Run1)	PSMs (Run2)	PSMs (Run3)	PSMs (Run4)	Distinct Peptides (Run1)	Distinct Peptides (Run2)	Distinct Peptides (Run3)	Distinct Peptides (Run4)	Unique Peptides (Run1)	Unique Peptides (Run2)	Unique Peptides (Run3)	Unique Peptides (Run4)	Ambiguity Members	Modifications	Protein Coverage (fraction)	Abundance (t=0_114)	Abundance (t=1_115)	Abundance (t=2_116)	Abundance (t=3_117)	Abundance (t=1_114)	Abundance (t=0_115)	Abundance (t=3_116)	Abundance (t=2_117)	Abundance (t=2_114)	Abundance (t=3_115)	Abundance (t=0_116)	Abundance (t=1_117)	Abundance (t=3_114)	Abundance (t=2_115)	Abundance (t=1_116)	Abundance (t=0_117)	Abundance (t=0)	Standard Deviation (t=0)	Standard Error (t=0)	Abundance (t=1)	Standard Deviation (t=1)	Standard Error (t=1)	Abundance (t=2)	Standard Deviation (t=2)	Standard Error (t=2)	Abundance (t=3)	Standard Deviation (t=3)	Standard Error (t=3)
PRT	P63017	Heat shock cognate 71 kDa protein	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	46	26	5	-3	1	1	1	1	1	1	1	1	1	1	1	1	null	0	0.34	17.3	45.4	234.3	26.7	45.82252834	20.34470807	31.44819833	271.9769613	286.1141095	31.6970172	16.51509633	68.77899312	30.84439527	283.0325889	54.48687448	16.18702563	17.58670751	1.897035247	1.095253811	53.62209899	10.93793842	6.315021691	268.8559149	23.82348667	13.75449644	30.1724027	2.342453014	1.352415878
PRT	P14602	Heat shock protein beta-1	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	14	100	-6	20	2	3	1	2	2	3	1	2	2	2	1	2	Q340U4,Q5K0U2,P8L901	0	0.12	343.54	643.4	2353	0.82	893.3569959	509.3562129	0.8703280099	2369.88544	2659.337273	1.013956057	477.9996074	1031.539362	1.135622127	3564.408615	1157.260777	211.1783086	385.5185322	136.6811754	78.91291343	931.3892837	220.1758783	127.1186026	2736.657832	569.4632404	328.7797552	0.9599765485	0.1430564589	0.08259368505
PRT	Q8K0U4	Heat shock 70 kDa protein 12A	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	120	120	21	39	1	1	1	1	1	1	1	1	1	1	1	1	null	0	0.14	6453.4	21.3	234.3	25.5	24.57083049	9725.949445	27.79987837	342.4708472	339.1391885	27.61650271	7419.782843	24.96532034	27.59736459	334.3343844	31.14666952	9003.08556	8150.554462	1485.826263	857.842193	25.49570509	4.109908402	2.372856722	312.561105	52.28085519	30.18436582	27.12843642	1.089455797	0.6289975977
PRT	Q61699	Heat shock protein 105 kDa	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	31	30	12	31	-2	1	1	1	1	1	1	1	1	1	1	1	1	null	0	0.08	346534.4	1232.5	54.2	1654.6	1715.04997	363520.6744	2977.115494	70.14750308	70.517314	1748.483298	181075.9557	1318.408071	2144.86276	62.7205008	1322.671448	352462.8461	310898.4691	86834.13016	50133.70842	1397.157372	215.960922	124.6850964	64.39632947	7.687998998	4.438668291	2131.265388	602.599732	347.9111175
PRT	P07901	Heat shock protein HSP 90-alpha	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	45	14	19	35	4	3	4	3	4	3	4	3	4	3	4	3	null	0-UNIMOD:35, 12-UNIMOD:35, 98-UNIMOD:35,727-UNIMOD:35	0.21	6434.6	123.3	234.3	110.2	158.9389923	6475.205011	177.6867536	266.8830076	251.6690321	121.5177559	5413.126383	125.5419061	116.8450207	242.6034553	136.2170466	5812.480954	6033.853087	513.061683	296.2163008	135.9994863	16.29832511	9.409842388	248.8638737	13.95059563	8.05438014	131.5623826	31.09825426	17.95458547

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	spectra_ref	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	pre	post	start	end 
COM	Sequence	PSM identifier	accession	Unqiue	Database	Database Version	Search Engine	Mascot score	Modifications	Spectra Reference	Retention Time	Charge	Experimental m/z	Calculated m/z	Pre	Post	Start	end 
PSM	QTQTFTTYSDNQPGVL	1	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	0-UNIMOD:214	ms_run[1]:scan=1296	1336.62	3	600.6218923	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	2	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,20-UNIMOD:214	ms_run[1]:scan=1300	1327.08	2	956.9534673	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	3	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[1]:scan=845	885.62	3	527.6131796	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	4	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	0-UNIMOD:214	ms_run[1]:scan=1155	1195.62	3	550.9127359	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	5	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	0-UNIMOD:214,12-UNIMOD:214,14-UNIMOD:214	ms_run[1]:scan=1064	1104.62	3	542.9870349	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	6	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	0-UNIMOD:214,3-UNIMOD:214	ms_run[1]:scan=2849	2876.08	2	1974.399785	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	7	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,23-UNIMOD:35	ms_run[1]:scan=2584	2611.08	2	1788.285564	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	8	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	76	0-UNIMOD:214,11-UNIMOD:35,13-UNIMOD:214	ms_run[1]:scan=1092	1132.62	3	450.5826992	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	9	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	87	0-UNIMOD:214,41-UNIMOD:214	ms_run[1]:scan=3157	3184.08	2	2405.61108	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	10	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	0-UNIMOD:214	ms_run[2]:scan=1530	1336.62	3	600.6243317	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	11	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,20-UNIMOD:214	ms_run[2]:scan=1534	1327.08	2	956.9697829	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	12	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	27	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[2]:scan=1079	885.62	3	527.6286638	527.6362	R	K	262	274
PSM	DWYPAHSR	13	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	13	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	13	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	14	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	33	0-UNIMOD:214	ms_run[2]:scan=1389	1195.62	3	550.9230929	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	15	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	0-UNIMOD:214,12-UNIMOD:214,14-UNIMOD:214	ms_run[2]:scan=1298	1104.62	3	542.9873474	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	16	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	0-UNIMOD:214,3-UNIMOD:214	ms_run[2]:scan=3083	2876.08	2	1974.379186	1974.3984	R	M	692	728
PSM	TLTIVDTGIGMTK	17	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	0-UNIMOD:214,13-UNIMOD:214	ms_run[2]:scan=1326	1132.62	3	450.5778381	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	18	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	25	0-UNIMOD:214,41-UNIMOD:214	ms_run[2]:scan=3391	3184.08	2	2405.596577	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	19	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	3	0-UNIMOD:214	ms_run[3]:scan=1062	1336.62	3	600.6219355	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	20	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	24	0-UNIMOD:214,20-UNIMOD:214	ms_run[3]:scan=1066	1327.08	2	956.9720273	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	21	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-4	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[3]:scan=611	885.62	3	527.651719	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	22	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	29	0-UNIMOD:214	ms_run[3]:scan=921	1195.62	3	550.9221479	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	23	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	39	0-UNIMOD:214,3-UNIMOD:214	ms_run[3]:scan=2615	2876.08	2	1974.419328	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	24	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	34	0-UNIMOD:214	ms_run[3]:scan=2350	2611.08	2	1788.298223	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	25	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	0-UNIMOD:214,13-UNIMOD:214	ms_run[3]:scan=858	1132.62	3	450.5423216	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	26	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	20	0-UNIMOD:214,41-UNIMOD:214	ms_run[3]:scan=2923	3184.08	2	2405.620339	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	27	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[4]:scan=2731	1336.62	3	600.6123009	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	28	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	39	null	ms_run[4]:scan=2735	1327.08	2	956.9765302	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	29	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	9-UNIMOD:4	ms_run[4]:scan=2280	885.62	3	527.6343404	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	30	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	20	null	ms_run[4]:scan=2590	1195.62	3	550.9284574	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	31	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[4]:scan=2499	1104.62	3	542.9715699	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	32	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	35	null	ms_run[4]:scan=4284	2876.08	2	1974.40429	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	33	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	0	23-UNIMOD:35	ms_run[4]:scan=4019	2611.08	2	1788.289062	1788.2886	K	M	695	728
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	34	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	11	0-UNIMOD:35	ms_run[4]:scan=4592	3184.08	2	2405.57421	2405.6084	-	E	1	41"#;

const ITRAQ_SQI: &str = r#"COM	This	line	serves	as	a	size	and	separator	hint	for	spreadsheet	applications.	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
COM	Report of a minimal Summary Quantification report" iTRAQ experiment, quantification on 4 study variables (t=0, t=1, t=2, t=3) reported, identifications reported.
MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Quantification
MTD	description	mzTab example file for reporting a summary report of quantification data quantified on the protein level
MTD	protein_search_engine_score[1]	[MS, MS:1001171, Mascot:score,]
MTD	psm_search_engine_score[1]	[MS, MS:1001171, Mascot:score,]
MTD	ms_run[1]-location	file://C:/path/to/my/file1.mzML
MTD	ms_run[2]-location	file://C:/path/to/my/file2.mzML
MTD	ms_run[3]-location	file://C:/path/to/my/file3.mzML
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000393, Relative quantification unit,]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	fixed_mod[2]	[UNIMOD, UNIMOD:214, iTRAQ4plex, ]
MTD	variable_mod[1]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	study_variable[1]-description	sample incubation at t=0
MTD	study_variable[2]-description	sample incubation at t=1
MTD	study_variable[3]-description	sample incubation at t=2
MTD	study_variable[4]-description	sample incubation at t=3

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	ambiguity_members	modifications	protein_abundance_study_variable[1]	protein_abundance_stdev_study_variable[1]	protein_abundance_std_error_study_variable[1]	protein_abundance_study_variable[2]	protein_abundance_stdev_study_variable[2]	protein_abundance_std_error_study_variable[2]	protein_abundance_study_variable[3]	protein_abundance_stdev_study_variable[3]	protein_abundance_std_error_study_variable[3]	protein_abundance_study_variable[4]	protein_abundance_stdev_study_variable[4]	protein_abundance_std_error_study_variable[4]
COM	Accession	Description	Taxonomie ID	Species	Database	Version	Search Engine	best Mascot score	Ambiguity Members	Modifications	Abundance (t=0)	Standard Deviation (t=0)	Standard Error (t=0)	Abundance (t=1)	Standard Deviation (t=1)	Standard Error (t=1)	Abundance (t=2)	Standard Deviation (t=2)	Standard Error (t=2)	Abundance (t=3)	Standard Deviation (t=3)	Standard Error (t=3)
PRT	P63017	Heat shock cognate 71 kDa protein	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	null	0	14.361	3.6605	2.1134	53.695	11.73	6.7722	276.59	45.934	26.52	29.767	3.1318	1.8082
PRT	P14602	Heat shock protein beta-1	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	Q340U4,Q5K0U2,P8L901	0	261.98	162.3	93.702	777.09	217.99	125.85	2865.4	514.48	297.04	0.8993	0.0893	0.0516
PRT	Q8K0U4	Heat shock 70 kDa protein 12A	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	null	0	7707.3	1310.5	756.62	23.923	5.0532	2.9175	270.59	33.407	19.287	28.595	2.9196	1.6856
PRT	Q61699	Heat shock protein 105 kDa	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	null	0	381591	92914	53644	1424	320.51	185.04	63.366	6.6935	3.8645	1737.5	70.203	40.532
PRT	P07901	Heat shock protein HSP 90-alpha	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	0-UNIMOD:35, 12-UNIMOD:35, 98-UNIMOD:35,727-UNIMOD:35	10051	2453.1	1416.3	130.13	7.2703	4.1975	278.63	36.593	21.127	126.18	11.102	6.4098

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	spectra_ref	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	pre	post	start	end
COM	Sequence	PSM identifier	accession	Unqiue	Database	Database Version	Search Engine	Mascot score	Modifications	Spectra Reference	Retention Time	Charge	Experimental m/z	Calculated m/z	Pre	Post	Start	end
PSM	QTQTFTTYSDNQPGVL	1	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	0-UNIMOD:214	ms_run[1]:scan=1296	1336.62	3	600.62	600.62	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	2	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,20-UNIMOD:214	ms_run[1]:scan=1300	1327.08	2	956.98	956.974	K	E	261	281
PSM	ALLRLHQECEKLK	3	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[1]:scan=845	885.62	3	527.635	527.636	R	K	262	274
PSM	MNQSNASPTLDGLFR	4	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	0-UNIMOD:214	ms_run[1]:scan=1155	1195.62	3	550.945	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	5	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	0-UNIMOD:214,12-UNIMOD:214,14-UNIMOD:214	ms_run[1]:scan=1064	1104.62	3	542.972	542.972	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	6	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	0-UNIMOD:214,3-UNIMOD:214	ms_run[1]:scan=2849	2876.08	2	1974.37	1974.4	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	7	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,23-UNIMOD:35	ms_run[1]:scan=2584	2611.08	2	1788.29	1788.29	K	M	695	728
PSM	TLTIVDTGIGMTK	8	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	76	0-UNIMOD:214,11-UNIMOD:35,13-UNIMOD:214	ms_run[1]:scan=1092	1132.62	3	450.551	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	9	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	87	0-UNIMOD:214,41-UNIMOD:214	ms_run[1]:scan=3157	3184.08	2	2405.58	2405.61	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	10	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	0-UNIMOD:214	ms_run[2]:scan=1530	1336.62	3	600.618	600.62	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	11	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	0-UNIMOD:214,20-UNIMOD:214	ms_run[2]:scan=1534	1327.08	2	956.98	956.974	K	E	261	281
PSM	ALLRLHQECEKLK	12	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	12	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[2]:scan=1079	885.62	3	527.625	527.636	R	K	262	274
PSM	DWYPAHSR	13	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.238	R	L	21	28
PSM	DWYPAHSR	13	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.238	K	E	143	150
PSM	DWYPAHSR	13	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	0-UNIMOD:214	ms_run[2]:scan=778	571.08	2	516.21	516.238	R	M	240	247
PSM	MNQSNASPTLDGLFR	14	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-4	0-UNIMOD:214	ms_run[2]:scan=1389	1195.62	3	550.935	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	15	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	0-UNIMOD:214,12-UNIMOD:214,14-UNIMOD:214	ms_run[2]:scan=1298	1104.62	3	542.994	542.972	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	16	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-4	0-UNIMOD:214,3-UNIMOD:214	ms_run[2]:scan=3083	2876.08	2	1974.39	1974.4	R	M	692	728
PSM	TLTIVDTGIGMTK	17	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-7	0-UNIMOD:214,13-UNIMOD:214	ms_run[2]:scan=1326	1132.62	3	450.584	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	18	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	28	0-UNIMOD:214,41-UNIMOD:214	ms_run[2]:scan=3391	3184.08	2	2405.62	2405.61	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	19	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	0-UNIMOD:214	ms_run[3]:scan=1062	1336.62	3	600.615	600.62	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	20	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	28	0-UNIMOD:214,20-UNIMOD:214	ms_run[3]:scan=1066	1327.08	2	956.979	956.974	K	E	261	281
PSM	ALLRLHQECEKLK	21	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	0-UNIMOD:214,9-UNIMOD:4,11-UNIMOD:214,13-UNIMOD:214	ms_run[3]:scan=611	885.62	3	527.631	527.636	R	K	262	274
PSM	MNQSNASPTLDGLFR	22	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	3	0-UNIMOD:214	ms_run[3]:scan=921	1195.62	3	550.924	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	23	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	25	0-UNIMOD:214,3-UNIMOD:214	ms_run[3]:scan=2615	2876.08	2	1974.42	1974.4	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	24	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-10	0-UNIMOD:214	ms_run[3]:scan=2350	2611.08	2	1788.29	1788.29	K	M	695	728
PSM	TLTIVDTGIGMTK	25	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	0-UNIMOD:214,13-UNIMOD:214	ms_run[3]:scan=858	1132.62	3	450.599	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	26	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-7	0-UNIMOD:214,41-UNIMOD:214	ms_run[3]:scan=2923	3184.08	2	2405.62	2405.61	-	E	1	41"#;

const LABELFREE_CQI: &str = r#"COM	This	line	serves	as	a	size	and	separator	hint	for	spreadsheet	applications.	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
COM	Report of a minimal "Complete Quantification report" label free experiment, quantification on 2 study variables (control/treatment), 3+3 assays (replicates) reported,identifications reported.
MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Complete
MTD	mzTab-type	Quantification
MTD	description	mzTab example file for reporting a summary report of quantification data quantified on the protein level
MTD	protein_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	psm_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	ms_run[1]-location	file://C:/path/to/my/file1.mzML
MTD	ms_run[2]-location	file://C:/path/to/my/file2.mzML
MTD	ms_run[3]-location	file://C:/path/to/my/file3.mzML
MTD	ms_run[4]-location	file://C:/path/to/my/file4.mzML
MTD	ms_run[5]-location	file://C:/path/to/my/file5.mzML
MTD	ms_run[6]-location	file://C:/path/to/my/file6.mzML
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000393, Relative quantification unit,]
MTD	software[1]	[MS, MS:1000752, TOPP software,]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	variable_mod[1]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	quantification_method	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[1]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[2]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[3]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[4]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[5]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[6]-quantification_reagent	[MS, MS:1002038, unlabeled sample, ]
MTD	assay[1]-ms_run_ref	ms_run[1]
MTD	assay[2]-ms_run_ref	ms_run[2]
MTD	assay[3]-ms_run_ref	ms_run[3]
MTD	assay[4]-ms_run_ref	ms_run[4]
MTD	assay[5]-ms_run_ref	ms_run[5]
MTD	assay[6]-ms_run_ref	ms_run[6]
MTD	study_variable[1]-assay_refs	assay[1], assay[2], assay[3]
MTD	study_variable[2]-assay_refs	assay[4], assay[5], assay[6]
MTD	study_variable[1]-description	heat shock response of control
MTD	study_variable[2]-description	heat shock response of treatment

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	search_engine_score[1]_ms_run[1]	search_engine_score[1]_ms_run[2]	search_engine_score[1]_ms_run[3]	search_engine_score[1]_ms_run[4]	search_engine_score[1]_ms_run[5]	search_engine_score[1]_ms_run[6]	num_psms_ms_run[1]	num_psms_ms_run[2]	num_psms_ms_run[3]	num_psms_ms_run[4]	num_psms_ms_run[5]	num_psms_ms_run[6]	num_peptides_distinct_ms_run[1]	num_peptides_distinct_ms_run[2]	num_peptides_distinct_ms_run[3]	num_peptides_distinct_ms_run[4]	num_peptides_distinct_ms_run[5]	num_peptides_distinct_ms_run[6]	num_peptides_unique_ms_run[1]	num_peptides_unique_ms_run[2]	num_peptides_unique_ms_run[3]	num_peptides_unique_ms_run[4]	num_peptides_unique_ms_run[5]	num_peptides_unique_ms_run[6]	ambiguity_members	modifications	protein_coverage	protein_abundance_assay[1]	protein_abundance_assay[2]	protein_abundance_assay[3]	protein_abundance_assay[4]	protein_abundance_assay[5]	protein_abundance_assay[6]	protein_abundance_study_variable[1]	protein_abundance_stdev_study_variable[1]	protein_abundance_std_error_study_variable[1]	protein_abundance_study_variable[2]	protein_abundance_stdev_study_variable[2]	protein_abundance_std_error_study_variable[2]
COM	Accession	Description	Taxonomie ID	Species	Database	Version	Search Engine	best Mascot score	Mascot score (HSPControlRep1)	Mascot score (HSPControlRep2)	Mascot score (HSPControlRep3)	Mascot score (HSPTreatmentRep1)	Mascot score (HSPTreatmentRep2)	Mascot score (HSPTreatmentRep3)	PSMs (HSPControlRep1)	PSMs (HSPControlRep2)	PSMs (HSPControlRep3)	PSMs (HSPTreatmentRep4)	PSMs (HSPTreatmentRep5)	PSMs (HSPTreatmentRep6)	Distinct Peptides (HSPControlRep1)	Distinct Peptides (HSPControlRep2)	Distinct Peptides (HSPControlRep3)	Distinct Peptides (HSPTreatmentRep4)	Distinct Peptides (HSPTreatmentRep5)	Distinct Peptides (HSPTreatmentRep6)	Unique Peptides (HSPControlRep1)	Unique Peptides (HSPControlRep2)	Unique Peptides (HSPControlRep3)	Unique Peptides (HSPTreatmentRep4)	Unique Peptides (HSPTreatmentRep5)	Unique Peptides (HSPTreatmentRep6)	Ambiguity Members	Modifications	Protein Coverage (fraction)	Abundance (HSPTreatmentRep1)	Abundance (HSPTreatmentRep2)	Abundance (HSPTreatmentRep3)	Abundance (HSPTreatmentRep4)	Abundance (HSPTreatmentRep5)	Abundance (HSPTreatmentRep6)	Abundance (HSPControl)	Standard Deviation (HSPControl)	Standard Error (HSPControl)	Abundance (HSPTreatment)	Standard Deviation (HSPTreatment)	Standard Error (HSPTreatment)
PRT	P63017	Heat shock cognate 71 kDa protein	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	46	26	36	-3	-1	null	1	1	1	1	1	0	1	1	1	1	1	0	1	1	1	1	1	0	null	0	0.34	34.3	40.43507695	41.12124635	266.9554147	234.4	271.0324163	38.61877444	3.755870949	2.168453103	257.4626103	20.07656548	11.59121048
PRT	P14602	Heat shock protein beta-1	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	100	100	-9	20	100	4	3	3	1	2	3	2	3	3	1	2	3	2	2	2	1	2	2	1	Q340U4,Q5K0U2,P8L901	0	0.12	98588.4	114212.9033	100070.7061	4709.411242	4345.7	6704.588342	104290.6698	8624.809914	4979.536326	5253.233195	1269.998146	733.2337713
PRT	Q8K0U4	Heat shock 70 kDa protein 12A	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	120	null	-2	39	36	-7	1	0	1	1	1	1	1	0	1	1	1	1	1	0	1	1	1	1	null	0	0.14	43.4	86.09123822	54.98032306	459.4934179	375.5	609.3477328	61.49052043	22.0776461	12.74653492	481.4470502	118.4595375	68.39264584
PRT	Q61699	Heat shock protein 105 kDa	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	30	31	36	-2	31	24	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	1	null	0	0.08	3432.54	3847.349077	3838.448278	11372.30364	9587.5	10303.56594	3706.112452	236.9624883	136.8103564	10421.12319	898.190283	518.5704017
PRT	P07901	Heat shock protein HSP 90-alpha	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	45	6	-2	35	8	3	4	3	4	3	2	4	4	3	4	3	2	4	4	3	4	3	2	4	null	12-UNIMOD:35, 98-UNIMOD:35,727-UNIMOD:35	0.21	3242354.3	3284123.069	3404460.592	633072.591	552426.4	618457.6276	3310312.654	84166.6994	48593.66656	601318.8729	42968.06623	24807.6246

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	spectra_ref	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	pre	post	start	end
COM	Sequence	PSM identifier	accession	Unqiue	Database	Database Version	Search Engine	Mascot score	Modifications	Spectra Reference	Retention Time	Charge	Experimental m/z	Calculated m/z	Pre	Post	Start	End
PSM	QTQTFTTYSDNQPGVL	1	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	null	ms_run[1]:scan=1296	1336.62	3	600.6006697	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	2	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	null	ms_run[1]:scan=1300	1327.08	2	956.9464833	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	3	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	9-UNIMOD:4	ms_run[1]:scan=845	885.62	3	527.6406579	527.6362	R	K	262	274
PSM	DWYPAHSR	4	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	4	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	4	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	5	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	null	ms_run[1]:scan=1155	1195.62	3	550.9282794	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	6	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=1064	1104.62	3	542.9688356	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	7	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=2849	2876.08	2	1974.400793	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	8	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	23-UNIMOD:35	ms_run[1]:scan=2584	2611.08	2	1788.281997	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	9	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	76	11-UNIMOD:35	ms_run[1]:scan=1092	1132.62	3	450.5920214	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	10	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	87	0-UNIMOD:35	ms_run[1]:scan=3157	3184.08	2	2405.587318	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	11	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	null	ms_run[2]:scan=1530	1336.62	3	600.6265518	600.6197	K	I	424	439
PSM	ALLRLHQECEKLK	12	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	9-UNIMOD:4	ms_run[2]:scan=1079	885.62	3	527.6362432	527.6362	R	K	262	274
PSM	DWYPAHSR	13	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	13	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	13	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	14	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	null	ms_run[2]:scan=1389	1195.62	3	550.9468571	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	15	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	16	null	ms_run[2]:scan=1298	1104.62	3	542.9666503	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	16	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	6	null	ms_run[2]:scan=3083	2876.08	2	1974.399035	1974.3984	R	M	692	728
PSM	TLTIVDTGIGMTK	17	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	21	null	ms_run[2]:scan=1326	1132.62	3	450.5400013	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	18	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-5	null	ms_run[2]:scan=3391	3184.08	2	2405.599817	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	19	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	null	ms_run[3]:scan=1062	1336.62	3	600.6484541	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	20	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	null	ms_run[3]:scan=1066	1327.08	2	956.9766608	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	21	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	9-UNIMOD:4	ms_run[3]:scan=611	885.62	3	527.6486368	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	22	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-9	null	ms_run[3]:scan=921	1195.62	3	550.9336303	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	23	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	null	ms_run[3]:scan=2615	2876.08	2	1974.392219	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	24	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	23-UNIMOD:35	ms_run[3]:scan=2350	2611.08	2	1788.28771	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	25	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	37	null	ms_run[3]:scan=858	1132.62	3	450.5960917	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	26	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	6	12-UNIMOD:35	ms_run[3]:scan=2923	3184.08	2	2405.604605	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	27	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[4]:scan=2731	1336.62	3	600.6123009	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	28	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	39	null	ms_run[4]:scan=2735	1327.08	2	956.9765302	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	29	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	9-UNIMOD:4	ms_run[4]:scan=2280	885.62	3	527.6343404	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	30	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	20	null	ms_run[4]:scan=2590	1195.62	3	550.9284574	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	31	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[4]:scan=2499	1104.62	3	542.9715699	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	32	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	35	null	ms_run[4]:scan=4284	2876.08	2	1974.40429	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	33	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	0	23-UNIMOD:35	ms_run[4]:scan=4019	2611.08	2	1788.289062	1788.2886	K	M	695	728
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	34	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	11	0-UNIMOD:35	ms_run[4]:scan=4592	3184.08	2	2405.57421	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	35	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-1	null	ms_run[5]:scan=2031	1336.62	3	600.5900228	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	36	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	null	ms_run[5]:scan=2035	1327.08	2	956.9477197	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	37	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	31	9-UNIMOD:4	ms_run[5]:scan=1580	885.62	3	527.6254449	527.6362	R	K	262	274
PSM	DWYPAHSR	38	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	38	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	38	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	39	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	32	null	ms_run[5]:scan=1890	1195.62	3	550.9120992	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	40	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	21	null	ms_run[5]:scan=1799	1104.62	3	542.9599424	542.9716	K	V	91	104
PSM	TLTIVDTGIGMTK	41	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	11-UNIMOD:35	ms_run[5]:scan=1827	1132.62	3	450.5534561	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	42	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-5	0-UNIMOD:35	ms_run[5]:scan=3892	3184.08	2	2405.594573	2405.6084	-	E	1	41
PSM	AVVNGYSASDTVGAGFAQAK	43	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-7	null	ms_run[6]:scan=1331	1327.08	2	956.9880766	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	44	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	24	9-UNIMOD:4	ms_run[6]:scan=876	885.62	3	527.64539	527.6362	R	K	262	274
PSM	DWYPAHSR	45	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	4	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	45	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	45	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	9	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	46	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	32	null	ms_run[6]:scan=1186	1195.62	3	550.9319012	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	47	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	3	null	ms_run[6]:scan=2880	2876.08	2	1974.377816	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	48	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	29	23-UNIMOD:35	ms_run[6]:scan=2615	2611.08	2	1788.294771	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	49	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	39	11-UNIMOD:35	ms_run[6]:scan=1123	1132.62	3	450.6038036	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	50	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	33	0-UNIMOD:35	ms_run[6]:scan=3188	3184.08	2	2405.605739	2405.6084	-	E	1	41"#;

const LABELFREE_SQI: &str = r#"COM	This	line	serves	as	a	size	and	separator	hint	for	spreadsheet	applications.	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-	-
COM	Report of a minimal "Summary Quantification report" label free experiment, quantification on 2 study variables (control/treatment) reported,identifications reported.
MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Quantification
MTD	description	mzTab example file for reporting a summary report of quantification data quantified on the protein level
MTD	protein_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	psm_search_engine_score[1]	[MS,MS:1001171,Mascot:score,]
MTD	ms_run[1]-location	file://C:/path/to/my/file1.mzML
MTD	ms_run[2]-location	file://C:/path/to/my/file2.mzML
MTD	ms_run[3]-location	file://C:/path/to/my/file3.mzML
MTD	ms_run[4]-location	file://C:/path/to/my/file4.mzML
MTD	ms_run[5]-location	file://C:/path/to/my/file5.mzML
MTD	ms_run[6]-location	file://C:/path/to/my/file6.mzML
MTD	protein-quantification_unit	[PRIDE, PRIDE:0000393, Relative quantification unit,]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	variable_mod[1]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	study_variable[1]-description	heat shock response of control
MTD	study_variable[2]-description	heat shock response of treatment

PRH	accession	description	taxid	species	database	database_version	search_engine	best_search_engine_score[1]	ambiguity_members	modifications	protein_abundance_study_variable[1]	protein_abundance_stdev_study_variable[1]	protein_abundance_std_error_study_variable[1]	protein_abundance_study_variable[2]	protein_abundance_stdev_study_variable[2]	protein_abundance_std_error_study_variable[2]
COM	Accession	Description	Taxonomie ID	Species	Database	Version	Search Engine	best Mascot score	Ambiguity Members	Modifications	Abundance (HSPControl)	Standard Deviation (HSPControl)	Standard Error (HSPControl)	Abundance (HSPTreatment)	Standard Deviation (HSPTreatment)	Standard Error (HSPTreatment)
PRT	P63017	Heat shock cognate 71 kDa protein	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	null	0	40.63833264	5.655317647	3.265099166	274.6104756	35.15499814	20.29674764
PRT	P14602	Heat shock protein beta-1	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	Q340U4,Q5K0U2,P8L901	0	117235.3899	17568.59128	10143.23091	5135.628975	1014.552511	585.7521654
PRT	Q8K0U4	Heat shock 70 kDa protein 12A	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	null	0	48.30102441	5.476244671	3.161711335	470.9562475	85.04568626	49.10114986
PRT	Q61699	Heat shock protein 105 kDa	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	37	null	0	4471.580171	1756.397919	1014.056812	9823.949819	267.8320859	154.6329269
PRT	P07901	Heat shock protein HSP 90-alpha	10090	Mus musculus	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	12-UNIMOD:35, 98-UNIMOD:35,727-UNIMOD:35	3468571.127	332005.0444	191683.2018	734432.1524	246781.4666	142479.3462

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	spectra_ref	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	pre	post	start	end
COM	Sequence	PSM identifier	accession	Unqiue	Database	Database Version	Search Engine	Mascot score	Modifications	Spectra Reference	Retention Time	Charge	Experimental m/z	Calculated m/z	Pre	Post	Start	End
PSM	QTQTFTTYSDNQPGVL	1	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	46	null	ms_run[1]:scan=1296	1336.62	3	600.6569942	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	2	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	null	ms_run[1]:scan=1300	1327.08	2	956.9749847	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	3	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	30	9-UNIMOD:4	ms_run[1]:scan=845	885.62	3	527.5989454	527.6362	R	K	262	274
PSM	DWYPAHSR	4	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	4	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	4	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[1]:scan=544	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	5	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	14	null	ms_run[1]:scan=1155	1195.62	3	550.9282581	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	6	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=1064	1104.62	3	542.9814129	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	7	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	45	null	ms_run[1]:scan=2849	2876.08	2	1974.419097	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	8	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	120	23-UNIMOD:35	ms_run[1]:scan=2584	2611.08	2	1788.273327	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	9	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	76	11-UNIMOD:35	ms_run[1]:scan=1092	1132.62	3	450.5972474	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	10	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	87	0-UNIMOD:35	ms_run[1]:scan=3157	3184.08	2	2405.590711	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	11	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	null	ms_run[2]:scan=1530	1336.62	3	600.6369267	600.6197	K	I	424	439
PSM	ALLRLHQECEKLK	12	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-8	9-UNIMOD:4	ms_run[2]:scan=1079	885.62	3	527.6681284	527.6362	R	K	262	274
PSM	DWYPAHSR	13	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	13	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	13	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[2]:scan=778	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	14	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	31	null	ms_run[2]:scan=1389	1195.62	3	550.9319625	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	15	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	20	null	ms_run[2]:scan=1298	1104.62	3	542.970027	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	16	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	2	null	ms_run[2]:scan=3083	2876.08	2	1974.385729	1974.3984	R	M	692	728
PSM	TLTIVDTGIGMTK	17	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	23	null	ms_run[2]:scan=1326	1132.62	3	450.5760358	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	18	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	18	null	ms_run[2]:scan=3391	3184.08	2	2405.607754	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	19	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	null	ms_run[3]:scan=1062	1336.62	3	600.6002396	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	20	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	null	ms_run[3]:scan=1066	1327.08	2	956.9780096	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	21	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	9-UNIMOD:4	ms_run[3]:scan=611	885.62	3	527.6504265	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	22	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	38	null	ms_run[3]:scan=921	1195.62	3	550.9318951	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	23	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	null	ms_run[3]:scan=2615	2876.08	2	1974.385759	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	24	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	34	23-UNIMOD:35	ms_run[3]:scan=2350	2611.08	2	1788.29715	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	25	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-8	null	ms_run[3]:scan=858	1132.62	3	450.5828076	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	26	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	12-UNIMOD:35	ms_run[3]:scan=2923	3184.08	2	2405.609401	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	27	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	22	null	ms_run[4]:scan=2731	1336.62	3	600.6215365	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	28	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	24	null	ms_run[4]:scan=2735	1327.08	2	956.9769799	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	29	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	19	9-UNIMOD:4	ms_run[4]:scan=2280	885.62	3	527.6318602	527.6362	R	K	262	274
PSM	MNQSNASPTLDGLFR	30	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	40	null	ms_run[4]:scan=2590	1195.62	3	550.9331388	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	31	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[4]:scan=2499	1104.62	3	542.9431778	542.9716	K	V	91	104
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	32	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	38	null	ms_run[4]:scan=4284	2876.08	2	1974.432596	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	33	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	10	23-UNIMOD:35	ms_run[4]:scan=4019	2611.08	2	1788.308805	1788.2886	K	M	695	728
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	34	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	11	0-UNIMOD:35	ms_run[4]:scan=4592	3184.08	2	2405.604604	2405.6084	-	E	1	41
PSM	QTQTFTTYSDNQPGVL	35	P63017	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	26	null	ms_run[5]:scan=2031	1336.62	3	600.631332	600.6197	K	I	424	439
PSM	AVVNGYSASDTVGAGFAQAK	36	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	32	null	ms_run[5]:scan=2035	1327.08	2	956.9393913	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	37	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	37	9-UNIMOD:4	ms_run[5]:scan=1580	885.62	3	527.6306181	527.6362	R	K	262	274
PSM	DWYPAHSR	38	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	38	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	38	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	100	null	ms_run[5]:scan=1279	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	39	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-3	null	ms_run[5]:scan=1890	1195.62	3	550.9468041	550.935	-	R	1	15
PSM	LWPFQVINEAGKPK	40	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	-2	null	ms_run[5]:scan=1799	1104.62	3	542.9779284	542.9716	K	V	91	104
PSM	TLTIVDTGIGMTK	41	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	2	11-UNIMOD:35	ms_run[5]:scan=1827	1132.62	3	450.5803135	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	42	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	6	0-UNIMOD:35	ms_run[5]:scan=3892	3184.08	2	2405.625379	2405.6084	-	E	1	41
PSM	AVVNGYSASDTVGAGFAQAK	43	Q8K0U4	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	25	null	ms_run[6]:scan=1331	1327.08	2	956.9900237	956.9736	K	E	261	281
PSM	ALLRLHQECEKLK	44	Q61699	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	6	9-UNIMOD:4	ms_run[6]:scan=876	885.62	3	527.6434782	527.6362	R	K	262	274
PSM	DWYPAHSR	45	P14602	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	17	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	R	L	21	28
PSM	DWYPAHSR	45	Q340U4	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	2	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	K	E	143	150
PSM	DWYPAHSR	45	P16627	0	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	8	null	ms_run[6]:scan=575	571.08	2	516.21	516.2383	R	M	240	247
PSM	MNQSNASPTLDGLFR	46	P14602	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	12	null	ms_run[6]:scan=1186	1195.62	3	550.9368104	550.935	-	R	1	15
PSM	MIKLGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	47	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	29	null	ms_run[6]:scan=2880	2876.08	2	1974.40836	1974.3984	R	M	692	728
PSM	LGLGIDEDDPTVDDTSAAVTEEMPPLEGDDDTSR	48	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	3	23-UNIMOD:35	ms_run[6]:scan=2615	2611.08	2	1788.288824	1788.2886	K	M	695	728
PSM	TLTIVDTGIGMTK	49	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	36	11-UNIMOD:35	ms_run[6]:scan=1123	1132.62	3	450.5838508	450.583	R	A	88	100
PSM	MPEETQTQDQPMEEEEVETFAFQAEIAQLMSLIINTFYSNK	50	P07901	1	UniProtKB	2013_08	[MS,MS:1001207,Mascot,]	11	0-UNIMOD:35	ms_run[6]:scan=3188	3184.08	2	2405.619464	2405.6084	-	E	1	41"#;

const CASANOVO_V3_2_0_A: &str = r"MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Identification
MTD	description	Casanovo identification file MM_peptides_casanovo_AspN
MTD	software[1]	[MS, MS:1003281, Casanovo, 3.2.0]
MTD	psm_search_engine_score[1]	[MS, MS:1001143, search engine specific score for PSMs, ]
MTD	ms_run[1]-location	file:///storage-hdd/douwe/casanovo/MM_Peng0013/2232_23614_AspN.mgf
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	fixed_mod[1]-site	C
MTD	variable_mod[1]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[1]-site	N
MTD	variable_mod[2]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[2]-site	Q
MTD	variable_mod[3]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[3]-site	M
MTD	variable_mod[4]	[UNIMOD, UNIMOD:385, Ammonia-loss, ]
MTD	variable_mod[4]-site	N-term
MTD	variable_mod[5]	[UNIMOD, UNIMOD:5, Carbamyl, ]
MTD	variable_mod[5]-site	N-term
MTD	variable_mod[6]	[UNIMOD, UNIMOD:1, Acetyl, ]
MTD	variable_mod[6]-site	N-term
MTD	software[1]-setting[1]	model = /storage-hdd/douwe/.cache/casanovo/casanovo_massivekb_v3_0_0.ckpt
MTD	software[1]-setting[2]	config_filename = /storage-hdd/douwe/miniconda3/envs/casanovo_env/lib/python3.10/site-packages/casanovo/config.yaml
MTD	software[1]-setting[3]	random_seed = 454
MTD	software[1]-setting[4]	n_peaks = 150
MTD	software[1]-setting[5]	min_mz = 50.0
MTD	software[1]-setting[6]	max_mz = 2500.0
MTD	software[1]-setting[7]	min_intensity = 0.01
MTD	software[1]-setting[8]	remove_precursor_tol = 2.0
MTD	software[1]-setting[9]	max_charge = 10
MTD	software[1]-setting[10]	precursor_mass_tol = 50.0
MTD	software[1]-setting[11]	isotope_error_range = (0, 1)
MTD	software[1]-setting[12]	dim_model = 512
MTD	software[1]-setting[13]	n_head = 8
MTD	software[1]-setting[14]	dim_feedforward = 1024
MTD	software[1]-setting[15]	n_layers = 9
MTD	software[1]-setting[16]	dropout = 0.0
MTD	software[1]-setting[17]	dim_intensity = None
MTD	software[1]-setting[18]	custom_encoder = None
MTD	software[1]-setting[19]	max_length = 100
MTD	software[1]-setting[21]	n_log = 1
MTD	software[1]-setting[22]	tb_summarywriter = None
MTD	software[1]-setting[23]	warmup_iters = 100000
MTD	software[1]-setting[24]	max_iters = 600000
MTD	software[1]-setting[25]	learning_rate = 0.0005
MTD	software[1]-setting[26]	weight_decay = 1e-05
MTD	software[1]-setting[27]	train_batch_size = 32
MTD	software[1]-setting[28]	predict_batch_size = 1024
MTD	software[1]-setting[29]	n_beams = 5
MTD	software[1]-setting[30]	logger = None
MTD	software[1]-setting[31]	max_epochs = 30
MTD	software[1]-setting[32]	num_sanity_val_steps = 0
MTD	software[1]-setting[33]	train_from_scratch = True
MTD	software[1]-setting[34]	save_model = True
MTD	software[1]-setting[35]	model_save_folder_path = 
MTD	software[1]-setting[36]	save_weights_only = True
MTD	software[1]-setting[37]	every_n_train_steps = 50000
MTD	software[1]-setting[38]	n_workers = 32
PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end	opt_ms_run[1]_aa_scores
PSM	GSGGR	0	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.735283625125885	null	null	1	413.26611328125	433.21537216688	ms_run[1]:index=0	null	null	null	null	0.25406,0.17583,0.22368,0.18075,0.48925
PSM	SSSSSR	1	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7253233976662159	null	null	1	511.4718017578125	610.2790921668799	ms_run[1]:index=1	null	null	null	null	0.35708,0.27899,0.19515,0.17837,0.12239,0.51608
PSM	-17.027QVAPETPEPTKEELLEDVK	2	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.3708054140210152	null	null	3	714.0145874023438	712.3633197002134	ms_run[1]:index=2	null	null	null	null	0.99922,0.93201,0.96866,0.04530,0.99950,0.96677,0.99117,0.40099,0.13996,0.99918,0.85622,0.40775,0.11147,0.18180,0.36664,0.41324,0.12600,0.68006,0.99869,0.99927
PSM	LAALLL	3	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6445258930325508	null	null	1	497.4560852050781	613.42832516688	ms_run[1]:index=3	null	null	null	null	0.34676,0.41201,0.23252,0.58530,0.31114,0.24512
PSM	DLLD	4	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.8191229030489922	null	null	1	405.2995910644531	475.23985516688	ms_run[1]:index=4	null	null	null	null	0.11139,0.15085,0.38070,0.08057
PSM	VVPELLHLAEEFLGGGK	5	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	0.47498726253123846	null	null	3	603.673828125	603.3366053668799	ms_run[1]:index=5	null	null	null	null	0.93643,0.04722,0.89374,0.57967,0.26205,0.18705,0.42346,0.28914,0.97375,0.51307,0.16243,0.08823,0.36511,0.98292,0.22967,0.27461,0.86622
PSM	FSSFKPEEEEFFR	6	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.592856129201559	null	null	3	559.9384155273438	560.2631053668799	ms_run[1]:index=6	null	null	null	null	0.99435,0.60932,0.26332,0.20480,0.11726,0.98989,0.52577,0.23362,0.25185,0.25347,0.16739,0.11135,0.57049
PSM	SPEELER	7	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.643719999917916	null	null	1	815.2142944335938	859.4155871668801	ms_run[1]:index=7	null	null	null	null	0.93320,0.11841,0.21546,0.18948,0.28040,0.12495,0.63206
PSM	KEELER	8	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7448172097404797	null	null	1	741.1951293945312	803.42575816688	ms_run[1]:index=8	null	null	null	null	0.24676,0.13971,0.14702,0.23043,0.13259,0.63459
PSM	AAGAGYK	9	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7750454034124101	null	null	1	577.1265258789062	637.33040316688	ms_run[1]:index=9	null	null	null	null	0.16576,0.32812,0.14013,0.23805,0.35547,0.18960,0.15755
PSM	ELEKL	10	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.8146590083837509	null	null	1	607.2094116210938	631.36611816688	ms_run[1]:index=10	null	null	null	null	0.11403,0.25258,0.14389,0.30997,0.10624
PSM	SLRELER	11	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6358094683715275	null	null	1	889.2332763671875	902.50540516688	ms_run[1]:index=11	null	null	null	null	0.94682,0.16216,0.34944,0.21870,0.16618,0.20836,0.49767
PSM	EEEER	12	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7128205999732018	null	null	1	667.17578125	691.28932416688	ms_run[1]:index=12	null	null	null	null	0.23836,0.23336,0.10908,0.15487,0.70023
PSM	TLVLRP	13	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.558969931056102	null	null	1	610.184326171875	698.45592816688	ms_run[1]:index=13	null	null	null	null	0.14757,0.15985,0.46905,0.86864,0.96035,0.04072
PSM	ALLLER	14	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6632984789709251	null	null	1	681.2283935546875	714.4508511668799	ms_run[1]:index=14	null	null	null	null	0.48801,0.20084,0.25201,0.36695,0.08368,0.62872
PSM	SLSSLTVTLKPQQEAVGGK	15	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5129245866678263	null	null	3	649.6820068359375	648.36510636688	ms_run[1]:index=15	null	null	null	null	0.99338,0.59771,0.22010,0.44879,0.18548,0.22883,0.48091,0.95719,0.32937,0.38608,0.99748,0.56994,0.20906,0.26193,0.15275,0.23017,0.98838,0.02397,0.99291
PSM	ELLR	16	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6342670135200024	null	null	1	413.26611328125	530.32967316688	ms_run[1]:index=16	null	null	null	null	0.17438,0.53024,0.10526,0.65305
PSM	ELLR	17	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7994027156382799	null	null	1	511.47210693359375	530.32967316688	ms_run[1]:index=17	null	null	null	null	0.19489,0.14747,0.09985,0.36018
PSM	PETPTELQGPDDEEEPEVK	18	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	0.5162208829271165	null	null	3	714.0145263671875	713.6585723668801	ms_run[1]:index=18	null	null	null	null	0.39791,0.27597,0.52044,0.97308,0.38256,0.21317,0.13434,0.22852,0.97315,0.50977,0.54235,0.27661,0.42599,0.50252,0.34173,0.99481,0.99594,0.12229,0.99706
PSM	LLLLT	19	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7800469532608986	null	null	1	497.4565124511719	572.40176716688	ms_run[1]:index=19	null	null	null	null	0.16727,0.21436,0.29215,0.35245,0.07354
PSM	YHLD	20	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7207246962934732	null	null	1	405.2998962402344	547.25108916688	ms_run[1]:index=20	null	null	null	null	0.33060,0.19889,0.51925,0.06835
PSM	-17.027QYLLLGEDEDDDEDR	21	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6228780914098024	null	null	3	603.3400268554688	603.25500436688	ms_run[1]:index=21	null	null	null	null	0.99333,0.33241,0.39100,0.14232,0.17767,0.14522,0.97413,0.38020,0.43828,0.28830,0.41420,0.23734,0.25463,0.12410,0.17882,0.56202
PSM	+42.011REELER	22	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6575511566230229	null	null	1	815.2136840820312	873.44247116688	ms_run[1]:index=22	null	null	null	null	0.59242,0.23652,0.21066,0.18524,0.26954,0.12411,0.77867
PSM	C+57.021LLLLR	23	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6975801015893619	null	null	1	741.1948852539062	787.48585716688	ms_run[1]:index=23	null	null	null	null	0.34809,0.14675,0.12477,0.22432,0.20907,0.76152
PSM	GSGGGYK	24	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.720152270581041	null	null	1	577.126220703125	625.29401716688	ms_run[1]:index=24	null	null	null	null	0.44734,0.09965,0.18184,0.40142,0.18695,0.14544,0.49630
PSM	GHNPEKPK	25	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.593485345132649	null	null	1	889.2333984375	906.47919116688	ms_run[1]:index=25	null	null	null	null	0.99997,0.59514,0.24540,0.38361,0.27849,0.22307,0.11219,0.41425
PSM	LGLLLK	26	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6089287549257278	null	null	1	607.20947265625	656.47052416688	ms_run[1]:index=26	null	null	null	null	0.17084,0.95180,0.45231,0.27613,0.16371,0.33165
PSM	LEDLAK	27	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7859395916263262	null	null	1	667.176513671875	688.3875821668801	ms_run[1]:index=27	null	null	null	null	0.15912,0.12875,0.13305,0.20670,0.09751,0.55924
PSM	EELEK	28	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7058303788304329	null	null	1	610.1837768554688	647.32464716688	ms_run[1]:index=28	null	null	null	null	0.20940,0.16574,0.44619,0.12282,0.52669
PSM	FPELREELEELLR	29	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6951866172827207	null	null	3	560.2728881835938	558.3016717002133	ms_run[1]:index=29	null	null	null	null	0.72327,0.14837,0.17906,0.23369,0.14013,0.20696,0.18871,0.14257,0.18361,0.15823,0.39111,0.49464,0.77224
PSM	GSGGR	30	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7082008183002472	null	null	1	413.26611328125	433.21537216688	ms_run[1]:index=30	null	null	null	null	0.26761,0.20572,0.26775,0.12307,0.59484
PSM	ELLR	31	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7764322571456432	null	null	1	511.47210693359375	530.32967316688	ms_run[1]:index=31	null	null	null	null	0.16861,0.17017,0.14896,0.40652
PSM	HYYY	32	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7147024851292372	null	null	1	497.45660400390625	645.26674016688	ms_run[1]:index=32	null	null	null	null	0.29846,0.42501,0.31307,0.10465
PSM	RPLPLPSSSSSSSSSLLPLDK	33	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.630400477420716	null	null	3	714.0147705078125	719.05834636688	ms_run[1]:index=33	null	null	null	null	0.94600,0.80461,0.39638,0.46708,0.21855,0.21456,0.20268,0.17582,0.20519,0.16950,0.19152,0.17588,0.17793,0.13533,0.12823,0.17359,0.29554,0.99159,0.91215,0.44311,0.33636
PSM	AAAAK	34	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.8378844454884529	null	null	1	405.2995910644531	431.26126016688	ms_run[1]:index=34	null	null	null	null	0.20668,0.17113,0.15158,0.10015,0.18104
PSM	FLSATLK	35	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6625255367585591	null	null	1	741.1953125	779.46615816688	ms_run[1]:index=35	null	null	null	null	0.21596,0.22240,0.70341,0.23933,0.15220,0.09006,0.73896
PSM	SPEVVYR	36	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.4961845608694213	null	null	1	815.2138061523438	849.4464941668799	ms_run[1]:index=36	null	null	null	null	0.50407,0.23244,0.13364,0.96076,0.82636,0.08735,0.78208
PSM	C+57.021LGSNK	37	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.4578092197577158	null	null	1	577.1262817382812	678.3239361668799	ms_run[1]:index=37	null	null	null	null	0.16027,0.46032,0.98156,0.62097,0.17419,0.85582
PSM	VKDDDYR	38	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.603866423879351	null	null	1	889.2324829101562	910.4264871668801	ms_run[1]:index=38	null	null	null	null	0.64687,0.54058,0.22728,0.43317,0.29016,0.11410,0.52077
PSM	GLSSHGR	39	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5621655200208937	null	null	1	667.176513671875	713.36891216688	ms_run[1]:index=39	null	null	null	null	0.49872,0.23509,0.96360,0.66934,0.21727,0.12017,0.36065
PSM	LYRPK	40	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7804442211985588	null	null	1	607.20947265625	676.41407216688	ms_run[1]:index=40	null	null	null	null	0.18311,0.12588,0.11594,0.19753,0.47532";

const CASANOVO_V3_2_0_B: &str = r"MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Identification
MTD	description	Casanovo identification file EH3420
MTD	software[1]	[MS, MS:1003281, Casanovo, 3.2.0]
MTD	psm_search_engine_score[1]	[MS, MS:1001143, search engine specific score for PSMs, ]
MTD	ms_run[1]-location	file:///net/noble/vol1/home/rowanon/proj/2022_rowanon_casanovo-ne/results/rowanon/2023-01-04_greiff-predictions/run1_tryp/EH3420.mzML
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	fixed_mod[1]-site	C
MTD	variable_mod[1]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[1]-site	N
MTD	variable_mod[2]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[2]-site	Q
MTD	variable_mod[3]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[3]-site	M
MTD	variable_mod[4]	[UNIMOD, UNIMOD:5, Carbamyl, ]
MTD	variable_mod[4]-site	N-term
MTD	variable_mod[5]	[UNIMOD, UNIMOD:385, Ammonia-loss, ]
MTD	variable_mod[5]-site	N-term
MTD	variable_mod[6]	[UNIMOD, UNIMOD:1, Acetyl, ]
MTD	variable_mod[6]-site	N-term
MTD	software[1]-setting[1]	model = /net/noble/vol1/home/rowanon/proj/casanovo-tensorboard/casanovo_massivekb.ckpt
MTD	software[1]-setting[2]	config_filename = ../config.yaml
MTD	software[1]-setting[3]	random_seed = 454
MTD	software[1]-setting[4]	n_peaks = 150
MTD	software[1]-setting[5]	min_mz = 50.0
MTD	software[1]-setting[6]	max_mz = 2500.0
MTD	software[1]-setting[7]	min_intensity = 0.01
MTD	software[1]-setting[8]	remove_precursor_tol = 2.0
MTD	software[1]-setting[9]	max_charge = 10
MTD	software[1]-setting[10]	precursor_mass_tol = 50.0
MTD	software[1]-setting[11]	isotope_error_range = (0, 1)
MTD	software[1]-setting[12]	dim_model = 512
MTD	software[1]-setting[13]	n_head = 8
MTD	software[1]-setting[14]	dim_feedforward = 1024
MTD	software[1]-setting[15]	n_layers = 9
MTD	software[1]-setting[16]	dropout = 0.0
MTD	software[1]-setting[17]	dim_intensity = None
MTD	software[1]-setting[18]	custom_encoder = None
MTD	software[1]-setting[19]	max_length = 100
MTD	software[1]-setting[21]	n_log = 1
MTD	software[1]-setting[22]	tb_summarywriter = None
MTD	software[1]-setting[23]	warmup_iters = 100000
MTD	software[1]-setting[24]	max_iters = 600000
MTD	software[1]-setting[25]	learning_rate = 0.0005
MTD	software[1]-setting[26]	weight_decay = 1e-05
MTD	software[1]-setting[27]	train_batch_size = 32
MTD	software[1]-setting[28]	predict_batch_size = 500
MTD	software[1]-setting[29]	n_beams = 5
MTD	software[1]-setting[30]	logger = None
MTD	software[1]-setting[31]	max_epochs = 30
MTD	software[1]-setting[32]	num_sanity_val_steps = 0
MTD	software[1]-setting[33]	train_from_scratch = True
MTD	software[1]-setting[34]	save_model = True
MTD	software[1]-setting[35]	model_save_folder_path = 
MTD	software[1]-setting[36]	save_weights_only = True
MTD	software[1]-setting[37]	every_n_train_steps = 50000
MTD	software[1]-setting[38]	n_workers = 2
PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end	opt_ms_run[1]_aa_scores
PSM	GRPPLPHLDH	0	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.4461659178137779	null	null	3	383.86480712890625	380.20871870021335	ms_run[1]:index=0	null	null	null	null	0.99904,0.98964,0.84902,0.83884,0.15015,0.29074,0.37499,0.27303,0.24778,0.52510
PSM	MSSQVGSK	1	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5847856290638447	null	null	2	414.07562255859375	412.20255281688003	ms_run[1]:index=1	null	null	null	null	0.66640,0.17903,0.68226,0.17542,0.09976,0.49653,0.16292,0.85939
PSM	AEKEEEDL	2	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.62057808379177	null	null	2	482.4651794433594	481.71928681688	ms_run[1]:index=2	null	null	null	null	0.78819,0.11135,0.54621,0.33442,0.07876,0.21017,0.95678,0.00949
PSM	SLVVRLR	3	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6672497125608581	null	null	2	422.064208984375	421.7821618168801	ms_run[1]:index=3	null	null	null	null	0.55333,0.13755,0.25205,0.09221,0.24880,0.30100,0.74432
PSM	RFVC+57.021LGR	4	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.49826620944908684	null	null	2	456.0058898925781	454.2501723168801	ms_run[1]:index=4	null	null	null	null	0.64091,0.11705,0.78480,0.14245,0.30823,0.88537,0.63334
PSM	M+15.995GRFELSLGSK	5	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6490602075031017	null	null	2	548.1034545898438	620.82135531688	ms_run[1]:index=5	null	null	null	null	0.00007,0.99833,0.33735,0.11539,0.11253,0.34670,0.10920,0.14156,0.74239,0.13897,0.81784
PSM	C+57.021EREEELR	6	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.7006897050887346	null	null	2	556.0919799804688	560.75621231688	ms_run[1]:index=6	null	null	null	null	0.36012,0.10643,0.24204,0.12737,0.26283,0.09381,0.36487,0.83702
PSM	GRVRRVK	7	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6030784249305725	null	null	2	429.5518798828125	435.79085281688003	ms_run[1]:index=7	null	null	null	null	0.74036,0.14125,0.25231,0.36535,0.23227,0.15053,0.89638
PSM	GLRPLRPPSRPGR	8	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6124410377098963	null	null	2	724.0575561523438	729.94185081688	ms_run[1]:index=8	null	null	null	null	0.94275,0.25993,0.28488,0.46152,0.19087,0.15771,0.20314,0.93311,0.29306,0.45651,0.28668,0.14996,0.41815
PSM	ERGRDSGGRSR	9	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6129438640041784	null	null	2	616.4884033203125	616.81377281688	ms_run[1]:index=9	null	null	null	null	0.58939,0.24340,0.23923,0.27543,0.16281,0.92394,0.42089,0.26815,0.50779,0.11710,0.50949
PSM	PRPGGGGGGEPGGGGGDAERSR	10	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6521947007964958	null	null	2	966.6580810546875	968.9520588168799	ms_run[1]:index=10	null	null	null	null	0.85867,0.23986,0.17449,0.38569,0.33624,0.22499,0.21487,0.14844,0.16734,0.16227,0.87289,0.26140,0.27287,0.27244,0.12829,0.13885,0.85297,0.56897,0.47480,0.24787,0.11504,0.53247
PSM	SRPLEEAGRPR	11	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5561954067511992	null	null	2	635.0299682617188	634.34691731688	ms_run[1]:index=11	null	null	null	null	0.75938,0.39603,0.79880,0.19599,0.30921,0.25349,0.79919,0.08429,0.50902,0.32400,0.45244
PSM	GC+57.021SC+57.021GGLQGSC+57.021C+57.021DK	12	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6422870276229722	null	null	2	774.0467529296875	773.2810868168801	ms_run[1]:index=12	null	null	null	null	0.93265,0.14167,0.32720,0.26710,0.43062,0.20436,0.42087,0.11947,0.81677,0.11050,0.38653,0.21537,0.11314,0.52174
PSM	GC+57.021PSC+57.021C+57.021EC+57.021ELLC+57.021GC+57.021C+57.021R	13	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6873424681834877	null	null	2	1042.597900390625	1039.36090281688	ms_run[1]:index=13	null	null	null	null	0.73166,0.22968,0.25070,0.08074,0.48584,0.17011,0.13177,0.65195,0.12170,0.21424,0.18802,0.92038,0.16399,0.15734,0.14085,0.36356
PSM	RPGPSSSSSSSSSSAPGK	14	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6451315006448163	null	null	2	824.1489868164062	825.38990281688	ms_run[1]:index=14	null	null	null	null	0.89683,0.28996,0.57223,0.14994,0.39359,0.38325,0.46525,0.23898,0.28539,0.22619,0.19727,0.28087,0.20874,0.23145,0.12492,0.87455,0.47887,0.08935
PSM	GPPRPSPGPPRPPGR	15	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.579569403330485	null	null	2	758.5711059570312	761.42349131688	ms_run[1]:index=15	null	null	null	null	0.99910,0.22756,0.27339,0.18721,0.90413,0.25386,0.23608,0.34753,0.36826,0.22318,0.36666,0.89724,0.22042,0.21014,0.59169
PSM	-17.027QFSVSGSEDGERSR	16	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6134571285675823	null	null	2	765.5602416992188	762.33968281688	ms_run[1]:index=16	null	null	null	null	0.00020,0.58778,0.37558,0.50266,0.11268,0.95565,0.20904,0.16704,0.15134,0.59550,0.62753,0.41759,0.20506,0.19376,0.69672
PSM	GSPRPPPPPRPPSR	17	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6202651421938623	null	null	2	750.5152587890625	747.91804131688	ms_run[1]:index=17	null	null	null	null	0.73851,0.09462,0.31551,0.48421,0.18692,0.18814,0.39923,0.48567,0.32907,0.50670,0.88701,0.23143,0.13997,0.32929
PSM	RPPQSLPSPEDAERSR	18	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5631687184795737	null	null	2	908.072509765625	911.4637378168801	ms_run[1]:index=18	null	null	null	null	0.83167,0.74420,0.60327,0.13304,0.27519,0.25185,0.26741,0.30498,0.15340,0.21123,0.87107,0.89104,0.36671,0.20911,0.13862,0.73650
PSM	+42.011EAEEEEEEEEEEKEEEK	19	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6765503593617015	null	null	2	1101.1851806640625	1097.9295123168804	ms_run[1]:index=19	null	null	null	null	0.93779,0.59237,0.19215,0.25614,0.20161,0.24077,0.33458,0.32610,0.38515,0.38685,0.37736,0.26192,0.18614,0.11733,0.21038,0.15161,0.20326,0.46061
PSM	RPRPDTEEEEDSRLSR	20	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6700926953926682	null	null	2	992.6118774414062	986.4775688168801	ms_run[1]:index=20	null	null	null	null	0.53384,0.37535,0.35582,0.13048,0.24458,0.26124,0.21336,0.23119,0.14443,0.13065,0.63876,0.73796,0.23265,0.10991,0.18494,0.75336
PSM	SSRSSHSSESSSSSDVRSSR	21	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5779158882796764	null	null	2	1051.5819091796875	1049.9708383168804	ms_run[1]:index=21	null	null	null	null	0.97392,0.65507,0.27572,0.43297,0.24515,0.22152,0.70110,0.41574,0.16471,0.38502,0.28616,0.22942,0.25594,0.15603,0.72250,0.84691,0.39921,0.54224,0.23392,0.29843
PSM	VPGGGGGGGGGGGGGGGSGGGGGGGGR	22	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6125287766809817	null	null	2	903.584716796875	885.3865533168798	ms_run[1]:index=22	null	null	null	null	0.51293,0.49640,0.33280,0.39725,0.30999,0.29846,0.33957,0.35192,0.46107,0.56359,0.51661,0.49513,0.44691,0.37404,0.42918,0.38970,0.54826,0.12636,0.41805,0.44943,0.50954,0.31963,0.27578,0.26108,0.19641,0.14786,0.49376
PSM	SEPEEPQRDAAFRR	23	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5877793074718544	null	null	2	843.60888671875	844.41097431688	ms_run[1]:index=23	null	null	null	null	0.75880,0.17664,0.33714,0.12124,0.38651,0.96035,0.34342,0.29106,0.34307,0.59936,0.22635,0.28584,0.21184,0.72946
PSM	DGRPLPGRPDLC+57.021RR	24	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6785634181329182	null	null	2	832.1341552734375	832.9417223168799	ms_run[1]:index=24	null	null	null	null	0.74849,0.55994,0.19498,0.20355,0.14082,0.13893,0.36600,0.15867,0.08365,0.74193,0.67841,0.14991,0.16886,0.16598
PSM	C+57.021PPPRPPRPSSPLK	25	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6233795301190445	null	null	2	785.0274658203125	793.43520981688	ms_run[1]:index=25	null	null	null	null	0.67031,0.51539,0.20601,0.43983,0.16074,0.09547,0.22632,0.28171,0.10565,0.19079,0.12818,0.97483,0.69951,0.57795
PSM	+42.011EAAGGGGGRGGGGGRGGGGRSR	26	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6432399587786716	null	null	2	918.0604858398438	921.4447358168798	ms_run[1]:index=26	null	null	null	null	0.71619,0.65161,0.42305,0.20894,0.39550,0.39213,0.31993,0.30224,0.40998,0.12583,0.35375,0.34001,0.33388,0.26046,0.30194,0.35603,0.18964,0.22208,0.33585,0.34460,0.33047,0.13062,0.76075
PSM	ERPDTPGRDGERSR	27	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6385833248496056	null	null	2	816.657470703125	814.39839381688	ms_run[1]:index=27	null	null	null	null	0.45552,0.18270,0.11395,0.39086,0.50339,0.21874,0.18221,0.11569,0.92970,0.46744,0.29390,0.35873,0.13508,0.71191
PSM	-17.027+42.011LAEEEEC+57.021RRK	28	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.712846930164839	null	null	2	683.128662109375	672.8142588168799	ms_run[1]:index=28	null	null	null	null	0.00009,0.64063,0.17719,0.19013,0.23281,0.13668,0.12174,0.50690,0.41178,0.27486,0.21216,0.54086
PSM	+42.011TASSEMTC+57.021EEEEK	29	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6095318155629295	null	null	2	784.0296630859375	786.8056273168801	ms_run[1]:index=29	null	null	null	null	0.99064,0.82696,0.34518,0.13102,0.19718,0.15842,0.43968,0.78204,0.13375,0.32085,0.38014,0.30150,0.18321,0.27599
PSM	PRPPLSEA	30	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5530537907034159	null	null	2	436.2544860839844	433.74015981688	ms_run[1]:index=30	null	null	null	null	0.62877,0.76531,0.45465,0.15062,0.24482,0.73938,0.27920,0.31283
PSM	+43.006-17.027+42.011PGYDDDHHR	31	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.4876084069419526	null	null	3	383.8638000488281	393.8201817002133	ms_run[1]:index=31	null	null	null	null	0.00000,0.95002,0.33121,0.89904,0.49555,0.29615,0.29130,0.20755,0.20814,0.98874,0.96862
PSM	GGRPGGGGGSGGGGGRPS	32	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6072375426689784	null	null	2	698.10791015625	692.3272458168801	ms_run[1]:index=32	null	null	null	null	0.97156,0.61504,0.35207,0.26685,0.31678,0.39027,0.38986,0.26909,0.46568,0.20504,0.16174,0.85534,0.40802,0.22217,0.34207,0.26627,0.44149,0.13039
PSM	-17.027QPVPVPSSSLK	33	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5739845062295597	null	null	2	566.0758056640625	561.3136888168799	ms_run[1]:index=33	null	null	null	null	0.92381,0.22186,0.54376,0.18260,0.71731,0.13466,0.96534,0.62857,0.17496,0.14423,0.15169,0.32341
PSM	Q+0.984GRPPPPPRPR	34	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6462019762653984	null	null	2	563.580078125	628.35454631688	ms_run[1]:index=34	null	null	null	null	0.00017,0.99755,0.31729,0.25215,0.19422,0.53948,0.48055,0.26469,0.26080,0.18458,0.40031
PSM	KLAQLAGR	35	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6669316366314888	null	null	2	422.06390380859375	428.77179481688	ms_run[1]:index=35	null	null	null	null	0.47135,0.35369,0.29972,0.16174,0.14064,0.24095,0.18199,0.81447
PSM	GRLLLEK	36	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6349193624087742	null	null	2	414.07513427734375	414.76872031688004	ms_run[1]:index=36	null	null	null	null	0.98834,0.22413,0.12946,0.13231,0.14462,0.13583,0.80088
PSM	Q+0.984GGRGRGSGRSR	37	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.6022851674351841	null	null	2	548.1034545898438	616.32176581688	ms_run[1]:index=37	null	null	null	null	0.00953,0.94340,0.29242,0.19026,0.35413,0.39287,0.94287,0.38538,0.21414,0.31340,0.17550,0.55868
PSM	GRALRRK	38	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5563970633915493	null	null	2	429.5519104003906	428.78302781688006	ms_run[1]:index=38	null	null	null	null	0.94940,0.20312,0.16045,0.16421,0.29413,0.39738,0.93653
PSM	SLDSERGR	39	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5517001785337925	null	null	2	456.0073547363281	460.23322981688005	ms_run[1]:index=39	null	null	null	null	0.86180,0.33256,0.34989,0.51239,0.29720,0.32683,0.18481,0.72092
PSM	GRPLSGSGSGLL	40	null	null	null	null	[MS, MS:1003281, Casanovo, 3.2.0]	-0.5062050955990951	null	null	2	556.092529296875	550.80656231688	ms_run[1]:index=40	null	null	null	null	0.99992,0.72701,0.68849,0.10854,0.25012,0.46942,0.37918,0.22213,0.23691,0.27216,0.78839,0.78325";

const CASANOVO_V4_2_1: &str = r"MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Identification
MTD	description	Casanovo identification file 20191211_F1_Ag5_peng0013_SA_her_tryp_c421
MTD	software[1]	[MS, MS:1003281, Casanovo, 4.2.1]
MTD	psm_search_engine_score[1]	[MS, MS:1001143, search engine specific score for PSMs, ]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	fixed_mod[1]-site	C
MTD	variable_mod[1]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[1]-site	N
MTD	variable_mod[2]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[2]-site	Q
MTD	variable_mod[3]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[3]-site	M
MTD	variable_mod[4]	[UNIMOD, UNIMOD:5, Carbamyl, ]
MTD	variable_mod[4]-site	N-term
MTD	variable_mod[5]	[UNIMOD, UNIMOD:1, Acetyl, ]
MTD	variable_mod[5]-site	N-term
MTD	variable_mod[6]	[UNIMOD, UNIMOD:385, Ammonia-loss, ]
MTD	variable_mod[6]-site	N-term
MTD	software[1]-setting[1]	model = /home/douwe/.cache/casanovo/casanovo_v4_2_0_v4_2_0.ckpt
MTD	software[1]-setting[2]	config_filename = default
MTD	software[1]-setting[3]	precursor_mass_tol = 50.0
MTD	software[1]-setting[4]	isotope_error_range = (0, 1)
MTD	software[1]-setting[5]	min_peptide_len = 6
MTD	software[1]-setting[6]	predict_batch_size = 1024
MTD	software[1]-setting[7]	n_beams = 1
MTD	software[1]-setting[8]	top_match = 1
MTD	software[1]-setting[9]	accelerator = auto
MTD	software[1]-setting[10]	devices = None
MTD	software[1]-setting[11]	random_seed = 454
MTD	software[1]-setting[12]	n_log = 1
MTD	software[1]-setting[13]	tb_summarywriter = None
MTD	software[1]-setting[14]	save_top_k = 5
MTD	software[1]-setting[15]	model_save_folder_path = 
MTD	software[1]-setting[16]	val_check_interval = 50000
MTD	software[1]-setting[17]	n_peaks = 150
MTD	software[1]-setting[18]	min_mz = 50.0
MTD	software[1]-setting[19]	max_mz = 2500.0
MTD	software[1]-setting[20]	min_intensity = 0.01
MTD	software[1]-setting[21]	remove_precursor_tol = 2.0
MTD	software[1]-setting[22]	max_charge = 10
MTD	software[1]-setting[23]	dim_model = 512
MTD	software[1]-setting[24]	n_head = 8
MTD	software[1]-setting[25]	dim_feedforward = 1024
MTD	software[1]-setting[26]	n_layers = 9
MTD	software[1]-setting[27]	dropout = 0.0
MTD	software[1]-setting[28]	dim_intensity = None
MTD	software[1]-setting[29]	max_length = 100
MTD	software[1]-setting[30]	warmup_iters = 100000
MTD	software[1]-setting[31]	cosine_schedule_period_iters = 600000
MTD	software[1]-setting[32]	learning_rate = 0.0005
MTD	software[1]-setting[33]	weight_decay = 1e-05
MTD	software[1]-setting[34]	train_label_smoothing = 0.01
MTD	software[1]-setting[35]	train_batch_size = 32
MTD	software[1]-setting[36]	max_epochs = 30
MTD	software[1]-setting[37]	num_sanity_val_steps = 0
MTD	software[1]-setting[38]	calculate_precision = False
MTD	software[1]-setting[40]	n_workers = 8
MTD	ms_run[1]-location	file:///home/douwe/casanovo/20191211_F1_Ag5_peng0013_SA_her_tryp.mgf
PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end	opt_ms_run[1]_aa_scores
PSM	KGWASDEEAEK	1	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.4668271640936533	null	null	3.0	414.85345	417.1946917002133	ms_run[1]:index=0	null	null	null	null	0.71591,0.34916,0.33638,0.36535,0.42044,0.60458,0.75123,0.50584,0.69560,0.38215,0.50953
PSM	+42.011M+15.995EDENRFLR	2	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.30652115832675586	null	null	2.0	635.74146	634.2904193168799	ms_run[1]:index=1	null	null	null	null	0.82141,0.60502,0.83577,0.70414,0.70981,0.81483,0.66209,0.46569,0.54220,0.62549
PSM	SPC+57.021DARC+57.021DVGRK	3	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5030156717850611	null	null	3.0	476.1977	474.2188490335467	ms_run[1]:index=2	null	null	null	null	0.56810,0.52984,0.65495,0.63893,0.41976,0.49837,0.36805,0.40090,0.48802,0.33700,0.35295,0.46221
PSM	LPVKQAD	4	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.9869352504611015	null	null	2.0	385.72415	385.72397881688	ms_run[1]:index=3	null	null	null	null	0.98968,0.98792,0.98733,0.98845,0.98846,0.98750,0.97786
PSM	DEEREEEALRR	5	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5342201627790928	null	null	2.0	713.79224	716.3447683168799	ms_run[1]:index=4	null	null	null	null	0.72203,0.45166,0.33947,0.30247,0.33964,0.30754,0.32689,0.44858,0.67535,0.57022,0.37796
PSM	DHYQDKTPLGDGPV	6	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.12772533694903054	null	null	2.0	770.8556	771.3651663168799	ms_run[1]:index=5	null	null	null	null	0.83481,0.92799,0.92870,0.83146,0.92028,0.93129,0.93092,0.93126,0.52231,0.75507,0.86222,0.93383,0.92821,0.91446
PSM	KPEEDDKLLGLLLDEKLDEDLLRLEFEELLHHLLK	7	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.8000616195301214	null	null	5.0	1975.1912	843.0628148068802	ms_run[1]:index=6	null	null	null	null	0.45919,0.22945,0.27336,0.18997,0.19180,0.22367,0.22266,0.22601,0.24597,0.20008,0.20599,0.16544,0.17550,0.17820,0.19395,0.17535,0.17482,0.16696,0.18901,0.18706,0.16894,0.16576,0.17201,0.18574,0.16877,0.17020,0.16283,0.18874,0.17663,0.17180,0.19102,0.20026,0.18106,0.19693,0.23231
PSM	TLSFKSDYEKLRDLDKFE	8	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.6495024808927586	null	null	4.0	559.5382	559.28894014188	ms_run[1]:index=7	null	null	null	null	0.81375,0.81804,0.81889,0.53001,0.41595,0.42108,0.45895,0.41373,0.37215,0.44860,0.43550,0.69372,0.79552,0.82069,0.80306,0.82243,0.81808,0.82004
PSM	TRSLFKLGYKEEEEKEKFE	9	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.6104551821947097	null	null	4.0	598.5629	598.3114096418801	ms_run[1]:index=8	null	null	null	null	0.79291,0.79848,0.41749,0.46241,0.39550,0.44279,0.47228,0.51870,0.54209,0.47879,0.40346,0.61118,0.56477,0.79214,0.79903,0.52375,0.78824,0.80227,0.80234
PSM	+42.011TELRKLKVKKPLLENLYFQA	10	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.5719454837116328	null	null	4.0	619.3616	619.1189786418798	ms_run[1]:index=9	null	null	null	null	0.77429,0.77845,0.78180,0.78385,0.78019,0.66530,0.77463,0.54082,0.47907,0.38980,0.55945,0.39603,0.50663,0.39947,0.39379,0.36030,0.39004,0.44699,0.38259,0.51537,0.70288
PSM	SEEEEEEEEEEEEEEEEEEEEK	11	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5807497157998707	null	null	4.0	701.8371	704.5046303918804	ms_run[1]:index=10	null	null	null	null	0.53498,0.32173,0.32228,0.35423,0.53905,0.42973,0.41647,0.41693,0.40582,0.35936,0.37880,0.40567,0.42720,0.42694,0.41314,0.40319,0.40896,0.40266,0.39474,0.32338,0.31822,0.53773
PSM	+42.011GEVGPAAPGRPR	12	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.4188014236944062	null	null	3.0	399.56613	402.55103336688	ms_run[1]:index=11	null	null	null	null	0.48574,0.38720,0.52595,0.47738,0.44245,0.65413,0.65194,0.75592,0.71475,0.45976,0.47650,0.54988,0.76942
PSM	EFGKDDFFLK	13	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.6841302947564558	null	null	3.0	417.21204	415.8765230335466	ms_run[1]:index=12	null	null	null	null	0.47184,0.22433,0.23184,0.26870,0.37998,0.22394,0.24172,0.23062,0.28609,0.26409
PSM	GVQAEKPSPGQLR	14	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5900209349180972	null	null	3.0	455.2458	456.25276436687994	ms_run[1]:index=13	null	null	null	null	0.47264,0.32825,0.32768,0.69512,0.42923,0.28495,0.25458,0.32470,0.30028,0.47913,0.39899,0.35678,0.38782
PSM	EKLKQELLEC+57.021DK	15	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.33597011405688065	null	null	3.0	511.5984	511.60447470021336	ms_run[1]:index=14	null	null	null	null	0.43114,0.26920,0.25682,0.23363,0.26460,0.23104,0.22980,0.23213,0.44671,0.45800,0.34734,0.30437
PSM	VNPRPPLPT	16	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.37625644803047176	null	null	2.0	497.7893	495.7901798168801	ms_run[1]:index=15	null	null	null	null	0.78954,0.80705,0.46197,0.55043,0.51690,0.68502,0.47104,0.50272,0.64631
PSM	RLSLLADKEDEEEEEEEEEREESEK	17	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5005050049378321	null	null	3.0	1026.4467	1027.128499700214	ms_run[1]:index=16	null	null	null	null	0.48983,0.35205,0.30345,0.35904,0.41366,0.40397,0.34164,0.38032,0.56106,0.50701,0.48330,0.43197,0.45054,0.53029,0.57068,0.73196,0.43954,0.50300,0.43214,0.47644,0.58200,0.65782,0.41789,0.68160,0.74074
PSM	YQQKPGKAPK	18	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.9906805482777682	null	null	3.0	382.2206	382.22062470021336	ms_run[1]:index=17	null	null	null	null	0.99106,0.99036,0.99120,0.99078,0.99045,0.99148,0.99104,0.99076,0.99007,0.99039
PSM	C+57.021N+0.984NVNHKPSNTK	19	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.11052017945509696	null	null	3.0	476.90186	471.8894937002133	ms_run[1]:index=18	null	null	null	null	0.94051,0.94121,0.94146,0.94086,0.94049,0.94115,0.94070,0.93987,0.94063,0.93688,0.64209,0.57897
PSM	QQKPGKAPK	20	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.9901066958904267	null	null	2.0	491.29605	491.29563431688007	ms_run[1]:index=19	null	null	null	null	0.99137,0.99006,0.99039,0.99133,0.99002,0.99128,0.98936,0.98879,0.98814
PSM	YQQKPGKAPK	21	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.9904172420501709	null	null	2.0	572.82764	572.82729881688	ms_run[1]:index=20	null	null	null	null	0.99111,0.98924,0.99100,0.99079,0.99087,0.99159,0.99202,0.98938,0.98912,0.98976
PSM	VYQQKPGKAPK	22	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.05620688696702325	null	null	2.0	621.81116	622.36150581688	ms_run[1]:index=21	null	null	null	null	0.96814,0.96651,0.96427,0.96865,0.96909,0.96734,0.96793,0.96582,0.96274,0.96367,0.69509
PSM	LNGHYEGVVR	23	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.2613199366764589	null	null	3.0	386.19342	381.86836203354665	ms_run[1]:index=22	null	null	null	null	0.67371,0.65330,0.63003,0.86526,0.86538,0.86516,0.67550,0.85691,0.48104,0.69511
PSM	YQGEKPGKAPK	24	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.9531161785125732	null	null	3.0	401.55554	401.55578436688	ms_run[1]:index=23	null	null	null	null	0.96872,0.97300,0.97126,0.97169,0.95169,0.96628,0.95899,0.82297,0.94970,0.96106,0.97073
PSM	+43.006KNNVNHKPSNTKVDK	25	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.09528706967830658	null	null	3.0	590.96466	589.3149710335467	ms_run[1]:index=24	null	null	null	null	0.94769,0.94860,0.94759,0.94794,0.94646,0.94418,0.94911,0.94785,0.94705,0.94709,0.94554,0.94825,0.94559,0.57480,0.68221,0.86243
PSM	+42.011VTLQKPK	26	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.2976924777030945	null	null	2.0	427.26663	428.26854931688007	ms_run[1]:index=25	null	null	null	null	0.83636,0.70101,0.59224,0.64751,0.58062,0.63313,0.76378,0.72141
PSM	KGRPPWPK	27	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.26536740528212654	null	null	2.0	482.78302	483.2876118168801	ms_run[1]:index=26	null	null	null	null	0.83896,0.82870,0.57694,0.75531,0.77734,0.52745,0.66851,0.77593
PSM	PPTQTLLSSSK	28	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5439809982975323	null	null	2.0	578.78656	579.82186931688	ms_run[1]:index=27	null	null	null	null	0.51458,0.41971,0.37353,0.41510,0.32155,0.36392,0.69286,0.36295,0.31507,0.45638,0.51349
PSM	QQPPSREGSSSR	29	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5449346384176841	null	null	2.0	648.7706	658.32109631688	ms_run[1]:index=28	null	null	null	null	0.55401,0.45291,0.54894,0.53092,0.71566,0.33439,0.29271,0.31285,0.36908,0.55939,0.41774,0.59972
PSM	QAAHVHLNLK	30	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.564722177657214	null	null	3.0	376.54584	377.55248536688003	ms_run[1]:index=29	null	null	null	null	0.58587,0.28512,0.35858,0.28698,0.34009,0.33840,0.48669,0.46075,0.45128,0.49001
PSM	YYM+15.995EKESPK	31	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.49976964890956876	null	null	3.0	394.87015	397.51811870021334	ms_run[1]:index=30	null	null	null	null	0.73877,0.72003,0.38079,0.36266,0.34880,0.35428,0.37396,0.54213,0.43632
PSM	TLLPTKGDETHK	32	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.3270964244237313	null	null	3.0	449.5596	447.2453780335467	ms_run[1]:index=31	null	null	null	null	0.82856,0.50023,0.44879,0.48431,0.47000,0.49420,0.77733,0.82265,0.82804,0.60577,0.82536,0.83141
PSM	+43.006GPAQVNHKPSNTK	33	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.13784322341283162	null	null	2.0	714.84985	710.8682098168799	ms_run[1]:index=32	null	null	null	null	0.92786,0.92476,0.92761,0.92677,0.92637,0.92712,0.92772,0.92176,0.91945,0.62280,0.92321,0.91419,0.57050,0.64587
PSM	YTVHGLRNTPSELLRNPK	34	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5817308190621828	null	null	4.0	523.2892	524.5406273918799	ms_run[1]:index=33	null	null	null	null	0.67938,0.40029,0.27915,0.37005,0.28451,0.31680,0.32418,0.37489,0.66163,0.31667,0.38008,0.32834,0.47863,0.49250,0.53643,0.34708,0.34254,0.32973
PSM	TLSKGEEYKHK	35	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	0.8675841962297758	null	null	3.0	440.56982	440.5699787002133	ms_run[1]:index=34	null	null	null	null	0.92861,0.92590,0.91094,0.57876,0.81192,0.82651,0.92971,0.79138,0.92464,0.92545,0.92812
PSM	LQEHPNLFNLER	36	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5211820252812825	null	null	3.0	503.59927	503.9318017002133	ms_run[1]:index=35	null	null	null	null	0.65360,0.30553,0.37917,0.30400,0.29878,0.34979,0.35247,0.66186,0.60661,0.48835,0.40691,0.68408
PSM	ATKPAEPAAPA	37	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.29470343639453256	null	null	2.0	508.78036	512.27710281688	ms_run[1]:index=36	null	null	null	null	0.49155,0.80568,0.59178,0.83962,0.84927,0.70907,0.56691,0.83208,0.81822,0.50034,0.61119
PSM	TPTTFTLSSLK	38	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.6585477354625862	null	null	2.0	601.83105	598.33206131688	ms_run[1]:index=37	null	null	null	null	0.49606,0.24439,0.25802,0.25530,0.24457,0.25336,0.24752,0.42064,0.27782,0.35260,0.38190
PSM	LNGGNNHTGEK	39	null	null	null	null	[MS, MS:1003281, Casanovo, 4.2.1]	-0.5278554943700631	null	null	2.0	579.2828	570.77324631688	ms_run[1]:index=38	null	null	null	null	0.59246,0.34809,0.47594,0.72999,0.35148,0.35545,0.30616,0.42239,0.51386,0.50807,0.34118";

const CONTRANOVO_V1_0_0: &str = r"MTD	mzTab-version	1.0.0
MTD	mzTab-mode	Summary
MTD	mzTab-type	Identification
MTD	description	Casanovo identification file 20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp_formatted4_pointnovo
MTD	software[1]	[MS, MS:1003281, Casanovo, 0.1]
MTD	psm_search_engine_score[1]	[MS, MS:1001143, search engine specific score for PSMs, ]
MTD	fixed_mod[1]	[UNIMOD, UNIMOD:4, Carbamidomethyl, ]
MTD	fixed_mod[1]-site	C
MTD	variable_mod[1]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[1]-site	N
MTD	variable_mod[2]	[UNIMOD, UNIMOD:7, Deamidated, ]
MTD	variable_mod[2]-site	Q
MTD	variable_mod[3]	[UNIMOD, UNIMOD:35, Oxidation, ]
MTD	variable_mod[3]-site	M
MTD	variable_mod[4]	[UNIMOD, UNIMOD:5, Carbamyl, ]
MTD	variable_mod[4]-site	N-term
MTD	variable_mod[5]	[UNIMOD, UNIMOD:1, Acetyl, ]
MTD	variable_mod[5]-site	N-term
MTD	variable_mod[6]	[UNIMOD, UNIMOD:385, Ammonia-loss, ]
MTD	variable_mod[6]-site	N-term
MTD	software[1]-setting[1]	peak_path = ../test_data/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp_formatted4_pointnovo.mgf
MTD	software[1]-setting[2]	model = ContraNovo/ControNovo.ckpt
MTD	software[1]-setting[3]	config_filename = /home/auke/ContraNovo/ContraNovo/config.yaml
MTD	software[1]-setting[4]	random_seed = 200
MTD	software[1]-setting[5]	n_peaks = 300
MTD	software[1]-setting[6]	min_mz = 50.5
MTD	software[1]-setting[7]	max_mz = 4500.0
MTD	software[1]-setting[8]	min_intensity = 0.0
MTD	software[1]-setting[9]	remove_precursor_tol = 2.0
MTD	software[1]-setting[10]	max_charge = 10
MTD	software[1]-setting[11]	precursor_mass_tol = 50.0
MTD	software[1]-setting[12]	isotope_error_range = (0, 1)
MTD	software[1]-setting[13]	dim_model = 512
MTD	software[1]-setting[14]	n_head = 8
MTD	software[1]-setting[15]	dim_feedforward = 1024
MTD	software[1]-setting[16]	n_layers = 9
MTD	software[1]-setting[17]	dropout = 0.18
MTD	software[1]-setting[18]	dim_intensity = None
MTD	software[1]-setting[19]	custom_encoder = None
MTD	software[1]-setting[20]	max_length = 100
MTD	software[1]-setting[22]	n_log = 1
MTD	software[1]-setting[23]	tb_summarywriter = None
MTD	software[1]-setting[24]	enable_neptune = True
MTD	software[1]-setting[25]	neptune_project = DeNovo/clip
MTD	software[1]-setting[26]	neptune_api_token = None
MTD	software[1]-setting[27]	tags = ['9-speice', 'bacillus', 'Lr = 0.0002 dp 0.15,0.4']
MTD	software[1]-setting[28]	n_nodes = 1
MTD	software[1]-setting[29]	train_from_resume = False
MTD	software[1]-setting[30]	warmup_iters = None
MTD	software[1]-setting[31]	max_iters = None
MTD	software[1]-setting[32]	max_epochs = 150
MTD	software[1]-setting[33]	warm_up_epochs = 1
MTD	software[1]-setting[34]	learning_rate = 0.0004
MTD	software[1]-setting[35]	weight_decay = 1e-05
MTD	software[1]-setting[36]	gradient_clip_val = 1.5
MTD	software[1]-setting[37]	gradient_clip_algorithm = norm
MTD	software[1]-setting[38]	accumulate_grad_batches = 1
MTD	software[1]-setting[39]	sync_batchnorm = False
MTD	software[1]-setting[40]	SWA = False
MTD	software[1]-setting[41]	train_batch_size = 16
MTD	software[1]-setting[42]	predict_batch_size = 512
MTD	software[1]-setting[43]	n_beams = 5
MTD	software[1]-setting[44]	logger = None
MTD	software[1]-setting[45]	num_sanity_val_steps = 0
MTD	software[1]-setting[46]	train_from_scratch = True
MTD	software[1]-setting[47]	save_model = True
MTD	software[1]-setting[48]	model_save_folder_path = ./clipcasa
MTD	software[1]-setting[49]	save_weights_only = True
MTD	software[1]-setting[50]	every_n_train_steps = 2500
MTD	software[1]-setting[51]	n_workers = 8
MTD	ms_run[1]-location	file:///home/auke/test_data/20230408_F1_UM4_Peng0013_SA_EXT00_her_01_tryp_formatted4_pointnovo.mgf
PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end	opt_ms_run[1]_aa_scores
PSM	KLEEEELQKTEEQQLEDKKEEEEEEEEWNKFDKDC+57.021LYSLSTGST	1	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.44305550374768	null	null	4	1353.1163330078125	1353.1176041418798	ms_run[1]:0	null	null	null	null	0.97496,0.77788,0.30081,0.36956,0.15416,0.31847,0.47165,0.30486,0.09371,0.41130,0.22300,0.17485,0.59561,0.12926,0.24904,0.20860,0.24107,0.28298,0.20957,0.33537,0.34288,0.12794,0.23933,0.39911,0.27495,0.69186,0.34929,0.33876,0.18698,0.40428,0.76069,0.93643,0.17847,0.98660,0.38581,0.28683,0.14022,0.99872,0.75276,0.99989,0.34070,0.54571,0.99998,0.99954
PSM	TTTGQEEEDKLTVKWEYEEEEKKEEEEEEEKEERPEEESLSSAST	2	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.514036755470766	null	null	4	1349.11279296875	1348.85249139188	ms_run[1]:0	null	null	null	null	0.99973,0.98592,0.91572,0.47818,0.48555,0.26641,0.29031,0.75039,0.11121,0.37211,0.53405,0.56592,0.11919,0.57706,0.53705,0.20620,0.55043,0.54644,0.47721,0.69136,0.37119,0.13717,0.25137,0.42755,0.34534,0.36220,0.45057,0.64610,0.69246,0.27098,0.42418,0.79126,0.58493,0.40822,0.11786,0.60609,0.90366,0.02768,0.98699,0.31950,0.99640,0.96313,0.32538,0.99977,0.26121
PSM	LSTTSDTLKQEEWEYAFKEDKLELEELEEDSKDFNKDSDNYTSGST	3	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.41807935386896133	null	null	4	1349.3675537109375	1349.36403189188	ms_run[1]:0	null	null	null	null	0.99998,0.96658,0.85805,0.43231,0.27559,0.40721,0.16042,0.19993,0.29198,0.12498,0.20795,0.31189,0.09612,0.18092,0.24717,0.39385,0.16277,0.42527,0.17207,0.24282,0.23854,0.15551,0.15816,0.17264,0.40371,0.20304,0.31455,0.16245,0.30171,0.49564,0.67912,0.16103,0.44414,0.06670,0.53131,0.68985,0.24201,0.84669,0.43852,0.05140,0.77743,0.99855,0.95179,0.85403,0.99993,0.73530
PSM	+43.006LKTEDLKEEEEEEEEEEEDEEKLELEEEKLDKLEDHLDSLTSGTS	4	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.3557239900464597	null	null	4	1348.869384765625	1348.8640658918803	ms_run[1]:0	null	null	null	null	1.00000,0.43703,0.65362,0.24261,0.14579,0.17645,0.30894,0.34765,0.15184,0.19118,0.19583,0.21353,0.21952,0.16750,0.28247,0.20801,0.23009,0.36673,0.18248,0.18102,0.26396,0.20715,0.28872,0.11996,0.21875,0.16509,0.28270,0.21085,0.29798,0.16432,0.17202,0.17211,0.28874,0.21747,0.30020,0.40727,0.19853,0.56933,0.29309,0.99702,0.32247,0.99880,0.98003,0.98018,0.99990,0.34439
PSM	VSSSANTSALELLFRPLLEDLDDDEDVDKEDSDKDC+57.021LDDDDFSLGTSST	5	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.4131583256685004	null	null	4	1353.1104736328125	1353.1046931418791	ms_run[1]:0	null	null	null	null	0.99983,0.84582,0.52171,0.43607,0.18439,0.22451,0.69929,0.72149,0.36816,0.32263,0.06194,0.36924,0.36457,0.07038,0.31396,0.64794,0.26952,0.31114,0.14229,0.21552,0.21006,0.20938,0.19417,0.12712,0.26782,0.26856,0.12600,0.30162,0.07219,0.25216,0.34055,0.39511,0.19573,0.46949,0.36790,0.58072,0.10060,0.54010,0.67349,0.64274,0.14173,0.22838,0.99958,0.77735,0.99987,0.30467,0.60988,0.99996,0.75742
PSM	ATTATAEALDKLEEWYAYLHSDLKDDEEEEEEKGDEDDKLDSLSSAST	6	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.43965837576737005	null	null	4	1348.856201171875	1348.8563078918794	ms_run[1]:0	null	null	null	null	0.99999,0.94427,0.66264,0.21490,0.39670,0.31395,0.32434,0.18998,0.35301,0.33828,0.52177,0.23028,0.20585,0.34579,0.15691,0.22099,0.55226,0.70039,0.27695,0.25140,0.35208,0.22060,0.21861,0.53947,0.55232,0.21715,0.30818,0.29232,0.38844,0.29181,0.22936,0.12010,0.62431,0.18308,0.43247,0.28199,0.26540,0.84728,0.14235,0.52586,0.72477,0.99544,0.53985,0.99978,0.96912,0.06582,0.99984,0.57516
PSM	LGTSQETEWKEEQFLAYAEEDKDAEDWKLEEKLEFDHEDTVTSGTS	7	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.4712260632735232	null	null	4	1349.1146240234375	1349.1127463918804	ms_run[1]:0	null	null	null	null	0.99996,0.87432,0.58201,0.98789,0.36327,0.20180,0.70200,0.27643,0.42053,0.24658,0.20502,0.19292,0.23883,0.27728,0.19913,0.64066,0.97691,0.49822,0.22642,0.43439,0.36102,0.35314,0.32060,0.11910,0.70667,0.19191,0.28517,0.23340,0.37241,0.26007,0.19294,0.59855,0.29867,0.13327,0.63639,0.48699,0.08726,0.27886,0.71112,0.98700,0.09813,0.97997,0.99585,0.96247,0.99988,0.48104
PSM	+43.006GSYGYFC+57.021MSFTSPRPPQSSSSSSSYYR	8	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.7111643430377755	null	null	4	779.8430786132812	779.8341831418801	ms_run[1]:0	null	null	null	null	1.00000,0.99928,0.61962,0.52775,0.51840,0.86637,0.53336,0.87908,0.85622,0.53942,0.99859,0.90330,0.44729,0.84114,0.91095,0.12435,0.98211,0.25006,0.51117,0.97233,0.45908,0.74452,0.60175,0.83218,0.74001,0.25748,0.99967,0.99712
PSM	LC+57.021N+0.984VNHKPSNTKVDK	9	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.9884099205334981	null	null	4	439.4770812988281	439.47680364188005	ms_run[1]:0	null	null	null	null	0.99999,0.99998,0.99994,0.99995,0.83759,0.99999,0.99999,0.99999,0.99998,0.99024,0.99996,0.99995,0.99992,0.99901,0.99967
PSM	HSKC+57.021YYGFPHQRYEEQYNSTYR	10	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.8234122578393329	null	null	4	729.0728759765625	729.07540264188	ms_run[1]:0	null	null	null	null	1.00000,0.83961,0.19054,0.95173,0.88183,0.34916,0.98129,0.98188,0.91361,0.60749,0.87774,0.91961,0.77453,0.27814,0.95546,0.94664,0.70434,0.96627,0.99985,0.99800,0.99998,0.99738
PSM	+43.006HEEQC+57.021YYSDTKPNREYYWNSTYR	11	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.7043743760635456	null	null	4	783.8417358398438	783.8348756418801	ms_run[1]:0	null	null	null	null	1.00000,0.99939,0.52497,0.23755,0.28598,0.87075,0.91967,0.82282,0.99843,0.99613,0.75445,0.63387,0.90550,0.15781,0.65341,0.14815,0.18379,0.33193,0.59467,0.89373,0.99833,0.99970,0.99994,0.99403
PSM	+43.006EDKYMC+57.021SDC+57.021GKPPSGC+57.021TC+57.021TVDHKPSNTK	12	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.7504192369649837	null	null	5	661.2783813476562	661.2786024068799	ms_run[1]:0	null	null	null	null	1.00000,0.77037,0.41518,0.34094,0.09274,0.68083,0.84236,0.97499,0.10502,0.97049,0.92313,0.03512,0.94022,0.77193,0.52368,0.81130,0.87513,0.81807,0.99240,0.00015,0.99995,0.99804,0.99999,0.99989,0.99999,0.99991,0.97226,0.99996,0.90811
PSM	LTC+57.021YYGEPEQVREEQNMNSTYR	13	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.8204205296933651	null	null	4	692.5576171875	692.5576173918801	ms_run[1]:0	null	null	null	null	0.99992,0.88430,0.79301,0.81179,0.94428,0.95777,0.47206,0.98215,0.29375,0.66168,0.98606,0.89750,0.98240,0.75751,0.70679,0.83775,0.10310,0.98099,0.99978,0.99953,0.99990,0.99722
PSM	LLLEPEEEEEEEEEEAGEEEEEEEEEEEEEEEEEEELQ+0.984RTEEEEAKEEALKR	14	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.7150767401147348	null	null	5	1280.1556396484375	1279.9337194068805	ms_run[1]:0	null	null	null	null	0.99969,0.57868,0.81944,0.57340,0.29217,0.73978,0.19699,0.91966,0.77690,0.58262,0.79238,0.87379,0.93634,0.93448,0.16881,0.30979,0.87217,0.84084,0.82761,0.65280,0.98648,0.97120,0.97184,0.98567,0.35210,0.98990,0.91811,0.99589,0.97367,0.90031,0.67746,0.91392,0.76999,0.40198,0.77688,0.78105,0.82760,0.60323,0.79419,0.08626,0.84870,0.60835,0.92646,0.68036,0.92797,0.70585,0.54158,0.83825,0.51298,0.25737,0.09268,0.87740
PSM	+43.006C+57.021YFC+57.021TDFDSTLRGEFTNMVSTYR	15	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.6927950376023849	null	null	3	971.7594604492188	971.7473990335467	ms_run[1]:0	null	null	null	null	1.00000,0.92055,0.99963,0.95899,0.78109,0.72074,0.97728,0.99953,0.54296,0.49947,0.29890,0.50482,0.97972,0.69091,0.37223,0.44133,0.18766,0.29558,0.72051,0.16928,0.98987,0.58448,0.99818,0.99336
PSM	AEGTGTGYSC+57.021DFDLNTVVEYEC+57.021RPGYR	16	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.8605344596284407	null	null	3	1039.453857421875	1039.4501903668802	ms_run[1]:0	null	null	null	null	0.99999,0.86460,0.99762,0.20947,0.99990,0.89011,0.99965,0.97949,0.84940,0.86107,0.99703,0.99240,0.98224,0.97638,0.03128,0.96294,0.99747,0.99902,0.99543,0.94746,0.93533,0.99426,0.83987,0.72686,0.31452,0.89586,0.99478
PSM	LC+57.021N+0.984VNHKPSNTKVDKK	17	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.9946432411670685	null	null	5	377.40289306640625	377.40189080688003	ms_run[1]:0	null	null	null	null	0.99999,0.99977,0.99950,0.99997,0.93791,0.99997,1.00000,1.00000,0.99999,0.99745,0.99997,0.99999,0.99970,0.99577,0.98451,0.99980
PSM	QGFTHGSSSSSSSYGGMDDYRDSSSSSSYR	18	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.7401801645755768	null	null	5	634.659423828125	634.6592528068801	ms_run[1]:0	null	null	null	null	0.99973,0.96678,0.96881,0.89851,0.13645,0.91631,0.21930,0.95760,0.79096,0.69497,0.62939,0.94909,0.35126,0.95382,0.57183,0.99624,0.79961,0.53239,0.22228,0.73127,0.95536,0.75359,0.87812,0.89091,0.87671,0.92214,0.94474,0.26590,0.48543,0.94589
PSM	-17.027QANC+57.021WGYTR	322	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	0.8814870677888393	null	null	2	569.7229614257812	569.7403618168801	ms_run[1]:0	null	null	null	null	1.00000,1.00000,0.99930,0.91362,0.97486,0.11735,0.99734,0.99371,0.82866,0.99004
PSM		1894	null	null	null	null	[MS, MS:1003281, Casanovo, 0.1]	nan	null	null	8	1906.1671142578125	3.25859705438	ms_run[1]:0	null	null	null	null	";
