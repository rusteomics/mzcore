macro_rules! test {
    ($case:literal, $name:ident) => {
        #[test]
        fn $name() {
            let mut peptides = 0;
            for read in crate::MZTabData::parse_reader(
                $case.as_bytes(),
                &mzcore::ontology::STATIC_ONTOLOGIES,
            ) {
                if let Ok(_) = read {
                    peptides += 1;
                }
            }
            println!("{peptides}");
        }
    };
}

test!(
    r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	ms_run[0]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]",
    invalid_ms_run_1
);

test!(
    r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[0]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]",
    invalid_ms_run_2
);

test!(
    r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[0]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]",
    invalid_ms_run_3
);

test!(
    r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end
PSM	LFDQAFGLPR	17699	IPI00025512	null	IPI human	2.31	null	null	null	null	null	null	null	ms_run[0]:spectrum=17699	null	null	28	37",
    invalid_ms_run_4
);

test!(
    r"MTD	mzTab-version	1.0 rc5
MTD	mzTab-mode	Complete
MTD	mzTab-type	Identification
MTD	mzTab-ID	1643
MTD	ms_run[1]-format	[MS, MS:1000564, PSI mzData file, ]
MTD	ms_run[1]-location	ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2005/12/PRD000001/PRIDE_Exp_Complete_Ac_1643.xml.gz
MTD	ms_run[1]-id_format	[MS, MS:1000777, spectrum identifier nativeID format, ]

PSH	sequence	PSM_ID	accession	unique	database	database_version	search_engine	search_engine_score[1]	modifications	retention_time	charge	exp_mass_to_charge	calc_mass_to_charge	spectra_ref	pre	post	start	end
PSM	QQVLDR	1661	223462890	null	NCBInr_2010_10	nr_101020.fasta	[MS, MS:1001207, Mascot, ]	37.76	8-MOD:01499	null	1	902.482117	902.518133	ms_run[1]:spectrum=1661	R	Y	20	25",
    invalid_mod_location
);

