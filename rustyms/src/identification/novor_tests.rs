#![allow(clippy::missing_panics_doc)]
use std::io::BufReader;

use crate::identification::{test_format, NovorData, NovorVersion};

#[test]
fn novor_old_denovo() {
    match test_format::<NovorData>(
        BufReader::new(DATA_OLD_DENOVO.as_bytes()),
        None,
        false,
        false,
        Some(NovorVersion::OldDenovo),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn novor_new_denovo() {
    match test_format::<NovorData>(
        BufReader::new(DATA_NEW_DENOVO.as_bytes()),
        None,
        false,
        false,
        Some(NovorVersion::NewDenovo),
    ) {
        Ok(n) => assert_eq!(n, 20),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

#[test]
fn novor_new_psm() {
    match test_format::<NovorData>(
        BufReader::new(DATA_NEW_PSM.as_bytes()),
        None,
        false,
        false,
        Some(NovorVersion::NewPSM),
    ) {
        Ok(n) => assert_eq!(n, 19),
        Err(e) => {
            println!("{e}");
            panic!("Failed identified peptides test");
        }
    }
}

const DATA_OLD_DENOVO: &str = r"Fraction,Scan #,m/z,z,Score,Peptide Mass,Error (ppm),Length,De Novo Peptide,DB Sequence
F1,18017,553.821533203125,2,97.9,1105.6284,0.1,9,LLLFWASTR,
F1,16407,561.81787109375,2,97.8,1121.6233,-1.9,9,LLLYWASTR,LLIYWASTR
F1,15003,454.247955322266,2,97.5,906.4811,0.3,8,DGFLQSLK,N(Deamidated)GFIQSLK
F1,16683,561.81884765625,2,97.4,1121.6233,-0.2,9,LLLYWASTR,LLIYWASTR
F1,7242,380.184997558594,2,96.9,758.3559,-0.6,7,ESGVPDR,ESGVPDR
F1,7490,380.184967041016,2,96.9,758.3559,-0.7,7,ESGVPDR,ESGVPDR
F1,9191,389.220947265625,2,96.9,776.428,-0.8,7,ELSSLTK,
F1,6345,391.741821289063,2,96.8,781.4698,-0.9,7,VDKPVPK,
F1,7346,399.852203369141,3,96.7,1196.5357,-0.7,9,SFQRSEC(Cam)QR,
F1,13671,565.310729980469,2,96.6,1128.6027,3.8,9,LVELEEELR,IVELEEELR
F1,11334,470.265930175781,2,96.6,938.5185,-1.3,8,QLPSPLER,
F1,11217,453.71923828125,2,96.4,905.4243,-0.4,8,ESGVPDFR,
F1,15375,577.813537597656,2,96.3,1153.6165,-3.5,10,LLLYSMTGTR,
F1,9423,438.734649658203,2,96.1,875.4535,1.4,8,MVTVPSSR,
F1,14966,619.295288085938,2,95.9,1236.575,0.8,9,TSYWMHWVK,TSYWMHWVK
F1,11965,451.720794677734,2,95.8,901.4269,0.2,6,WM(O)HWVK,
F1,7806,464.2369,2,95.8,926.457,2.4,8,PNLNADQR,
F1,13719,451.720458984375,2,95.7,901.4269,-0.6,6,WM(O)HWVK,
F1,9345,496.262512207031,2,95.6,990.5094,1.1,9,SSQSLLNSR,
F1,9309,409.229461669922,2,95.6,816.4454,-1.2,7,QSLLNSR,";

const DATA_NEW_DENOVO: &str = r"# id, scanNum, RT, mz(data), z, pepMass(denovo), err(data-denovo), ppm(1e6*err/(mz*z)), score, peptide, aaScore 
0, 10234, 3058.1, 706.3292, 4, 2821.2907, -0.0032, -1.1, 24.8, MM(O)EKGPHTLNSPGRGMFM(O)TTLLDPE, 11-15-39-35-2-4-3-1-5-44-2-1-6-3-2-34-46-6-1-3-37-57-68-79-91
1, 8236, 2472.1, 543.2781, 3, 1626.8154, -0.0029, -1.8, 73.1, DLAPGTYLHWVREA, 25-53-25-3-43-52-94-96-99-94-98-97-90-90
2, 1846, 587.9, 446.7351, 2, 891.4491, 0.0065, 7.3, 40.9, FLWANEL, 69-42-63-79-8-24-1
3, 6786, 2051.4, 394.6718, 2, 787.3244, 0.0046, 5.8, 58.1, TYM(O)NML, 90-94-27-13-32-99
4, 13791, 4112.9, 671.0255, 3, 2010.0496, 0.0049, 2.4, 37.6, VSYASVVALLTTLHCAYTV, 18-1-28-2-43-44-25-1-3-71-63-35-57-76-3-3-62-60-70
5, 1193, 387.7, 433.2311, 2, 864.4415, 0.0060, 7.0, 50.1, TKMPEFL, 73-80-88-15-2-6-99
6, 4089, 1266.8, 550.6080, 3, 1648.8110, -0.0088, -5.3, 47.8, LPHPWATHDVDHPK, 49-52-37-2-1-1-1-6-80-89-98-99-85-82
7, 13600, 4054.7, 489.2448, 2, 976.4866, -0.0115, -11.7, 48.6, FADPQLSSL, 57-33-1-4-25-96-62-58-99
8, 7476, 2251.5, 507.2335, 3, 1518.6885, -0.0100, -6.6, 59.5, SVVGTHMHEVTHAD, 55-58-58-46-1-1-32-30-96-93-98-99-99-93
9, 6088, 1841.9, 574.9354, 3, 1721.7728, 0.0116, 6.7, 38.5, TAALGMDFVM(O)CWHKP, 78-85-72-99-73-21-1-1-1-1-2-74-99-8-2
10, 8861, 2652.9, 498.9156, 3, 1493.7184, 0.0064, 4.3, 74.2, DKTAYLEMAPERA, 36-79-99-99-99-99-99-85-88-1-2-92-90
11, 10093, 3017.3, 581.8321, 2, 1161.6434, 0.0062, 5.4, 43.9, PGSVFLFPVSL, 4-6-29-62-85-98-92-44-1-1-1
12, 8021, 2410.0, 411.5299, 3, 1231.5656, 0.0022, 1.8, 45.4, ESYLHCPVRE, 10-10-80-98-80-1-1-26-34-82
13, 7196, 2170.5, 537.5913, 3, 1609.7559, -0.0039, -2.4, 56.3, SLRLSECHSYYVPG, 62-85-95-97-76-11-1-1-32-67-87-87-51-3
14, 14034, 4187.9, 531.2891, 2, 1060.5627, 0.0010, 0.9, 59.4, LASVMLFPPS, 25-26-25-15-22-91-99-89-95-99
15, 14076, 4201.0, 759.3683, 3, 2275.0575, 0.0256, 11.2, 26.9, M(O)VDM(O)AGVLVCKGFYPSDLAVE, 12-32-1-1-4-1-1-6-2-36-32-17-28-54-3-3-42-80-48-47-87
16, 2935, 918.8, 572.8251, 2, 1143.6288, 0.0069, 6.0, 42.0, SGKWLPETKV, 15-7-49-37-31-89-2-1-91-91
17, 8699, 2605.9, 498.6289, 5, 2488.1052, 0.0030, 1.2, 45.4, DVFSPCANHVVWMHEALHNDPA, 43-45-77-3-4-21-25-47-50-34-55-52-55-63-59-50-63-59-31-49-53-16
18, 7042, 2125.7, 497.7650, 4, 1987.0064, 0.0245, 12.3, 46.8, SFFLQHLNKAGPASRWE, 40-86-98-93-22-1-1-11-1-2-7-1-15-79-70-88-89
19, 11745, 3499.8, 558.0223, 4, 2228.0651, -0.0050, -2.2, 21.4, VKFDWYVNPNGHVVSGGEGPA, 18-4-44-53-44-28-25-24-3-1-6-1-3-43-4-2-6-46-54-3-3";

const DATA_NEW_PSM: &str = r"#id, spectraId, scanNum, RT, mz, z, pepMass, err, ppm, score, protein, start, length, origin, peptide, noPTMPeptide, aac, allProteins
1, 14, 14034, 4187.9, 531.2891, 2, 1060.5593, 0.0044, 4.1, 2.3662, 6, 5, 10, APSVFIFPPS, APSVFIFPPS, APSVFIFPPS, 2-5-17-18-23-52-52-52-38-37, 5@6;112@7
2, 20, 9784, 2928.1, 429.2599, 3, 1284.7554, 0.0024, 1.9, 5.3597, 1, 183, 11, YRVVSVLTVLH, YRVVSVLTVLH, YRVVSVLTVLH, 36-37-52-52-43-52-52-52-42-44-44, 183@1;302@2;180@4
3, 21, 4244, 1311.6, 452.5754, 3, 1354.6980, 0.0063, 4.6, 6.1235, 6, 70, 12, STLTLSKADYEK, STLTLSKADYEK, STLTLSKADYEK, 36-44-33-52-48-52-35-43-38-39-36-28, 70@6;177@7
4, 24, 5873, 1780.7, 446.9401, 3, 1337.7919, 0.0065, 4.9, 6.0435, 1, 210, 13, ALPAPIEKTISKA, ALPAPIEKTISKA, ALPAPIEKTISKA, 33-37-32-33-17-24-52-48-52-37-33-52-29, 210@1;329@2;203@8
5, 32, 12831, 3828.9, 591.7968, 2, 1181.5757, 0.0034, 2.8, 3.3047, 1, 154, 9, PEVKFNWYV, PEVKFN(Deamidated)WYV, PEVKFNWYV, 44-40-52-40-38-40-52-52-44, 154@1;273@2
6, 34, 302, 100.1, 401.9044, 3, 1202.6871, 0.0042, 3.5, 2.8449, 4, 212, 11, IEKTISKAKGQ, IEKTISKAKGQ(Deamidated), IEKTISKAKGQ, 22-20-22-52-3-10-7-45-52-48-29, 212@4;215@1;334@2
7, 36, 5985, 1813.0, 485.7864, 2, 969.5535, 0.0047, 4.8, 2.8170, 13, 288, 8, FLYSKLTV, FLYSKLTV, FLYSKLTV, 44-44-32-52-18-3-10-16, 288@13;335@5;288@1;284@3;407@2;287@43
8, 42, 4010, 1243.8, 520.5062, 4, 2077.9891, 0.0067, 3.2, 3.5192, 6, 79, 18, YEKHKVYACEVTHQGLSS, YEKHKVYACEVTHQGLSS, YEKHKVYACEVTHQGLSS, 21-11-32-42-32-32-24-21-9-7-9-14-5-6-3-5-4-18, 79@6;186@7
9, 43, 10443, 3118.7, 503.2653, 3, 1506.7606, 0.0135, 8.9, 5.5589, 6, 7, 13, SVFIFPPSDEQLK, SVFIFPPSDEQ(Deamidated)LK, SVFIFPPSDEQLK, 36-44-52-52-52-52-26-30-35-52-52-29-30, 7@6;114@7
10, 47, 2780, 872.3, 618.3141, 2, 1234.6081, 0.0055, 4.4, 4.2817, 1, 274, 11, YKTTPPVLDSD, YKTTPPVLDSD, YKTTPPVLDSD, 32-32-2-2-38-27-43-52-35-32-39, 274@1;393@2;271@4
11, 54, 10063, 3008.6, 528.7844, 4, 2111.0952, 0.0133, 6.3, 9.0608, 2, 32, 18, YTIHWVRQAPGKGLEWVA, YTIHWVRQ(Deamidated)APGKGLEWVA, YTIHWVRQAPGKGLEWVA, 19-38-25-29-39-33-52-52-27-32-33-52-36-39-52-34-31-37, 32@2
12, 57, 4199, 1298.6, 434.7617, 2, 867.5066, 0.0023, 2.6, 2.9065, 1, 212, 8, PAPIEKTI, PAPIEKTI, PAPIEKTI, 44-44-33-32-27-52-28-32, 212@1;331@2;205@8;210@13;259@5;208@3;206@76;206@77
13, 58, 3296, 1033.3, 414.5737, 3, 1240.6928, 0.0064, 5.1, 4.2274, 6, 34, 10, PREAKVQWKV, PREAKVQ(Deamidated)WKV, PREAKVQWKV, 37-37-40-52-40-45-52-52-40-25, 34@6;141@7
14, 66, 3836, 1193.3, 472.6027, 3, 1414.7820, 0.0043, 3.0, 6.7772, 1, 3, 14, TKGPSVFPLAPSSK, TKGPSVFPLAPSSK, TKGPSVFPLAPSSK, 44-31-15-22-10-25-52-52-38-36-44-37-45-44, 3@1;122@2
15, 79, 8566, 2567.4, 504.2584, 2, 1006.4971, 0.0051, 5.1, 2.6684, 3, 52, 10, TFPAVLQSSG, TFPAVLQ(Deamidated)SSG, TFPAVLQSSG, 29-37-52-34-36-38-35-1-2-1, 52@3;52@4;52@1;171@2;52@5
16, 83, 2920, 914.4, 585.5712, 4, 2338.2533, 0.0025, 1.1, 4.3492, 1, 212, 21, PAPIEKTISKAKGQPREPQVY, PAPIEKTISKAKGQ(Deamidated)PREPQ(Deamidated)VY, PAPIEKTISKAKGQPREPQVY, 44-13-6-20-33-6-5-32-32-19-24-23-25-4-20-27-6-17-28-28-27, 212@1;331@2
17, 90, 1299, 420.4, 394.2071, 3, 1179.5884, 0.0110, 9.3, 4.4532, 6, 95, 11, SSPVTKSFNRG, SSPVTKSFN(Deamidated)RG, SSPVTKSFNRG, 36-44-25-28-18-33-52-31-52-33-22, 95@6;202@7
18, 100, 7456, 2245.7, 444.5684, 3, 1330.6769, 0.0065, 4.9, 4.1649, 6, 38, 11, KVQWKVDNALQ, KVQ(Deamidated)WKVDN(Deamidated)ALQ(Deamidated), KVQWKVDNALQ, 44-37-52-32-52-40-39-52-52-24-39, 38@6;145@7
19, 105, 4125, 1277.3, 523.2606, 4, 2089.0004, 0.0130, 6.2, 1.8430, 74, 543, 18, EKVTDILDVWFDSGVTHE, EKVTDILDVWFDSGVTHE, EKVTDILDVWFDSGVTHE, 36-44-4-1-7-8-3-6-5-7-8-7-3-2-33-30-26-27, 543@74;543@75";
