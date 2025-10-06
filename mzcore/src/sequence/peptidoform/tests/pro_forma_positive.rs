use crate::parse_test;

parse_test!("AA", positive_example_1);
parse_test!("A[+1]", positive_example_2);
parse_test!("AA[+1]", positive_example_3);
//parse_test!("A(AAAA)[+1][+2]", positive_example_4);
parse_test!("UWAKJDNLASNOIJPojkjjdakjn[U:Oxidation]", positive_example_5);
parse_test!("[+1]-A[+1]-[+1]", positive_example_6);
parse_test!("AA+AA", positive_example_7);
parse_test!(
    "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
    positive_example_8
);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_9
);
parse_test!("EM[Oxidation]EVEES[Phospho]PEK", positive_example_10);
parse_test!(
    "EM[R:L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    positive_example_11
);
parse_test!("EMEVTK[X:DSS#XL1]SESPEK", positive_example_12);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK", positive_example_13);
parse_test!("EM[+15.9949]EVEES[+79.9663]PEK", positive_example_14);
parse_test!("EM[U:+15.995]EVEES[U:+79.966]PEK", positive_example_15);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_16);
parse_test!("RTAAX[+367.0537]WT", positive_example_17);
parse_test!(
    "{Glycan:Hex}EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_18
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
    positive_example_19
);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_20
);
parse_test!(
    "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
    positive_example_21
);
//parse_test!("<[MOD:01090]@C>ATPEILTCNSIGCLK", positive_example_22);
//parse_test!("[Phospho]?EM[Oxidation]EVTSESPEK", positive_example_23);
//parse_test!(
//    "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
//    positive_example_24
//);
parse_test!(
    "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
    positive_example_25
);
parse_test!(
    "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
    positive_example_26
);
// parse_test!(
//     "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.90)]ES[#s1(0.90)]PEK",
//     positive_example_27
// );
// parse_test!("PROT(EOSFORMS)[+19.0523]ISK", positive_example_28);
// parse_test!(
//     "PROT(EOC[Carbamidomethyl]FORMS)[+19.0523]ISK",
//     positive_example_29
// );
parse_test!("SEQUEN[Formula:C12H20O2]CE", positive_example_30);
parse_test!("SEQUEN[Formula:HN-1O2]CE", positive_example_31);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", positive_example_32);
parse_test!("SEQUEN[Glycan:HexNAc]CE", positive_example_33);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", positive_example_34);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", positive_example_35);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_36
);
// parse_test!(
//     "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
//     positive_example_37
// ); This modification cannot be placed on N termini (in the database), so invalid
parse_test!(
    "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_38
);
parse_test!("ELVIS[Phospho|+79.966331]K", positive_example_39);
parse_test!("ELVIS[Phospho|Obs:+79.978]K", positive_example_40);
parse_test!("ELV[INFO:xxxxx]IS", positive_example_41);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
    positive_example_42
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    positive_example_43
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
    positive_example_44
);
parse_test!("<13C>ATPEILTVNSIGQLK", positive_example_45);
parse_test!("EMEVEESPEK/2", positive_example_46);
parse_test!("EMEVEESPEK+ELVISLIVER", positive_example_47);
parse_test!("EMEVEESPEK/2+ELVISLIVER/3", positive_example_48);
parse_test!(
    "A[X:DSS#XL1]//B[#XL1]+C[X:DSS#XL1]//D[#XL1]",
    positive_example_49
);
parse_test!("<[Carbamidomethyl]@C>ATPEILTCNSIGCLK", positive_example_50);
parse_test!("<[Oxidation]@C,M>MTPEILTCNSIGCLK", positive_example_51);
parse_test!("<[TMT6plex]@K,N-term>ATPEILTCNSIGCLK", positive_example_52);
parse_test!(
    "<[TMT6plex]@K,N-term:A>ATPEILTCNSIGCLK",
    positive_example_53
);
parse_test!(
    "<[TMT6plex]@K,N-term:A,N-term:B>ATPEILTCNSIGCLK",
    positive_example_54
);
parse_test!("EM[Oxidation]EVEES[Phospho]PEK", positive_example_55);
parse_test!(
    "EM[L-methionine sulfoxide]EVEES[O-phospho-L-serine]PEK",
    positive_example_56
);
parse_test!(
    "EM[R:L-methionine sulfone]EVEES[O-phospho-L-serine]PEK",
    positive_example_57
);
parse_test!("EMEVTK[X:DSS#XL1]SESPEK", positive_example_58);
parse_test!("NEEYN[GNO:G59626AS]K", positive_example_59);
parse_test!("NEEYN[G:G59626AS]K", positive_example_60);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK", positive_example_61);
parse_test!(
    "EM[M:L-methionine sulfoxide]EVEES[M:O-phospho-L-serine]PEK",
    positive_example_62
);
parse_test!(
    "EM[U:Oxidation]EVEES[M:O-phospho-L-serine]PEK",
    positive_example_63
);
parse_test!(
    "EM[Oxidation]EVEES[O-phospho-L-serine]PEK",
    positive_example_64
);
parse_test!(
    "EM[Oxidation]EVE[Cation:Mg[II]]ES[Phospho]PEK",
    positive_example_65
);
parse_test!("EM[MOD:00719]EVEES[MOD:00046]PEK", positive_example_66);
parse_test!("EM[UNIMOD:35]EVEES[UNIMOD:56]PEK", positive_example_67);
parse_test!(
    "EM[RESID:AA0581]EVEES[RESID:AA0037]PEK",
    positive_example_68
);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK[#XL1]", positive_example_69);
parse_test!(
    "EMK[XLMOD:02000#XL1]EVTKSE[XLMOD:02010#XL2]SK[#XL1]PEK[#XL2]AR",
    positive_example_70
);
parse_test!("EMEVTK[XLMOD:02001#XL1]SESPEK", positive_example_71);
parse_test!("EMEVTK[XLMOD:02001]SESPEK", positive_example_72);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[XLMOD:02001#XL1]SESPEK",
    positive_example_73
);
parse_test!(
    "SEK[XLMOD:02001#XL1]UENCE//EMEVTK[#XL1]SESPEK",
    positive_example_74
);
parse_test!("EVTSEKC[MOD:00034#XL1]LEMSC[#XL1]EFD", positive_example_75);
parse_test!(
    "EVTSEKC[L-cystine (cross-link)#XL1]LEMSC[#XL1]EFD",
    positive_example_76
);
parse_test!(
    "FVNQHLC[MOD:00034#XL1]GSHLVEALYLVC[MOD:00034#XL2]GERGFFYTPKA//GIVEQC[MOD:00034#XL3]C[#XL1]TSIC[#XL3]SLYQLENYC[#XL2]N",
    positive_example_77
);
parse_test!(
    "EVTSEKC[XLMOD:02009#XL1]LEMSC[#XL1]EFD",
    positive_example_79
);
parse_test!(
    "EVTSEKC[X:Disulfide#XL1]LEMSC[#XL1]EFD",
    positive_example_80
);
parse_test!(
    "EVTSEKC[half cystine]LEMSC[half cystine]EFD",
    positive_example_81
);
parse_test!(
    "EVTSEKC[MOD:00798]LEMSC[MOD:00798]EFDEVTSEKC[MOD:00798]LEMSC[MOD:00798]EFD",
    positive_example_82
);
parse_test!("EVTSEKC[UNIMOD:374#XL1]LEMSC[#XL1]EFD", positive_example_83);
parse_test!("EVTSEKC[Dehydro#XL1]LEMSC[#XL1]EFD", positive_example_84);
// parse_test!(
//     "ETFGD[MOD:00093#BRANCH]//R[#BRANCH]ATER",
//     positive_example_85
// ); This modification cannot be placed on N termini (in the database), so invalid
parse_test!(
    "AVTKYTSSK-[MOD:00134#BRANCH]//AGKQLEDGRTLSDYNIQKESTLHLVLRLRG-[#BRANCH]",
    positive_example_86
);
parse_test!("NEEYN[GNO:G59626AS]K", positive_example_87);
parse_test!(
    "YPVLN[GNO:G62765YT]VTMPN[GNO:G02815KT]NSNGKFDK",
    positive_example_88
);
parse_test!("EM[+15.9949]EVEES[+79.9663]PEK", positive_example_89);
parse_test!("EM[+15.995]EVEES[-18.01]PEK", positive_example_90);
parse_test!("EM[U:+15.9949]EVEES[U:+79.9663]PEK", positive_example_91);
parse_test!("EM[U:+15.995]EVEES[U:+79.966]PEK", positive_example_92);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_93);
parse_test!("EM[U:+15.995]EVEES[Obs:+79.978]PEK", positive_example_94);
parse_test!("RTAAX[+367.0537]WT", positive_example_95);
parse_test!("SEQUEN[Formula:C12H20O2]CE", positive_example_96);
parse_test!("SEQUEN[Formula:[13C2]CH6N]CE", positive_example_97);
parse_test!("SEQUEN[Formula:[13C2][12C-2]H2N]CE", positive_example_98);
parse_test!("SEQUEN[Glycan:HexNAc1Hex2]CE", positive_example_99);
parse_test!(
    "[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK",
    positive_example_100
);
parse_test!(
    "[iTRAQ4plex]-EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_101
);
parse_test!(
    "{Glycan:Hex}EM[U:Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_102
);
parse_test!(
    "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]",
    positive_example_103
);
parse_test!(
    "{Glycan:Hex}[iTRAQ4plex]-EM[Oxidation]EVNES[Phospho]PEK[iTRAQ4plex]-[Methyl]",
    positive_example_104
);
parse_test!("{Glycan:Hex}{Glycan:NeuAc}EMEVNESPEK", positive_example_105);
// parse_test!("[Phospho]?EM[Oxidation]EVTSESPEK", positive_example_106);
// parse_test!(
//     "[Phospho][Phospho]?[Acetyl]-EM[Oxidation]EVTSESPEK",
//     positive_example_107
// );
// parse_test!(
//     "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
//     positive_example_108
// );
// parse_test!(
//     "[Phospho]^2?[Acetyl]-EM[Oxidation]EVTSESPEK",
//     positive_example_109
// );
parse_test!(
    "EM[Oxidation]EVT[#g1]S[#g1]ES[Phospho#g1]PEK",
    positive_example_110
);
// parse_test!("PRT(ESFRMS)[+19.0523]ISK", positive_example_111);
// parse_test!(
//     "PRT(EC[Carbamidomethyl]FRMS)[+19.0523]ISK",
//     positive_example_112
// );
parse_test!(
    "EM[Oxidation]EVT[#g1(0.01)]S[#g1(0.09)]ES[Phospho#g1(0.90)]PEK",
    positive_example_113
);
// parse_test!(
//     "[Phospho#s1]?EM[Oxidation]EVT[#s1(0.01)]S[#s1(0.09)]ES[#s1(0.90)]PEK",
//     positive_example_114
// );
// parse_test!("MPGLVDSNPAPPESQEKKPLK(PCCACPETKKARDACIIEKGEEHCGHLIEAHKECMRALGFKI)[Oxidation][Oxidation][half cystine][half cystine]", positive_example_115);
parse_test!("<13C>ATPEILTVNSIGQLK", positive_example_116);
parse_test!("<15N>ATPEILTVNSIGQLK", positive_example_117);
parse_test!("<D>ATPEILTVNSIGQLK", positive_example_118);
parse_test!("<13C><15N>ATPEILTVNSIGQLK", positive_example_119);
parse_test!(
    "<[S-carboxamidomethyl-L-cysteine]@C>ATPEILTCNSIGCLK",
    positive_example_120
);
// parse_test!("<[MOD:01090]@C>ATPEILTCNSIGCLK", positive_example_121);
parse_test!("<[Oxidation]@C,M>MTPEILTCNSIGCLK", positive_example_122);
// parse_test!(
//     "<[MOD:01090]@C>[Phospho]?EM[Oxidation]EVTSECSPEK",
//     positive_example_123
// );
// parse_test!(
//     "<[MOD:01090]@C>[Acetyl]-EM[Oxidation]EVTSECSPEK",
//     positive_example_124
// );
parse_test!(
    "(?DQ)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_125
);
parse_test!(
    "(?N)NGTWEM[Oxidation]ESNENFEGYM[Oxidation]K",
    positive_example_126
);
parse_test!("ELV[INFO:AnyString]IS", positive_example_127);
parse_test!("ELV[info:AnyString]IS", positive_example_128);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered]K",
    positive_example_129
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:really awesome]K",
    positive_example_130
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created on 2021-06]K",
    positive_example_131
);
parse_test!(
    "ELVIS[Phospho|INFO:newly discovered|INFO:Created by software Tool1]K",
    positive_example_132
);
parse_test!("ELVIS[U:Phospho|+79.966331]K", positive_example_133);
parse_test!("ELVIS[U:Phospho|Obs:+79.978]K", positive_example_134);
parse_test!("ELVIS[Phospho|O-phospho-L-serine]K", positive_example_135);
parse_test!("ELVIS[UNIMOD:21|MOD:00046]K", positive_example_136);
parse_test!("ELVIS[UNIMOD:21|Phospho]K", positive_example_137);
parse_test!(
    "ELVIS[Phospho|O-phospho-L-serine|Obs:+79.966]K",
    positive_example_138
);
parse_test!("ELVIS[Obs:+79.966|Phospho|Sulfo]K", positive_example_139);
parse_test!("EMEVEESPEK/2", positive_example_140);
parse_test!("EM[U:Oxidation]EVEES[U:Phospho]PEK/3", positive_example_141);
parse_test!(
    "[U:iTRAQ4plex]-EM[U:Oxidation]EVNES[U:Phospho]PEK[U:iTRAQ4plex]-[U:Methyl]/3",
    positive_example_142
);
parse_test!("EMEVEESPEK/3[+2Na+,+H+]", positive_example_143);
parse_test!("EMEVEESPEK/1[+2Na+,-H+]", positive_example_144);
parse_test!("EMEVEESPEK/-2[2I-]", positive_example_145);
// parse_test!("EMEVEESPEK/-1[+e-]", positive_example_146);
parse_test!("EMEVEESPEK/2+ELVISLIVER/3", positive_example_147);
parse_test!("AA(?AA)", positive_example_148);
parse_test!("AA(?AA)AA", positive_example_149);
// parse_test!("[dehydro]^3?[gln->pyro-glu]-QSC", positive_example_150);
parse_test!("[deamidated#1]-FEEAQ[#1]A", positive_example_151);
parse_test!("[#1]-FEEAQ[deamidated#1]A", positive_example_152);
parse_test!("AHAM[oxidation#1]TEG-[#1]", positive_example_153);
parse_test!("AHAM[#1]TEG-[oxidation#1]", positive_example_154);
