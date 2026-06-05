use std::{error::Error, fmt::Display};

use context_error::*;
use mzcv::curie;
use serde::{Deserialize, Serialize};

/// Error type to indicate that the given name was not a recognised species from the IMGT database
#[derive(Clone, Copy, Debug, Default)]
pub struct NotASpecies;

impl Display for NotASpecies {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Not a species in the IMGT database")
    }
}

impl Error for NotASpecies {}

macro_rules! species {
    ($($identifier:ident, $common:expr, $imgt:expr, $scientific:expr, $cv:ident:$curie:literal)*) => {
        /// All species available in the IMGT dataset. Look at the main documentation to see which actually have data provided.
        #[derive(Clone, Copy, Debug, bincode::Decode, Deserialize, bincode::Encode, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
        #[non_exhaustive]
        pub enum Species {
            $(
            #[doc = $imgt]
            $identifier,
            )*
        }

        impl Display for Species {
            fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                write!(f, "{}", self.common_name())
            }
        }

        impl Species {
            /// The common name for this species, e.g. `Human`
            pub const fn common_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $common,)*
                }
            }

            /// The name IMGT uses to identify this species, e.g. `Homo sapiens (human)`
            pub const fn imgt_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $imgt,)*
                }
            }

            /// The common name for this species, e.g. `Homo sapiens`
            pub const fn scientific_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $scientific,)*
                }
            }

            /// The enum name for this species, e.g. `HomoSapiens`
            pub const fn ident(&self) -> &'static str {
                match self {
                    $(Self::$identifier => stringify!($identifier),)*
                }
            }

            /// The CURIE to identify this species within the `NCBITaxonomy` Ontology
            pub const fn curie(&self) -> mzcv::Curie {
                match self {
                    $(Self::$identifier => curie!($cv:$curie),)*
                }
            }

            /// Get the species name from IMGT name tag.
            /// # Errors
            /// `Err` when the name could not be recognised (it is case-sensitive).
            /// `Ok(None)` when it is recognised as a species used by IMGT, but it is not a proper species (vector/plasmid etc.).
            pub fn from_imgt(s: &str) -> Result<Option<Self>, NotASpecies> {
                match s {
                    $($imgt => Ok(Some(Self::$identifier)),)*
                    "synthetic construct" | "synthetic construct (synthetic construct)" |
                    "unidentified" | "unclassified sequences" | "unidentified cloning vector" |
                    "Cloning vector AbVec-hIgG1" |
                    "Cloning vector AbVec-hIgKappa" |
                    "Cloning vector pASK88-huHRS3-VH-EP3/1" |
                    "Cloning vector pchiIGHG1" |
                    "Cloning vector pchiIGKC" |
                    "Cloning vector pCL" |
                    "Cloning vector pCLZip" |
                    "Cloning vector pMAB136" |
                    "Cloning vector pUR4546" |
                    "Cloning vector pUR4585" |
                    "Cloning vector pZGT5" |
                    "Expression vector p28BIOH-LIC4" |
                    "Expression vector pFUSE-HEAVY" |
                    "Expression vector pFUSE-hFc2-adapt-scFv" |
                    "Expression vector pFUSE-LIGHT" |
                    "Expression vector pFUSE-mFc2-adapt-scFv" |
                    "Expression vector pFUSE-rFc2-adapt-scFv" |
                    "Expression vector pHIN-PEP" |
                    "Expression vector pHIN-TRI" |
                    "Expression vector pSFV4" |
                    "Expression vector pTH-HIN" |
                    "Hepacivirus hominis" |
                    "Phagemid vector pGALD7" |
                    "Phagemid vector pGALD7DL" |
                    "Phagemid vector pGALD7DLFN" |
                    "Phagemid vector pGALD9" |
                    "Phagemid vector pGALD9DL" |
                    "Phagemid vector pGALD9DLFN" |
                    "Phagemid vector pMID21" |
                    "Enterobacteria phage M13 vector DY3F63"
                        => Ok(None),
                    _ => Err(NotASpecies),
                }
            }
        }

        /// The list with all usable names
        const SPECIES_PARSE_LIST: &[(&str, Species)] = &[
            $(($common, Species::$identifier),)*
            $(($scientific, Species::$identifier),)*
        ];

        impl std::str::FromStr for Species {
            type Err = BoxedError<'static, BasicKind>;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let s = s.trim().to_lowercase();
                for (name, species) in SPECIES_PARSE_LIST {
                    if name.eq_ignore_ascii_case(&s) {
                        return Ok(*species);
                    }
                }

                let mut options: Vec<(usize, Self)> = Vec::with_capacity(5);
                for (name, species) in SPECIES_PARSE_LIST {
                    let d = mzcv::text::levenshtein_distance(&name.to_lowercase(), &s);
                    let index = options.binary_search_by(|a| a.0.cmp(&d)).unwrap_or_else(|i| i);
                    if index < 5 {
                        if options.len() == 5 {
                            options.remove(4);
                        }
                        options.insert(index, (d, *species));
                    }
                }

                Err(BoxedError::new(BasicKind::Error,
                        "Unknown species name",
                        "The provided name could not be recognised as a species name.",
                        Context::default().lines(0, s.as_str()).to_owned()
                    ).suggestions(options.into_iter().map(|(_, s)| s.to_string()))
                )
            }
        }

        impl std::convert::TryFrom<mzcv::Curie> for Species {
            type Error = ();
            fn try_from(value: mzcv::Curie) -> Result<Self, Self::Error> {
                match value {
                    $(curie!($cv:$curie) => Ok(Self::$identifier),)*
                    _ => Err(())
                }
            }
        }
    };
}

species!(
    AcanthopagrusSchlegelii, "Black porgy", "Acanthopagrus schlegelii (black porgy)", "Acanthopagrus schlegelii",NCBITaxon:72011
    AcipenserBaerii, "Siberian sturgeon", "Acipenser baerii (Siberian sturgeon)", "Acipenser baerii",NCBITaxon:27689
    AcipenserGueldenstaedtii, "Russian sturgeon", "Acipenser gueldenstaedtii (Russian sturgeon)", "Acipenser gueldenstaedtii",NCBITaxon:7902
    AcipenserRuthenus, "Sterlet", "Acipenser ruthenus (sterlet)", "Acipenser ruthenus",NCBITaxon:7906
    AcipenserSchrenckii, "Amur sturgeon", "Acipenser schrenckii (Amur sturgeon)", "Acipenser schrenckii",NCBITaxon:111304
    AcipenserSinensis, "Chinese sturgeon", "Acipenser sinensis (Chinese sturgeon)", "Acipenser sinensis",NCBITaxon:61970
    AiluropodaMelanoleuca, "Giant panda", "Ailuropoda melanoleuca (giant panda)", "Ailuropoda melanoleuca",NCBITaxon:9646
    AlligatorSinensis, "Chinese alligator", "Alligator sinensis (Chinese alligator)", "Alligator sinensis",NCBITaxon:38654
    AmblyrajaGeorgiana, "Antarctic starry skate", "Amblyraja georgiana (Antarctic starry skate)", "Amblyraja georgiana",NCBITaxon:1051904
    AmblyrajaHyperborea, "Arctic skate", "Amblyraja hyperborea (Arctic skate)", "Amblyraja hyperborea",NCBITaxon:315322
    AmbystomaMexicanum, "Axolotl", "Ambystoma mexicanum (axolotl)", "Ambystoma mexicanum",NCBITaxon:8296
    AmeivaAmeiva, "Jungle runners", "Ameiva ameiva", "Ameiva ameiva",NCBITaxon:8535
    AmiaCalva, "Bowfin", "Amia calva (bowfin)", "Amia calva",NCBITaxon:7924
    AmphiprionClarkii, "Yellowtail clownfish", "Amphiprion clarkii (yellowtail clownfish)", "Amphiprion clarkii",NCBITaxon:80970
    AnarhichasMinor, "Spotted wolffish", "Anarhichas minor (spotted wolffish)", "Anarhichas minor",NCBITaxon:65739
    AnasPlatyrhynchos, "Mallard", "Anas platyrhynchos (mallard)", "Anas platyrhynchos",NCBITaxon:8839
    AnguillaAnguilla, "European eel", "Anguilla anguilla (European eel)", "Anguilla anguilla",NCBITaxon:7936
    AnguillaJaponica, "Japanese eel", "Anguilla japonica (Japanese eel)", "Anguilla japonica",NCBITaxon:7937
    AnolisCarolinensis, "Green anole", "Anolis carolinensis (green anole)", "Anolis carolinensis",NCBITaxon:28377
    AnoplopomaFimbria, "Sablefish", "Anoplopoma fimbria (sablefish)", "Anoplopoma fimbria",NCBITaxon:229290
    AnserAnser, "Domestic goose", "Anser anser (Domestic goose)", "Anser anser",NCBITaxon:8843
    AnserAnserDomesticus, "Domestic goose", "Anser anser domesticus", "Anser anser domesticus",NCBITaxon:8843 // TODO: maybe merge
    AnserCaerulescens, "Snow goose", "Anser caerulescens (Snow goose)", "Anser caerulescens",NCBITaxon:8849
    AnserSpGIGHV2011, "Geese GIGHV2011", "Anser sp. GIGHV2011", "Anser sp. GIGHV2011",NCBITaxon:980882
    AnserSpGIGLV2009, "Geese GIGLV2009", "Anser sp. GIGLV2009", "Anser sp. GIGLV2009",NCBITaxon:713105
    AotusAzarai, "Azara's night monkey", "Aotus azarai (Azara's night monkey)", "Aotus azarai",NCBITaxon:30591
    AotusNancymaae, "Ma's night monkey", "Aotus nancymaae (Ma's night monkey)", "Aotus nancymaae",NCBITaxon:37293
    AotusTrivirgatus, "Douroucouli", "Aotus trivirgatus (douroucouli)", "Aotus trivirgatus",NCBITaxon:9505
    ApisCerana, "Asiatic honeybee", "Apis cerana (Asiatic honeybee)", "Apis cerana",NCBITaxon:7461
    ArgyropelecusHemigymnus, "Half-naked hatchetfish", "Argyropelecus hemigymnus (half-naked hatchetfish)", "Argyropelecus hemigymnus",NCBITaxon:68504
    AtelesBelzebuth, "White-bellied spider monkey", "Ateles belzebuth (white-bellied spider monkey)", "Ateles belzebuth",NCBITaxon:9507
    AtelesGeoffroyi, "Black-handed spider monkey", "Ateles geoffroyi (black-handed spider monkey)", "Ateles geoffroyi",NCBITaxon:9509
    AtherinaBoyeri, "Big-scale sand smelt", "Atherina boyeri (big-scale sand smelt)", "Atherina boyeri",NCBITaxon:87785
    BalaenopteraAcutorostrata, "Minke whale", "Balaenoptera acutorostrata (minke whale)", "Balaenoptera acutorostrata",NCBITaxon:9767
    BalaenopteraOmurai, "Omura's baleen whale", "Balaenoptera omurai (Omura's baleen whale)", "Balaenoptera omurai",NCBITaxon:255217
    BathyrajaAlbomaculata, "White-dotted skate", "Bathyraja albomaculata (white-dotted skate)", "Bathyraja albomaculata",NCBITaxon:458555
    BathyrajaBrachyurops, "Broadnose skate", "Bathyraja brachyurops (broadnose skate)", "Bathyraja brachyurops",NCBITaxon:458556
    BathyrajaEatonii, "Eaton's skate", "Bathyraja eatonii (Eaton's skate)", "Bathyraja eatonii",NCBITaxon:298348
    BosGaurus, "Gaur", "Bos gaurus (gaur)", "Bos gaurus",NCBITaxon:9904
    BosIndicus, "Domestic zebu", "Bos indicus (zebu cattle)", "Bos indicus",NCBITaxon:9915
    BosJavanicus, "Banteng", "Bos javanicus (banteng)", "Bos javanicus",NCBITaxon:9906
    BosTaurus, "Domestic bovine", "Bos taurus (bovine)", "Bos taurus",NCBITaxon:9913
    BosTaurusXBosIndicus, "Bos taurus and Bos indicus cross", "Bos taurus x Bos indicus", "Bos taurus x Bos indicus",NCBITaxon:30522
    BosIndicisXBosTaurus, "Bos indicus and Bos taurus cross", "Bos indicus x Bos taurus (hybrid cattle)", "Bos indicus x Bos taurus",NCBITaxon:30522 // TODO: maybe merge
    BovichtusDiacanthus, "Tristan clipfish", "Bovichtus diacanthus", "Bovichtus diacanthus",NCBITaxon:483258
    BubalusBubalis, "Water buffalo", "Bubalus bubalis (water buffalo)", "Bubalus bubalis",NCBITaxon:89462
    BuergeriaBuergeri, "Buerger's frog", "Buergeria buergeri (Buerger's frog)", "Buergeria buergeri",NCBITaxon:191197
    CaimanCrocodilus, "Spectacled caiman", "Caiman crocodilus (spectacled caiman)", "Caiman crocodilus",NCBITaxon:8499
    CairinaMoschata, "Muscovy duck", "Cairina moschata (Muscovy duck)", "Cairina moschata",NCBITaxon:8855
    CallithrixJacchus, "white-tufted-ear marmoset", "Callithrix jacchus (white-tufted-ear marmoset)", "Callithrix jacchus",NCBITaxon:9483
    CallorhinchusMilii, "Elephant shark", "Callorhinchus milii (elephant shark)", "Callorhinchus milii",NCBITaxon:7868
    Camelidae, "Camels", "Camelidae", "Camelidae",NCBITaxon:9835
    CamelusBactrianus, "Bactrian camel", "Camelus bactrianus (Bactrian camel)", "Camelus bactrianus",NCBITaxon:9837
    CamelusDromedarius, "Arabian camel", "Camelus dromedarius (Arabian camel)", "Camelus dromedarius",NCBITaxon:9838
    CanisLupus, "Gray wolf", "Canis lupus (gray wolf)", "Canis lupus",NCBITaxon:9612
    CanisLupusFamiliaris, "Domestic dog", "Canis lupus familiaris (dog)", "Canis lupus familiaris",NCBITaxon:9615
    CanisSp, "Dogs", "Canis sp.", "Canis sp.",NCBITaxon:9616
    CapraHircus, "Domestic goat", "Capra hircus (goat)", "Capra hircus",NCBITaxon:9925
    CarassiusAuratus, "Goldfish", "Carassius auratus (goldfish)", "Carassius auratus",NCBITaxon:7957
    CarassiusLangsdorfii, "Japanese silver crucian carp", "Carassius langsdorfii (Japanese silver crucian carp)", "Carassius langsdorfii",NCBITaxon:138676
    CarcharhinusLeucas, "Bull shark", "Carcharhinus leucas (bull shark)", "Carcharhinus leucas",NCBITaxon:46612
    CarcharhinusPlumbeus, "Sandbar shark", "Carcharhinus plumbeus (sandbar shark)", "Carcharhinus plumbeus",NCBITaxon:7808
    CarlitoSyrichta, "Philippine tarsier", "Carlito syrichta (Philippine tarsier)", "Carlito syrichta",NCBITaxon:1868482
    CarolliaPerspicillata, "Seba's short-tailed bat", "Carollia perspicillata (Seba's short-tailed bat)", "Carollia perspicillata",NCBITaxon:40233
    CaviaPorcellus, "Domestic guinea pig", "Cavia porcellus (domestic guinea pig)", "Cavia porcellus",NCBITaxon:10141
    CephalopachusBancanus, "Horsfield's tarsier", "Cephalopachus bancanus (Horsfield's tarsier)", "Cephalopachus bancanus",NCBITaxon:9477
    CeratotheriumSimum, "White rhinoceros", "Ceratotherium simum (white rhinoceros)", "Ceratotherium simum",NCBITaxon:9807
    CercocebusAtys, "Sooty mangabey", "Cercocebus atys (sooty mangabey)", "Cercocebus atys",NCBITaxon:9531
    CercocebusTorquatus, "Collared mangabey", "Cercocebus torquatus (collared mangabey)", "Cercocebus torquatus",NCBITaxon:9530
    CervusElaphusHispanicus, "Spanish red deer", "Cervus elaphus hispanicus", "Cervus elaphus hispanicus",NCBITaxon:141655
    ChaenocephalusAceratus, "Blackfin icefish", "Chaenocephalus aceratus (blackfin icefish)", "Chaenocephalus aceratus",NCBITaxon:36190
    ChampsocephalusEsox, "Pike icefish", "Champsocephalus esox (pike icefish)", "Champsocephalus esox",NCBITaxon:159716
    ChannaArgus, "Northern snakehead", "Channa argus (northern snakehead)", "Channa argus",NCBITaxon:215402
    ChannaStriata, "Snakehead murrel", "Channa striata (snakehead murrel)", "Channa striata",NCBITaxon:64152
    ChaunaTorquata, "Southern screamer", "Chauna torquata (southern screamer)", "Chauna torquata",NCBITaxon:30388
    ChelonAuratus, "Golden grey mullet", "Chelon auratus (golden grey mullet)", "Chelon auratus",NCBITaxon:48191
    ChionodracoHamatus, "Antarctic icefish", "Chionodraco hamatus (Antarctic icefish)", "Chionodraco hamatus",NCBITaxon:36188
    ChionodracoRastrospinosus, "Ocellated icefish", "Chionodraco rastrospinosus (ocellated icefish)", "Chionodraco rastrospinosus",NCBITaxon:34790
    ChlorocebusAethiops, "Grivet", "Chlorocebus aethiops (grivet)", "Chlorocebus aethiops",NCBITaxon:9534
    ClupeaPallasii, "Pacific herring", "Clupea pallasii (Pacific herring)", "Clupea pallasii",NCBITaxon:30724
    ColobusGuereza, "Mantled guereza", "Colobus guereza (mantled guereza)", "Colobus guereza",NCBITaxon:33548
    ColobusPolykomos, "King colobus", "Colobus polykomos (king colobus)", "Colobus polykomos",NCBITaxon:9572
    CricetinaeSp, "Hamster", "Cricetinae gen. sp. (Hamster)", "Cricetinae sp.",NCBITaxon:36483
    CricetulusMigratorius, "Armenian hamster", "Cricetulus migratorius (Armenian hamster)", "Cricetulus migratorius",NCBITaxon:3122392
    CrocodylusSiamensis, "Siamese crocodile", "Crocodylus siamensis (Siamese crocodile)", "Crocodylus siamensis",NCBITaxon:68455
    CtenopharyngodonIdella, "Grass carp", "Ctenopharyngodon idella (grass carp)", "Ctenopharyngodon idella",NCBITaxon:7959
    CygnodracoMawsoni, "Mawson's dragonfish", "Cygnodraco mawsoni (Mawson's dragonfish)", "Cygnodraco mawsoni",NCBITaxon:8216
    CynoglossusSemilaevis, "Tongue sole", "Cynoglossus semilaevis (tongue sole)", "Cynoglossus semilaevis",NCBITaxon:244447
    CynopterusSphinx, "Indian short-nosed fruit bat", "Cynopterus sphinx (Indian short-nosed fruit bat)", "Cynopterus sphinx",NCBITaxon:9400
    CyprinusCarpio, "Common carp", "Cyprinus carpio (common carp)", "Cyprinus carpio",NCBITaxon:7962
    DanioRerio, "Zebrafish", "Danio rerio (zebrafish)", "Danio rerio",NCBITaxon:7955
    DaubentoniaMadagascariensis, "Aye-aye", "Daubentonia madagascariensis (aye-aye)", "Daubentonia madagascariensis",NCBITaxon:31869
    DelphinapterusLeucas, "Beluga whale", "Delphinapterus leucas (beluga whale)", "Delphinapterus leucas",NCBITaxon:9749
    DelphinusCapensis, "Long-beaked common dolphin", "Delphinus capensis (long-beaked common dolphin)", "Delphinus capensis",NCBITaxon:103584
    DicentrarchusLabrax, "European seabass", "Dicentrarchus labrax (European seabass)", "Dicentrarchus labrax",NCBITaxon:13489
    DissostichusMawsoni, "Antarctic toothfish", "Dissostichus mawsoni (Antarctic toothfish)", "Dissostichus mawsoni",NCBITaxon:36200
    DrosophilaMelanogaster, "Fruit fly", "Drosophila melanogaster (fruit fly)", "Drosophila melanogaster",NCBITaxon:7227
    ElapheTaeniura, "Beauty snake", "Elaphe taeniura (beauty snake)", "Elaphe taeniura",NCBITaxon:74398
    EleginopsMaclovinus, "Patagonian blennie", "Eleginops maclovinus (Patagonian blennie)", "Eleginops maclovinus",NCBITaxon:56733
    ElopsAaurus, "Ladyfish", "Elops saurus (ladyfish)", "Elops saurus",NCBITaxon:7928
    EpinephelusAkaara, "Hong Kong grouper", "Epinephelus akaara (Hong Kong grouper)", "Epinephelus akaara",NCBITaxon:215347
    EpinephelusCoioides, "Orange-spotted grouper", "Epinephelus coioides (orange-spotted grouper)", "Epinephelus coioides",NCBITaxon:94232
    EptatretusBurgeri, "Inshore hagfish", "Eptatretus burgeri (inshore hagfish)", "Eptatretus burgeri",NCBITaxon:7764
    EptesicusFuscus, "Big brown bat", "Eptesicus fuscus (big brown bat)", "Eptesicus fuscus",NCBITaxon:29078
    EquusAsinus, "Ass", "Equus asinus (ass)", "Equus asinus",NCBITaxon:9793
    EquusBurchelliiAntiquorum, "Burchell's zebra", "Equus burchellii antiquorum", "Equus burchellii antiquorum",NCBITaxon:89252
    EquusCaballus, "Domestic horse", "Equus caballus (horse)", "Equus caballus",NCBITaxon:9796
    EquusQuaggaBurchellii, "Burchell's zebra", "Equus quagga burchellii (Burchell's zebra)", "Equus quagga burchellii",NCBITaxon:89252
    ErythrocebusPatas, "Red guenon", "Erythrocebus patas (red guenon)", "Erythrocebus patas",NCBITaxon:9538
    EscherichiaColi, "E. coli", "Escherichia coli (E. coli)", "Escherichia coli",NCBITaxon:562
    EsoxLucius, "Northern pike", "Esox lucius (northern pike)", "Esox lucius",NCBITaxon:8010
    EublepharisMacularius, "Leopard gecko", "Eublepharis macularius (Leopard gecko)", "Eublepharis macularius",NCBITaxon:481883
    EulemurFulvus, "Brown lemur", "Eulemur fulvus (brown lemur)", "Eulemur fulvus",NCBITaxon:13515
    FelineLeukemiaVirus, "Feline leukemia virus", "Feline leukemia virus", "Feline leukemia virus",NCBITaxon:11768
    FelisCatus, "Domestic cat", "Felis catus (domestic cat)", "Felis catus",NCBITaxon:9685
    FelisSp, "Cats", "Felis sp.", "Felis sp.",NCBITaxon:9687
    GadusMorhua, "Atlantic cod", "Gadus morhua (Atlantic cod)", "Gadus morhua",NCBITaxon:8049
    GalagoSenegalensis, "Senegal galago", "Galago senegalensis (Senegal galago)", "Galago senegalensis",NCBITaxon:9465
    GallusGallus, "Domestic chicken", "Gallus gallus (chicken)", "Gallus gallus",NCBITaxon:9031
    GasterosteusAculeatus, "Three-spined stickleback", "Gasterosteus aculeatus (three-spined stickleback)", "Gasterosteus aculeatus",NCBITaxon:69293
    GinglymostomaCirratum, "Nurse shark", "Ginglymostoma cirratum (nurse shark)", "Ginglymostoma cirratum",NCBITaxon:7801
    GobionotothenGibberifrons, "Humped rockcod", "Gobionotothen gibberifrons (humped rockcod)", "Gobionotothen gibberifrons",NCBITaxon:36202
    GorillaGorilla, "Western gorilla", "Gorilla gorilla (western gorilla)", "Gorilla gorilla",NCBITaxon:9593
    GorillaGorillaGorilla, "Western lowland gorilla", "Gorilla gorilla gorilla (western lowland gorilla)", "Gorilla gorilla gorilla",NCBITaxon:9595
    GrampusGriseus, "Risso's dolphin", "Grampus griseus (Risso's dolphin)", "Grampus griseus",NCBITaxon:83653
    GymnodracoAcuticeps, "Ploughfish", "Gymnodraco acuticeps", "Gymnodraco acuticeps",NCBITaxon:8218
    GymnogypsCalifornianus, "California condor", "Gymnogyps californianus (California condor)", "Gymnogyps californianus",NCBITaxon:33616
    HaemorhousMexicanus, "House finch", "Haemorhous mexicanus (house finch)", "Haemorhous mexicanus",NCBITaxon:30427
    HemibagrusMacropterus, "Largefin longbarbel catfish", "Hemibagrus macropterus", "Hemibagrus macropterus",NCBITaxon:175789
    HepacivirusC, "Hepacivirus C", "Hepacivirus C", "Hepacivirus C",NCBITaxon:3052230
    HeterocephalusGlaber, "Naked mole-rat", "Heterocephalus glaber (naked mole-rat)", "Heterocephalus glaber",NCBITaxon:10181
    HeterodontusFrancisci, "Horn shark", "Heterodontus francisci (horn shark)", "Heterodontus francisci",NCBITaxon:7792
    HippoglossusHippoglossus, "Atlantic halibut", "Hippoglossus hippoglossus (Atlantic halibut)", "Hippoglossus hippoglossus",NCBITaxon:8267
    HistiodracoVelifer, "Histiodraco velifer", "Histiodraco velifer", "Histiodraco velifer",NCBITaxon:36208
    HomoSapiens, "Human", "Homo sapiens (human)", "Homo sapiens",NCBITaxon:9606
    HoolockHoolock, "Hoolock gibbon", "Hoolock hoolock (hoolock gibbon)", "Hoolock hoolock",NCBITaxon:61851
    HusoHuso, "Beluga", "Huso huso (beluga)", "Huso huso",NCBITaxon:61971
    HydrolagusColliei, "Spotted ratfish", "Hydrolagus colliei (spotted ratfish)", "Hydrolagus colliei",NCBITaxon:7873
    HylobatesLar, "Common gibbon", "Hylobates lar (common gibbon)", "Hylobates lar",NCBITaxon:9580
    IctalurusPunctatus, "Channel catfish", "Ictalurus punctatus (channel catfish)", "Ictalurus punctatus",NCBITaxon:7998
    IsoodonMacrourus, "Northern brown bandicoot", "Isoodon macrourus (northern brown bandicoot)", "Isoodon macrourus",NCBITaxon:37698
    KogiaSima, "Dwarf sperm whale", "Kogia sima (dwarf sperm whale)", "Kogia sima",NCBITaxon:9752
    LabeobarbusIntermedius, "Labeobarbus intermedius", "Labeobarbus intermedius", "Labeobarbus intermedius",NCBITaxon:40831
    LabeoRohita, "Rohu", "Labeo rohita (rohu)", "Labeo rohita",NCBITaxon:84645
    LamaGlama, "Llama", "Lama glama (llama)", "Lama glama",NCBITaxon:9844
    LarimichthysCrocea, "Large yellow croaker", "Larimichthys crocea (large yellow croaker)", "Larimichthys crocea",NCBITaxon:215358
    LatimeriaChalumnae, "Coelacanth", "Latimeria chalumnae (coelacanth)", "Latimeria chalumnae",NCBITaxon:7897
    LatimeriaMenadoensis, "Menado coelacanth", "Latimeria menadoensis (Menado coelacanth)", "Latimeria menadoensis",NCBITaxon:106881
    LatrisLineata, "Striped trumpeter", "Latris lineata (striped trumpeter)", "Latris lineata",NCBITaxon:97162
    LemurCatta, "Ring-tailed lemur", "Lemur catta (Ring-tailed lemur)", "Lemur catta",NCBITaxon:9447
    LeontopithecusRosalia, "Golden lion tamarin", "Leontopithecus rosalia (golden lion tamarin)", "Leontopithecus rosalia",NCBITaxon:30588
    LepilemurRuficaudatus, "Red-tailed sportive lemur", "Lepilemur ruficaudatus (red-tailed sportive lemur)", "Lepilemur ruficaudatus",NCBITaxon:78866
    LepisosteusOsseus, "Longnose gar", "Lepisosteus osseus (longnose gar)", "Lepisosteus osseus",NCBITaxon:34771
    LepusAmericanus, "Snowshoe hare", "Lepus americanus (snowshoe hare)", "Lepus americanus",NCBITaxon:48086
    LepusCalifornicus, "Black-tailed jackrabbit", "Lepus californicus (black-tailed jackrabbit)", "Lepus californicus",NCBITaxon:48087
    LepusCallotis, "White-sided jackrabbit", "Lepus callotis (white-sided jackrabbit)", "Lepus callotis",NCBITaxon:62619
    LepusCapensis, "Brown hare", "Lepus capensis (brown hare)", "Lepus capensis",NCBITaxon:9981
    LepusCastroviejoi, "Broom Hare", "Lepus castroviejoi (Broom Hare)", "Lepus castroviejoi",NCBITaxon:187398
    LepusEuropaeus, "European hare", "Lepus europaeus (European hare)", "Lepus europaeus",NCBITaxon:9983
    LepusGranatensis, "Granada hare", "Lepus granatensis (Granada hare)", "Lepus granatensis",NCBITaxon:100182
    LepusSaxatilis, "Scrub hare", "Lepus saxatilis (scrub hare)", "Lepus saxatilis",NCBITaxon:63224
    LepusTimidus, "Mountain hare", "Lepus timidus (Mountain hare)", "Lepus timidus",NCBITaxon:62621
    LeucorajaErinacea, "Little skate", "Leucoraja erinacea (little skate)", "Leucoraja erinacea",NCBITaxon:7782
    LipotesVexillifer, "Yangtze River dolphin", "Lipotes vexillifer (Yangtze River dolphin)", "Lipotes vexillifer",NCBITaxon:118797
    LutjanusSanguineus, "Humphead snapper", "Lutjanus sanguineus (humphead snapper)", "Lutjanus sanguineus",NCBITaxon:264213
    MacacaArctoides, "Stump-tailed macaque", "Macaca arctoides (stump-tailed macaque)", "Macaca arctoides",NCBITaxon:9540
    MacacaAssamensis, "Assam macaque", "Macaca assamensis (Assam macaque)", "Macaca assamensis",NCBITaxon:9551
    MacacaCyclopis, "Taiwan macaque", "Macaca cyclopis (Taiwan macaque)", "Macaca cyclopis",NCBITaxon:78449
    MacacaFascicularis, "Crab-eating macaque", "Macaca fascicularis (crab-eating macaque)", "Macaca fascicularis",NCBITaxon:9541
    MacacaMulatta, "Rhesus monkey", "Macaca mulatta (Rhesus monkey)", "Macaca mulatta",NCBITaxon:9544
    MacacaNemestrina, "Pig-tailed macaque", "Macaca nemestrina (pig-tailed macaque)", "Macaca nemestrina",NCBITaxon:9545
    MacacaSilenus, "Liontail macaque", "Macaca silenus (liontail macaque)", "Macaca silenus",NCBITaxon:54601
    MacacaThibetana, "Pere David's macaque", "Macaca thibetana (Pere David's macaque)", "Macaca thibetana",NCBITaxon:54602
    MarecaStrepera, "Gadwall", "Mareca strepera (gadwall)", "Mareca strepera",NCBITaxon:75861
    MarmotaHimalayana, "Himalayan marmot", "Marmota himalayana (Himalayan marmot)", "Marmota himalayana",NCBITaxon:93163
    MarmotaMonax, "Woodchuck", "Marmota monax (woodchuck)", "Marmota monax",NCBITaxon:9995
    MauremysMutica, "Yellowpond turtle", "Mauremys mutica (yellowpond turtle)", "Mauremys mutica",NCBITaxon:74926
    MelanogrammusAeglefinus, "Haddock", "Melanogrammus aeglefinus (haddock)", "Melanogrammus aeglefinus",NCBITaxon:8056
    MeleagrisGallopavo, "Turkey", "Meleagris gallopavo (turkey)", "Meleagris gallopavo",NCBITaxon:9103
    MerionesUnguiculatus, "Mongolian gerbil", "Meriones unguiculatus (Mongolian gerbil)", "Meriones unguiculatus",NCBITaxon:10047
    MesocricetusAuratus, "Golden hamster", "Mesocricetus auratus (golden hamster)", "Mesocricetus auratus",NCBITaxon:10036
    MicrocebusMurinus, "Gray mouse lemur", "Microcebus murinus (gray mouse lemur)", "Microcebus murinus",NCBITaxon:30608
    MonodelphisDomestica, "Gray short-tailed opossum", "Monodelphis domestica (gray short-tailed opossum)", "Monodelphis domestica",NCBITaxon:13616
    MurinaeSp, "Old world rats and mice", "Murinae gen. sp.", "Murinae sp.",NCBITaxon:3148847
    Mus, "mouse", "Mus (mouse)", "Mus",NCBITaxon:10088
    MusCookii, "Cook's mouse", "Mus cookii (Cook's mouse)", "Mus cookii",NCBITaxon:10098
    MusMinutoides, "Southern African pygmy mouse", "Mus minutoides (Southern African pygmy mouse)", "Mus minutoides",NCBITaxon:10105
    MusMusculus, "House mouse", "Mus musculus (house mouse)", "Mus musculus",NCBITaxon:10090
    MusMusculusCastaneus, "Southeastern Asian house mouse", "Mus musculus castaneus (southeastern Asian house mouse)", "Mus musculus castaneus",NCBITaxon:10091
    MusMusculusDomesticus, "Western European house mouse", "Mus musculus domesticus (western European house mouse)", "Mus musculus domesticus",NCBITaxon:10092
    MusMusculusMolossinus, "Japanese wild mouse", "Mus musculus molossinus (Japanese wild mouse)", "Mus musculus molossinus",NCBITaxon:57486
    MusMusculusMusculus, "Eastern European house mouse", "Mus musculus musculus (eastern European house mouse)", "Mus musculus musculus",NCBITaxon:39442
    MusPahari, "Shrew mouse", "Mus pahari (shrew mouse)", "Mus pahari",NCBITaxon:10093
    MusSaxicola, "Spiny mouse", "Mus saxicola (spiny mouse)", "Mus saxicola",NCBITaxon:10094
    MusSp, "Mice", "Mus sp. (mice)", "Mus sp.",NCBITaxon:10095
    MusSpretus, "Western wild mouse", "Mus spretus (western wild mouse)", "Mus spretus",NCBITaxon:10096
    MustelaPutoriusFuro, "Domestic ferret", "Mustela putorius furo (domestic ferret)", "Mustela putorius furo",NCBITaxon:9669
    MustelaSp, "Ferret", "Mustela sp.", "Mustela sp.",NCBITaxon:30549
    MyotisLucifugus, "Little brown bat", "Myotis lucifugus (little brown bat)", "Myotis lucifugus",NCBITaxon:59463
    NeophocaenaPhocaenoides, "Indo-Pacific finless porpoise", "Neophocaena phocaenoides (Indo-Pacific finless porpoise)", "Neophocaena phocaenoides",NCBITaxon:34892
    NeogaleVison, "American mink", "Neogale vison (American mink)", "Neogale vison",NCBITaxon:452646
    NomascusConcolor, "Black crested gibbon", "Nomascus concolor (Black crested gibbon)", "Nomascus concolor",NCBITaxon:29089
    NotamacropusEugenii, "Tammar wallaby", "Notamacropus eugenii (tammar wallaby)", "Notamacropus eugenii",NCBITaxon:9315
    NothocricetulusMigratorius, "Armenian hamster", "Nothocricetulus migratorius (Armenian hamster)", "Nothocricetulus migratorius",NCBITaxon:3122392
    NototheniaCoriiceps, "Black rockcod", "Notothenia coriiceps (black rockcod)", "Notothenia coriiceps",NCBITaxon:8208
    NycticebusCoucang, "Slow loris", "Nycticebus coucang (slow loris)", "Nycticebus coucang",NCBITaxon:9470
    OncorhynchusGorbuscha, "Pink salmon", "Oncorhynchus gorbuscha (pink salmon)", "Oncorhynchus gorbuscha",NCBITaxon:8017
    OncorhynchusMykiss, "Rainbow trout", "Oncorhynchus mykiss (rainbow trout)", "Oncorhynchus mykiss",NCBITaxon:8022
    OncorhynchusTshawytscha, "Chinook salmon", "Oncorhynchus tshawytscha (Chinook salmon)", "Oncorhynchus tshawytscha",NCBITaxon:74940
    OrectolobusMaculatus, "Spotted wobbegong", "Orectolobus maculatus (spotted wobbegong)", "Orectolobus maculatus",NCBITaxon:168098
    OreochromisNiloticus, "Nile tilapia", "Oreochromis niloticus (Nile tilapia)", "Oreochromis niloticus",NCBITaxon:8128
    OrnithorhynchusAnatinus, "Platypus", "Ornithorhynchus anatinus (platypus)", "Ornithorhynchus anatinus",NCBITaxon:9258
    OryctolagusCuniculus, "Rabbit", "Oryctolagus cuniculus (rabbit)", "Oryctolagus cuniculus",NCBITaxon:9986
    OryctolagusCuniculusAlgirus, "European rabbit", "Oryctolagus cuniculus algirus", "Oryctolagus cuniculus algirus",NCBITaxon:230741
    OryctolagusCuniculusCuniculus, "Rabbit", "Oryctolagus cuniculus cuniculus", "Oryctolagus cuniculus cuniculus",NCBITaxon:568996
    OryziasLatipes, "Japanese medaka", "Oryzias latipes (Japanese medaka)", "Oryzias latipes",NCBITaxon:8090
    OryziasMelastigma, "Indian medaka", "Oryzias melastigma (Indian medaka)", "Oryzias melastigma",NCBITaxon:30732
    OtolemurCrassicaudatus, "Thick-tailed bush baby", "Otolemur crassicaudatus (thick-tailed bush baby)", "Otolemur crassicaudatus",NCBITaxon:9463
    OvisAries, "Domestic sheep", "Ovis aries (sheep)", "Ovis aries",NCBITaxon:9940
    OvisSp, "Sheep", "Ovis sp.", "Ovis sp.",NCBITaxon:9939
    PacifastacusLeniusculus, "Signal crayfish", "Pacifastacus leniusculus (signal crayfish)", "Pacifastacus leniusculus",NCBITaxon:6720
    PagetopsisMacropterus, "Pagetopsis macropterus", "Pagetopsis macropterus", "Pagetopsis macropterus",NCBITaxon:36194
    PagrusMajor, "Red seabream", "Pagrus major (red seabream)", "Pagrus major",NCBITaxon:143350
    PangasianodonHypophthalmus, "Striped catfish", "Pangasianodon hypophthalmus (striped catfish)", "Pangasianodon hypophthalmus",NCBITaxon:310915
    PanPaniscus, "Pygmy chimpanzee", "Pan paniscus (pygmy chimpanzee)", "Pan paniscus",NCBITaxon:9597
    PantheraPardus, "Leopard", "Panthera pardus (leopard)", "Panthera pardus",NCBITaxon:9691
    PanTroglodytes, "Chimpanzee", "Pan troglodytes (chimpanzee)", "Pan troglodytes",NCBITaxon:9598
    PanTroglodytesVerus, "Western chimpanzee", "Pan troglodytes verus", "Pan troglodytes verus",NCBITaxon:37012
    PapioAnubis, "Olive baboon", "Papio anubis (olive baboon)", "Papio anubis",NCBITaxon:9555
    PapioAnubisAnubis, "Olive baboon anubis", "Papio anubis anubis", "Papio anubis anubis",NCBITaxon:211508
    PapioHamadryas, "Hamadryas baboon", "Papio hamadryas (hamadryas baboon)", "Papio hamadryas",NCBITaxon:9557
    PapioPapio, "Guinea baboon", "Papio papio (Guinea baboon)", "Papio papio",NCBITaxon:100937
    ParalichthysOlivaceus, "Japanese flounder", "Paralichthys olivaceus (Japanese flounder)", "Paralichthys olivaceus",NCBITaxon:8255
    PelodiscusSinensis, "Chinese soft-shelled turtle", "Pelodiscus sinensis (Chinese soft-shelled turtle)", "Pelodiscus sinensis",NCBITaxon:13735
    PelteobagrusFulvidraco, "Yellow catfish", "Pelteobagrus fulvidraco (yellow catfish)", "Pelteobagrus fulvidraco",NCBITaxon:1234273
    PerdixPerdix, "Grey partridge", "Perdix perdix (grey partridge)", "Perdix perdix",NCBITaxon:9052
    PeromyscusManiculatus, "North American deer mouse", "Peromyscus maniculatus (North American deer mouse)", "Peromyscus maniculatus",NCBITaxon:10042
    PetromyzonMarinus, "Sea lamprey", "Petromyzon marinus (sea lamprey)", "Petromyzon marinus",NCBITaxon:7757
    PhascogaleCalura, "Red-tailed phascogale", "Phascogale calura (red-tailed phascogale)", "Phascogale calura",NCBITaxon:38776
    PhasianusColchicus, "Ring-necked pheasant", "Phasianus colchicus (Ring-necked pheasant)", "Phasianus colchicus",NCBITaxon:9054
    PhyseterCatodon, "Sperm whale", "Physeter catodon (sperm whale)", "Physeter catodon",NCBITaxon:9755
    PitheciaPithecia, "White-faced saki", "Pithecia pithecia (white-faced saki)", "Pithecia pithecia",NCBITaxon:43777
    PlataleaAjaja, "Roseate spoonbil", "Platalea ajaja", "Platalea ajaja",NCBITaxon:371920
    Platyrrhini, "New World monkeys", "Platyrrhini (New World monkeys)", "Platyrrhini",NCBITaxon:9479
    PlecoglossusAltivelisAltivelis, "Ayu sweetfish", "Plecoglossus altivelis altivelis", "Plecoglossus altivelis altivelis",NCBITaxon:281464
    PleurodelesWaltl, "Iberian ribbed newt", "Pleurodeles waltl (Iberian ribbed newt)", "Pleurodeles waltl",NCBITaxon:8319
    PogonophryneScotti, "Pogonophryne scotti", "Pogonophryne scotti", "Pogonophryne scotti",NCBITaxon:36210
    PolyprionOxygeneios, "Hāpuku", "Polyprion oxygeneios", "Polyprion oxygeneios",NCBITaxon:334912
    PongoAbelii, "Sumatran orangutan", "Pongo abelii (Sumatran orangutan)", "Pongo abelii",NCBITaxon:9601
    PongoPygmaeus, "Bornean orangutan", "Pongo pygmaeus (Bornean orangutan)", "Pongo pygmaeus",NCBITaxon:9600
    PresbytisComata, "Grizzled leaf monkey", "Presbytis comata (grizzled leaf monkey)", "Presbytis comata",NCBITaxon:78452
    PresbytisFemoralis, "Banded leaf monkey", "Presbytis femoralis (banded leaf monkey)", "Presbytis femoralis",NCBITaxon:78453
    PresbytisMelalophos, "Mitred leaf monkey", "Presbytis melalophos (mitred leaf monkey)", "Presbytis melalophos",NCBITaxon:78451
    PropithecusVerreauxi, "White sifaka", "Propithecus verreauxi (white sifaka)", "Propithecus verreauxi",NCBITaxon:34825
    ProtopterusAethiopicus, "Marbled lungfish", "Protopterus aethiopicus (marbled lungfish)", "Protopterus aethiopicus",NCBITaxon:7886
    PseudobatosProductus, "Shovelnose guitarfish", "Pseudobatos productus (shovelnose guitarfish)", "Pseudobatos productus",NCBITaxon:170824
    PteropusAlecto, "Black flying fox", "Pteropus alecto (black flying fox)", "Pteropus alecto",NCBITaxon:9402
    PythonBivittatus, "Burmese python", "Python bivittatus (Burmese python)", "Python bivittatus",NCBITaxon:176946
    RachycentronCanadum, "Cobia", "Rachycentron canadum (cobia)", "Rachycentron canadum",NCBITaxon:141264
    RajaEglanteria, "Clearnose skate", "Raja eglanteria (clearnose skate)", "Raja eglanteria",NCBITaxon:3360502
    RattusFuscipes, "Bush rat", "Rattus fuscipes (bush rat)", "Rattus fuscipes",NCBITaxon:10119
    RattusLeucopus, "Mottle-tailed rat", "Rattus leucopus (mottle-tailed rat)", "Rattus leucopus",NCBITaxon:10115
    RattusNorvegicus, "Norway rat", "Rattus norvegicus (Norway rat)", "Rattus norvegicus",NCBITaxon:10116
    RattusRattus, "Black rat", "Rattus rattus (black rat)", "Rattus rattus",NCBITaxon:10117
    RattusSordidus, "Australian dusky field rat", "Rattus sordidus (Australian dusky field rat)", "Rattus sordidus",NCBITaxon:10120
    RattusSp, "Rats", "Rattus sp. (rats)", "Rattus sp.",NCBITaxon:10118
    RattusTunneyi, "Tunney's rat", "Rattus tunneyi (Tunney's rat)", "Rattus tunneyi",NCBITaxon:10121
    RattusVillosissimus, "Long-haired rat", "Rattus villosissimus (long-haired rat)", "Rattus villosissimus",NCBITaxon:10122
    RhinocerosUnicornis, "Greater Indian rhinoceros", "Rhinoceros unicornis (greater Indian rhinoceros)", "Rhinoceros unicornis",NCBITaxon:9809
    RousettusLeschenaultii, "Leschenault's rousette", "Rousettus leschenaultii (Leschenault's rousette)", "Rousettus leschenaultii",NCBITaxon:9408
    SaccharomycesCerevisiae, "baker's yeast", "Saccharomyces cerevisiae (baker's yeast)", "Saccharomyces cerevisiae",NCBITaxon:4932
    SaguinusLabiatus, "Red-chested mustached tamarin", "Saguinus labiatus (red-chested mustached tamarin)", "Saguinus labiatus",NCBITaxon:78454
    SaguinusMidas, "Midas tamarin", "Saguinus midas (Midas tamarin)", "Saguinus midas",NCBITaxon:30586
    SaguinusOedipus, "Cotton-top tamarin", "Saguinus oedipus (cotton-top tamarin)", "Saguinus oedipus",NCBITaxon:9490
    SaimiriBoliviensisBoliviensis, "Bolivian squirrel monkey", "Saimiri boliviensis boliviensis (Bolivian squirrel monkey)", "Saimiri boliviensis boliviensis",NCBITaxon:39432
    SaimiriSciureus, "Common squirrel monkey", "Saimiri sciureus (common squirrel monkey)", "Saimiri sciureus",NCBITaxon:9521
    SalmoMarmoratus, "Salmo marmoratus", "Salmo marmoratus", "Salmo marmoratus",NCBITaxon:33518
    SalmoSalar, "Atlantic salmon", "Salmo salar (Atlantic salmon)", "Salmo salar",NCBITaxon:8030
    SalmoTrutta, "River trout", "Salmo trutta (river trout)", "Salmo trutta",NCBITaxon:8032
    SalvelinusAlpinus, "Arctic char", "Salvelinus alpinus (Arctic char)", "Salvelinus alpinus",NCBITaxon:8036
    SanderVitreus, "Walleye", "Sander vitreus (walleye)", "Sander vitreus",NCBITaxon:283036
    SapajusApella, "Tufted capuchin", "Sapajus apella (Tufted capuchin)", "Sapajus apella",NCBITaxon:9515
    SchizophyllumCommune, "Schizophyllum commune", "Schizophyllum commune", "Schizophyllum commune",NCBITaxon:5334
    SciaenopsOcellatus, "Red drum", "Sciaenops ocellatus (red drum)", "Sciaenops ocellatus",NCBITaxon:76340
    ScophthalmusMaximus, "Turbot", "Scophthalmus maximus (turbot)", "Scophthalmus maximus",NCBITaxon:52904
    ScyliorhinusCanicula, "Smaller spotted catshark", "Scyliorhinus canicula (smaller spotted catshark)", "Scyliorhinus canicula",NCBITaxon:7830
    SeriolaQuinqueradiata, "Japanese amberjack", "Seriola quinqueradiata (Japanese amberjack)", "Seriola quinqueradiata",NCBITaxon:8161
    SilurusAsotus, "Amur catfish", "Silurus asotus (Amur catfish)", "Silurus asotus",NCBITaxon:30991
    SilurusMeridionalis, "Silurus meridionalis", "Silurus meridionalis", "Silurus meridionalis",NCBITaxon:175797
    SinipercaChuatsi, "Mandarin fish", "Siniperca chuatsi (mandarin fish)", "Siniperca chuatsi",NCBITaxon:119488
    SousaChinensis, "Indo-pacific humpbacked dolphin", "Sousa chinensis (Indo-pacific humpbacked dolphin)", "Sousa chinensis",NCBITaxon:103600
    SparusAurata, "Gilthead seabream", "Sparus aurata (gilthead seabream)", "Sparus aurata",NCBITaxon:8175
    SphoeroidesNephelus, "Southern puffer", "Sphoeroides nephelus (southern puffer)", "Sphoeroides nephelus",NCBITaxon:39110
    SqualusAcanthias, "Spiny dogfish", "Squalus acanthias (spiny dogfish)", "Squalus acanthias",NCBITaxon:7797
    StegastesLeucostictus, "Beaugregory", "Stegastes leucostictus (beaugregory)", "Stegastes leucostictus",NCBITaxon:258452
    StegastesPartitus, "Bicolor damselfish", "Stegastes partitus (bicolor damselfish)", "Stegastes partitus",NCBITaxon:144197
    StenellaAttenuata, "Bridled dolphin", "Stenella attenuata (bridled dolphin)", "Stenella attenuata",NCBITaxon:9735
    StenellaCoeruleoalba, "Striped dolphin", "Stenella coeruleoalba (striped dolphin)", "Stenella coeruleoalba",NCBITaxon:9737
    StreptomycesAvidinii, "Streptomyces avidinii", "Streptomyces avidinii", "Streptomyces avidinii",NCBITaxon:1895
    StreptomycesViridochromogenes, "Streptomyces viridochromogenes", "Streptomyces viridochromogenes", "Streptomyces viridochromogenes",NCBITaxon:1938
    StruthioCamelus, "African ostrich", "Struthio camelus (African ostrich)", "Struthio camelus",NCBITaxon:8801
    SuncusMurinus, "House shrew", "Suncus murinus (house shrew)", "Suncus murinus",NCBITaxon:9378
    SusScrofa, "Domestic pig", "Sus scrofa (pig)", "Sus scrofa",NCBITaxon:9823
    SylvilagusCunicularis, "Mexican cottontail", "Sylvilagus cunicularis (Mexican cottontail)", "Sylvilagus cunicularis",NCBITaxon:3370854
    SylvilagusFloridanus, "Eastern cottontail", "Sylvilagus floridanus (eastern cottontail)", "Sylvilagus floridanus",NCBITaxon:9988
    SymphalangusSyndactylus, "Siamang", "Symphalangus syndactylus (siamang)", "Symphalangus syndactylus",NCBITaxon:9590
    TachyglossusAculeatus, "Australian echidna", "Tachyglossus aculeatus (Australian echidna)", "Tachyglossus aculeatus",NCBITaxon:9261
    TachysurusFulvidraco, "Yellow catfish", "Tachysurus fulvidraco (yellow catfish)", "Tachysurus fulvidraco",NCBITaxon:1234273
    TachysurusVachellii, "Tachysurus vachellii", "Tachysurus vachellii", "Tachysurus vachellii",NCBITaxon:175792
    TaeniopygiaGuttata, "Zebra finch", "Taeniopygia guttata (zebra finch)", "Taeniopygia guttata",NCBITaxon:59729
    TakifuguRubripes, "Torafugu", "Takifugu rubripes (torafugu)", "Takifugu rubripes",NCBITaxon:31033
    TarsiusDentatus, "Dian's tarsier", "Tarsius dentatus (Dian's tarsier)", "Tarsius dentatus",NCBITaxon:449501
    TarsiusLariang, "Lariang tarsier", "Tarsius lariang (Lariang tarsier)", "Tarsius lariang",NCBITaxon:630277
    TarsiusSyrichta, "Philippine tarsier", "Tarsius syrichta (Philippine tarsier)", "Tarsius syrichta",NCBITaxon:1868482
    TetraodonNigroviridis, "Spotted green pufferfish", "Tetraodon nigroviridis (spotted green pufferfish)", "Tetraodon nigroviridis",NCBITaxon:99883
    TrachemysScripta, "Red-eared slider turtle", "Trachemys scripta (red-eared slider turtle)", "Trachemys scripta",NCBITaxon:34903
    TrachemysScriptaElegans, "Red-eared slider turtle elegans", "Trachemys scripta elegans", "Trachemys scripta elegans",NCBITaxon:31138
    TrachypithecusCristatus, "Silvery lutung", "Trachypithecus cristatus (Silvery lutung)", "Trachypithecus cristatus",NCBITaxon:122765
    TrachypithecusObscurus, "Dusky leaf-monkey", "Trachypithecus obscurus (Dusky leaf-monkey)", "Trachypithecus obscurus",NCBITaxon:54181
    TrematomusBernacchii, "Emerald rockcod", "Trematomus bernacchii (emerald rockcod)", "Trematomus bernacchii",NCBITaxon:40690
    TrematomusHansoni, "Striped rockcod", "Trematomus hansoni (striped rockcod)", "Trematomus hansoni",NCBITaxon:36196
    TrematomusLoennbergii, "Deepwater notothen", "Trematomus loennbergii (deepwater notothen)", "Trematomus loennbergii",NCBITaxon:40692
    TrematomusNewnesi, "Dusky notothen", "Trematomus newnesi (dusky notothen)", "Trematomus newnesi",NCBITaxon:35730
    TrematomusPennellii, "Sharp-spined notothen", "Trematomus pennellii (sharp-spined notothen)", "Trematomus pennellii",NCBITaxon:36198
    TriakisScyllium, "Banded houndshark", "Triakis scyllium (banded houndshark)", "Triakis scyllium",NCBITaxon:30494
    TrichechusManatusLatirostris, "Florida manatee", "Trichechus manatus latirostris (Florida manatee)", "Trichechus manatus latirostris",NCBITaxon:127582
    TrichosurusVulpecula, "Common brushtail", "Trichosurus vulpecula (common brushtail)", "Trichosurus vulpecula",NCBITaxon:9337
    TursiopsAduncus, "Indo-pacific bottlenose dolphin", "Tursiops aduncus (Indo-pacific bottlenose dolphin)", "Tursiops aduncus",NCBITaxon:79784
    TursiopsTruncatus, "Common bottlenose dolphin", "Tursiops truncatus (common bottlenose dolphin)", "Tursiops truncatus",NCBITaxon:9739
    VicugnaPacos, "Alpaca", "Vicugna pacos (alpaca)", "Vicugna pacos",NCBITaxon:30538
    Xenopus, "Xenopus", "Xenopus", "Xenopus",NCBITaxon:262014
    XenopusLaevis, "African clawed frog", "Xenopus laevis (African clawed frog)", "Xenopus laevis",NCBITaxon:8355
    XenopusLaevisOrGilli, "African or Cape clawed frog", "Xenopus laevis/gilli", "Xenopus laevis/gilli",NCBITaxon:8359
    XenopusSp, "Clawed frog", "Xenopus sp. (clawed frog)", "Xenopus sp.",NCBITaxon:8356
    XenopusTropicalis, "Tropical clawed frog", "Xenopus tropicalis (tropical clawed frog)", "Xenopus tropicalis",NCBITaxon:8364
);
