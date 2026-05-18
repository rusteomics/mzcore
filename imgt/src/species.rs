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
    ($($identifier:ident, $common:expr, $imgt:expr, $scientific:expr, $curie:expr)*) => {
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
            /// The common name for this species, eg `Human`
            pub const fn common_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $common,)*
                }
            }

            /// The name IMGT uses to identify this species, eg `Homo sapiens (human)`
            pub const fn imgt_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $imgt,)*
                }
            }

            /// The common name for this species, eg `Homo sapiens`
            pub const fn scientific_name(&self) -> &'static str {
                match self {
                    $(Self::$identifier => $scientific,)*
                }
            }

            /// The enum name for this species, eg `HomoSapiens`
            pub const fn ident(&self) -> &'static str {
                match self {
                    $(Self::$identifier => stringify!($identifier),)*
                }
            }

            /// The CURIE to identify this species within the NCBITaxonomy Ontology
            pub const fn curie(&self) -> mzcv::Curie {
                match self {
                    $(Self::$identifier => $curie,)*
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
    };
}

species!(
    AcanthopagrusSchlegelii, "Black porgy", "Acanthopagrus schlegelii (black porgy)", "Acanthopagrus schlegelii", curie!(NCBITaxon:72011)
    AcipenserBaerii, "Siberian sturgeon", "Acipenser baerii (Siberian sturgeon)", "Acipenser baerii", curie!(NCBITaxon:27689)
    AcipenserGueldenstaedtii, "Russian sturgeon", "Acipenser gueldenstaedtii (Russian sturgeon)", "Acipenser gueldenstaedtii", curie!(NCBITaxon:7902)
    AcipenserRuthenus, "Sterlet", "Acipenser ruthenus (sterlet)", "Acipenser ruthenus", curie!(NCBITaxon:7906)
    AcipenserSchrenckii, "Amur sturgeon", "Acipenser schrenckii (Amur sturgeon)", "Acipenser schrenckii", curie!(NCBITaxon:111304)
    AcipenserSinensis, "Chinese sturgeon", "Acipenser sinensis (Chinese sturgeon)", "Acipenser sinensis", curie!(NCBITaxon:61970)
    AiluropodaMelanoleuca, "Giant panda", "Ailuropoda melanoleuca (giant panda)", "Ailuropoda melanoleuca", curie!(NCBITaxon:9646)
    AlligatorSinensis, "Chinese alligator", "Alligator sinensis (Chinese alligator)", "Alligator sinensis", curie!(NCBITaxon:38654)
    AmblyrajaGeorgiana, "Antarctic starry skate", "Amblyraja georgiana (Antarctic starry skate)", "Amblyraja georgiana", curie!(NCBITaxon:1051904)
    AmblyrajaHyperborea, "Arctic skate", "Amblyraja hyperborea (Arctic skate)", "Amblyraja hyperborea", curie!(NCBITaxon:315322)
    AmbystomaMexicanum, "Axolotl", "Ambystoma mexicanum (axolotl)", "Ambystoma mexicanum", curie!(NCBITaxon:8296)
    AmeivaAmeiva, "Jungle runners", "Ameiva ameiva", "Ameiva ameiva", curie!(NCBITaxon:8535)
    AmiaCalva, "Bowfin", "Amia calva (bowfin)", "Amia calva", curie!(NCBITaxon:7924)
    AmphiprionClarkii, "Yellowtail clownfish", "Amphiprion clarkii (yellowtail clownfish)", "Amphiprion clarkii", curie!(NCBITaxon:80970)
    AnarhichasMinor, "Spotted wolffish", "Anarhichas minor (spotted wolffish)", "Anarhichas minor", curie!(NCBITaxon:65739)
    AnasPlatyrhynchos, "Mallard", "Anas platyrhynchos (mallard)", "Anas platyrhynchos", curie!(NCBITaxon:8839)
    AnguillaAnguilla, "European eel", "Anguilla anguilla (European eel)", "Anguilla anguilla", curie!(NCBITaxon:7936)
    AnguillaJaponica, "Japanese eel", "Anguilla japonica (Japanese eel)", "Anguilla japonica", curie!(NCBITaxon:7937)
    AnolisCarolinensis, "Green anole", "Anolis carolinensis (green anole)", "Anolis carolinensis", curie!(NCBITaxon:28377)
    AnoplopomaFimbria, "Sablefish", "Anoplopoma fimbria (sablefish)", "Anoplopoma fimbria", curie!(NCBITaxon:229290)
    AnserAnser, "Domestic goose", "Anser anser (Domestic goose)", "Anser anser", curie!(NCBITaxon:8843)
    AnserAnserDomesticus, "Domestic goose", "Anser anser domesticus", "Anser anser domesticus", curie!(NCBITaxon:8843) // TODO: maybe merge?
    AnserCaerulescens, "Snow goose", "Anser caerulescens (Snow goose)", "Anser caerulescens", curie!(NCBITaxon:8849)
    AnserSpGIGHV2011, "Geese GIGHV2011", "Anser sp. GIGHV2011", "Anser sp. GIGHV2011", curie!(NCBITaxon:980882)
    AnserSpGIGLV2009, "Geese GIGLV2009", "Anser sp. GIGLV2009", "Anser sp. GIGLV2009", curie!(NCBITaxon:713105)
    AotusAzarai, "Azara's night monkey", "Aotus azarai (Azara's night monkey)", "Aotus azarai", curie!(NCBITaxon:30591)
    AotusNancymaae, "Ma's night monkey", "Aotus nancymaae (Ma's night monkey)", "Aotus nancymaae", curie!(NCBITaxon:37293)
    AotusTrivirgatus, "Douroucouli", "Aotus trivirgatus (douroucouli)", "Aotus trivirgatus", curie!(NCBITaxon:9505)
    ApisCerana, "Asiatic honeybee", "Apis cerana (Asiatic honeybee)", "Apis cerana", curie!(NCBITaxon:7461)
    ArgyropelecusHemigymnus, "Half-naked hatchetfish", "Argyropelecus hemigymnus (half-naked hatchetfish)", "Argyropelecus hemigymnus", curie!(NCBITaxon:68504)
    AtelesBelzebuth, "White-bellied spider monkey", "Ateles belzebuth (white-bellied spider monkey)", "Ateles belzebuth", curie!(NCBITaxon:9507)
    AtelesGeoffroyi, "Black-handed spider monkey", "Ateles geoffroyi (black-handed spider monkey)", "Ateles geoffroyi", curie!(NCBITaxon:9509)
    AtherinaBoyeri, "Big-scale sand smelt", "Atherina boyeri (big-scale sand smelt)", "Atherina boyeri", curie!(NCBITaxon:87785)
    BalaenopteraAcutorostrata, "Minke whale", "Balaenoptera acutorostrata (minke whale)", "Balaenoptera acutorostrata", curie!(NCBITaxon:9767)
    BalaenopteraOmurai, "Omura's baleen whale", "Balaenoptera omurai (Omura's baleen whale)", "Balaenoptera omurai", curie!(NCBITaxon:255217)
    BathyrajaAlbomaculata, "White-dotted skate", "Bathyraja albomaculata (white-dotted skate)", "Bathyraja albomaculata", curie!(NCBITaxon:458555)
    BathyrajaBrachyurops, "Broadnose skate", "Bathyraja brachyurops (broadnose skate)", "Bathyraja brachyurops", curie!(NCBITaxon:458556)
    BathyrajaEatonii, "Eaton's skate", "Bathyraja eatonii (Eaton's skate)", "Bathyraja eatonii", curie!(NCBITaxon:298348)
    BosGaurus, "Gaur", "Bos gaurus (gaur)", "Bos gaurus", curie!(NCBITaxon:9904)
    BosIndicus, "Domestic zebu", "Bos indicus (zebu cattle)", "Bos indicus", curie!(NCBITaxon:9915)
    BosJavanicus, "Banteng", "Bos javanicus (banteng)", "Bos javanicus", curie!(NCBITaxon:9906)
    BosTaurus, "Domestic bovine", "Bos taurus (bovine)", "Bos taurus", curie!(NCBITaxon:9913)
    BosTaurusXBosIndicus, "Bos taurus and Bos indicus cross", "Bos taurus x Bos indicus", "Bos taurus x Bos indicus", curie!(NCBITaxon:30522)
    BosIndicisXBosTaurus, "Bos indicus and Bos taurus cross", "Bos indicus x Bos taurus (hybrid cattle)", "Bos indicus x Bos taurus", curie!(NCBITaxon:30522) // TODO: maybe merge?
    BovichtusDiacanthus, "Tristan clipfish", "Bovichtus diacanthus", "Bovichtus diacanthus", curie!(NCBITaxon:483258)
    BubalusBubalis, "Water buffalo", "Bubalus bubalis (water buffalo)", "Bubalus bubalis", curie!(NCBITaxon:89462)
    BuergeriaBuergeri, "Buerger's frog", "Buergeria buergeri (Buerger's frog)", "Buergeria buergeri", curie!(NCBITaxon:191197)
    CaimanCrocodilus, "Spectacled caiman", "Caiman crocodilus (spectacled caiman)", "Caiman crocodilus", curie!(NCBITaxon:8499)
    CairinaMoschata, "Muscovy duck", "Cairina moschata (Muscovy duck)", "Cairina moschata", curie!(NCBITaxon:8855)
    CallithrixJacchus, "white-tufted-ear marmoset", "Callithrix jacchus (white-tufted-ear marmoset)", "Callithrix jacchus", curie!(NCBITaxon:9483)
    CallorhinchusMilii, "Elephant shark", "Callorhinchus milii (elephant shark)", "Callorhinchus milii", curie!(NCBITaxon:7868)
    Camelidae, "Camels", "Camelidae", "Camelidae", curie!(NCBITaxon:9835)
    CamelusBactrianus, "Bactrian camel", "Camelus bactrianus (Bactrian camel)", "Camelus bactrianus", curie!(NCBITaxon:9837)
    CamelusDromedarius, "Arabian camel", "Camelus dromedarius (Arabian camel)", "Camelus dromedarius", curie!(NCBITaxon:9838)
    CanisLupus, "Gray wolf", "Canis lupus (gray wolf)", "Canis lupus", curie!(NCBITaxon:9612)
    CanisLupusFamiliaris, "Domestic dog", "Canis lupus familiaris (dog)", "Canis lupus familiaris", curie!(NCBITaxon:9615)
    CanisSp, "Dogs", "Canis sp.", "Canis sp.", curie!(NCBITaxon:9616)
    CapraHircus, "Domestic goat", "Capra hircus (goat)", "Capra hircus", curie!(NCBITaxon:9925)
    CarassiusAuratus, "Goldfish", "Carassius auratus (goldfish)", "Carassius auratus", curie!(NCBITaxon:7957)
    CarassiusLangsdorfii, "Japanese silver crucian carp", "Carassius langsdorfii (Japanese silver crucian carp)", "Carassius langsdorfii", curie!(NCBITaxon:138676)
    CarcharhinusLeucas, "Bull shark", "Carcharhinus leucas (bull shark)", "Carcharhinus leucas", curie!(NCBITaxon:46612)
    CarcharhinusPlumbeus, "Sandbar shark", "Carcharhinus plumbeus (sandbar shark)", "Carcharhinus plumbeus", curie!(NCBITaxon:7808)
    CarlitoSyrichta, "Philippine tarsier", "Carlito syrichta (Philippine tarsier)", "Carlito syrichta", curie!(NCBITaxon:1868482)
    CarolliaPerspicillata, "Seba's short-tailed bat", "Carollia perspicillata (Seba's short-tailed bat)", "Carollia perspicillata", curie!(NCBITaxon:40233)
    CaviaPorcellus, "Domestic guinea pig", "Cavia porcellus (domestic guinea pig)", "Cavia porcellus", curie!(NCBITaxon:10141)
    CephalopachusBancanus, "Horsfield's tarsier", "Cephalopachus bancanus (Horsfield's tarsier)", "Cephalopachus bancanus", curie!(NCBITaxon:9477) // TODO: not done from this point
    CeratotheriumSimum, "White rhinoceros", "Ceratotherium simum (white rhinoceros)", "Ceratotherium simum", curie!(NCBITaxon:9807)
    CercocebusAtys, "Sooty mangabey", "Cercocebus atys (sooty mangabey)", "Cercocebus atys", curie!(NCBITaxon:9531)
    CercocebusTorquatus, "Collared mangabey", "Cercocebus torquatus (collared mangabey)", "Cercocebus torquatus", curie!(NCBITaxon:9530)
    CervusElaphusHispanicus, "Spanish red deer", "Cervus elaphus hispanicus", "Cervus elaphus hispanicus", curie!(NCBITaxon:141655)
    ChaenocephalusAceratus, "Blackfin icefish", "Chaenocephalus aceratus (blackfin icefish)", "Chaenocephalus aceratus", curie!(NCBITaxon:36190)
    ChampsocephalusEsox, "Pike icefish", "Champsocephalus esox (pike icefish)", "Champsocephalus esox", curie!(NCBITaxon:159716)
    ChannaArgus, "Northern snakehead", "Channa argus (northern snakehead)", "Channa argus", curie!(NCBITaxon:215402)
    ChannaStriata, "Snakehead murrel", "Channa striata (snakehead murrel)", "Channa striata", curie!(NCBITaxon:64152)
    ChaunaTorquata, "Southern screamer", "Chauna torquata (southern screamer)", "Chauna torquata", curie!(NCBITaxon:30388)
    ChelonAuratus, "Golden grey mullet", "Chelon auratus (golden grey mullet)", "Chelon auratus", curie!(NCBITaxon:48191)
    ChionodracoHamatus, "Antarctic icefish", "Chionodraco hamatus (Antarctic icefish)", "Chionodraco hamatus", curie!(NCBITaxon:36188)
    ChionodracoRastrospinosus, "Ocellated icefish", "Chionodraco rastrospinosus (ocellated icefish)", "Chionodraco rastrospinosus", curie!(NCBITaxon:34790)
    ChlorocebusAethiops, "Grivet", "Chlorocebus aethiops (grivet)", "Chlorocebus aethiops", curie!(NCBITaxon:9534)
    ClupeaPallasii, "Pacific herring", "Clupea pallasii (Pacific herring)", "Clupea pallasii", curie!(NCBITaxon:30724)
    ColobusGuereza, "Mantled guereza", "Colobus guereza (mantled guereza)", "Colobus guereza", curie!(NCBITaxon:33548)
    ColobusPolykomos, "King colobus", "Colobus polykomos (king colobus)", "Colobus polykomos", curie!(NCBITaxon:9572)
    CricetinaeSp, "Hamster", "Cricetinae gen. sp. (Hamster)", "Cricetinae sp.", curie!(NCBITaxon:36483)
    CricetulusMigratorius, "Armenian hamster", "Cricetulus migratorius (Armenian hamster)", "Cricetulus migratorius", curie!(NCBITaxon:3122392)
    CrocodylusSiamensis, "Siamese crocodile", "Crocodylus siamensis (Siamese crocodile)", "Crocodylus siamensis", curie!(NCBITaxon:68455)
    CtenopharyngodonIdella, "Grass carp", "Ctenopharyngodon idella (grass carp)", "Ctenopharyngodon idella", curie!(NCBITaxon:7959)
    CygnodracoMawsoni, "Mawson's dragonfish", "Cygnodraco mawsoni (Mawson's dragonfish)", "Cygnodraco mawsoni", curie!(NCBITaxon:8216)
    CynoglossusSemilaevis, "Tongue sole", "Cynoglossus semilaevis (tongue sole)", "Cynoglossus semilaevis", curie!(NCBITaxon:244447)
    CynopterusSphinx, "Indian short-nosed fruit bat", "Cynopterus sphinx (Indian short-nosed fruit bat)", "Cynopterus sphinx", curie!(NCBITaxon:9400)
    CyprinusCarpio, "Common carp", "Cyprinus carpio (common carp)", "Cyprinus carpio", curie!(NCBITaxon:7962)
    DanioRerio, "Zebrafish", "Danio rerio (zebrafish)", "Danio rerio", curie!(NCBITaxon:7955)
    DaubentoniaMadagascariensis, "Aye-aye", "Daubentonia madagascariensis (aye-aye)", "Daubentonia madagascariensis", curie!(NCBITaxon:31869)
    DelphinapterusLeucas, "Beluga whale", "Delphinapterus leucas (beluga whale)", "Delphinapterus leucas", curie!(NCBITaxon:9749)
    DelphinusCapensis, "Long-beaked common dolphin", "Delphinus capensis (long-beaked common dolphin)", "Delphinus capensis", curie!(NCBITaxon:103584)
    DicentrarchusLabrax, "European seabass", "Dicentrarchus labrax (European seabass)", "Dicentrarchus labrax", curie!(NCBITaxon:13489)
    DissostichusMawsoni, "Antarctic toothfish", "Dissostichus mawsoni (Antarctic toothfish)", "Dissostichus mawsoni", curie!(NCBITaxon:36200)
    DrosophilaMelanogaster, "Fruit fly", "Drosophila melanogaster (fruit fly)", "Drosophila melanogaster", curie!(NCBITaxon:7227)
    ElapheTaeniura, "Beauty snake", "Elaphe taeniura (beauty snake)", "Elaphe taeniura", curie!(NCBITaxon:74398)
    EleginopsMaclovinus, "Patagonian blennie", "Eleginops maclovinus (Patagonian blennie)", "Eleginops maclovinus", curie!(NCBITaxon:56733)
    ElopsAaurus, "Ladyfish", "Elops saurus (ladyfish)", "Elops saurus", curie!(NCBITaxon:7928)
    EpinephelusAkaara, "Hong Kong grouper", "Epinephelus akaara (Hong Kong grouper)", "Epinephelus akaara", curie!(NCBITaxon:215347)
    EpinephelusCoioides, "Orange-spotted grouper", "Epinephelus coioides (orange-spotted grouper)", "Epinephelus coioides", curie!(NCBITaxon:94232)
    EptatretusBurgeri, "Inshore hagfish", "Eptatretus burgeri (inshore hagfish)", "Eptatretus burgeri", curie!(NCBITaxon:7764)
    EptesicusFuscus, "Big brown bat", "Eptesicus fuscus (big brown bat)", "Eptesicus fuscus", curie!(NCBITaxon:29078)
    EquusAsinus, "Ass", "Equus asinus (ass)", "Equus asinus", curie!(NCBITaxon:9793)
    EquusBurchelliiAntiquorum, "Burchell's zebra", "Equus burchellii antiquorum", "Equus burchellii antiquorum", curie!(NCBITaxon:89252)
    EquusCaballus, "Domestic horse", "Equus caballus (horse)", "Equus caballus", curie!(NCBITaxon:9796)
    EquusQuaggaBurchellii, "Burchell's zebra", "Equus quagga burchellii (Burchell's zebra)", "Equus quagga burchellii", curie!(NCBITaxon:89252)
    ErythrocebusPatas, "Red guenon", "Erythrocebus patas (red guenon)", "Erythrocebus patas", curie!(NCBITaxon:9538)
    EscherichiaColi, "E. coli", "Escherichia coli (E. coli)", "Escherichia coli", curie!(NCBITaxon:562)
    EsoxLucius, "Northern pike", "Esox lucius (northern pike)", "Esox lucius", curie!(NCBITaxon:8010)
    EublepharisMacularius, "Leopard gecko", "Eublepharis macularius (Leopard gecko)", "Eublepharis macularius", curie!(NCBITaxon:481883)
    EulemurFulvus, "Brown lemur", "Eulemur fulvus (brown lemur)", "Eulemur fulvus", curie!(NCBITaxon:13515)
    FelineLeukemiaVirus, "Feline leukemia virus", "Feline leukemia virus", "Feline leukemia virus", curie!(NCBITaxon:11768)
    FelisCatus, "Domestic cat", "Felis catus (domestic cat)", "Felis catus", curie!(NCBITaxon:9685)
    FelisSp, "Cats", "Felis sp.", "Felis sp.", curie!(NCBITaxon:9687)
    GadusMorhua, "Atlantic cod", "Gadus morhua (Atlantic cod)", "Gadus morhua", curie!(NCBITaxon:8049)
    GalagoSenegalensis, "Senegal galago", "Galago senegalensis (Senegal galago)", "Galago senegalensis", curie!(NCBITaxon:9465)
    GallusGallus, "Domestic chicken", "Gallus gallus (chicken)", "Gallus gallus", curie!(NCBITaxon:9031)
    GasterosteusAculeatus, "Three-spined stickleback", "Gasterosteus aculeatus (three-spined stickleback)", "Gasterosteus aculeatus", curie!(NCBITaxon:69293)
    GinglymostomaCirratum, "Nurse shark", "Ginglymostoma cirratum (nurse shark)", "Ginglymostoma cirratum", curie!(NCBITaxon:7801)
    GobionotothenGibberifrons, "Humped rockcod", "Gobionotothen gibberifrons (humped rockcod)", "Gobionotothen gibberifrons", curie!(NCBITaxon:36202)
    GorillaGorilla, "Western gorilla", "Gorilla gorilla (western gorilla)", "Gorilla gorilla", curie!(NCBITaxon:9593)
    GorillaGorillaGorilla, "Western lowland gorilla", "Gorilla gorilla gorilla (western lowland gorilla)", "Gorilla gorilla gorilla", curie!(NCBITaxon:9595)
    GrampusGriseus, "Risso's dolphin", "Grampus griseus (Risso's dolphin)", "Grampus griseus", curie!(NCBITaxon:83653)
    GymnodracoAcuticeps, "Ploughfish", "Gymnodraco acuticeps", "Gymnodraco acuticeps", curie!(NCBITaxon:8218)
    GymnogypsCalifornianus, "California condor", "Gymnogyps californianus (California condor)", "Gymnogyps californianus", curie!(NCBITaxon:33616)
    HaemorhousMexicanus, "House finch", "Haemorhous mexicanus (house finch)", "Haemorhous mexicanus", curie!(NCBITaxon:30427)
    HemibagrusMacropterus, "Largefin longbarbel catfish", "Hemibagrus macropterus", "Hemibagrus macropterus", curie!(NCBITaxon:175789)
    HepacivirusC, "Hepacivirus C", "Hepacivirus C", "Hepacivirus C", curie!(NCBITaxon:3052230)
    HeterocephalusGlaber, "Naked mole-rat", "Heterocephalus glaber (naked mole-rat)", "Heterocephalus glaber", curie!(NCBITaxon:10181)
    HeterodontusFrancisci, "Horn shark", "Heterodontus francisci (horn shark)", "Heterodontus francisci", curie!(NCBITaxon:7792)
    HippoglossusHippoglossus, "Atlantic halibut", "Hippoglossus hippoglossus (Atlantic halibut)", "Hippoglossus hippoglossus", curie!(NCBITaxon:8267)
    HistiodracoVelifer, "Histiodraco velifer", "Histiodraco velifer", "Histiodraco velifer", curie!(NCBITaxon:36208)
    HomoSapiens, "Human", "Homo sapiens (human)", "Homo sapiens", curie!(NCBITaxon:9606)
    HoolockHoolock, "Hoolock gibbon", "Hoolock hoolock (hoolock gibbon)", "Hoolock hoolock", curie!(NCBITaxon:61851)
    HusoHuso, "Beluga", "Huso huso (beluga)", "Huso huso", curie!(NCBITaxon:61971)
    HydrolagusColliei, "Spotted ratfish", "Hydrolagus colliei (spotted ratfish)", "Hydrolagus colliei", curie!(NCBITaxon:7873)
    HylobatesLar, "Common gibbon", "Hylobates lar (common gibbon)", "Hylobates lar", curie!(NCBITaxon:9580)
    IctalurusPunctatus, "Channel catfish", "Ictalurus punctatus (channel catfish)", "Ictalurus punctatus", curie!(NCBITaxon:7998)
    IsoodonMacrourus, "Northern brown bandicoot", "Isoodon macrourus (northern brown bandicoot)", "Isoodon macrourus", curie!(NCBITaxon:37698)
    KogiaSima, "Dwarf sperm whale", "Kogia sima (dwarf sperm whale)", "Kogia sima", curie!(NCBITaxon:9752)
    LabeobarbusIntermedius, "Labeobarbus intermedius", "Labeobarbus intermedius", "Labeobarbus intermedius", curie!(NCBITaxon:40831)
    LabeoRohita, "Rohu", "Labeo rohita (rohu)", "Labeo rohita", curie!(NCBITaxon:84645)
    LamaGlama, "Llama", "Lama glama (llama)", "Lama glama", curie!(NCBITaxon:9844)
    LarimichthysCrocea, "Large yellow croaker", "Larimichthys crocea (large yellow croaker)", "Larimichthys crocea", curie!(NCBITaxon:215358)
    LatimeriaChalumnae, "Coelacanth", "Latimeria chalumnae (coelacanth)", "Latimeria chalumnae", curie!(NCBITaxon:7897)
    LatimeriaMenadoensis, "Menado coelacanth", "Latimeria menadoensis (Menado coelacanth)", "Latimeria menadoensis", curie!(NCBITaxon:106881)
    LatrisLineata, "Striped trumpeter", "Latris lineata (striped trumpeter)", "Latris lineata", curie!(NCBITaxon:97162)
    LemurCatta, "Ring-tailed lemur", "Lemur catta (Ring-tailed lemur)", "Lemur catta", curie!(NCBITaxon:9447)
    LeontopithecusRosalia, "Golden lion tamarin", "Leontopithecus rosalia (golden lion tamarin)", "Leontopithecus rosalia", curie!(NCBITaxon:30588)
    LepilemurRuficaudatus, "Red-tailed sportive lemur", "Lepilemur ruficaudatus (red-tailed sportive lemur)", "Lepilemur ruficaudatus", curie!(NCBITaxon:78866)
    LepisosteusOsseus, "Longnose gar", "Lepisosteus osseus (longnose gar)", "Lepisosteus osseus", curie!(NCBITaxon:34771)
    LepusAmericanus, "Snowshoe hare", "Lepus americanus (snowshoe hare)", "Lepus americanus", curie!(NCBITaxon:48086)
    LepusCalifornicus, "Black-tailed jackrabbit", "Lepus californicus (black-tailed jackrabbit)", "Lepus californicus", curie!(NCBITaxon:48087)
    LepusCallotis, "White-sided jackrabbit", "Lepus callotis (white-sided jackrabbit)", "Lepus callotis", curie!(NCBITaxon:62619)
    LepusCapensis, "Brown hare", "Lepus capensis (brown hare)", "Lepus capensis", curie!(NCBITaxon:9981)
    LepusCastroviejoi, "Broom Hare", "Lepus castroviejoi (Broom Hare)", "Lepus castroviejoi", curie!(NCBITaxon:187398)
    LepusEuropaeus, "European hare", "Lepus europaeus (European hare)", "Lepus europaeus", curie!(NCBITaxon:9983)
    LepusGranatensis, "Granada hare", "Lepus granatensis (Granada hare)", "Lepus granatensis", curie!(NCBITaxon:100182)
    LepusSaxatilis, "Scrub hare", "Lepus saxatilis (scrub hare)", "Lepus saxatilis", curie!(NCBITaxon:63224)
    LepusTimidus, "Mountain hare", "Lepus timidus (Mountain hare)", "Lepus timidus", curie!(NCBITaxon:62621)
    LeucorajaErinacea, "Little skate", "Leucoraja erinacea (little skate)", "Leucoraja erinacea", curie!(NCBITaxon:7782)
    LipotesVexillifer, "Yangtze River dolphin", "Lipotes vexillifer (Yangtze River dolphin)", "Lipotes vexillifer", curie!(NCBITaxon:118797)
    LutjanusSanguineus, "Humphead snapper", "Lutjanus sanguineus (humphead snapper)", "Lutjanus sanguineus", curie!(NCBITaxon:264213)
    MacacaArctoides, "Stump-tailed macaque", "Macaca arctoides (stump-tailed macaque)", "Macaca arctoides", curie!(NCBITaxon:9540)
    MacacaAssamensis, "Assam macaque", "Macaca assamensis (Assam macaque)", "Macaca assamensis", curie!(NCBITaxon:9551)
    MacacaCyclopis, "Taiwan macaque", "Macaca cyclopis (Taiwan macaque)", "Macaca cyclopis", curie!(NCBITaxon:78449)
    MacacaFascicularis, "Crab-eating macaque", "Macaca fascicularis (crab-eating macaque)", "Macaca fascicularis", curie!(NCBITaxon:9541)
    MacacaMulatta, "Rhesus monkey", "Macaca mulatta (Rhesus monkey)", "Macaca mulatta", curie!(NCBITaxon:9544)
    MacacaNemestrina, "Pig-tailed macaque", "Macaca nemestrina (pig-tailed macaque)", "Macaca nemestrina", curie!(NCBITaxon:9545)
    MacacaSilenus, "Liontail macaque", "Macaca silenus (liontail macaque)", "Macaca silenus", curie!(NCBITaxon:54601)
    MacacaThibetana, "Pere David's macaque", "Macaca thibetana (Pere David's macaque)", "Macaca thibetana", curie!(NCBITaxon:54602)
    MarecaStrepera, "Gadwall", "Mareca strepera (gadwall)", "Mareca strepera", curie!(NCBITaxon:75861)
    MarmotaHimalayana, "Himalayan marmot", "Marmota himalayana (Himalayan marmot)", "Marmota himalayana", curie!(NCBITaxon:93163)
    MarmotaMonax, "Woodchuck", "Marmota monax (woodchuck)", "Marmota monax", curie!(NCBITaxon:9995)
    MauremysMutica, "Yellowpond turtle", "Mauremys mutica (yellowpond turtle)", "Mauremys mutica", curie!(NCBITaxon:74926)
    MelanogrammusAeglefinus, "Haddock", "Melanogrammus aeglefinus (haddock)", "Melanogrammus aeglefinus", curie!(NCBITaxon:8056)
    MeleagrisGallopavo, "Turkey", "Meleagris gallopavo (turkey)", "Meleagris gallopavo", curie!(NCBITaxon:9103)
    MerionesUnguiculatus, "Mongolian gerbil", "Meriones unguiculatus (Mongolian gerbil)", "Meriones unguiculatus", curie!(NCBITaxon:10047)
    MesocricetusAuratus, "Golden hamster", "Mesocricetus auratus (golden hamster)", "Mesocricetus auratus", curie!(NCBITaxon:10036)
    MicrocebusMurinus, "Gray mouse lemur", "Microcebus murinus (gray mouse lemur)", "Microcebus murinus", curie!(NCBITaxon:30608)
    MonodelphisDomestica, "Gray short-tailed opossum", "Monodelphis domestica (gray short-tailed opossum)", "Monodelphis domestica", curie!(NCBITaxon:13616)
    MurinaeSp, "Old world rats and mice", "Murinae gen. sp.", "Murinae sp.", curie!(NCBITaxon:3148847)
    Mus, "mouse", "Mus (mouse)", "Mus", curie!(NCBITaxon:10088)
    MusCookii, "Cook's mouse", "Mus cookii (Cook's mouse)", "Mus cookii", curie!(NCBITaxon:10098)
    MusMinutoides, "Southern African pygmy mouse", "Mus minutoides (Southern African pygmy mouse)", "Mus minutoides", curie!(NCBITaxon:10105)
    MusMusculus, "House mouse", "Mus musculus (house mouse)", "Mus musculus", curie!(NCBITaxon:10090)
    MusMusculusCastaneus, "Southeastern Asian house mouse", "Mus musculus castaneus (southeastern Asian house mouse)", "Mus musculus castaneus", curie!(NCBITaxon:10091)
    MusMusculusDomesticus, "Western European house mouse", "Mus musculus domesticus (western European house mouse)", "Mus musculus domesticus", curie!(NCBITaxon:10092)
    MusMusculusMolossinus, "Japanese wild mouse", "Mus musculus molossinus (Japanese wild mouse)", "Mus musculus molossinus", curie!(NCBITaxon:57486)
    MusMusculusMusculus, "Eastern European house mouse", "Mus musculus musculus (eastern European house mouse)", "Mus musculus musculus", curie!(NCBITaxon:39442)
    MusPahari, "Shrew mouse", "Mus pahari (shrew mouse)", "Mus pahari", curie!(NCBITaxon:10093)
    MusSaxicola, "Spiny mouse", "Mus saxicola (spiny mouse)", "Mus saxicola", curie!(NCBITaxon:10094)
    MusSp, "Mice", "Mus sp. (mice)", "Mus sp.", curie!(NCBITaxon:10095)
    MusSpretus, "Western wild mouse", "Mus spretus (western wild mouse)", "Mus spretus", curie!(NCBITaxon:10096)
    MustelaPutoriusFuro, "Domestic ferret", "Mustela putorius furo (domestic ferret)", "Mustela putorius furo", curie!(NCBITaxon:9669)
    MustelaSp, "Ferret", "Mustela sp.", "Mustela sp.", curie!(NCBITaxon:30549)
    MyotisLucifugus, "Little brown bat", "Myotis lucifugus (little brown bat)", "Myotis lucifugus", curie!(NCBITaxon:59463)
    NeophocaenaPhocaenoides, "Indo-Pacific finless porpoise", "Neophocaena phocaenoides (Indo-Pacific finless porpoise)", "Neophocaena phocaenoides", curie!(NCBITaxon:34892)
    NeogaleVison, "American mink", "Neogale vison (American mink)", "Neogale vison", curie!(NCBITaxon:452646)
    NomascusConcolor, "Black crested gibbon", "Nomascus concolor (Black crested gibbon)", "Nomascus concolor", curie!(NCBITaxon:29089)
    NotamacropusEugenii, "Tammar wallaby", "Notamacropus eugenii (tammar wallaby)", "Notamacropus eugenii", curie!(NCBITaxon:9315)
    NothocricetulusMigratorius, "Armenian hamster", "Nothocricetulus migratorius (Armenian hamster)", "Nothocricetulus migratorius", curie!(NCBITaxon:3122392)
    NototheniaCoriiceps, "Black rockcod", "Notothenia coriiceps (black rockcod)", "Notothenia coriiceps", curie!(NCBITaxon:8208)
    NycticebusCoucang, "Slow loris", "Nycticebus coucang (slow loris)", "Nycticebus coucang", curie!(NCBITaxon:9470)
    OncorhynchusGorbuscha, "Pink salmon", "Oncorhynchus gorbuscha (pink salmon)", "Oncorhynchus gorbuscha", curie!(NCBITaxon:8017)
    OncorhynchusMykiss, "Rainbow trout", "Oncorhynchus mykiss (rainbow trout)", "Oncorhynchus mykiss", curie!(NCBITaxon:8022)
    OncorhynchusTshawytscha, "Chinook salmon", "Oncorhynchus tshawytscha (Chinook salmon)", "Oncorhynchus tshawytscha", curie!(NCBITaxon:74940)
    OrectolobusMaculatus, "Spotted wobbegong", "Orectolobus maculatus (spotted wobbegong)", "Orectolobus maculatus", curie!(NCBITaxon:168098)
    OreochromisNiloticus, "Nile tilapia", "Oreochromis niloticus (Nile tilapia)", "Oreochromis niloticus", curie!(NCBITaxon:8128)
    OrnithorhynchusAnatinus, "Platypus", "Ornithorhynchus anatinus (platypus)", "Ornithorhynchus anatinus", curie!(NCBITaxon:9258)
    OryctolagusCuniculus, "Rabbit", "Oryctolagus cuniculus (rabbit)", "Oryctolagus cuniculus", curie!(NCBITaxon:9986)
    OryctolagusCuniculusAlgirus, "European rabbit", "Oryctolagus cuniculus algirus", "Oryctolagus cuniculus algirus", curie!(NCBITaxon:230741)
    OryctolagusCuniculusCuniculus, "Rabbit", "Oryctolagus cuniculus cuniculus", "Oryctolagus cuniculus cuniculus", curie!(NCBITaxon:568996)
    OryziasLatipes, "Japanese medaka", "Oryzias latipes (Japanese medaka)", "Oryzias latipes", curie!(NCBITaxon:8090)
    OryziasMelastigma, "Indian medaka", "Oryzias melastigma (Indian medaka)", "Oryzias melastigma", curie!(NCBITaxon:30732)
    OtolemurCrassicaudatus, "Thick-tailed bush baby", "Otolemur crassicaudatus (thick-tailed bush baby)", "Otolemur crassicaudatus", curie!(NCBITaxon:9463)
    OvisAries, "Domestic sheep", "Ovis aries (sheep)", "Ovis aries", curie!(NCBITaxon:9940)
    OvisSp, "Sheep", "Ovis sp.", "Ovis sp.", curie!(NCBITaxon:9939)
    PacifastacusLeniusculus, "Signal crayfish", "Pacifastacus leniusculus (signal crayfish)", "Pacifastacus leniusculus", curie!(NCBITaxon:6720)
    PagetopsisMacropterus, "Pagetopsis macropterus", "Pagetopsis macropterus", "Pagetopsis macropterus", curie!(NCBITaxon:36194)
    PagrusMajor, "Red seabream", "Pagrus major (red seabream)", "Pagrus major", curie!(NCBITaxon:143350)
    PangasianodonHypophthalmus, "Striped catfish", "Pangasianodon hypophthalmus (striped catfish)", "Pangasianodon hypophthalmus", curie!(NCBITaxon:310915)
    PanPaniscus, "Pygmy chimpanzee", "Pan paniscus (pygmy chimpanzee)", "Pan paniscus", curie!(NCBITaxon:9597)
    PantheraPardus, "Leopard", "Panthera pardus (leopard)", "Panthera pardus", curie!(NCBITaxon:9691)
    PanTroglodytes, "Chimpanzee", "Pan troglodytes (chimpanzee)", "Pan troglodytes", curie!(NCBITaxon:9598)
    PanTroglodytesVerus, "Western chimpanzee", "Pan troglodytes verus", "Pan troglodytes verus", curie!(NCBITaxon:37012)
    PapioAnubis, "Olive baboon", "Papio anubis (olive baboon)", "Papio anubis", curie!(NCBITaxon:9555)
    PapioAnubisAnubis, "Olive baboon anubis", "Papio anubis anubis", "Papio anubis anubis", curie!(NCBITaxon:211508)
    PapioHamadryas, "Hamadryas baboon", "Papio hamadryas (hamadryas baboon)", "Papio hamadryas", curie!(NCBITaxon:9557)
    PapioPapio, "Guinea baboon", "Papio papio (Guinea baboon)", "Papio papio", curie!(NCBITaxon:100937)
    ParalichthysOlivaceus, "Japanese flounder", "Paralichthys olivaceus (Japanese flounder)", "Paralichthys olivaceus", curie!(NCBITaxon:8255)
    PelodiscusSinensis, "Chinese soft-shelled turtle", "Pelodiscus sinensis (Chinese soft-shelled turtle)", "Pelodiscus sinensis", curie!(NCBITaxon:13735)
    PelteobagrusFulvidraco, "Yellow catfish", "Pelteobagrus fulvidraco (yellow catfish)", "Pelteobagrus fulvidraco", curie!(NCBITaxon:1234273)
    PerdixPerdix, "Grey partridge", "Perdix perdix (grey partridge)", "Perdix perdix", curie!(NCBITaxon:9052)
    PeromyscusManiculatus, "North American deer mouse", "Peromyscus maniculatus (North American deer mouse)", "Peromyscus maniculatus", curie!(NCBITaxon:10042)
    PetromyzonMarinus, "Sea lamprey", "Petromyzon marinus (sea lamprey)", "Petromyzon marinus", curie!(NCBITaxon:7757)
    PhascogaleCalura, "Red-tailed phascogale", "Phascogale calura (red-tailed phascogale)", "Phascogale calura", curie!(NCBITaxon:38776)
    PhasianusColchicus, "Ring-necked pheasant", "Phasianus colchicus (Ring-necked pheasant)", "Phasianus colchicus", curie!(NCBITaxon:9054)
    PhyseterCatodon, "Sperm whale", "Physeter catodon (sperm whale)", "Physeter catodon", curie!(NCBITaxon:9755)
    PitheciaPithecia, "White-faced saki", "Pithecia pithecia (white-faced saki)", "Pithecia pithecia", curie!(NCBITaxon:43777)
    PlataleaAjaja, "Roseate spoonbil", "Platalea ajaja", "Platalea ajaja", curie!(NCBITaxon:371920)
    Platyrrhini, "New World monkeys", "Platyrrhini (New World monkeys)", "Platyrrhini", curie!(NCBITaxon:9479)
    PlecoglossusAltivelisAltivelis, "Ayu sweetfish", "Plecoglossus altivelis altivelis", "Plecoglossus altivelis altivelis", curie!(NCBITaxon:281464)
    PleurodelesWaltl, "Iberian ribbed newt", "Pleurodeles waltl (Iberian ribbed newt)", "Pleurodeles waltl", curie!(NCBITaxon:8319)
    PogonophryneScotti, "Pogonophryne scotti", "Pogonophryne scotti", "Pogonophryne scotti", curie!(NCBITaxon:36210)
    PolyprionOxygeneios, "Hāpuku", "Polyprion oxygeneios", "Polyprion oxygeneios", curie!(NCBITaxon:334912)
    PongoAbelii, "Sumatran orangutan", "Pongo abelii (Sumatran orangutan)", "Pongo abelii", curie!(NCBITaxon:9601)
    PongoPygmaeus, "Bornean orangutan", "Pongo pygmaeus (Bornean orangutan)", "Pongo pygmaeus", curie!(NCBITaxon:9600)
    PresbytisComata, "Grizzled leaf monkey", "Presbytis comata (grizzled leaf monkey)", "Presbytis comata", curie!(NCBITaxon:78452)
    PresbytisFemoralis, "Banded leaf monkey", "Presbytis femoralis (banded leaf monkey)", "Presbytis femoralis", curie!(NCBITaxon:78453)
    PresbytisMelalophos, "Mitred leaf monkey", "Presbytis melalophos (mitred leaf monkey)", "Presbytis melalophos", curie!(NCBITaxon:78451)
    PropithecusVerreauxi, "White sifaka", "Propithecus verreauxi (white sifaka)", "Propithecus verreauxi", curie!(NCBITaxon:34825)
    ProtopterusAethiopicus, "Marbled lungfish", "Protopterus aethiopicus (marbled lungfish)", "Protopterus aethiopicus", curie!(NCBITaxon:7886)
    PseudobatosProductus, "Shovelnose guitarfish", "Pseudobatos productus (shovelnose guitarfish)", "Pseudobatos productus", curie!(NCBITaxon:170824)
    PteropusAlecto, "Black flying fox", "Pteropus alecto (black flying fox)", "Pteropus alecto", curie!(NCBITaxon:9402)
    PythonBivittatus, "Burmese python", "Python bivittatus (Burmese python)", "Python bivittatus", curie!(NCBITaxon:176946)
    RachycentronCanadum, "Cobia", "Rachycentron canadum (cobia)", "Rachycentron canadum", curie!(NCBITaxon:141264)
    RajaEglanteria, "Clearnose skate", "Raja eglanteria (clearnose skate)", "Raja eglanteria", curie!(NCBITaxon:3360502)
    RattusFuscipes, "Bush rat", "Rattus fuscipes (bush rat)", "Rattus fuscipes", curie!(NCBITaxon:10119)
    RattusLeucopus, "Mottle-tailed rat", "Rattus leucopus (mottle-tailed rat)", "Rattus leucopus", curie!(NCBITaxon:10115)
    RattusNorvegicus, "Norway rat", "Rattus norvegicus (Norway rat)", "Rattus norvegicus", curie!(NCBITaxon:10116)
    RattusRattus, "Black rat", "Rattus rattus (black rat)", "Rattus rattus", curie!(NCBITaxon:10117)
    RattusSordidus, "Australian dusky field rat", "Rattus sordidus (Australian dusky field rat)", "Rattus sordidus", curie!(NCBITaxon:10120)
    RattusSp, "Rats", "Rattus sp. (rats)", "Rattus sp.", curie!(NCBITaxon:10118)
    RattusTunneyi, "Tunney's rat", "Rattus tunneyi (Tunney's rat)", "Rattus tunneyi", curie!(NCBITaxon:10121)
    RattusVillosissimus, "Long-haired rat", "Rattus villosissimus (long-haired rat)", "Rattus villosissimus", curie!(NCBITaxon:10122)
    RhinocerosUnicornis, "Greater Indian rhinoceros", "Rhinoceros unicornis (greater Indian rhinoceros)", "Rhinoceros unicornis", curie!(NCBITaxon:9809)
    RousettusLeschenaultii, "Leschenault's rousette", "Rousettus leschenaultii (Leschenault's rousette)", "Rousettus leschenaultii", curie!(NCBITaxon:9408)
    SaccharomycesCerevisiae, "baker's yeast", "Saccharomyces cerevisiae (baker's yeast)", "Saccharomyces cerevisiae", curie!(NCBITaxon:4932)
    SaguinusLabiatus, "Red-chested mustached tamarin", "Saguinus labiatus (red-chested mustached tamarin)", "Saguinus labiatus", curie!(NCBITaxon:78454)
    SaguinusMidas, "Midas tamarin", "Saguinus midas (Midas tamarin)", "Saguinus midas", curie!(NCBITaxon:30586)
    SaguinusOedipus, "Cotton-top tamarin", "Saguinus oedipus (cotton-top tamarin)", "Saguinus oedipus", curie!(NCBITaxon:9490)
    SaimiriBoliviensisBoliviensis, "Bolivian squirrel monkey", "Saimiri boliviensis boliviensis (Bolivian squirrel monkey)", "Saimiri boliviensis boliviensis", curie!(NCBITaxon:39432)
    SaimiriSciureus, "Common squirrel monkey", "Saimiri sciureus (common squirrel monkey)", "Saimiri sciureus", curie!(NCBITaxon:9521)
    SalmoMarmoratus, "Salmo marmoratus", "Salmo marmoratus", "Salmo marmoratus", curie!(NCBITaxon:33518)
    SalmoSalar, "Atlantic salmon", "Salmo salar (Atlantic salmon)", "Salmo salar", curie!(NCBITaxon:8030)
    SalmoTrutta, "River trout", "Salmo trutta (river trout)", "Salmo trutta", curie!(NCBITaxon:8032)
    SalvelinusAlpinus, "Arctic char", "Salvelinus alpinus (Arctic char)", "Salvelinus alpinus", curie!(NCBITaxon:8036)
    SanderVitreus, "Walleye", "Sander vitreus (walleye)", "Sander vitreus", curie!(NCBITaxon:283036)
    SapajusApella, "Tufted capuchin", "Sapajus apella (Tufted capuchin)", "Sapajus apella", curie!(NCBITaxon:9515)
    SchizophyllumCommune, "Schizophyllum commune", "Schizophyllum commune", "Schizophyllum commune", curie!(NCBITaxon:5334)
    SciaenopsOcellatus, "Red drum", "Sciaenops ocellatus (red drum)", "Sciaenops ocellatus", curie!(NCBITaxon:76340)
    ScophthalmusMaximus, "Turbot", "Scophthalmus maximus (turbot)", "Scophthalmus maximus", curie!(NCBITaxon:52904)
    ScyliorhinusCanicula, "Smaller spotted catshark", "Scyliorhinus canicula (smaller spotted catshark)", "Scyliorhinus canicula", curie!(NCBITaxon:7830)
    SeriolaQuinqueradiata, "Japanese amberjack", "Seriola quinqueradiata (Japanese amberjack)", "Seriola quinqueradiata", curie!(NCBITaxon:8161)
    SilurusAsotus, "Amur catfish", "Silurus asotus (Amur catfish)", "Silurus asotus", curie!(NCBITaxon:30991)
    SilurusMeridionalis, "Silurus meridionalis", "Silurus meridionalis", "Silurus meridionalis", curie!(NCBITaxon:175797)
    SinipercaChuatsi, "Mandarin fish", "Siniperca chuatsi (mandarin fish)", "Siniperca chuatsi", curie!(NCBITaxon:119488)
    SousaChinensis, "Indo-pacific humpbacked dolphin", "Sousa chinensis (Indo-pacific humpbacked dolphin)", "Sousa chinensis", curie!(NCBITaxon:103600)
    SparusAurata, "Gilthead seabream", "Sparus aurata (gilthead seabream)", "Sparus aurata", curie!(NCBITaxon:8175)
    SphoeroidesNephelus, "Southern puffer", "Sphoeroides nephelus (southern puffer)", "Sphoeroides nephelus", curie!(NCBITaxon:39110)
    SqualusAcanthias, "Spiny dogfish", "Squalus acanthias (spiny dogfish)", "Squalus acanthias", curie!(NCBITaxon:7797)
    StegastesLeucostictus, "Beaugregory", "Stegastes leucostictus (beaugregory)", "Stegastes leucostictus", curie!(NCBITaxon:258452)
    StegastesPartitus, "Bicolor damselfish", "Stegastes partitus (bicolor damselfish)", "Stegastes partitus", curie!(NCBITaxon:144197)
    StenellaAttenuata, "Bridled dolphin", "Stenella attenuata (bridled dolphin)", "Stenella attenuata", curie!(NCBITaxon:9735)
    StenellaCoeruleoalba, "Striped dolphin", "Stenella coeruleoalba (striped dolphin)", "Stenella coeruleoalba", curie!(NCBITaxon:9737)
    StreptomycesAvidinii, "Streptomyces avidinii", "Streptomyces avidinii", "Streptomyces avidinii", curie!(NCBITaxon:1895)
    StreptomycesViridochromogenes, "Streptomyces viridochromogenes", "Streptomyces viridochromogenes", "Streptomyces viridochromogenes", curie!(NCBITaxon:1938)
    StruthioCamelus, "African ostrich", "Struthio camelus (African ostrich)", "Struthio camelus", curie!(NCBITaxon:8801)
    SuncusMurinus, "House shrew", "Suncus murinus (house shrew)", "Suncus murinus", curie!(NCBITaxon:9378)
    SusScrofa, "Domestic pig", "Sus scrofa (pig)", "Sus scrofa", curie!(NCBITaxon:9823)
    SylvilagusCunicularis, "Mexican cottontail", "Sylvilagus cunicularis (Mexican cottontail)", "Sylvilagus cunicularis", curie!(NCBITaxon:3370854)
    SylvilagusFloridanus, "Eastern cottontail", "Sylvilagus floridanus (eastern cottontail)", "Sylvilagus floridanus", curie!(NCBITaxon:9988)
    SymphalangusSyndactylus, "Siamang", "Symphalangus syndactylus (siamang)", "Symphalangus syndactylus", curie!(NCBITaxon:9590)
    TachyglossusAculeatus, "Australian echidna", "Tachyglossus aculeatus (Australian echidna)", "Tachyglossus aculeatus", curie!(NCBITaxon:9261)
    TachysurusFulvidraco, "Yellow catfish", "Tachysurus fulvidraco (yellow catfish)", "Tachysurus fulvidraco", curie!(NCBITaxon:1234273)
    TachysurusVachellii, "Tachysurus vachellii", "Tachysurus vachellii", "Tachysurus vachellii", curie!(NCBITaxon:175792)
    TaeniopygiaGuttata, "Zebra finch", "Taeniopygia guttata (zebra finch)", "Taeniopygia guttata", curie!(NCBITaxon:59729)
    TakifuguRubripes, "Torafugu", "Takifugu rubripes (torafugu)", "Takifugu rubripes", curie!(NCBITaxon:31033)
    TarsiusDentatus, "Dian's tarsier", "Tarsius dentatus (Dian's tarsier)", "Tarsius dentatus", curie!(NCBITaxon:449501)
    TarsiusLariang, "Lariang tarsier", "Tarsius lariang (Lariang tarsier)", "Tarsius lariang", curie!(NCBITaxon:630277)
    TarsiusSyrichta, "Philippine tarsier", "Tarsius syrichta (Philippine tarsier)", "Tarsius syrichta", curie!(NCBITaxon:1868482)
    TetraodonNigroviridis, "Spotted green pufferfish", "Tetraodon nigroviridis (spotted green pufferfish)", "Tetraodon nigroviridis", curie!(NCBITaxon:99883)
    TrachemysScripta, "Red-eared slider turtle", "Trachemys scripta (red-eared slider turtle)", "Trachemys scripta", curie!(NCBITaxon:34903)
    TrachemysScriptaElegans, "Red-eared slider turtle elegans", "Trachemys scripta elegans", "Trachemys scripta elegans", curie!(NCBITaxon:31138)
    TrachypithecusCristatus, "Silvery lutung", "Trachypithecus cristatus (Silvery lutung)", "Trachypithecus cristatus", curie!(NCBITaxon:122765)
    TrachypithecusObscurus, "Dusky leaf-monkey", "Trachypithecus obscurus (Dusky leaf-monkey)", "Trachypithecus obscurus", curie!(NCBITaxon:54181)
    TrematomusBernacchii, "Emerald rockcod", "Trematomus bernacchii (emerald rockcod)", "Trematomus bernacchii", curie!(NCBITaxon:40690)
    TrematomusHansoni, "Striped rockcod", "Trematomus hansoni (striped rockcod)", "Trematomus hansoni", curie!(NCBITaxon:36196)
    TrematomusLoennbergii, "Deepwater notothen", "Trematomus loennbergii (deepwater notothen)", "Trematomus loennbergii", curie!(NCBITaxon:40692)
    TrematomusNewnesi, "Dusky notothen", "Trematomus newnesi (dusky notothen)", "Trematomus newnesi", curie!(NCBITaxon:35730)
    TrematomusPennellii, "Sharp-spined notothen", "Trematomus pennellii (sharp-spined notothen)", "Trematomus pennellii", curie!(NCBITaxon:36198)
    TriakisScyllium, "Banded houndshark", "Triakis scyllium (banded houndshark)", "Triakis scyllium", curie!(NCBITaxon:30494)
    TrichechusManatusLatirostris, "Florida manatee", "Trichechus manatus latirostris (Florida manatee)", "Trichechus manatus latirostris", curie!(NCBITaxon:127582)
    TrichosurusVulpecula, "Common brushtail", "Trichosurus vulpecula (common brushtail)", "Trichosurus vulpecula", curie!(NCBITaxon:9337)
    TursiopsAduncus, "Indo-pacific bottlenose dolphin", "Tursiops aduncus (Indo-pacific bottlenose dolphin)", "Tursiops aduncus", curie!(NCBITaxon:79784)
    TursiopsTruncatus, "Common bottlenose dolphin", "Tursiops truncatus (common bottlenose dolphin)", "Tursiops truncatus", curie!(NCBITaxon:9739)
    VicugnaPacos, "Alpaca", "Vicugna pacos (alpaca)", "Vicugna pacos", curie!(NCBITaxon:30538)
    Xenopus, "Xenopus", "Xenopus", "Xenopus", curie!(NCBITaxon:262014)
    XenopusLaevis, "African clawed frog", "Xenopus laevis (African clawed frog)", "Xenopus laevis", curie!(NCBITaxon:8355)
    XenopusLaevisOrGilli, "African or Cape clawed frog", "Xenopus laevis/gilli", "Xenopus laevis/gilli", curie!(NCBITaxon:8359)
    XenopusSp, "Clawed frog", "Xenopus sp. (clawed frog)", "Xenopus sp.", curie!(NCBITaxon:8356)
    XenopusTropicalis, "Tropical clawed frog", "Xenopus tropicalis (tropical clawed frog)", "Xenopus tropicalis", curie!(NCBITaxon:8364)
);
