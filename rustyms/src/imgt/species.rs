#![allow(dead_code, unreachable_patterns)]
use serde::{Deserialize, Serialize};
use std::fmt::Display;

macro_rules! species {
    ($($identifier:ident, $common:expr, $imgt:expr, $scientific:expr)*) => {
        /// All species available in the IMGT dataset. Look at the main documentation to see which actually have data provided.
        #[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Ord, PartialOrd)]
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

            /// Get the species name from IMGT name tag.
            /// # Errors
            /// `Err` when the name could not be recognised (it is case sensitive).
            /// `Ok(None)` when it is recognised as a species used by IMGT, but it is not a proper species (vector/plasmid etc).
            pub fn from_imgt(s: &str) -> Result<Option<Self>, ()> {
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
                    _ => Err(()),
                }
            }
        }

        /// The list with all usable names
        const SPECIES_PARSE_LIST: &[(&str, Species)] = &[
            $(($common, Species::$identifier),)*
            $(($scientific, Species::$identifier),)*
        ];

        impl std::str::FromStr for Species {
            type Err = crate::error::CustomError;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                let s = s.trim().to_lowercase();
                for (name, species) in SPECIES_PARSE_LIST {
                    if name.to_lowercase() == s {
                        return Ok(*species);
                    }
                }

                let options: Vec<String> = SPECIES_PARSE_LIST.iter()
                    .map(|option| option.0.to_lowercase())
                    .collect();
                let options: Vec<&str> = options.iter()
                    .map(|option| option.as_str())
                    .collect();

                Err(crate::error::CustomError::error("Unknown species name", "The provided name could not be recognised as a species name.", crate::error::Context::show(s.as_str())).with_suggestions(similar::get_close_matches(s.as_str(), &options, 3, 0.75)))
            }
        }
    };
}

species!(
    AcanthopagrusSchlegelii, "Black porgy", "Acanthopagrus schlegelii (black porgy)", "Acanthopagrus schlegelii"
    AcipenserBaerii, "Siberian sturgeon", "Acipenser baerii (Siberian sturgeon)", "Acipenser baerii"
    AcipenserGueldenstaedtii, "Russian sturgeon", "Acipenser gueldenstaedtii (Russian sturgeon)", "Acipenser gueldenstaedtii"
    AcipenserRuthenus, "Sterlet", "Acipenser ruthenus (sterlet)", "Acipenser ruthenus"
    AcipenserSchrenckii, "Amur sturgeon", "Acipenser schrenckii (Amur sturgeon)", "Acipenser schrenckii"
    AcipenserSinensis, "Chinese sturgeon", "Acipenser sinensis (Chinese sturgeon)", "Acipenser sinensis"
    AiluropodaMelanoleuca, "Giant panda", "Ailuropoda melanoleuca (giant panda)", "Ailuropoda melanoleuca"
    AlligatorSinensis, "Chinese alligator", "Alligator sinensis (Chinese alligator)", "Alligator sinensis"
    AmblyrajaGeorgiana, "Antarctic starry skate", "Amblyraja georgiana (Antarctic starry skate)", "Amblyraja georgiana"
    AmblyrajaHyperborea, "Arctic skate", "Amblyraja hyperborea (Arctic skate)", "Amblyraja hyperborea"
    AmbystomaMexicanum, "Axolotl", "Ambystoma mexicanum (axolotl)", "Ambystoma mexicanum"
    AmeivaAmeiva, "Jungle runners", "Ameiva ameiva", "Ameiva ameiva"
    AmiaCalva, "Bowfin", "Amia calva (bowfin)", "Amia calva"
    AmphiprionClarkii, "Yellowtail clownfish", "Amphiprion clarkii (yellowtail clownfish)", "Amphiprion clarkii"
    AnarhichasMinor, "Spotted wolffish", "Anarhichas minor (spotted wolffish)", "Anarhichas minor"
    AnasPlatyrhynchos, "Mallard", "Anas platyrhynchos (mallard)", "Anas platyrhynchos"
    AnguillaAnguilla, "European eel", "Anguilla anguilla (European eel)", "Anguilla anguilla"
    AnguillaJaponica, "Japanese eel", "Anguilla japonica (Japanese eel)", "Anguilla japonica"
    AnolisCarolinensis, "Green anole", "Anolis carolinensis (green anole)", "Anolis carolinensis"
    AnoplopomaFimbria, "Sablefish", "Anoplopoma fimbria (sablefish)", "Anoplopoma fimbria"
    AnserAnser, "Domestic goose", "Anser anser (Domestic goose)", "Anser anser"
    AnserAnserDomesticus, "Domestic goose", "Anser anser domesticus", "Anser anser domesticus"
    AnserCaerulescens, "Snow goose", "Anser caerulescens (Snow goose)", "Anser caerulescens"
    AnserSpGIGHV2011, "Geese GIGHV2011", "Anser sp. GIGHV2011", "Anser sp. GIGHV2011"
    AnserSpGIGLV2009, "Geese GIGLV2009", "Anser sp. GIGLV2009", "Anser sp. GIGLV2009"
    AotusAzarai, "Azara's night monkey", "Aotus azarai (Azara's night monkey)", "Aotus azarai"
    AotusNancymaae, "Ma's night monkey", "Aotus nancymaae (Ma's night monkey)", "Aotus nancymaae"
    AotusTrivirgatus, "Douroucouli", "Aotus trivirgatus (douroucouli)", "Aotus trivirgatus"
    ApisCerana, "Asiatic honeybee", "Apis cerana (Asiatic honeybee)", "Apis cerana"
    ArgyropelecusHemigymnus, "Half-naked hatchetfish", "Argyropelecus hemigymnus (half-naked hatchetfish)", "Argyropelecus hemigymnus"
    AtelesBelzebuth, "White-bellied spider monkey", "Ateles belzebuth (white-bellied spider monkey)", "Ateles belzebuth"
    AtelesGeoffroyi, "Black-handed spider monkey", "Ateles geoffroyi (black-handed spider monkey)", "Ateles geoffroyi"
    AtherinaBoyeri, "Big-scale sand smelt", "Atherina boyeri (big-scale sand smelt)", "Atherina boyeri"
    BalaenopteraAcutorostrata, "Minke whale", "Balaenoptera acutorostrata (minke whale)", "Balaenoptera acutorostrata"
    BalaenopteraOmurai, "Omura's baleen whale", "Balaenoptera omurai (Omura's baleen whale)", "Balaenoptera omurai"
    BathyrajaAlbomaculata, "White-dotted skate", "Bathyraja albomaculata (white-dotted skate)", "Bathyraja albomaculata"
    BathyrajaBrachyurops, "Broadnose skate", "Bathyraja brachyurops (broadnose skate)", "Bathyraja brachyurops"
    BathyrajaEatonii, "Eaton's skate", "Bathyraja eatonii (Eaton's skate)", "Bathyraja eatonii"
    BosGaurus, "Gaur", "Bos gaurus (gaur)", "Bos gaurus"
    BosIndicus, "Domestic zebu", "Bos indicus (zebu cattle)", "Bos indicus"
    BosJavanicus, "Banteng", "Bos javanicus (banteng)", "Bos javanicus"
    BosTaurus, "Domestic bovine", "Bos taurus (bovine)", "Bos taurus"
    BosTaurusXBosIndicus, "Bos taurus and Bos indicus cross", "Bos taurus x Bos indicus", "Bos taurus x Bos indicus"
    BosIndicisXBosTaurus, "Bos indicus and Bos taurus cross", "Bos indicus x Bos taurus (hybrid cattle)", "Bos indicus x Bos taurus"
    BovichtusDiacanthus, "Tristan clipfish", "Bovichtus diacanthus", "Bovichtus diacanthus"
    BubalusBubalis, "Water buffalo", "Bubalus bubalis (water buffalo)", "Bubalus bubalis"
    BuergeriaBuergeri, "Buerger's frog", "Buergeria buergeri (Buerger's frog)", "Buergeria buergeri"
    CaimanCrocodilus, "Spectacled caiman", "Caiman crocodilus (spectacled caiman)", "Caiman crocodilus"
    CairinaMoschata, "Muscovy duck", "Cairina moschata (Muscovy duck)", "Cairina moschata"
    CallithrixJacchus, "white-tufted-ear marmoset", "Callithrix jacchus (white-tufted-ear marmoset)", "Callithrix jacchus"
    CallorhinchusMilii, "Elephant shark", "Callorhinchus milii (elephant shark)", "Callorhinchus milii"
    Camelidae, "Camels", "Camelidae", "Camelidae"
    CamelusBactrianus, "Bactrian camel", "Camelus bactrianus (Bactrian camel)", "Camelus bactrianus"
    CamelusDromedarius, "Arabian camel", "Camelus dromedarius (Arabian camel)", "Camelus dromedarius"
    CanisLupus, "Gray wolf", "Canis lupus (gray wolf)", "Canis lupus"
    CanisLupusFamiliaris, "Domestic dog", "Canis lupus familiaris (dog)", "Canis lupus familiaris"
    CanisSp, "Dogs", "Canis sp.", "Canis sp."
    CapraHircus, "Domestic goat", "Capra hircus (goat)", "Capra hircus"
    CarassiusAuratus, "Goldfish", "Carassius auratus (goldfish)", "Carassius auratus"
    CarassiusLangsdorfii, "Japanese silver crucian carp", "Carassius langsdorfii (Japanese silver crucian carp)", "Carassius langsdorfii"
    CarcharhinusLeucas, "Bull shark", "Carcharhinus leucas (bull shark)", "Carcharhinus leucas"
    CarcharhinusPlumbeus, "Sandbar shark", "Carcharhinus plumbeus (sandbar shark)", "Carcharhinus plumbeus"
    CarlitoSyrichta, "Philippine tarsier", "Carlito syrichta (Philippine tarsier)", "Carlito syrichta"
    CarolliaPerspicillata, "Seba's short-tailed bat", "Carollia perspicillata (Seba's short-tailed bat)", "Carollia perspicillata"
    CaviaPorcellus, "Domestic guinea pig", "Cavia porcellus (domestic guinea pig)", "Cavia porcellus"
    CephalopachusBancanus, "Horsfield's tarsier", "Cephalopachus bancanus (Horsfield's tarsier)", "Cephalopachus bancanus"
    CeratotheriumSimum, "White rhinoceros", "Ceratotherium simum (white rhinoceros)", "Ceratotherium simum"
    CercocebusAtys, "Sooty mangabey", "Cercocebus atys (sooty mangabey)", "Cercocebus atys"
    CercocebusTorquatus, "Collared mangabey", "Cercocebus torquatus (collared mangabey)", "Cercocebus torquatus"
    CervusElaphusHispanicus, "Spanish red deer", "Cervus elaphus hispanicus", "Cervus elaphus hispanicus"
    ChaenocephalusAceratus, "Blackfin icefish", "Chaenocephalus aceratus (blackfin icefish)", "Chaenocephalus aceratus"
    ChampsocephalusEsox, "Pike icefish", "Champsocephalus esox (pike icefish)", "Champsocephalus esox"
    ChannaArgus, "Northern snakehead", "Channa argus (northern snakehead)", "Channa argus"
    ChannaStriata, "Snakehead murrel", "Channa striata (snakehead murrel)", "Channa striata"
    ChaunaTorquata, "Southern screamer", "Chauna torquata (southern screamer)", "Chauna torquata"
    ChelonAuratus, "Golden grey mullet", "Chelon auratus (golden grey mullet)", "Chelon auratus"
    ChionodracoHamatus, "Antarctic icefish", "Chionodraco hamatus (Antarctic icefish)", "Chionodraco hamatus"
    ChionodracoRastrospinosus, "Ocellated icefish", "Chionodraco rastrospinosus (ocellated icefish)", "Chionodraco rastrospinosus"
    ChlorocebusAethiops, "Grivet", "Chlorocebus aethiops (grivet)", "Chlorocebus aethiops"
    ClupeaPallasii, "Pacific herring", "Clupea pallasii (Pacific herring)", "Clupea pallasii"
    ColobusGuereza, "Mantled guereza", "Colobus guereza (mantled guereza)", "Colobus guereza"
    ColobusPolykomos, "King colobus", "Colobus polykomos (king colobus)", "Colobus polykomos"
    CricetinaeSp, "Hamster", "Cricetinae gen. sp. (Hamster)", "Cricetinae sp."
    CricetulusMigratorius, "Armenian hamster", "Cricetulus migratorius (Armenian hamster)", "Cricetulus migratorius"
    CrocodylusSiamensis, "Siamese crocodile", "Crocodylus siamensis (Siamese crocodile)", "Crocodylus siamensis"
    CtenopharyngodonIdella, "Grass carp", "Ctenopharyngodon idella (grass carp)", "Ctenopharyngodon idella"
    CygnodracoMawsoni, "Mawson's dragonfish", "Cygnodraco mawsoni (Mawson's dragonfish)", "Cygnodraco mawsoni"
    CynoglossusSemilaevis, "Tongue sole", "Cynoglossus semilaevis (tongue sole)", "Cynoglossus semilaevis"
    CynopterusSphinx, "Indian short-nosed fruit bat", "Cynopterus sphinx (Indian short-nosed fruit bat)", "Cynopterus sphinx"
    CyprinusCarpio, "Common carp", "Cyprinus carpio (common carp)", "Cyprinus carpio"
    DanioRerio, "Zebrafish", "Danio rerio (zebrafish)", "Danio rerio"
    DaubentoniaMadagascariensis, "Aye-aye", "Daubentonia madagascariensis (aye-aye)", "Daubentonia madagascariensis"
    DelphinapterusLeucas, "Beluga whale", "Delphinapterus leucas (beluga whale)", "Delphinapterus leucas"
    DelphinusCapensis, "Long-beaked common dolphin", "Delphinus capensis (long-beaked common dolphin)", "Delphinus capensis"
    DicentrarchusLabrax, "European seabass", "Dicentrarchus labrax (European seabass)", "Dicentrarchus labrax"
    DissostichusMawsoni, "Antarctic toothfish", "Dissostichus mawsoni (Antarctic toothfish)", "Dissostichus mawsoni"
    DrosophilaMelanogaster, "Fruit fly", "Drosophila melanogaster (fruit fly)", "Drosophila melanogaster"
    ElapheTaeniura, "Beauty snake", "Elaphe taeniura (beauty snake)", "Elaphe taeniura"
    EleginopsMaclovinus, "Patagonian blennie", "Eleginops maclovinus (Patagonian blennie)", "Eleginops maclovinus"
    ElopsAaurus, "Ladyfish", "Elops saurus (ladyfish)", "Elops saurus"
    EpinephelusAkaara, "Hong Kong grouper", "Epinephelus akaara (Hong Kong grouper)", "Epinephelus akaara"
    EpinephelusCoioides, "Orange-spotted grouper", "Epinephelus coioides (orange-spotted grouper)", "Epinephelus coioides"
    EptatretusBurgeri, "Inshore hagfish", "Eptatretus burgeri (inshore hagfish)", "Eptatretus burgeri"
    EptesicusFuscus, "Big brown bat", "Eptesicus fuscus (big brown bat)", "Eptesicus fuscus"
    EquusAsinus, "Ass", "Equus asinus (ass)", "Equus asinus"
    EquusBurchelliiAntiquorum, "Burchell's zebra", "Equus burchellii antiquorum", "Equus burchellii antiquorum"
    EquusCaballus, "Domestic horse", "Equus caballus (horse)", "Equus caballus"
    EquusQuaggaBurchellii, "Burchell's zebra", "Equus quagga burchellii (Burchell's zebra)", "Equus quagga burchellii"
    ErythrocebusPatas, "Red guenon", "Erythrocebus patas (red guenon)", "Erythrocebus patas"
    EscherichiaColi, "E. coli", "Escherichia coli (E. coli)", "Escherichia coli"
    EsoxLucius, "Northern pike", "Esox lucius (northern pike)", "Esox lucius"
    EublepharisMacularius, "Leopard gecko", "Eublepharis macularius (Leopard gecko)", "Eublepharis macularius"
    EulemurFulvus, "Brown lemur", "Eulemur fulvus (brown lemur)", "Eulemur fulvus"
    FelineLeukemiaVirus, "Feline leukemia virus", "Feline leukemia virus", "Feline leukemia virus"
    FelisCatus, "Domestic cat", "Felis catus (domestic cat)", "Felis catus"
    FelisSp, "Cats", "Felis sp.", "Felis sp."
    GadusMorhua, "Atlantic cod", "Gadus morhua (Atlantic cod)", "Gadus morhua"
    GalagoSenegalensis, "Senegal galago", "Galago senegalensis (Senegal galago)", "Galago senegalensis"
    GallusGallus, "Domestic chicken", "Gallus gallus (chicken)", "Gallus gallus"
    GasterosteusAculeatus, "Three-spined stickleback", "Gasterosteus aculeatus (three-spined stickleback)", "Gasterosteus aculeatus"
    GinglymostomaCirratum, "Nurse shark", "Ginglymostoma cirratum (nurse shark)", "Ginglymostoma cirratum"
    GobionotothenGibberifrons, "Humped rockcod", "Gobionotothen gibberifrons (humped rockcod)", "Gobionotothen gibberifrons"
    GorillaGorilla, "Western gorilla", "Gorilla gorilla (western gorilla)", "Gorilla gorilla"
    GorillaGorillaGorilla, "Western lowland gorilla", "Gorilla gorilla gorilla (western lowland gorilla)", "Gorilla gorilla gorilla"
    GrampusGriseus, "Risso's dolphin", "Grampus griseus (Risso's dolphin)", "Grampus griseus"
    GymnodracoAcuticeps, "Ploughfish", "Gymnodraco acuticeps", "Gymnodraco acuticeps"
    GymnogypsCalifornianus, "California condor", "Gymnogyps californianus (California condor)", "Gymnogyps californianus"
    HaemorhousMexicanus, "House finch", "Haemorhous mexicanus (house finch)", "Haemorhous mexicanus"
    HemibagrusMacropterus, "Largefin longbarbel catfish", "Hemibagrus macropterus", "Hemibagrus macropterus"
    HepacivirusC, "Hepacivirus C", "Hepacivirus C", "Hepacivirus C"
    HeterocephalusGlaber, "Naked mole-rat", "Heterocephalus glaber (naked mole-rat)", "Heterocephalus glaber"
    HeterodontusFrancisci, "Horn shark", "Heterodontus francisci (horn shark)", "Heterodontus francisci"
    HippoglossusHippoglossus, "Atlantic halibut", "Hippoglossus hippoglossus (Atlantic halibut)", "Hippoglossus hippoglossus"
    HistiodracoVelifer, "Histiodraco velifer", "Histiodraco velifer", "Histiodraco velifer"
    HomoSapiens, "Human", "Homo sapiens (human)", "Homo sapiens"
    HoolockHoolock, "Hoolock gibbon", "Hoolock hoolock (hoolock gibbon)", "Hoolock hoolock"
    HusoHuso, "Beluga", "Huso huso (beluga)", "Huso huso"
    HydrolagusColliei, "Spotted ratfish", "Hydrolagus colliei (spotted ratfish)", "Hydrolagus colliei"
    HylobatesLar, "Common gibbon", "Hylobates lar (common gibbon)", "Hylobates lar"
    IctalurusPunctatus, "Channel catfish", "Ictalurus punctatus (channel catfish)", "Ictalurus punctatus"
    IsoodonMacrourus, "Northern brown bandicoot", "Isoodon macrourus (northern brown bandicoot)", "Isoodon macrourus"
    KogiaSima, "Dwarf sperm whale", "Kogia sima (dwarf sperm whale)", "Kogia sima"
    LabeobarbusIntermedius, "Labeobarbus intermedius", "Labeobarbus intermedius", "Labeobarbus intermedius"
    LabeoRohita, "Rohu", "Labeo rohita (rohu)", "Labeo rohita"
    LamaGlama, "Llama", "Lama glama (llama)", "Lama glama"
    LarimichthysCrocea, "Large yellow croaker", "Larimichthys crocea (large yellow croaker)", "Larimichthys crocea"
    LatimeriaChalumnae, "Coelacanth", "Latimeria chalumnae (coelacanth)", "Latimeria chalumnae"
    LatimeriaMenadoensis, "Menado coelacanth", "Latimeria menadoensis (Menado coelacanth)", "Latimeria menadoensis"
    LatrisLineata, "Striped trumpeter", "Latris lineata (striped trumpeter)", "Latris lineata"
    LemurCatta, "Ring-tailed lemur", "Lemur catta (Ring-tailed lemur)", "Lemur catta"
    LeontopithecusRosalia, "Golden lion tamarin", "Leontopithecus rosalia (golden lion tamarin)", "Leontopithecus rosalia"
    LepilemurRuficaudatus, "Red-tailed sportive lemur", "Lepilemur ruficaudatus (red-tailed sportive lemur)", "Lepilemur ruficaudatus"
    LepisosteusOsseus, "Longnose gar", "Lepisosteus osseus (longnose gar)", "Lepisosteus osseus"
    LepusAmericanus, "Snowshoe hare", "Lepus americanus (snowshoe hare)", "Lepus americanus"
    LepusCalifornicus, "Black-tailed jackrabbit", "Lepus californicus (black-tailed jackrabbit)", "Lepus californicus"
    LepusCallotis, "White-sided jackrabbit", "Lepus callotis (white-sided jackrabbit)", "Lepus callotis"
    LepusCapensis, "Brown hare", "Lepus capensis (brown hare)", "Lepus capensis"
    LepusCastroviejoi, "Broom Hare", "Lepus castroviejoi (Broom Hare)", "Lepus castroviejoi"
    LepusEuropaeus, "European hare", "Lepus europaeus (European hare)", "Lepus europaeus"
    LepusGranatensis, "Granada hare", "Lepus granatensis (Granada hare)", "Lepus granatensis"
    LepusSaxatilis, "Scrub hare", "Lepus saxatilis (scrub hare)", "Lepus saxatilis"
    LepusTimidus, "Mountain hare", "Lepus timidus (Mountain hare)", "Lepus timidus"
    LeucorajaErinacea, "Little skate", "Leucoraja erinacea (little skate)", "Leucoraja erinacea"
    LipotesVexillifer, "Yangtze River dolphin", "Lipotes vexillifer (Yangtze River dolphin)", "Lipotes vexillifer"
    LutjanusSanguineus, "Humphead snapper", "Lutjanus sanguineus (humphead snapper)", "Lutjanus sanguineus"
    MacacaArctoides, "Stump-tailed macaque", "Macaca arctoides (stump-tailed macaque)", "Macaca arctoides"
    MacacaAssamensis, "Assam macaque", "Macaca assamensis (Assam macaque)", "Macaca assamensis"
    MacacaCyclopis, "Taiwan macaque", "Macaca cyclopis (Taiwan macaque)", "Macaca cyclopis"
    MacacaFascicularis, "Crab-eating macaque", "Macaca fascicularis (crab-eating macaque)", "Macaca fascicularis"
    MacacaMulatta, "Rhesus monkey", "Macaca mulatta (Rhesus monkey)", "Macaca mulatta"
    MacacaNemestrina, "Pig-tailed macaque", "Macaca nemestrina (pig-tailed macaque)", "Macaca nemestrina"
    MacacaSilenus, "Liontail macaque", "Macaca silenus (liontail macaque)", "Macaca silenus"
    MacacaThibetana, "Pere David's macaque", "Macaca thibetana (Pere David's macaque)", "Macaca thibetana"
    MarecaStrepera, "Gadwall", "Mareca strepera (gadwall)", "Mareca strepera"
    MarmotaHimalayana, "Himalayan marmot", "Marmota himalayana (Himalayan marmot)", "Marmota himalayana"
    MarmotaMonax, "Woodchuck", "Marmota monax (woodchuck)", "Marmota monax"
    MauremysMutica, "Yellowpond turtle", "Mauremys mutica (yellowpond turtle)", "Mauremys mutica"
    MelanogrammusAeglefinus, "Haddock", "Melanogrammus aeglefinus (haddock)", "Melanogrammus aeglefinus"
    MeleagrisGallopavo, "Turkey", "Meleagris gallopavo (turkey)", "Meleagris gallopavo"
    MerionesUnguiculatus, "Mongolian gerbil", "Meriones unguiculatus (Mongolian gerbil)", "Meriones unguiculatus"
    MesocricetusAuratus, "Golden hamster", "Mesocricetus auratus (golden hamster)", "Mesocricetus auratus"
    MicrocebusMurinus, "Gray mouse lemur", "Microcebus murinus (gray mouse lemur)", "Microcebus murinus"
    MonodelphisDomestica, "Gray short-tailed opossum", "Monodelphis domestica (gray short-tailed opossum)", "Monodelphis domestica"
    MurinaeSp, "Old world rats and mice", "Murinae gen. sp.", "Murinae sp."
    Mus, "mouse", "Mus (mouse)", "Mus"
    MusCookii, "Cook's mouse", "Mus cookii (Cook's mouse)", "Mus cookii"
    MusMinutoides, "Southern African pygmy mouse", "Mus minutoides (Southern African pygmy mouse)", "Mus minutoides"
    MusMusculus, "House mouse", "Mus musculus (house mouse)", "Mus musculus"
    MusMusculusCastaneus, "Southeastern Asian house mouse", "Mus musculus castaneus (southeastern Asian house mouse)", "Mus musculus castaneus"
    MusMusculusDomesticus, "Western European house mouse", "Mus musculus domesticus (western European house mouse)", "Mus musculus domesticus"
    MusMusculusMolossinus, "Japanese wild mouse", "Mus musculus molossinus (Japanese wild mouse)", "Mus musculus molossinus"
    MusMusculusMusculus, "Eastern European house mouse", "Mus musculus musculus (eastern European house mouse)", "Mus musculus musculus"
    MusPahari, "Shrew mouse", "Mus pahari (shrew mouse)", "Mus pahari"
    MusSaxicola, "Spiny mouse", "Mus saxicola (spiny mouse)", "Mus saxicola"
    MusSp, "Mice", "Mus sp. (mice)", "Mus sp."
    MusSpretus, "Western wild mouse", "Mus spretus (western wild mouse)", "Mus spretus"
    MustelaPutoriusFuro, "Domestic ferret", "Mustela putorius furo (domestic ferret)", "Mustela putorius furo"
    MustelaSp, "Ferret", "Mustela sp.", "Mustela sp."
    MyotisLucifugus,"Little brown bat", "Myotis lucifugus (little brown bat)", "Myotis lucifugus"
    NeophocaenaPhocaenoides, "Indo-Pacific finless porpoise", "Neophocaena phocaenoides (Indo-Pacific finless porpoise)", "Neophocaena phocaenoides"
    NeogaleVison, "American mink", "Neogale vison (American mink)", "Neogale vison"
    NomascusConcolor, "Black crested gibbon", "Nomascus concolor (Black crested gibbon)", "Nomascus concolor"
    NotamacropusEugenii, "Tammar wallaby", "Notamacropus eugenii (tammar wallaby)", "Notamacropus eugenii"
    NothocricetulusMigratorius, "Armenian hamster", "Nothocricetulus migratorius (Armenian hamster)", "Nothocricetulus migratorius"
    NototheniaCoriiceps, "Black rockcod", "Notothenia coriiceps (black rockcod)", "Notothenia coriiceps"
    NycticebusCoucang, "Slow loris", "Nycticebus coucang (slow loris)", "Nycticebus coucang"
    OncorhynchusGorbuscha, "Pink salmon", "Oncorhynchus gorbuscha (pink salmon)", "Oncorhynchus gorbuscha"
    OncorhynchusMykiss, "Rainbow trout", "Oncorhynchus mykiss (rainbow trout)", "Oncorhynchus mykiss"
    OncorhynchusTshawytscha, "Chinook salmon", "Oncorhynchus tshawytscha (Chinook salmon)", "Oncorhynchus tshawytscha"
    OrectolobusMaculatus, "Spotted wobbegong", "Orectolobus maculatus (spotted wobbegong)", "Orectolobus maculatus"
    OreochromisNiloticus, "Nile tilapia", "Oreochromis niloticus (Nile tilapia)", "Oreochromis niloticus"
    OrnithorhynchusAnatinus, "Platypus", "Ornithorhynchus anatinus (platypus)", "Ornithorhynchus anatinus"
    OryctolagusCuniculus, "Rabbit", "Oryctolagus cuniculus (rabbit)", "Oryctolagus cuniculus"
    OryctolagusCuniculusAlgirus, "European rabbit", "Oryctolagus cuniculus algirus", "Oryctolagus cuniculus algirus"
    OryctolagusCuniculusCuniculus, "Rabbit", "Oryctolagus cuniculus cuniculus", "Oryctolagus cuniculus cuniculus"
    OryziasLatipes, "Japanese medaka", "Oryzias latipes (Japanese medaka)", "Oryzias latipes"
    OryziasMelastigma, "Indian medaka", "Oryzias melastigma (Indian medaka)", "Oryzias melastigma"
    OtolemurCrassicaudatus, "Thick-tailed bush baby", "Otolemur crassicaudatus (thick-tailed bush baby)", "Otolemur crassicaudatus"
    OvisAries, "Domestic sheep", "Ovis aries (sheep)", "Ovis aries"
    OvisSp, "Sheep", "Ovis sp.", "Ovis sp."
    PacifastacusLeniusculus, "Signal crayfish", "Pacifastacus leniusculus (signal crayfish)", "Pacifastacus leniusculus"
    PagetopsisMacropterus, "Pagetopsis macropterus", "Pagetopsis macropterus", "Pagetopsis macropterus"
    PagrusMajor, "Red seabream", "Pagrus major (red seabream)", "Pagrus major"
    PangasianodonHypophthalmus, "Striped catfish", "Pangasianodon hypophthalmus (striped catfish)", "Pangasianodon hypophthalmus"
    PanPaniscus, "Pygmy chimpanzee", "Pan paniscus (pygmy chimpanzee)", "Pan paniscus"
    PantheraPardus, "Leopard", "Panthera pardus (leopard)", "Panthera pardus"
    PanTroglodytes, "Chimpanzee", "Pan troglodytes (chimpanzee)", "Pan troglodytes"
    PanTroglodytesVerus, "Western chimpanzee", "Pan troglodytes verus", "Pan troglodytes verus"
    PapioAnubis, "Olive baboon", "Papio anubis (olive baboon)", "Papio anubis"
    PapioAnubisAnubis, "Olive baboon anubis", "Papio anubis anubis", "Papio anubis anubis"
    PapioHamadryas, "Hamadryas baboon", "Papio hamadryas (hamadryas baboon)", "Papio hamadryas"
    PapioPapio, "Guinea baboon", "Papio papio (Guinea baboon)", "Papio papio"
    ParalichthysOlivaceus, "Japanese flounder", "Paralichthys olivaceus (Japanese flounder)", "Paralichthys olivaceus"
    PelodiscusSinensis, "Chinese soft-shelled turtle", "Pelodiscus sinensis (Chinese soft-shelled turtle)", "Pelodiscus sinensis"
    PelteobagrusFulvidraco, "Yellow catfish", "Pelteobagrus fulvidraco (yellow catfish)", "Pelteobagrus fulvidraco"
    PerdixPerdix, "Grey partridge", "Perdix perdix (grey partridge)", "Perdix perdix"
    PeromyscusManiculatus, "North American deer mouse", "Peromyscus maniculatus (North American deer mouse)", "Peromyscus maniculatus"
    PetromyzonMarinus, "Sea lamprey", "Petromyzon marinus (sea lamprey)", "Petromyzon marinus"
    PhascogaleCalura, "Red-tailed phascogale", "Phascogale calura (red-tailed phascogale)", "Phascogale calura"
    PhasianusColchicus, "Ring-necked pheasant", "Phasianus colchicus (Ring-necked pheasant)", "Phasianus colchicus"
    PhyseterCatodon, "Sperm whale", "Physeter catodon (sperm whale)", "Physeter catodon"
    PitheciaPithecia, "White-faced saki", "Pithecia pithecia (white-faced saki)", "Pithecia pithecia"
    PlataleaAjaja, "Roseate spoonbil", "Platalea ajaja", "Platalea ajaja"
    Platyrrhini, "New World monkeys", "Platyrrhini (New World monkeys)", "Platyrrhini"
    PlecoglossusAltivelisAltivelis, "Ayu sweetfish", "Plecoglossus altivelis altivelis", "Plecoglossus altivelis altivelis"
    PleurodelesWaltl, "Iberian ribbed newt", "Pleurodeles waltl (Iberian ribbed newt)", "Pleurodeles waltl"
    PogonophryneScotti, "Pogonophryne scotti", "Pogonophryne scotti", "Pogonophryne scotti"
    PolyprionOxygeneios, "Hāpuku", "Polyprion oxygeneios", "Polyprion oxygeneios"
    PongoAbelii, "Sumatran orangutan", "Pongo abelii (Sumatran orangutan)", "Pongo abelii"
    PongoPygmaeus, "Bornean orangutan", "Pongo pygmaeus (Bornean orangutan)", "Pongo pygmaeus"
    PresbytisComata, "Grizzled leaf monkey", "Presbytis comata (grizzled leaf monkey)", "Presbytis comata"
    PresbytisFemoralis, "Banded leaf monkey", "Presbytis femoralis (banded leaf monkey)", "Presbytis femoralis"
    PresbytisMelalophos, "Mitred leaf monkey", "Presbytis melalophos (mitred leaf monkey)", "Presbytis melalophos"
    PropithecusVerreauxi, "White sifaka", "Propithecus verreauxi (white sifaka)", "Propithecus verreauxi"
    ProtopterusAethiopicus, "Marbled lungfish", "Protopterus aethiopicus (marbled lungfish)", "Protopterus aethiopicus"
    PseudobatosProductus, "Shovelnose guitarfish", "Pseudobatos productus (shovelnose guitarfish)", "Pseudobatos productus"
    PteropusAlecto, "Black flying fox", "Pteropus alecto (black flying fox)", "Pteropus alecto"
    PythonBivittatus, "Burmese python", "Python bivittatus (Burmese python)", "Python bivittatus"
    RachycentronCanadum, "Cobia", "Rachycentron canadum (cobia)", "Rachycentron canadum"
    RajaEglanteria, "Clearnose skate", "Raja eglanteria (clearnose skate)", "Raja eglanteria"
    RattusFuscipes, "Bush rat", "Rattus fuscipes (bush rat)", "Rattus fuscipes"
    RattusLeucopus, "Mottle-tailed rat", "Rattus leucopus (mottle-tailed rat)", "Rattus leucopus"
    RattusNorvegicus, "Norway rat", "Rattus norvegicus (Norway rat)", "Rattus norvegicus"
    RattusRattus, "Black rat", "Rattus rattus (black rat)", "Rattus rattus"
    RattusSordidus, "Australian dusky field rat", "Rattus sordidus (Australian dusky field rat)", "Rattus sordidus"
    RattusSp, "Rats", "Rattus sp. (rats)", "Rattus sp."
    RattusTunneyi, "Tunney's rat", "Rattus tunneyi (Tunney's rat)", "Rattus tunneyi"
    RattusVillosissimus, "Long-haired rat", "Rattus villosissimus (long-haired rat)", "Rattus villosissimus"
    RhinocerosUnicornis, "Greater Indian rhinoceros", "Rhinoceros unicornis (greater Indian rhinoceros)", "Rhinoceros unicornis"
    RousettusLeschenaultii, "Leschenault's rousette", "Rousettus leschenaultii (Leschenault's rousette)", "Rousettus leschenaultii"
    SaccharomycesCerevisiae, "baker's yeast", "Saccharomyces cerevisiae (baker's yeast)", "Saccharomyces cerevisiae"
    SaguinusLabiatus, "Red-chested mustached tamarin", "Saguinus labiatus (red-chested mustached tamarin)", "Saguinus labiatus"
    SaguinusMidas, "Midas tamarin", "Saguinus midas (Midas tamarin)", "Saguinus midas"
    SaguinusOedipus, "Cotton-top tamarin", "Saguinus oedipus (cotton-top tamarin)", "Saguinus oedipus"
    SaimiriBoliviensisBoliviensis, "Bolivian squirrel monkey", "Saimiri boliviensis boliviensis (Bolivian squirrel monkey)", "Saimiri boliviensis boliviensis"
    SaimiriSciureus, "Common squirrel monkey", "Saimiri sciureus (common squirrel monkey)", "Saimiri sciureus"
    SalmoMarmoratus, "Salmo marmoratus", "Salmo marmoratus", "Salmo marmoratus"
    SalmoSalar, "Atlantic salmon", "Salmo salar (Atlantic salmon)", "Salmo salar"
    SalmoTrutta, "River trout", "Salmo trutta (river trout)", "Salmo trutta"
    SalvelinusAlpinus, "Arctic char", "Salvelinus alpinus (Arctic char)", "Salvelinus alpinus"
    SanderVitreus, "Walleye", "Sander vitreus (walleye)", "Sander vitreus"
    SapajusApella, "Tufted capuchin", "Sapajus apella (Tufted capuchin)", "Sapajus apella"
    SchizophyllumCommune, "Schizophyllum commune", "Schizophyllum commune", "Schizophyllum commune"
    SciaenopsOcellatus, "Red drum", "Sciaenops ocellatus (red drum)", "Sciaenops ocellatus"
    ScophthalmusMaximus, "Turbot", "Scophthalmus maximus (turbot)", "Scophthalmus maximus"
    ScyliorhinusCanicula, "Smaller spotted catshark", "Scyliorhinus canicula (smaller spotted catshark)", "Scyliorhinus canicula"
    SeriolaQuinqueradiata, "Japanese amberjack", "Seriola quinqueradiata (Japanese amberjack)", "Seriola quinqueradiata"
    SilurusAsotus, "Amur catfish", "Silurus asotus (Amur catfish)", "Silurus asotus"
    SilurusMeridionalis, "Silurus meridionalis", "Silurus meridionalis", "Silurus meridionalis"
    SinipercaChuatsi, "Mandarin fish", "Siniperca chuatsi (mandarin fish)", "Siniperca chuatsi"
    SousaChinensis, "Indo-pacific humpbacked dolphin", "Sousa chinensis (Indo-pacific humpbacked dolphin)", "Sousa chinensis"
    SparusAurata, "Gilthead seabream", "Sparus aurata (gilthead seabream)", "Sparus aurata"
    SphoeroidesNephelus, "Southern puffer", "Sphoeroides nephelus (southern puffer)", "Sphoeroides nephelus"
    SqualusAcanthias, "Spiny dogfish", "Squalus acanthias (spiny dogfish)", "Squalus acanthias"
    StegastesLeucostictus, "Beaugregory", "Stegastes leucostictus (beaugregory)", "Stegastes leucostictus"
    StegastesPartitus, "Bicolor damselfish", "Stegastes partitus (bicolor damselfish)", "Stegastes partitus"
    StenellaAttenuata, "Bridled dolphin", "Stenella attenuata (bridled dolphin)", "Stenella attenuata"
    StenellaCoeruleoalba, "Striped dolphin", "Stenella coeruleoalba (striped dolphin)", "Stenella coeruleoalba"
    StreptomycesAvidinii, "Streptomyces avidinii", "Streptomyces avidinii", "Streptomyces avidinii"
    StreptomycesViridochromogenes, "Streptomyces viridochromogenes", "Streptomyces viridochromogenes", "Streptomyces viridochromogenes"
    StruthioCamelus, "African ostrich", "Struthio camelus (African ostrich)", "Struthio camelus"
    SuncusMurinus, "House shrew", "Suncus murinus (house shrew)", "Suncus murinus"
    SusScrofa, "Domestic pig", "Sus scrofa (pig)", "Sus scrofa"
    SylvilagusCunicularis, "Mexican cottontail", "Sylvilagus cunicularis (Mexican cottontail)", "Sylvilagus cunicularis"
    SylvilagusFloridanus, "Eastern cottontail", "Sylvilagus floridanus (eastern cottontail)", "Sylvilagus floridanus"
    SymphalangusSyndactylus, "Siamang", "Symphalangus syndactylus (siamang)", "Symphalangus syndactylus"
    TachyglossusAculeatus, "Australian echidna", "Tachyglossus aculeatus (Australian echidna)", "Tachyglossus aculeatus"
    TachysurusFulvidraco, "Yellow catfish", "Tachysurus fulvidraco (yellow catfish)", "Tachysurus fulvidraco"
    TachysurusVachellii, "Tachysurus vachellii", "Tachysurus vachellii", "Tachysurus vachellii"
    TaeniopygiaGuttata, "Zebra finch", "Taeniopygia guttata (zebra finch)", "Taeniopygia guttata"
    TakifuguRubripes, "Torafugu", "Takifugu rubripes (torafugu)", "Takifugu rubripes"
    TarsiusDentatus, "Dian's tarsier", "Tarsius dentatus (Dian's tarsier)", "Tarsius dentatus"
    TarsiusLariang, "Lariang tarsier", "Tarsius lariang (Lariang tarsier)", "Tarsius lariang"
    TarsiusSyrichta, "Philippine tarsier", "Tarsius syrichta (Philippine tarsier)", "Tarsius syrichta"
    TetraodonNigroviridis, "Spotted green pufferfish", "Tetraodon nigroviridis (spotted green pufferfish)", "Tetraodon nigroviridis"
    TrachemysScripta, "Red-eared slider turtle", "Trachemys scripta (red-eared slider turtle)", "Trachemys scripta"
    TrachemysScriptaElegans, "Red-eared slider turtle elegans", "Trachemys scripta elegans", "Trachemys scripta elegans"
    TrachypithecusCristatus, "Silvery lutung", "Trachypithecus cristatus (Silvery lutung)", "Trachypithecus cristatus"
    TrachypithecusObscurus, "Dusky leaf-monkey", "Trachypithecus obscurus (Dusky leaf-monkey)", "Trachypithecus obscurus"
    TrematomusBernacchii, "Emerald rockcod", "Trematomus bernacchii (emerald rockcod)", "Trematomus bernacchii"
    TrematomusHansoni, "Striped rockcod", "Trematomus hansoni (striped rockcod)", "Trematomus hansoni"
    TrematomusLoennbergii, "Deepwater notothen", "Trematomus loennbergii (deepwater notothen)", "Trematomus loennbergii"
    TrematomusNewnesi, "Dusky notothen", "Trematomus newnesi (dusky notothen)", "Trematomus newnesi"
    TrematomusPennellii, "Sharp-spined notothen", "Trematomus pennellii (sharp-spined notothen)", "Trematomus pennellii"
    TriakisScyllium, "Banded houndshark", "Triakis scyllium (banded houndshark)", "Triakis scyllium"
    TrichechusManatusLatirostris, "Florida manatee", "Trichechus manatus latirostris (Florida manatee)", "Trichechus manatus latirostris"
    TrichosurusVulpecula, "Common brushtail", "Trichosurus vulpecula (common brushtail)", "Trichosurus vulpecula"
    TursiopsAduncus, "Indo-pacific bottlenose dolphin", "Tursiops aduncus (Indo-pacific bottlenose dolphin)", "Tursiops aduncus"
    TursiopsTruncatus, "Common bottlenose dolphin", "Tursiops truncatus (common bottlenose dolphin)", "Tursiops truncatus"
    VicugnaPacos, "Alpaca", "Vicugna pacos (alpaca)", "Vicugna pacos"
    Xenopus, "Xenopus", "Xenopus", "Xenopus"
    XenopusLaevis, "African clawed frog", "Xenopus laevis (African clawed frog)", "Xenopus laevis"
    XenopusLaevisOrGilli, "African or Cape clawed frog", "Xenopus laevis/gilli", "Xenopus laevis/gilli"
    XenopusSp, "Clawed frog", "Xenopus sp. (clawed frog)", "Xenopus sp."
    XenopusTropicalis, "Tropical clawed frog", "Xenopus tropicalis (tropical clawed frog)", "Xenopus tropicalis"
);
