// @generated
#![allow(non_snake_case,non_upper_case_globals)]
use std::sync::LazyLock;
use bincode::config::Configuration;
use super::{Germlines, Species};
/// Get the germlines for any of the available species. See the main documentation for which species have which data available.
pub(super) fn germlines(species: Species) -> Option<&'static Germlines> {match species {
Species::AnarhichasMinor => Some(&ANARHICHASMINOR),
Species::BosTaurus => Some(&BOSTAURUS),
Species::CamelusDromedarius => Some(&CAMELUSDROMEDARIUS),
Species::CanisLupusFamiliaris => Some(&CANISLUPUSFAMILIARIS),
Species::CapraHircus => Some(&CAPRAHIRCUS),
Species::CarcharhinusPlumbeus => Some(&CARCHARHINUSPLUMBEUS),
Species::CercocebusAtys => Some(&CERCOCEBUSATYS),
Species::ChaenocephalusAceratus => Some(&CHAENOCEPHALUSACERATUS),
Species::CyprinusCarpio => Some(&CYPRINUSCARPIO),
Species::DanioRerio => Some(&DANIORERIO),
Species::DicentrarchusLabrax => Some(&DICENTRARCHUSLABRAX),
Species::EquusCaballus => Some(&EQUUSCABALLUS),
Species::FelisCatus => Some(&FELISCATUS),
Species::GadusMorhua => Some(&GADUSMORHUA),
Species::GallusGallus => Some(&GALLUSGALLUS),
Species::GasterosteusAculeatus => Some(&GASTEROSTEUSACULEATUS),
Species::GinglymostomaCirratum => Some(&GINGLYMOSTOMACIRRATUM),
Species::GorillaGorilla => Some(&GORILLAGORILLA),
Species::GorillaGorillaGorilla => Some(&GORILLAGORILLAGORILLA),
Species::HeterodontusFrancisci => Some(&HETERODONTUSFRANCISCI),
Species::HomoSapiens => Some(&HOMOSAPIENS),
Species::HydrolagusColliei => Some(&HYDROLAGUSCOLLIEI),
Species::HylobatesLar => Some(&HYLOBATESLAR),
Species::IctalurusPunctatus => Some(&ICTALURUSPUNCTATUS),
Species::LemurCatta => Some(&LEMURCATTA),
Species::LeucorajaErinacea => Some(&LEUCORAJAERINACEA),
Species::MacacaArctoides => Some(&MACACAARCTOIDES),
Species::MacacaCyclopis => Some(&MACACACYCLOPIS),
Species::MacacaFascicularis => Some(&MACACAFASCICULARIS),
Species::MacacaMulatta => Some(&MACACAMULATTA),
Species::MacacaNemestrina => Some(&MACACANEMESTRINA),
Species::MacacaSilenus => Some(&MACACASILENUS),
Species::MacacaThibetana => Some(&MACACATHIBETANA),
Species::MesocricetusAuratus => Some(&MESOCRICETUSAURATUS),
Species::MonodelphisDomestica => Some(&MONODELPHISDOMESTICA),
Species::MusCookii => Some(&MUSCOOKII),
Species::MusMinutoides => Some(&MUSMINUTOIDES),
Species::MusMusculus => Some(&MUSMUSCULUS),
Species::MusMusculusCastaneus => Some(&MUSMUSCULUSCASTANEUS),
Species::MusMusculusDomesticus => Some(&MUSMUSCULUSDOMESTICUS),
Species::MusMusculusMolossinus => Some(&MUSMUSCULUSMOLOSSINUS),
Species::MusMusculusMusculus => Some(&MUSMUSCULUSMUSCULUS),
Species::MusPahari => Some(&MUSPAHARI),
Species::MusSaxicola => Some(&MUSSAXICOLA),
Species::MusSp => Some(&MUSSP),
Species::MusSpretus => Some(&MUSSPRETUS),
Species::MustelaPutoriusFuro => Some(&MUSTELAPUTORIUSFURO),
Species::NeogaleVison => Some(&NEOGALEVISON),
Species::NototheniaCoriiceps => Some(&NOTOTHENIACORIICEPS),
Species::OncorhynchusMykiss => Some(&ONCORHYNCHUSMYKISS),
Species::OrnithorhynchusAnatinus => Some(&ORNITHORHYNCHUSANATINUS),
Species::OryctolagusCuniculus => Some(&ORYCTOLAGUSCUNICULUS),
Species::OryctolagusCuniculusAlgirus => Some(&ORYCTOLAGUSCUNICULUSALGIRUS),
Species::OryctolagusCuniculusCuniculus => Some(&ORYCTOLAGUSCUNICULUSCUNICULUS),
Species::OvisAries => Some(&OVISARIES),
Species::PanTroglodytes => Some(&PANTROGLODYTES),
Species::PapioAnubisAnubis => Some(&PAPIOANUBISANUBIS),
Species::PongoAbelii => Some(&PONGOABELII),
Species::PongoPygmaeus => Some(&PONGOPYGMAEUS),
Species::ProtopterusAethiopicus => Some(&PROTOPTERUSAETHIOPICUS),
Species::RajaEglanteria => Some(&RAJAEGLANTERIA),
Species::RattusNorvegicus => Some(&RATTUSNORVEGICUS),
Species::RattusRattus => Some(&RATTUSRATTUS),
Species::SalmoSalar => Some(&SALMOSALAR),
Species::SalmoTrutta => Some(&SALMOTRUTTA),
Species::SeriolaQuinqueradiata => Some(&SERIOLAQUINQUERADIATA),
Species::SinipercaChuatsi => Some(&SINIPERCACHUATSI),
Species::SusScrofa => Some(&SUSSCROFA),
Species::TrematomusBernacchii => Some(&TREMATOMUSBERNACCHII),
Species::VicugnaPacos => Some(&VICUGNAPACOS),
Species::XenopusLaevisOrGilli => Some(&XENOPUSLAEVISORGILLI),
_=>None}}
/// Get all germlines in one iterator, see the main documentation for more information about the available germlines
pub(super) fn all_germlines() -> impl Iterator<Item = &'static Germlines> {
[
&*ANARHICHASMINOR,
&*BOSTAURUS,
&*CAMELUSDROMEDARIUS,
&*CANISLUPUSFAMILIARIS,
&*CAPRAHIRCUS,
&*CARCHARHINUSPLUMBEUS,
&*CERCOCEBUSATYS,
&*CHAENOCEPHALUSACERATUS,
&*CYPRINUSCARPIO,
&*DANIORERIO,
&*DICENTRARCHUSLABRAX,
&*EQUUSCABALLUS,
&*FELISCATUS,
&*GADUSMORHUA,
&*GALLUSGALLUS,
&*GASTEROSTEUSACULEATUS,
&*GINGLYMOSTOMACIRRATUM,
&*GORILLAGORILLA,
&*GORILLAGORILLAGORILLA,
&*HETERODONTUSFRANCISCI,
&*HOMOSAPIENS,
&*HYDROLAGUSCOLLIEI,
&*HYLOBATESLAR,
&*ICTALURUSPUNCTATUS,
&*LEMURCATTA,
&*LEUCORAJAERINACEA,
&*MACACAARCTOIDES,
&*MACACACYCLOPIS,
&*MACACAFASCICULARIS,
&*MACACAMULATTA,
&*MACACANEMESTRINA,
&*MACACASILENUS,
&*MACACATHIBETANA,
&*MESOCRICETUSAURATUS,
&*MONODELPHISDOMESTICA,
&*MUSCOOKII,
&*MUSMINUTOIDES,
&*MUSMUSCULUS,
&*MUSMUSCULUSCASTANEUS,
&*MUSMUSCULUSDOMESTICUS,
&*MUSMUSCULUSMOLOSSINUS,
&*MUSMUSCULUSMUSCULUS,
&*MUSPAHARI,
&*MUSSAXICOLA,
&*MUSSP,
&*MUSSPRETUS,
&*MUSTELAPUTORIUSFURO,
&*NEOGALEVISON,
&*NOTOTHENIACORIICEPS,
&*ONCORHYNCHUSMYKISS,
&*ORNITHORHYNCHUSANATINUS,
&*ORYCTOLAGUSCUNICULUS,
&*ORYCTOLAGUSCUNICULUSALGIRUS,
&*ORYCTOLAGUSCUNICULUSCUNICULUS,
&*OVISARIES,
&*PANTROGLODYTES,
&*PAPIOANUBISANUBIS,
&*PONGOABELII,
&*PONGOPYGMAEUS,
&*PROTOPTERUSAETHIOPICUS,
&*RAJAEGLANTERIA,
&*RATTUSNORVEGICUS,
&*RATTUSRATTUS,
&*SALMOSALAR,
&*SALMOTRUTTA,
&*SERIOLAQUINQUERADIATA,
&*SINIPERCACHUATSI,
&*SUSSCROFA,
&*TREMATOMUSBERNACCHII,
&*VICUGNAPACOS,
&*XENOPUSLAEVISORGILLI,
].into_iter()
}
/// Get all germlines in one parallel iterator, see the main documentation for more information about the available germlines
#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
pub(super) fn par_germlines() -> impl ParallelIterator<Item = &'static Germlines> {
[
&*ANARHICHASMINOR,
&*BOSTAURUS,
&*CAMELUSDROMEDARIUS,
&*CANISLUPUSFAMILIARIS,
&*CAPRAHIRCUS,
&*CARCHARHINUSPLUMBEUS,
&*CERCOCEBUSATYS,
&*CHAENOCEPHALUSACERATUS,
&*CYPRINUSCARPIO,
&*DANIORERIO,
&*DICENTRARCHUSLABRAX,
&*EQUUSCABALLUS,
&*FELISCATUS,
&*GADUSMORHUA,
&*GALLUSGALLUS,
&*GASTEROSTEUSACULEATUS,
&*GINGLYMOSTOMACIRRATUM,
&*GORILLAGORILLA,
&*GORILLAGORILLAGORILLA,
&*HETERODONTUSFRANCISCI,
&*HOMOSAPIENS,
&*HYDROLAGUSCOLLIEI,
&*HYLOBATESLAR,
&*ICTALURUSPUNCTATUS,
&*LEMURCATTA,
&*LEUCORAJAERINACEA,
&*MACACAARCTOIDES,
&*MACACACYCLOPIS,
&*MACACAFASCICULARIS,
&*MACACAMULATTA,
&*MACACANEMESTRINA,
&*MACACASILENUS,
&*MACACATHIBETANA,
&*MESOCRICETUSAURATUS,
&*MONODELPHISDOMESTICA,
&*MUSCOOKII,
&*MUSMINUTOIDES,
&*MUSMUSCULUS,
&*MUSMUSCULUSCASTANEUS,
&*MUSMUSCULUSDOMESTICUS,
&*MUSMUSCULUSMOLOSSINUS,
&*MUSMUSCULUSMUSCULUS,
&*MUSPAHARI,
&*MUSSAXICOLA,
&*MUSSP,
&*MUSSPRETUS,
&*MUSTELAPUTORIUSFURO,
&*NEOGALEVISON,
&*NOTOTHENIACORIICEPS,
&*ONCORHYNCHUSMYKISS,
&*ORNITHORHYNCHUSANATINUS,
&*ORYCTOLAGUSCUNICULUS,
&*ORYCTOLAGUSCUNICULUSALGIRUS,
&*ORYCTOLAGUSCUNICULUSCUNICULUS,
&*OVISARIES,
&*PANTROGLODYTES,
&*PAPIOANUBISANUBIS,
&*PONGOABELII,
&*PONGOPYGMAEUS,
&*PROTOPTERUSAETHIOPICUS,
&*RAJAEGLANTERIA,
&*RATTUSNORVEGICUS,
&*RATTUSRATTUS,
&*SALMOSALAR,
&*SALMOTRUTTA,
&*SERIOLAQUINQUERADIATA,
&*SINIPERCACHUATSI,
&*SUSSCROFA,
&*TREMATOMUSBERNACCHII,
&*VICUGNAPACOS,
&*XENOPUSLAEVISORGILLI,
].into_par_iter()
}
static ANARHICHASMINOR: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Spotted wolffish.bin"),Configuration::default(),).unwrap().0});
static BOSTAURUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic bovine.bin"),Configuration::default(),).unwrap().0});
static CAMELUSDROMEDARIUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Arabian camel.bin"),Configuration::default(),).unwrap().0});
static CANISLUPUSFAMILIARIS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic dog.bin"),Configuration::default(),).unwrap().0});
static CAPRAHIRCUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic goat.bin"),Configuration::default(),).unwrap().0});
static CARCHARHINUSPLUMBEUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Sandbar shark.bin"),Configuration::default(),).unwrap().0});
static CERCOCEBUSATYS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Sooty mangabey.bin"),Configuration::default(),).unwrap().0});
static CHAENOCEPHALUSACERATUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Blackfin icefish.bin"),Configuration::default(),).unwrap().0});
static CYPRINUSCARPIO: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Common carp.bin"),Configuration::default(),).unwrap().0});
static DANIORERIO: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Zebrafish.bin"),Configuration::default(),).unwrap().0});
static DICENTRARCHUSLABRAX: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("European seabass.bin"),Configuration::default(),).unwrap().0});
static EQUUSCABALLUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic horse.bin"),Configuration::default(),).unwrap().0});
static FELISCATUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic cat.bin"),Configuration::default(),).unwrap().0});
static GADUSMORHUA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Atlantic cod.bin"),Configuration::default(),).unwrap().0});
static GALLUSGALLUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic chicken.bin"),Configuration::default(),).unwrap().0});
static GASTEROSTEUSACULEATUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Three-spined stickleback.bin"),Configuration::default(),).unwrap().0});
static GINGLYMOSTOMACIRRATUM: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Nurse shark.bin"),Configuration::default(),).unwrap().0});
static GORILLAGORILLA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Western gorilla.bin"),Configuration::default(),).unwrap().0});
static GORILLAGORILLAGORILLA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Western lowland gorilla.bin"),Configuration::default(),).unwrap().0});
static HETERODONTUSFRANCISCI: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Horn shark.bin"),Configuration::default(),).unwrap().0});
static HOMOSAPIENS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Human.bin"),Configuration::default(),).unwrap().0});
static HYDROLAGUSCOLLIEI: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Spotted ratfish.bin"),Configuration::default(),).unwrap().0});
static HYLOBATESLAR: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Common gibbon.bin"),Configuration::default(),).unwrap().0});
static ICTALURUSPUNCTATUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Channel catfish.bin"),Configuration::default(),).unwrap().0});
static LEMURCATTA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Ring-tailed lemur.bin"),Configuration::default(),).unwrap().0});
static LEUCORAJAERINACEA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Little skate.bin"),Configuration::default(),).unwrap().0});
static MACACAARCTOIDES: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Stump-tailed macaque.bin"),Configuration::default(),).unwrap().0});
static MACACACYCLOPIS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Taiwan macaque.bin"),Configuration::default(),).unwrap().0});
static MACACAFASCICULARIS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Crab-eating macaque.bin"),Configuration::default(),).unwrap().0});
static MACACAMULATTA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Rhesus monkey.bin"),Configuration::default(),).unwrap().0});
static MACACANEMESTRINA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Pig-tailed macaque.bin"),Configuration::default(),).unwrap().0});
static MACACASILENUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Liontail macaque.bin"),Configuration::default(),).unwrap().0});
static MACACATHIBETANA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Pere David's macaque.bin"),Configuration::default(),).unwrap().0});
static MESOCRICETUSAURATUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Golden hamster.bin"),Configuration::default(),).unwrap().0});
static MONODELPHISDOMESTICA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Gray short-tailed opossum.bin"),Configuration::default(),).unwrap().0});
static MUSCOOKII: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Cook's mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMINUTOIDES: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Southern African pygmy mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMUSCULUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("House mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMUSCULUSCASTANEUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Southeastern Asian house mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMUSCULUSDOMESTICUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Western European house mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMUSCULUSMOLOSSINUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Japanese wild mouse.bin"),Configuration::default(),).unwrap().0});
static MUSMUSCULUSMUSCULUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Eastern European house mouse.bin"),Configuration::default(),).unwrap().0});
static MUSPAHARI: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Shrew mouse.bin"),Configuration::default(),).unwrap().0});
static MUSSAXICOLA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Spiny mouse.bin"),Configuration::default(),).unwrap().0});
static MUSSP: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Mice.bin"),Configuration::default(),).unwrap().0});
static MUSSPRETUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Western wild mouse.bin"),Configuration::default(),).unwrap().0});
static MUSTELAPUTORIUSFURO: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic ferret.bin"),Configuration::default(),).unwrap().0});
static NEOGALEVISON: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("American mink.bin"),Configuration::default(),).unwrap().0});
static NOTOTHENIACORIICEPS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Black rockcod.bin"),Configuration::default(),).unwrap().0});
static ONCORHYNCHUSMYKISS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Rainbow trout.bin"),Configuration::default(),).unwrap().0});
static ORNITHORHYNCHUSANATINUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Platypus.bin"),Configuration::default(),).unwrap().0});
static ORYCTOLAGUSCUNICULUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Rabbit.bin"),Configuration::default(),).unwrap().0});
static ORYCTOLAGUSCUNICULUSALGIRUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("European rabbit.bin"),Configuration::default(),).unwrap().0});
static ORYCTOLAGUSCUNICULUSCUNICULUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Rabbit.bin"),Configuration::default(),).unwrap().0});
static OVISARIES: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic sheep.bin"),Configuration::default(),).unwrap().0});
static PANTROGLODYTES: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Chimpanzee.bin"),Configuration::default(),).unwrap().0});
static PAPIOANUBISANUBIS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Olive baboon anubis.bin"),Configuration::default(),).unwrap().0});
static PONGOABELII: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Sumatran orangutan.bin"),Configuration::default(),).unwrap().0});
static PONGOPYGMAEUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Bornean orangutan.bin"),Configuration::default(),).unwrap().0});
static PROTOPTERUSAETHIOPICUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Marbled lungfish.bin"),Configuration::default(),).unwrap().0});
static RAJAEGLANTERIA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Clearnose skate.bin"),Configuration::default(),).unwrap().0});
static RATTUSNORVEGICUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Norway rat.bin"),Configuration::default(),).unwrap().0});
static RATTUSRATTUS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Black rat.bin"),Configuration::default(),).unwrap().0});
static SALMOSALAR: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Atlantic salmon.bin"),Configuration::default(),).unwrap().0});
static SALMOTRUTTA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("River trout.bin"),Configuration::default(),).unwrap().0});
static SERIOLAQUINQUERADIATA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Japanese amberjack.bin"),Configuration::default(),).unwrap().0});
static SINIPERCACHUATSI: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Mandarin fish.bin"),Configuration::default(),).unwrap().0});
static SUSSCROFA: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Domestic pig.bin"),Configuration::default(),).unwrap().0});
static TREMATOMUSBERNACCHII: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Emerald rockcod.bin"),Configuration::default(),).unwrap().0});
static VICUGNAPACOS: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("Alpaca.bin"),Configuration::default(),).unwrap().0});
static XENOPUSLAEVISORGILLI: LazyLock<Germlines> = LazyLock::new(|| {bincode::serde::decode_from_slice::<Germlines, Configuration>(include_bytes!("African or Cape clawed frog.bin"),Configuration::default(),).unwrap().0});
