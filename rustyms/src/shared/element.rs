use serde::{Deserialize, Serialize};

use crate::system::f64::Mass;

/// The elements (and electrons)
#[derive(
    Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Default, Serialize, Deserialize,
)]
pub enum Element {
    /// Not necessarily an element but handy to have: electron
    Electron = 0,
    /// Element Hydrogen (H) atomic number: 1
    #[default]
    H = 1,
    /// Element Helium (He) atomic number: 2
    He,
    /// Element Lithium (Li) atomic number: 3
    Li,
    /// Element Beryllium (Be) atomic number: 4
    Be,
    /// Element Boron (B) atomic number: 5
    B,
    /// Element Carbon (C) atomic number: 6
    C,
    /// Element Nitrogen (N) atomic number: 7
    N,
    /// Element Oxygen (O) atomic number: 8
    O,
    /// Element Fluorine (F) atomic number: 9
    F,
    /// Element Neon (Ne) atomic number: 10
    Ne,
    /// Element Sodium (Na) atomic number: 11
    Na,
    /// Element Magnesium (Mg) atomic number: 12
    Mg,
    /// Element Aluminium (Al) atomic number: 13
    Al,
    /// Element Silicon (Si) atomic number: 14
    Si,
    /// Element Phosphorus (P) atomic number: 15
    P,
    /// Element Sulphur (S) atomic number: 16
    S,
    /// Element Chlorine (Cl) atomic number: 17
    Cl,
    /// Element Argon (Ar) atomic number: 18
    Ar,
    /// Element Potassium (K) atomic number: 19
    K,
    /// Element Calcium (Ca) atomic number: 20
    Ca,
    /// Element Scandium (Sc) atomic number: 21
    Sc,
    /// Element Titanium (Ti) atomic number: 22
    Ti,
    /// Element Vanadium (V) atomic number: 23
    V,
    /// Element Chromium (Cr) atomic number: 24
    Cr,
    /// Element Manganese (Mn) atomic number: 25
    Mn,
    /// Element Iron (Fe) atomic number: 26
    Fe,
    /// Element Cobalt (Co) atomic number: 27
    Co,
    /// Element Nickel (Ni) atomic number: 28
    Ni,
    /// Element Copper (Cu) atomic number: 29
    Cu,
    /// Element Zinc (Zn) atomic number: 30
    Zn,
    /// Element Gallium (Ga) atomic number: 31
    Ga,
    /// Element Germanium (Ge) atomic number: 32
    Ge,
    /// Element Arsenic (As) atomic number: 33
    As,
    /// Element Selenium (Se) atomic number: 34
    Se,
    /// Element Bromine (Br) atomic number: 35
    Br,
    /// Element Krypton (Kr) atomic number: 36
    Kr,
    /// Element Rubidium (Rb) atomic number: 37
    Rb,
    /// Element Strontium (Sr) atomic number: 38
    Sr,
    /// Element Yttrium (Y) atomic number: 39
    Y,
    /// Element Zirconium (Zr) atomic number: 40
    Zr,
    /// Element Niobium (Nb) atomic number: 41
    Nb,
    /// Element Molybdenum (Mo) atomic number: 42
    Mo,
    /// Element Technetium (Tc) atomic number: 43
    Tc,
    /// Element Ruthenium (Ru) atomic number: 44
    Ru,
    /// Element Rhodium (Rh) atomic number: 45
    Rh,
    /// Element Palladium (Pd) atomic number: 46
    Pd,
    /// Element Silver (Ag) atomic number: 47
    Ag,
    /// Element Cadmium (Cd) atomic number: 48
    Cd,
    /// Element Indium (In) atomic number: 49
    In,
    /// Element Tin (Sn) atomic number: 50
    Sn,
    /// Element Antimony (Sb) atomic number: 51
    Sb,
    /// Element Tellurium (Te) atomic number: 52
    Te,
    /// Element Iodine (I) atomic number: 53
    I,
    /// Element Xenon (Xe) atomic number: 54
    Xe,
    /// Element Caesium (Cs) atomic number: 55
    Cs,
    /// Element Barium (Ba) atomic number: 56
    Ba,
    /// Element Lanthanum (La) atomic number: 57
    La,
    /// Element Cerium (Ce) atomic number: 58
    Ce,
    /// Element Praseodymium (Pr) atomic number: 59
    Pr,
    /// Element Neodymium (Nd) atomic number: 60
    Nd,
    /// Element Promethium (Pm) atomic number: 61
    Pm,
    /// Element Samarium (Sm) atomic number: 62
    Sm,
    /// Element Europium (Eu) atomic number: 63
    Eu,
    /// Element Gadolinium (Gd) atomic number: 64
    Gd,
    /// Element Terbium (Tb) atomic number: 65
    Tb,
    /// Element Dysprosium (Dy) atomic number: 66
    Dy,
    /// Element Holmium (Ho) atomic number: 67
    Ho,
    /// Element Erbium (Er) atomic number: 68
    Er,
    /// Element Thulium (Tm) atomic number: 69
    Tm,
    /// Element Ytterbium (Yb) atomic number: 70
    Yb,
    /// Element Lutetium (Lu) atomic number: 71
    Lu,
    /// Element Hafnium (Hf) atomic number: 72
    Hf,
    /// Element Tantalum (Ta) atomic number: 73
    Ta,
    /// Element Tungsten (W) atomic number: 74
    W,
    /// Element Rhenium (Re) atomic number: 75
    Re,
    /// Element Osmium (Os) atomic number: 76
    Os,
    /// Element Iridium (Ir) atomic number: 77
    Ir,
    /// Element Platinum (Pt) atomic number: 78
    Pt,
    /// Element Gold (Au) atomic number: 79
    Au,
    /// Element Mercury (Hg) atomic number: 80
    Hg,
    /// Element Thallium (Tl) atomic number: 81
    Tl,
    /// Element Lead (Pb) atomic number: 82
    Pb,
    /// Element Bismuth (Bi) atomic number: 83
    Bi,
    /// Element Polonium (Po) atomic number: 84
    Po,
    /// Element Astatine (At) atomic number: 85
    At,
    /// Element Radon (Rn) atomic number: 86
    Rn,
    /// Element Francium (Fr) atomic number: 87
    Fr,
    /// Element Radium (Ra) atomic number: 88
    Ra,
    /// Element Actinium (Ac) atomic number: 89
    Ac,
    /// Element Thorium (Th) atomic number: 90
    Th,
    /// Element Protactinium (Pa) atomic number: 91
    Pa,
    /// Element Uranium (U) atomic number: 92
    U,
    /// Element Neptunium (Np) atomic number: 93
    Np,
    /// Element Plutonium (Pu) atomic number: 94
    Pu,
    /// Element Americium (Am) atomic number: 95
    Am,
    /// Element Curium (Cm) atomic number: 96
    Cm,
    /// Element Berkelium (Bk) atomic number: 97
    Bk,
    /// Element Californium (Cf) atomic number: 98
    Cf,
    /// Element Einsteinium (Es) atomic number: 99
    Es,
    /// Element Fermium (Fm) atomic number: 100
    Fm,
    /// Element Mendelevium (Md) atomic number: 101
    Md,
    /// Element Nobelium (No) atomic number: 102
    No,
    /// Element Lawrencium (Lr) atomic number: 103
    Lr,
    /// Element Rutherfordium (Rf) atomic number: 104
    Rf,
    /// Element Dubnium (Db) atomic number: 105
    Db,
    /// Element Seaborgium (Sg) atomic number: 106
    Sg,
    /// Element Bohrium (Bh) atomic number: 107
    Bh,
    /// Element Hassium (Hs) atomic number: 108
    Hs,
    /// Element Meitnerium (Mt) atomic number: 109
    Mt,
    /// Element Darmstadtium (Ds) atomic number: 110
    Ds,
    /// Element Roentgenium (Rg) atomic number: 111
    Rg,
    /// Element Copernicium (Cn) atomic number: 112
    Cn,
    /// Element Nihonium (Nh) atomic number: 113
    Nh,
    /// Element Flerovium (Fl) atomic number: 114
    Fl,
    /// Element Moscovium (Mc) atomic number: 115
    Mc,
    /// Element Livermorium (Lv) atomic number: 116
    Lv,
    /// Element Tennessine (Ts) atomic number: 117
    Ts,
    /// Element Oganesson (Og) atomic number: 118
    Og,
}

impl TryFrom<&str> for Element {
    type Error = ();
    #[expect(clippy::too_many_lines)]
    fn try_from(value: &str) -> Result<Self, Self::Error> {
        match value.to_ascii_lowercase().as_str() {
            "h" => Ok(Self::H),
            "he" => Ok(Self::He),
            "li" => Ok(Self::Li),
            "be" => Ok(Self::Be),
            "b" => Ok(Self::B),
            "c" => Ok(Self::C),
            "n" => Ok(Self::N),
            "o" => Ok(Self::O),
            "f" => Ok(Self::F),
            "ne" => Ok(Self::Ne),
            "na" => Ok(Self::Na),
            "mg" => Ok(Self::Mg),
            "al" => Ok(Self::Al),
            "si" => Ok(Self::Si),
            "p" => Ok(Self::P),
            "s" => Ok(Self::S),
            "cl" => Ok(Self::Cl),
            "ar" => Ok(Self::Ar),
            "k" => Ok(Self::K),
            "ca" => Ok(Self::Ca),
            "sc" => Ok(Self::Sc),
            "ti" => Ok(Self::Ti),
            "v" => Ok(Self::V),
            "cr" => Ok(Self::Cr),
            "mn" => Ok(Self::Mn),
            "fe" => Ok(Self::Fe),
            "co" => Ok(Self::Co),
            "ni" => Ok(Self::Ni),
            "cu" => Ok(Self::Cu),
            "zn" => Ok(Self::Zn),
            "ga" => Ok(Self::Ga),
            "ge" => Ok(Self::Ge),
            "as" => Ok(Self::As),
            "se" => Ok(Self::Se),
            "br" => Ok(Self::Br),
            "kr" => Ok(Self::Kr),
            "rb" => Ok(Self::Rb),
            "sr" => Ok(Self::Sr),
            "y" => Ok(Self::Y),
            "zr" => Ok(Self::Zr),
            "nb" => Ok(Self::Nb),
            "mo" => Ok(Self::Mo),
            "tc" => Ok(Self::Tc),
            "ru" => Ok(Self::Ru),
            "rh" => Ok(Self::Rh),
            "pd" => Ok(Self::Pd),
            "ag" => Ok(Self::Ag),
            "cd" => Ok(Self::Cd),
            "in" => Ok(Self::In),
            "sn" => Ok(Self::Sn),
            "sb" => Ok(Self::Sb),
            "te" => Ok(Self::Te),
            "i" => Ok(Self::I),
            "xe" => Ok(Self::Xe),
            "cs" => Ok(Self::Cs),
            "ba" => Ok(Self::Ba),
            "la" => Ok(Self::La),
            "ce" => Ok(Self::Ce),
            "pr" => Ok(Self::Pr),
            "nd" => Ok(Self::Nd),
            "pm" => Ok(Self::Pm),
            "sm" => Ok(Self::Sm),
            "eu" => Ok(Self::Eu),
            "gd" => Ok(Self::Gd),
            "tb" => Ok(Self::Tb),
            "dy" => Ok(Self::Dy),
            "ho" => Ok(Self::Ho),
            "er" => Ok(Self::Er),
            "tm" => Ok(Self::Tm),
            "yb" => Ok(Self::Yb),
            "lu" => Ok(Self::Lu),
            "hf" => Ok(Self::Hf),
            "ta" => Ok(Self::Ta),
            "w" => Ok(Self::W),
            "re" => Ok(Self::Re),
            "os" => Ok(Self::Os),
            "ir" => Ok(Self::Ir),
            "pt" => Ok(Self::Pt),
            "au" => Ok(Self::Au),
            "hg" => Ok(Self::Hg),
            "tl" => Ok(Self::Tl),
            "pb" => Ok(Self::Pb),
            "bi" => Ok(Self::Bi),
            "po" => Ok(Self::Po),
            "at" => Ok(Self::At),
            "rn" => Ok(Self::Rn),
            "fr" => Ok(Self::Fr),
            "ra" => Ok(Self::Ra),
            "ac" => Ok(Self::Ac),
            "th" => Ok(Self::Th),
            "pa" => Ok(Self::Pa),
            "u" => Ok(Self::U),
            "np" => Ok(Self::Np),
            "pu" => Ok(Self::Pu),
            "am" => Ok(Self::Am),
            "cm" => Ok(Self::Cm),
            "bk" => Ok(Self::Bk),
            "cf" => Ok(Self::Cf),
            "es" => Ok(Self::Es),
            "fm" => Ok(Self::Fm),
            "md" => Ok(Self::Md),
            "no" => Ok(Self::No),
            "lr" => Ok(Self::Lr),
            "rf" => Ok(Self::Rf),
            "db" => Ok(Self::Db),
            "sg" => Ok(Self::Sg),
            "bh" => Ok(Self::Bh),
            "hs" => Ok(Self::Hs),
            "mt" => Ok(Self::Mt),
            "ds" => Ok(Self::Ds),
            "rg" => Ok(Self::Rg),
            "cn" => Ok(Self::Cn),
            "nh" => Ok(Self::Nh),
            "fl" => Ok(Self::Fl),
            "mc" => Ok(Self::Mc),
            "lv" => Ok(Self::Lv),
            "ts" => Ok(Self::Ts),
            "og" => Ok(Self::Og),
            _ => Err(()),
        }
    }
}

impl TryFrom<usize> for Element {
    type Error = ();
    #[expect(clippy::too_many_lines)]
    fn try_from(value: usize) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(Self::Electron),
            1 => Ok(Self::H),
            2 => Ok(Self::He),
            3 => Ok(Self::Li),
            4 => Ok(Self::Be),
            5 => Ok(Self::B),
            6 => Ok(Self::C),
            7 => Ok(Self::N),
            8 => Ok(Self::O),
            9 => Ok(Self::F),
            10 => Ok(Self::Ne),
            11 => Ok(Self::Na),
            12 => Ok(Self::Mg),
            13 => Ok(Self::Al),
            14 => Ok(Self::Si),
            15 => Ok(Self::P),
            16 => Ok(Self::S),
            17 => Ok(Self::Cl),
            18 => Ok(Self::Ar),
            19 => Ok(Self::K),
            20 => Ok(Self::Ca),
            21 => Ok(Self::Sc),
            22 => Ok(Self::Ti),
            23 => Ok(Self::V),
            24 => Ok(Self::Cr),
            25 => Ok(Self::Mn),
            26 => Ok(Self::Fe),
            27 => Ok(Self::Co),
            28 => Ok(Self::Ni),
            29 => Ok(Self::Cu),
            30 => Ok(Self::Zn),
            31 => Ok(Self::Ga),
            32 => Ok(Self::Ge),
            33 => Ok(Self::As),
            34 => Ok(Self::Se),
            35 => Ok(Self::Br),
            36 => Ok(Self::Kr),
            37 => Ok(Self::Rb),
            38 => Ok(Self::Sr),
            39 => Ok(Self::Y),
            40 => Ok(Self::Zr),
            41 => Ok(Self::Nb),
            42 => Ok(Self::Mo),
            43 => Ok(Self::Tc),
            44 => Ok(Self::Ru),
            45 => Ok(Self::Rh),
            46 => Ok(Self::Pd),
            47 => Ok(Self::Ag),
            48 => Ok(Self::Cd),
            49 => Ok(Self::In),
            50 => Ok(Self::Sn),
            51 => Ok(Self::Sb),
            52 => Ok(Self::Te),
            53 => Ok(Self::I),
            54 => Ok(Self::Xe),
            55 => Ok(Self::Cs),
            56 => Ok(Self::Ba),
            57 => Ok(Self::La),
            58 => Ok(Self::Ce),
            59 => Ok(Self::Pr),
            60 => Ok(Self::Nd),
            61 => Ok(Self::Pm),
            62 => Ok(Self::Sm),
            63 => Ok(Self::Eu),
            64 => Ok(Self::Gd),
            65 => Ok(Self::Tb),
            66 => Ok(Self::Dy),
            67 => Ok(Self::Ho),
            68 => Ok(Self::Er),
            69 => Ok(Self::Tm),
            70 => Ok(Self::Yb),
            71 => Ok(Self::Lu),
            72 => Ok(Self::Hf),
            73 => Ok(Self::Ta),
            74 => Ok(Self::W),
            75 => Ok(Self::Re),
            76 => Ok(Self::Os),
            77 => Ok(Self::Ir),
            78 => Ok(Self::Pt),
            79 => Ok(Self::Au),
            80 => Ok(Self::Hg),
            81 => Ok(Self::Tl),
            82 => Ok(Self::Pb),
            83 => Ok(Self::Bi),
            84 => Ok(Self::Po),
            85 => Ok(Self::At),
            86 => Ok(Self::Rn),
            87 => Ok(Self::Fr),
            88 => Ok(Self::Ra),
            89 => Ok(Self::Ac),
            90 => Ok(Self::Th),
            91 => Ok(Self::Pa),
            92 => Ok(Self::U),
            93 => Ok(Self::Np),
            94 => Ok(Self::Pu),
            95 => Ok(Self::Am),
            96 => Ok(Self::Cm),
            97 => Ok(Self::Bk),
            98 => Ok(Self::Cf),
            99 => Ok(Self::Es),
            100 => Ok(Self::Fm),
            101 => Ok(Self::Md),
            102 => Ok(Self::No),
            103 => Ok(Self::Lr),
            104 => Ok(Self::Rf),
            105 => Ok(Self::Db),
            106 => Ok(Self::Sg),
            107 => Ok(Self::Bh),
            108 => Ok(Self::Hs),
            109 => Ok(Self::Mt),
            110 => Ok(Self::Ds),
            111 => Ok(Self::Rg),
            112 => Ok(Self::Cn),
            113 => Ok(Self::Nh),
            114 => Ok(Self::Fl),
            115 => Ok(Self::Mc),
            116 => Ok(Self::Lv),
            117 => Ok(Self::Ts),
            118 => Ok(Self::Og),
            _ => Err(()),
        }
    }
}

impl Element {
    pub fn symbol(&self) -> &'static str {
        match self {
            Self::H => "H",
            Self::He => "He",
            Self::Li => "Li",
            Self::Be => "Be",
            Self::B => "B",
            Self::C => "C",
            Self::N => "N",
            Self::O => "O",
            Self::F => "F",
            Self::Ne => "Ne",
            Self::Na => "Na",
            Self::Mg => "Mg",
            Self::Al => "Al",
            Self::Si => "Si",
            Self::P => "P",
            Self::S => "S",
            Self::Cl => "Cl",
            Self::Ar => "Ar",
            Self::K => "K",
            Self::Ca => "Ca",
            Self::Sc => "Sc",
            Self::Ti => "Ti",
            Self::V => "V",
            Self::Cr => "Cr",
            Self::Mn => "Mn",
            Self::Fe => "Fe",
            Self::Co => "Co",
            Self::Ni => "Ni",
            Self::Cu => "Cu",
            Self::Zn => "Zn",
            Self::Ga => "Ga",
            Self::Ge => "Ge",
            Self::As => "As",
            Self::Se => "Se",
            Self::Br => "Br",
            Self::Kr => "Kr",
            Self::Rb => "Rb",
            Self::Sr => "Sr",
            Self::Y => "Y",
            Self::Zr => "Zr",
            Self::Nb => "Nb",
            Self::Mo => "Mo",
            Self::Tc => "Tc",
            Self::Ru => "Ru",
            Self::Rh => "Rh",
            Self::Pd => "Pd",
            Self::Ag => "Ag",
            Self::Cd => "Cd",
            Self::In => "In",
            Self::Sn => "Sn",
            Self::Sb => "Sb",
            Self::Te => "Te",
            Self::I => "I",
            Self::Xe => "Xe",
            Self::Cs => "Cs",
            Self::Ba => "Ba",
            Self::La => "La",
            Self::Ce => "Ce",
            Self::Pr => "Pr",
            Self::Nd => "Nd",
            Self::Pm => "Pm",
            Self::Sm => "Sm",
            Self::Eu => "Eu",
            Self::Gd => "Gd",
            Self::Tb => "Tb",
            Self::Dy => "Dy",
            Self::Ho => "Ho",
            Self::Er => "Er",
            Self::Tm => "Tm",
            Self::Yb => "Yb",
            Self::Lu => "Lu",
            Self::Hf => "Hf",
            Self::Ta => "Ta",
            Self::W => "W",
            Self::Re => "Re",
            Self::Os => "Os",
            Self::Ir => "Ir",
            Self::Pt => "Pt",
            Self::Au => "Au",
            Self::Hg => "Hg",
            Self::Tl => "Tl",
            Self::Pb => "Pb",
            Self::Bi => "Bi",
            Self::Po => "Po",
            Self::At => "At",
            Self::Rn => "Rn",
            Self::Fr => "Fr",
            Self::Ra => "Ra",
            Self::Ac => "Ac",
            Self::Th => "Th",
            Self::Pa => "Pa",
            Self::U => "U",
            Self::Np => "Np",
            Self::Pu => "Pu",
            Self::Am => "Am",
            Self::Cm => "Cm",
            Self::Bk => "Bk",
            Self::Cf => "Cf",
            Self::Es => "Es",
            Self::Fm => "Fm",
            Self::Md => "Md",
            Self::No => "No",
            Self::Lr => "Lr",
            Self::Rf => "Rf",
            Self::Db => "Db",
            Self::Sg => "Sg",
            Self::Bh => "Bh",
            Self::Hs => "Hs",
            Self::Mt => "Mt",
            Self::Ds => "Ds",
            Self::Rg => "Rg",
            Self::Cn => "Cn",
            Self::Nh => "Nh",
            Self::Fl => "Fl",
            Self::Mc => "Mc",
            Self::Lv => "Lv",
            Self::Ts => "Ts",
            Self::Og => "Og",
            // Self::Proton => "Proton",
            Self::Electron => "e",
        }
    }
}

impl std::fmt::Display for Element {
    #[expect(clippy::too_many_lines)]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.symbol())
    }
}

/// All elements sorted so that single characters come after two character element symbols (needed for greedy parsing)
pub const ELEMENT_PARSE_LIST: &[(&str, Element)] = &[
    ("he", Element::He),
    ("li", Element::Li),
    ("be", Element::Be),
    ("ne", Element::Ne),
    ("na", Element::Na),
    ("mg", Element::Mg),
    ("al", Element::Al),
    ("si", Element::Si),
    ("cl", Element::Cl),
    ("ar", Element::Ar),
    ("ca", Element::Ca),
    ("sc", Element::Sc),
    ("ti", Element::Ti),
    ("cr", Element::Cr),
    ("mn", Element::Mn),
    ("fe", Element::Fe),
    ("co", Element::Co),
    ("ni", Element::Ni),
    ("cu", Element::Cu),
    ("zn", Element::Zn),
    ("ga", Element::Ga),
    ("ge", Element::Ge),
    ("as", Element::As),
    ("se", Element::Se),
    ("br", Element::Br),
    ("kr", Element::Kr),
    ("rb", Element::Rb),
    ("sr", Element::Sr),
    ("zr", Element::Zr),
    ("nb", Element::Nb),
    ("mo", Element::Mo),
    ("tc", Element::Tc),
    ("ru", Element::Ru),
    ("rh", Element::Rh),
    ("pd", Element::Pd),
    ("ag", Element::Ag),
    ("cd", Element::Cd),
    ("in", Element::In),
    ("sn", Element::Sn),
    ("sb", Element::Sb),
    ("te", Element::Te),
    ("xe", Element::Xe),
    ("cs", Element::Cs),
    ("ba", Element::Ba),
    ("la", Element::La),
    ("ce", Element::Ce),
    ("pr", Element::Pr),
    ("nd", Element::Nd),
    ("pm", Element::Pm),
    ("sm", Element::Sm),
    ("eu", Element::Eu),
    ("gd", Element::Gd),
    ("tb", Element::Tb),
    ("dy", Element::Dy),
    ("ho", Element::Ho),
    ("er", Element::Er),
    ("tm", Element::Tm),
    ("yb", Element::Yb),
    ("lu", Element::Lu),
    ("hf", Element::Hf),
    ("ta", Element::Ta),
    ("re", Element::Re),
    ("os", Element::Os),
    ("ir", Element::Ir),
    ("pt", Element::Pt),
    ("au", Element::Au),
    ("hg", Element::Hg),
    ("tl", Element::Tl),
    ("pb", Element::Pb),
    ("bi", Element::Bi),
    ("po", Element::Po),
    ("at", Element::At),
    ("rn", Element::Rn),
    ("fr", Element::Fr),
    ("ra", Element::Ra),
    ("ac", Element::Ac),
    ("th", Element::Th),
    ("pa", Element::Pa),
    ("np", Element::Np),
    ("pu", Element::Pu),
    ("am", Element::Am),
    ("cm", Element::Cm),
    ("bk", Element::Bk),
    ("cf", Element::Cf),
    ("es", Element::Es),
    ("fm", Element::Fm),
    ("md", Element::Md),
    ("no", Element::No),
    ("lr", Element::Lr),
    ("rf", Element::Rf),
    ("db", Element::Db),
    ("sg", Element::Sg),
    ("bh", Element::Bh),
    ("hs", Element::Hs),
    ("mt", Element::Mt),
    ("ds", Element::Ds),
    ("rg", Element::Rg),
    ("cn", Element::Cn),
    ("nh", Element::Nh),
    ("fl", Element::Fl),
    ("mc", Element::Mc),
    ("lv", Element::Lv),
    ("ts", Element::Ts),
    ("og", Element::Og),
    ("u", Element::U),
    ("w", Element::W),
    ("i", Element::I),
    ("y", Element::Y),
    ("v", Element::V),
    ("k", Element::K),
    ("s", Element::S),
    ("b", Element::B),
    ("c", Element::C),
    ("n", Element::N),
    ("o", Element::O),
    ("f", Element::F),
    ("h", Element::H),
    ("p", Element::P),
];

/// Most common elements selected for abundance in biological systems as well as unambiguous parsing.
/// Sorted so that single characters come after two character element symbols (needed for greedy parsing).
pub const COMMON_ELEMENT_PARSE_LIST: &[(&str, Element)] = &[
    ("he", Element::He),
    ("li", Element::Li),
    ("na", Element::Na),
    ("mg", Element::Mg),
    ("al", Element::Al),
    ("cl", Element::Cl),
    ("ca", Element::Ca),
    ("fe", Element::Fe),
    ("cu", Element::Cu),
    ("zn", Element::Zn),
    ("se", Element::Se),
    ("ag", Element::Ag),
    ("au", Element::Au),
    ("i", Element::I),
    ("k", Element::K),
    ("s", Element::S),
    ("b", Element::B),
    ("c", Element::C),
    ("n", Element::N),
    ("o", Element::O),
    ("f", Element::F),
    ("h", Element::H),
    ("p", Element::P),
];

#[expect(clippy::redundant_pub_crate)]
/// The shared type to send the data from all the elements from build time to compile time
pub(crate) type ElementalData = Vec<(Option<Mass>, Option<Mass>, Vec<(u16, Mass, f64)>)>;
