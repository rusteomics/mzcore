use crate::chemistry::Element;

#[derive(Debug)]
pub struct Wurcs {
    pub residues: Vec<Residue>,
    pub sequence: Vec<u8>,
    pub linkage: Vec<Linkage>,
}

#[derive(Debug)]
pub struct Residue {
    pub backbone: BackBone,
    pub anomeric: Option<(Option<u8>, AnomericSymbol)>,
    pub mods: Vec<Mod>,
}

#[derive(Debug)]
pub enum BackBone {
    Defined(TerminalCarbon, Vec<Carbon>, TerminalCarbon),
    Repeating(Option<TerminalCarbon>, Vec<Carbon>, Option<TerminalCarbon>),
}

/// The carbon descriptors extended with the options in: https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.6b00650/suppl_file/ci6b00650_si_001.pdf section 2.8.
#[derive(Clone, Copy, Debug)]
pub enum Carbon {
    /// 'd' deoxy `H-C-H`
    Methylene,
    /// 'C' `X-C-X`
    Dual,
    /// '1' `X-C-H`
    HydroxyLeft,
    /// '2' `H-C-X`
    HydroxyRight,
    /// '3' `C(X)(H)`
    HydroxyOpposite,
    /// '4' `C(H)(X)`
    HydroxySame,
    /// 'x' one of 1 or 2
    HydroxyUnknown,
    /// '5' `X-C-Y`
    DualLeft,
    /// '6' `Y-C-X`
    DualRight,
    /// '7' `C(X)(Y)`
    DualOpposite,
    /// '8' `C(Y)(X)`
    DualSame,
    /// 'X' one of 5 or 6
    DualUnknown,
    /// 'O' `C=O`
    Ketone,
    /// `e`
    DoubleHydroxyEntgegen,
    /// `z`
    DoubleHydroxyZusammen,
    /// `n`
    DoubleHydroxyNoIsomer,
    /// `f`
    DoubleHydroxyUnknown,
    /// `E`
    DoubleEntgegen,
    /// `Z`
    DoubleZusammen,
    /// `N`
    DoubleNoIsomer,
    /// `F`
    DoubleUnknown,
    /// `K`
    DoubleBonded,
    /// `T`
    TripleBonded,
    /// 'a' anomeric
    Hemiketal,
    /// 'U'
    KetoneOrHemiketal,
    /// 'Q' can be any of the other carbon descriptors
    Unknown,
}

#[derive(Clone, Copy, Debug)]
#[allow(non_camel_case_types)]
pub enum TerminalCarbon {
    /// 'm'
    CHHH,
    /// 'M'
    CXXX,
    /// 'h'
    CHHX,
    /// 'c'
    CXXH,
    /// 'C'
    CXXY,
    /// '1'
    CXYH,
    /// '2'
    CYXH,
    /// '3'
    CXYH_opposite,
    /// '4'
    CYXH_same,
    /// 'x'
    CXYH_unknown,
    /// '5'
    CXYZ,
    /// '6'
    CYXZ,
    /// '7'
    CXYZ_opposite,
    /// '8'
    CYXZ_same,
    /// 'X'
    CXYZ_unknown,
    /// 'o' -C=XH
    CXH,
    /// 'A' -C=XY
    CXY,
    /// 'n' =CHH
    CHH,
    /// 'N' =CXX
    CXX,
    /// 'e'
    Double_CXH_entgegen,
    /// 'z'
    Double_CXH_zusammen,
    /// 'f'
    Double_CXH_unknown,
    /// 'E'
    Double_CXY_entgegen,
    /// 'Z'
    Double_CXY_zusammen,
    /// 'F'
    Double_CXY_unknown,
    /// 'T' ≡C-X or -C≡X
    Triple_CX,
    /// 'K'
    Double_CX,
    /// 't'
    Triple_CH,
    /// 'a'
    Hemiacetal,
    /// 'u'
    AldehydeOrHemiacetal,
    /// 'Q'
    Unknown,
}

#[derive(Clone, Copy, Debug)]
pub enum AnomericSymbol {
    /// 'a'
    Alpha,
    /// 'b'
    Beta,
    /// 'u'
    Up,
    /// 'd'
    Down,
    /// 'x'
    Unknown,
    /// 'o'
    None,
}

#[derive(Debug)]
pub struct Mod {
    pub lips: Vec<LIPOption>,
    pub modification: Vec<MAPSymbol>,
}

#[derive(Clone, Debug)]
pub enum LIPOption {
    Known(LIP),
    Statistic(bool, Probability, LIP),
    Alternative(Vec<LIP>),
}

#[derive(Clone, Copy, Debug)]
pub struct LIP {
    pub position: Option<u8>,
    pub direction: Direction,
    pub star_index: u8,
}

#[derive(Clone, Copy, Debug)]
pub enum Direction {
    /// 'u'
    Upside,
    /// 'd'
    Downside,
    /// 't'
    Tres,
    /// 'a'
    Same,
    /// 'b'
    Opposite,
    /// 'c'
    Third,
    /// 'x'
    Unknown,
    /// 'e'
    Entgegen,
    /// 'z'
    Zusammen,
    /// 'f'
    UnknownGeometricalIsomerism,
    /// 'n'
    Obvious,
}

#[derive(Clone, Copy, Debug)]
pub enum MAPSymbol {
    Element(Element),
    /// '*' or '*n'
    Star(Option<u8>),
    /// '/n'
    Branch(u8),
    /// '$n'
    Cyclic(u8),
    /// '='
    DoubleBond,
    /// '#'
    TripleBond,
    Chirality(Chirality),
    AromaticStart,
    AromaticEnd,
}

#[derive(Clone, Copy, Debug)]
pub enum Chirality {
    R,
    S,
    E,
    Z,
    Unknown,
}

#[derive(Debug)]
pub enum Linkage {
    Known(LIN),
    Repeated(Repeat, LIN),
}

#[derive(Debug)]
pub struct LIN {
    pub lips: Vec<GLIPOption>,
    pub modification: Vec<MAPSymbol>,
}

#[derive(Clone, Debug)]
pub enum GLIPOption {
    Known(GLIP),
    Statistic(bool, Probability, GLIP),
    Alternative(Vec<GLIP>),
    RESAlternative(bool, Vec<GLIP>),
}

#[derive(Clone, Copy, Debug)]
pub struct GLIP {
    /// Encoded using base 52
    pub res_index: u8,
    /// Position or unknown
    pub position: Option<u8>,
    pub direction: Direction,
    pub star_index: u8,
}

#[derive(Clone, Copy, Debug)]
pub enum Probability {
    Single(Option<f32>),
    Range(Option<f32>, Option<f32>),
}

#[derive(Clone, Copy, Debug)]
pub enum Repeat {
    Single(Option<u8>),
    /// Max, min
    Range(Option<u8>, Option<u8>),
}
