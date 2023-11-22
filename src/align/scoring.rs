/// The type of a single match step
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq, Hash)]
pub enum MatchType {
    /// Aminoacid + Mass identity
    FullIdentity,
    /// Aminoacid + Mass mismatch
    IdentityMassMismatch,
    /// Full mismatch
    #[default]
    Mismatch,
    /// Set of aminoacids + mods with the same mass but different sequence
    Isobaric,
    /// Set of aminoacids + mods in a different order in the two sequences
    Switched,
    /// A gap
    Gap,
}

pub const MISMATCH: i8 = -1;
pub const MASS_MISMATCH_PENALTY: i8 = -1;
pub const SWITCHED: i8 = 3;
pub const ISOMASS: i8 = 2;
pub const GAP_START_PENALTY: i8 = -5;
pub const GAP_EXTEND_PENALTY: i8 = -1;

/// The (slightly modified) blosum62 matrix for aminoacid homology scoring
pub const BLOSUM62: &[&[i8]] = include!("blosum62.txt");
