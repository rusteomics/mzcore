//! Code to make alignments of two peptides based on mass mistakes, and genetic information.
//!
//! A mass based alignment handles the case in which multiple amino acids are wrong, but the total mass
//! of this set of amino acids is equal to the mass of a set of different amino acids on the other peptide.
//! This is quite common in mass spectrometry where mistakes based on mass coincidences are very common.
//! For example `N` has the same mass as `GG`, so if we want to make a mass spectrometry faithful alignment
//! of `ANA` with `AGGA` the result should reflect this fact:
//!
//! ```text
//! Identity: 0.500 (2/4), Similarity: 0.750 (3/4), Gaps: 0.000 (0/4), Score: 0.706 (12/17),
//! Equal mass, Tolerance: 10 ppm, Alignment: global
//! Start: A 0 B 0, Path: 1=1:2i1=
//!
//! AN·A A
//! AGGA B
//!  ╶╴
//! ```
//! _Generated using this algorithm bound to a cli tool: <https://github.com/snijderlab/align-cli>_
//! ```rust
//! use rustyms::{*, align::*};
//! let a = LinearPeptide::pro_forma("ANA", None).unwrap().into_simple_linear().unwrap();
//! let b = LinearPeptide::pro_forma("AGGA", None).unwrap().into_simple_linear().unwrap();
//! let alignment = align::<4, SimpleLinear, SimpleLinear>(&a, &b, AlignScoring::default(), AlignType::GLOBAL);
//! assert_eq!(alignment.short(), "1=1:2i1=");
//! assert_eq!(alignment.ppm().value, 0.0);
//! ```

mod align_type;
mod alignment;
#[cfg(test)]
mod bad_alignments;
mod diagonal_array;
mod mass_alignment;
mod piece;
mod scoring;
#[cfg(test)]
mod test_alignments;

#[cfg(feature = "imgt")]
mod consecutive;
#[cfg(feature = "imgt")]
pub use consecutive::*;

pub use align_type::{AlignType, Side};
pub use alignment::{Alignment, Score, Stats};
pub use mass_alignment::align;
pub use piece::Piece;
pub use scoring::{AlignScoring, MatchType};

/// Different scoring matrices that can be used.
/// Matrices from: <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/util/tables/> and <https://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/>
///
/// The UO columns are added, for these the B/J/Z score is the rounded down average of the corresponding non ambiguous AAs. All UO scores are exactly the same for all matrices (except identity).
pub mod matrix {
    use super::scoring;
    pub use scoring::matrices::*;
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{peptide::SimpleLinear, LinearPeptide};

    use super::{scoring::AlignScoring, AlignType, Alignment};

    fn align<'a, const STEPS: u16>(
        a: &'a LinearPeptide<SimpleLinear>,
        b: &'a LinearPeptide<SimpleLinear>,
    ) -> Alignment<'a, SimpleLinear, SimpleLinear> {
        super::align::<STEPS, SimpleLinear, SimpleLinear>(
            a,
            b,
            AlignScoring::default(),
            AlignType::GLOBAL,
        )
    }

    fn linear(aa: &str) -> LinearPeptide<SimpleLinear> {
        LinearPeptide::pro_forma(aa, None)
            .unwrap()
            .into_simple_linear()
            .unwrap()
    }

    #[test]
    fn simple_1() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<1>(&a, &b));
        assert_eq!(c.short(), "1=1X1=1X2=");
        assert_eq!(
            Alignment::create_from_path(
                &a,
                &b,
                0,
                0,
                &c.short(),
                AlignScoring::default(),
                AlignType::GLOBAL,
                1
            )
            .unwrap(),
            c
        );
    }

    #[test]
    fn simple_4() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<4>(&a, &b));
        assert_eq!(c.short(), "1=1:2i2:1i2=");
        assert_eq!(
            Alignment::create_from_path(
                &a,
                &b,
                0,
                0,
                &c.short(),
                AlignScoring::default(),
                AlignType::GLOBAL,
                4
            )
            .unwrap(),
            c
        );
    }

    #[test]
    fn simple_unbounded() {
        let a = linear("ANGARS");
        let b = linear("AGGQRS");
        let c = dbg!(align::<{ u16::MAX }>(&a, &b));
        assert_eq!(c.short(), "1=1:2i2:1i2=");
        assert_eq!(
            Alignment::create_from_path(
                &a,
                &b,
                0,
                0,
                &c.short(),
                AlignScoring::default(),
                AlignType::GLOBAL,
                u16::MAX
            )
            .unwrap(),
            c
        );
    }
}
