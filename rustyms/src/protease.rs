use std::sync::LazyLock;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{AminoAcid, SequenceElement};

/// A protease defined by it ability to cut at any site identified by the right amino acids at the n and c terminal.
/// Each position is identified by an option, a none means that there is no specificity at this position. If there is
/// a specificity at a certain position any amino acid that is contained in the set is allowed (see
/// [`crate::CheckedAminoAcid::canonical_identical`]).
///
/// #TODO: Add the possibility to define the effiency of amino acids or certain motifs
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Protease {
    /// The amino acids n terminal of the cut site.
    pub before: Vec<Option<Vec<AminoAcid>>>,
    /// The amino acids c terminal of the cut site.
    pub after: Vec<Option<Vec<AminoAcid>>>,
}

impl Protease {
    /// Define a simple protease that cuts exactly between the specified sequences.
    pub fn new(n_term: Vec<AminoAcid>, c_term: Vec<AminoAcid>) -> Self {
        Self {
            before: vec![Some(n_term)],
            after: vec![Some(c_term)],
        }
    }

    /// Define a protease that cuts on the c terminal side of the provided amino acids.
    pub fn c_terminal_of(residues: Vec<AminoAcid>) -> Self {
        Self {
            before: vec![Some(residues)],
            after: Vec::new(),
        }
    }

    /// Define a protease that cuts on the n terminal side of the provided amino acids.
    pub fn n_terminal_of(residues: Vec<AminoAcid>) -> Self {
        Self {
            before: Vec::new(),
            after: vec![Some(residues)],
        }
    }

    /// Helper function to get a list of all amino acids except the ones given
    pub fn get_exclusive(exclude: &[AminoAcid]) -> Vec<AminoAcid> {
        AminoAcid::ALL_AMINO_ACIDS
            .iter()
            .copied()
            .filter(|aa| !exclude.contains(aa))
            .collect_vec()
    }

    /// All locations in the given sequence where this protease could cut
    /// Note if a cutsite is "after" the last element in the sequence, it will not rapport that, only the cutsites inside the sequence
    pub fn match_locations<T>(&self, sequence: &[SequenceElement<T>]) -> Vec<usize> {
        let upper = sequence
            .len()
            .saturating_sub(self.after.len())
            .min(sequence.len().saturating_sub(1));
        (self.before.len()..=upper)
            .filter(|i| self.matches_at(&sequence[i - self.before.len()..i + self.after.len()]))
            .collect_vec()
    }

    fn matches_at<T>(&self, slice: &[SequenceElement<T>]) -> bool {
        debug_assert!(slice.len() == self.before.len() + self.after.len());
        'positions: for (actual, pattern) in slice
            .iter()
            .zip(self.before.iter().chain(self.after.iter()))
        {
            if let Some(pattern) = pattern {
                for option in pattern {
                    if option.canonical_identical(actual.aminoacid.aminoacid()) {
                        continue 'positions;
                    }
                }
                return false;
            }
        }
        true
    }
}

/// Some well known and widely used proteases
pub mod known_proteases {
    use super::*;

    /// `Trypsin` cuts after Lysine (K) or Arginine (R), unless followed by Proline (P)
    pub static TRYPSIN: LazyLock<Protease> = LazyLock::new(|| {
        Protease::new(
            vec![AminoAcid::Lysine, AminoAcid::Arginine],
            Protease::get_exclusive(&[AminoAcid::Proline]),
        )
    });

    /// `Chymotrypsin` cuts after Phenylalanine (F), Tryptophan (W), Tyrosine (Y), unless followed by Proline (P)
    pub static CHYMOTRYPSIN: LazyLock<Protease> = LazyLock::new(|| {
        Protease::new(
            vec![
                AminoAcid::Phenylalanine,
                AminoAcid::Tryptophan,
                AminoAcid::Tyrosine,
            ],
            Protease::get_exclusive(&[AminoAcid::Proline]),
        )
    });

    /// `Pepsin` (pH > 2) cuts after Phenylalanine (F), Tryptophan (W), Tyrosine (Y), Leucine (L)
    pub static PEPSIN: LazyLock<Protease> = LazyLock::new(|| {
        Protease::c_terminal_of(vec![
            AminoAcid::Phenylalanine,
            AminoAcid::Tryptophan,
            AminoAcid::Tyrosine,
            AminoAcid::Leucine,
        ])
    });

    /// `AspN` cuts before Aspartic acid (D)
    pub static ASPN: LazyLock<Protease> =
        LazyLock::new(|| Protease::n_terminal_of(vec![AminoAcid::AsparticAcid]));

    /// `GluC` cuts after Glutamic acid (E)
    pub static GLUC: LazyLock<Protease> =
        LazyLock::new(|| Protease::c_terminal_of(vec![AminoAcid::GlutamicAcid]));

    /// `LysC` cuts after Lysine (K)
    pub static LYSC: LazyLock<Protease> =
        LazyLock::new(|| Protease::c_terminal_of(vec![AminoAcid::Lysine]));

    /// `ArgC` cuts after Arginine (R)
    pub static ARGC: LazyLock<Protease> =
        LazyLock::new(|| Protease::c_terminal_of(vec![AminoAcid::Arginine]));
}

#[cfg(test)]
#[allow(clippy::missing_panics_doc)]
mod tests {
    use crate::{Linear, Peptidoform};

    use super::*;

    pub struct ProteaseTestCase {
        pub sequence: Peptidoform<Linear>,
        pub expected_cut_sites: Vec<usize>,
        pub expected_peptides: Vec<Peptidoform<Linear>>,
    }

    /// Generic test function for all proteases
    pub fn test_protease(protease: &Protease, test_case: &ProteaseTestCase) {
        // Test cut sites
        let cut_sites = protease.match_locations(test_case.sequence.sequence());

        assert_eq!(
            cut_sites, test_case.expected_cut_sites,
            "Incorrect cut sites: found '{cut_sites:?}' expected '{:?}'",
            test_case.expected_cut_sites
        );

        // Test peptides
        let peptides = test_case.sequence.digest(protease, 0);

        if peptides.len() != test_case.expected_peptides.len() {
            for peptide in &peptides {
                println!("{peptide}");
            }
            panic!("Incorrect number of peptides")
        }

        for (i, peptide) in peptides.iter().enumerate() {
            assert_eq!(
                peptide, &test_case.expected_peptides[i],
                "Peptides don't match: found '{peptide}' expected '{}'",
                test_case.expected_peptides[i]
            );
        }
    }

    fn str_to_peptidoform(str_peptide: &str) -> Peptidoform<Linear> {
        Peptidoform::pro_forma(str_peptide, None)
            .unwrap()
            .into_linear()
            .unwrap()
    }

    #[test]
    fn trypsin() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("AKRPGKR"),
                expected_cut_sites: vec![2, 6],
                expected_peptides: vec![
                    str_to_peptidoform("AK"),
                    str_to_peptidoform("RPGK"),
                    str_to_peptidoform("R"),
                ],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("ARAKGCVLRPKDGR"),
                expected_cut_sites: vec![2, 4, 11],
                expected_peptides: vec![
                    str_to_peptidoform("AR"),
                    str_to_peptidoform("AK"),
                    str_to_peptidoform("GCVLRPK"),
                    str_to_peptidoform("DGR"),
                ],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::TRYPSIN, &test_case);
        }
    }

    #[test]
    fn chymotrypsin() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("AFWYPLGF"),
                expected_cut_sites: vec![2, 3],
                expected_peptides: vec![
                    str_to_peptidoform("AF"),
                    str_to_peptidoform("W"),
                    str_to_peptidoform("YPLGF"),
                ],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("AVFUDGWTYPMSR"),
                expected_cut_sites: vec![3, 7],
                expected_peptides: vec![
                    str_to_peptidoform("AVF"),
                    str_to_peptidoform("UDGW"),
                    str_to_peptidoform("TYPMSR"),
                ],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::CHYMOTRYPSIN, &test_case);
        }
    }

    #[test]
    fn pepsin() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("AACVFLPAKLURF"),
                expected_cut_sites: vec![5, 6, 10],
                expected_peptides: vec![
                    str_to_peptidoform("AACVF"),
                    str_to_peptidoform("L"),
                    str_to_peptidoform("PAKL"),
                    str_to_peptidoform("URF"),
                ],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("GFLPKDLVMSRG"),
                expected_cut_sites: vec![2, 3, 7],
                expected_peptides: vec![
                    str_to_peptidoform("GF"),
                    str_to_peptidoform("L"),
                    str_to_peptidoform("PKDL"),
                    str_to_peptidoform("VMSRG"),
                ],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::PEPSIN, &test_case);
        }
    }

    #[test]
    fn elastase() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FASGVPRKT"),
                expected_cut_sites: vec![2, 3, 4, 5],
                expected_peptides: vec![
                    str_to_peptidoform("FA"),
                    str_to_peptidoform("S"),
                    str_to_peptidoform("G"),
                    str_to_peptidoform("V"),
                    str_to_peptidoform("PRKT"),
                ],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PGAVSLFTGK"),
                expected_cut_sites: vec![2, 3, 4, 5, 6, 9],
                expected_peptides: vec![
                    str_to_peptidoform("PG"),
                    str_to_peptidoform("A"),
                    str_to_peptidoform("V"),
                    str_to_peptidoform("S"),
                    str_to_peptidoform("L"),
                    str_to_peptidoform("FTG"),
                    str_to_peptidoform("K"),
                ],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::ELASTASE, &test_case);
        }
    }

    #[test]
    fn aspn() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FARDKPGLFD"),
                expected_cut_sites: vec![3, 9],
                expected_peptides: vec![
                    str_to_peptidoform("FAR"),
                    str_to_peptidoform("DKPGLF"),
                    str_to_peptidoform("D"),
                ],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PFKDLTMSR"),
                expected_cut_sites: vec![3],
                expected_peptides: vec![str_to_peptidoform("PFK"), str_to_peptidoform("DLTMSR")],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::ASPN, &test_case);
        }
    }

    #[test]
    fn gluc() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FAREDKPGLF"),
                expected_cut_sites: vec![4],
                expected_peptides: vec![str_to_peptidoform("FARE"), str_to_peptidoform("DKPGLF")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PFKELGTMSR"),
                expected_cut_sites: vec![4],
                expected_peptides: vec![str_to_peptidoform("PFKE"), str_to_peptidoform("LGTMSR")],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::GLUC, &test_case);
        }
    }

    #[test]
    fn lysc() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FARKDPGLF"),
                expected_cut_sites: vec![4],
                expected_peptides: vec![str_to_peptidoform("FARK"), str_to_peptidoform("DPGLF")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PFKDLTKMSR"),
                expected_cut_sites: vec![3, 7],
                expected_peptides: vec![
                    str_to_peptidoform("PFK"),
                    str_to_peptidoform("DLTK"),
                    str_to_peptidoform("MSR"),
                ],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::LYSC, &test_case);
        }
    }

    #[test]
    fn argc() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FARKDPGLF"),
                expected_cut_sites: vec![3],
                expected_peptides: vec![str_to_peptidoform("FAR"), str_to_peptidoform("KDPGLF")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PFKDLRTMSR"),
                expected_cut_sites: vec![6],
                expected_peptides: vec![str_to_peptidoform("PFKDLR"), str_to_peptidoform("TMSR")],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::ARGC, &test_case);
        }
    }
}
