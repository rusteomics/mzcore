use std::sync::LazyLock;

use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::sequence::{AminoAcid, SequenceElement};

/// A protease defined by it ability to cut at any site identified by the right amino acids at the n and c terminal.
/// Each position is identified by an option, a none means that there is no specificity at this position. If there is
/// a specificity at a certain position any amino acid that is contained in the set is allowed (see
/// [`crate::CheckedAminoAcid::canonical_identical`]).
///
/// A standard set of proteases can be found here [`known_proteases`].
///
/// # Examples
///
/// ## Basic digestion with Trypsin
/// ```rust
/// # use rustyms::prelude::*;
/// # use rustyms::sequence::known_proteases;
/// // Define a protease and sequence to digest
/// let trypsin = &known_proteases::TRYPSIN;
/// let sequence = Peptidoform::pro_forma("SIADIRGGKSLAIEGCRTKM", None).unwrap().into_linear().unwrap();
///
/// // Run an in-silico digest with 0 missed cleavages and only allow peptides ranging from 4 to 40 amino acids long
/// let digest = sequence.digest(trypsin, 0, 4..40);
///
/// assert_eq!(digest.len(), 2);
/// assert_eq!(digest[0].to_string(), "SIADIR");
/// ```
///
/// ## Finding cut sites in a sequence
/// ```rust
/// # use rustyms::prelude::*;
/// # use rustyms::sequence::known_proteases;
/// // Define a protease and sequence to digest
/// let trypsin = &known_proteases::TRYPSIN;
/// let sequence = Peptidoform::pro_forma("SIADIRGRKM", None).unwrap().into_linear().unwrap();
///
/// // Get all locations where trypsin would cut
/// let cut_sites = trypsin.match_locations(sequence.sequence());
///
/// assert_eq!(cut_sites, vec![6, 8, 9]);
/// ```
///
/// ## Using different proteases
/// ```rust
/// # use rustyms::prelude::*;
/// # use rustyms::sequence::known_proteases;
/// // Define a sequence to digest
/// let sequence = Peptidoform::pro_forma("FAREDKPGLF", None).unwrap().into_linear().unwrap();
///
/// // Compare different digestion patterns
/// let gluc_digest = sequence.digest(&known_proteases::GLUC, 0, 4..40);
/// let lysc_digest = sequence.digest(&known_proteases::LYSC, 0, 4..40);
///
/// assert_eq!(gluc_digest[0].to_string(), "FARE");
/// assert_eq!(lysc_digest[0].to_string(), "FAREDK");
/// ```
///
/// ## Creating a custom protease
/// ```rust
/// # use rustyms::prelude::*;
/// // Define a custom protease that cuts after Histidine (H)
/// let his_protease = Protease::c_terminal_of(vec![AminoAcid::Histidine]);
///
/// // Test it on a sequence
/// let sequence = Peptidoform::pro_forma("AAHFGHKLM", None).unwrap().into_linear().unwrap();
/// let digest = sequence.digest(&his_protease, 0, 1..40);
///
/// assert_eq!(digest.len(), 3);
/// assert_eq!(digest[0].to_string(), "AAH");
/// assert_eq!(digest[1].to_string(), "FGH");
/// assert_eq!(digest[2].to_string(), "KLM");
/// ```
///
/// ## Digestion with missed cleavages
/// ```rust
/// # use rustyms::prelude::*;
/// # use rustyms::sequence::known_proteases;
/// // Define a protease and sequence to digest
/// let trypsin = &known_proteases::TRYPSIN;
/// let sequence = Peptidoform::pro_forma("AARKFGKPLM", None).unwrap().into_linear().unwrap();
///
/// // With 0 missed cleavages
/// let digest_0 = sequence.digest(trypsin, 0, 1..40);
/// assert_eq!(digest_0.len(), 3);
///
/// // With 1 missed cleavage
/// let digest_1 = sequence.digest(trypsin, 1, 1..40);
/// assert_eq!(digest_1.len(), 5); // Original 3 peptides + 2 peptides with 1 missed cleavage
/// ```
///
/// ## Digestion with a complex custom protease
/// ```rust
/// # use rustyms::prelude::*;
/// # use rustyms::sequence::known_proteases;
/// // A protease that:
/// // 1. Requires two phenylalanine residues followed by any amino acid before the cut site (FF*)
/// // 2. And requires any amino acid except proline followed by either histidine or tryptophan after the cut site (*H/W)
/// let custom_protease = Protease {
///     before: vec![
///         Some(vec![AminoAcid::Phenylalanine]),
///         Some(vec![AminoAcid::Phenylalanine]),
///         None
///     ],
///     after: vec![
///         Some(Protease::get_exclusive(&[AminoAcid::Proline])),
///         Some(vec![AminoAcid::Histidine, AminoAcid::Tryptophan])
///     ],
/// };
///
/// // A sequence with multiple potential cut sites:
/// // - Position 10: FF(M)-(A)H - valid cut site
/// // - Position 18: FF(G)-(V)H - valid cut site
/// // - Position 26: FF(A)-(P)H - invalid due to proline after cut site
/// // - Position 34: FF(T)-(S)W - valid cut site
/// let long_sequence = Peptidoform::pro_forma("SLDARETFFMAHKDGFFGVHPDTFFAPHPHYFFTSWNVPG", None)
///     .unwrap()
///     .into_linear()
///     .unwrap();
///
/// // Find all cut sites
/// let cut_sites = custom_protease.match_locations(long_sequence.sequence());
/// assert_eq!(cut_sites, vec![10, 18, 34]);
///
/// // Perform a digestion with 0 missed cleavages
/// let digest = long_sequence.digest(&custom_protease, 0, 1..40);
///
/// assert_eq!(digest.len(), 4);
/// assert_eq!(digest[0].to_string(), "SLDARETFFM");
/// assert_eq!(digest[1].to_string(), "AHKDGFFG");
/// assert_eq!(digest[2].to_string(), "VHPDTFFAPHPHYFFT");
/// assert_eq!(digest[3].to_string(), "SWNVPG");
///
/// // Perform a digestion with 1 missed cleavage
/// let digest_missed = long_sequence.digest(&custom_protease, 1, 1..40);
///
/// // Should include original fragments plus those with 1 missed cleavage
/// assert_eq!(digest_missed.len(), 7);
///
/// // Original fragments
/// assert!(digest_missed.iter().any(|p| p.to_string() == "SLDARETFFM"));
/// assert!(digest_missed.iter().any(|p| p.to_string() == "AHKDGFFG"));
/// assert!(digest_missed.iter().any(|p| p.to_string() == "VHPDTFFAPHPHYFFT"));
/// assert!(digest_missed.iter().any(|p| p.to_string() == "SWNVPG"));
///
/// // Fragments with 1 missed cleavage
/// assert!(digest_missed.iter().any(|p| p.to_string() == "SLDARETFFMAHKDGFFG"));
/// assert!(digest_missed.iter().any(|p| p.to_string() == "AHKDGFFGVHPDTFFAPHPHYFFT"));
/// assert!(digest_missed.iter().any(|p| p.to_string() == "VHPDTFFAPHPHYFFTSWNVPG"));
/// ```
#[derive(Clone, Debug, Deserialize, Eq, Hash, PartialEq, Serialize)]
pub struct Protease {
    /// The amino acids n terminal of the cut site.
    pub before: Vec<Option<Vec<AminoAcid>>>,
    /// The amino acids c terminal of the cut site.
    pub after: Vec<Option<Vec<AminoAcid>>>,
}

impl Protease {
    /// Define a protease that cuts exactly between the specified sequences.
    /// ```rust
    /// # use rustyms::prelude::*;
    /// # let cetuximab = "QVQLKQSGPGLVQPSQSLSITCTVSGFSLTNYGVHWVRQSPGKGLEWLGVIWSGGNTDYNTPFTSRLSINKDNSKSQVFFKMNSLQSNDTAIYYCARALTYYDYEFAYWGQGTLVTVSAASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK";
    /// // FabDELLO is a protease designed to create Fabs from antibodies
    /// // Define the cut site for FabDello "KSCDK / THTCPPCP"
    /// let fab_dello = Protease::between_stretches(
    ///     &"KSCDK".chars().map(|c| AminoAcid::try_from(c).unwrap()).collect::<Vec<_>>(),
    ///     &"THTCPPCP".chars().map(|c| AminoAcid::try_from(c).unwrap()).collect::<Vec<_>>());
    ///
    /// // Get the sequence of the antibody Cetuximab
    /// let sequence = Peptidoform::pro_forma(cetuximab, None).unwrap().into_linear().unwrap();
    ///
    /// // Run an in-silico digest with 0 missed cleavages
    /// let digest = sequence.digest(&fab_dello, 0, ..);
    ///
    /// assert_eq!(digest.len(), 2);
    /// ```
    pub fn between_stretches(before: &[AminoAcid], after: &[AminoAcid]) -> Self {
        Self {
            before: before.iter().map(|aa| Some(vec![*aa])).collect_vec(),
            after: after.iter().map(|aa| Some(vec![*aa])).collect_vec(),
        }
    }

    /// Define a protease that cuts exactly between the specified options before the site and the specified options after the site.
    /// ```rust
    /// # use rustyms::prelude::*;
    /// // Define a protease and sequence to digest
    /// let chymotrypsin = Protease::between_options(
    ///        vec![
    ///            AminoAcid::Phenylalanine,
    ///            AminoAcid::Tryptophan,
    ///            AminoAcid::Tyrosine,
    ///        ],
    ///        Protease::get_exclusive(&[AminoAcid::Proline]),
    ///    );
    /// let sequence = Peptidoform::pro_forma("AFWYPLGF", None).unwrap().into_linear().unwrap();
    ///
    /// // Run an in-silico digest with 0 missed cleavages and only allow peptides ranging from 4 to 40 amino acids long
    /// let digest = sequence.digest(&chymotrypsin, 0, 4..40);
    ///
    /// assert_eq!(digest.len(), 1);
    /// assert_eq!(digest[0].to_string(), "YPLGF");
    /// ```
    pub fn between_options(before: Vec<AminoAcid>, after: Vec<AminoAcid>) -> Self {
        Self {
            before: vec![Some(before)],
            after: vec![Some(after)],
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
        Protease::between_options(
            vec![AminoAcid::Lysine, AminoAcid::Arginine],
            Protease::get_exclusive(&[AminoAcid::Proline]),
        )
    });

    /// `Chymotrypsin` cuts after Phenylalanine (F), Tryptophan (W), Tyrosine (Y), unless followed by Proline (P)
    pub static CHYMOTRYPSIN: LazyLock<Protease> = LazyLock::new(|| {
        Protease::between_options(
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
    use crate::sequence::{Linear, Peptidoform};

    use super::*;

    pub(super) struct ProteaseTestCase {
        pub sequence: Peptidoform<Linear>,
        pub expected_cut_sites: Vec<usize>,
        pub expected_peptides: Vec<Peptidoform<Linear>>,
    }

    /// Generic test function for all proteases
    pub(super) fn test_protease(protease: &Protease, test_case: &ProteaseTestCase) {
        // Test cut sites
        let cut_sites = protease.match_locations(test_case.sequence.sequence());

        assert_eq!(
            cut_sites, test_case.expected_cut_sites,
            "Incorrect cut sites: found '{cut_sites:?}' expected '{:?}'",
            test_case.expected_cut_sites
        );

        // Test peptides
        let peptides = test_case.sequence.digest(protease, 0, 4..40);

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
                expected_peptides: vec![str_to_peptidoform("RPGK")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("ARAKGCVLRPKDGR"),
                expected_cut_sites: vec![2, 4, 11],
                expected_peptides: vec![str_to_peptidoform("GCVLRPK")],
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
                expected_peptides: vec![str_to_peptidoform("YPLGF")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("AVFUDGWTYPMSR"),
                expected_cut_sites: vec![3, 7],
                expected_peptides: vec![str_to_peptidoform("UDGW"), str_to_peptidoform("TYPMSR")],
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
                expected_peptides: vec![str_to_peptidoform("AACVF"), str_to_peptidoform("PAKL")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("GFLPKDLVMSRG"),
                expected_cut_sites: vec![2, 3, 7],
                expected_peptides: vec![str_to_peptidoform("PKDL"), str_to_peptidoform("VMSRG")],
            },
        ];

        for test_case in test_cases {
            test_protease(&known_proteases::PEPSIN, &test_case);
        }
    }

    #[test]
    fn aspn() {
        let test_cases = vec![
            ProteaseTestCase {
                sequence: str_to_peptidoform("FARDKPGLFD"),
                expected_cut_sites: vec![3, 9],
                expected_peptides: vec![str_to_peptidoform("DKPGLF")],
            },
            ProteaseTestCase {
                sequence: str_to_peptidoform("PFKDLTMSR"),
                expected_cut_sites: vec![3],
                expected_peptides: vec![str_to_peptidoform("DLTMSR")],
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
                expected_peptides: vec![str_to_peptidoform("DLTK")],
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
                expected_peptides: vec![str_to_peptidoform("KDPGLF")],
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
