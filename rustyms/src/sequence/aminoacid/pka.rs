//! Module used to store and calculate pKa and isoelectric point values for a given [`AminoAcid`] or [Peptidoform] respectively

use serde::{Deserialize, Serialize};

use crate::sequence::{
    AminoAcid, AtMax, Peptidoform, SemiAmbiguous, SimpleModificationInner, properties::ChargeClass,
};

use super::is_amino_acid::IsAminoAcid;

/// A source for pKa values, which can be used to calculate the pKa for peptidoforms.
pub trait PKaSource<AA: IsAminoAcid> {
    /// Get the pKa values for the given amino acid and modifications.
    #[allow(non_snake_case)]
    fn pKa(
        amino_acid: AA,
        side_chain_modifications: impl Iterator<Item = impl AsRef<SimpleModificationInner>>,
        n_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
        c_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
    ) -> Option<AminoAcidPKa>;
}

impl<Complexity: AtMax<SemiAmbiguous>> Peptidoform<Complexity> {
    /// Get the calculated isoelectric point (pI) for the peptidoform, or None if any sequence elements lack pKa values.
    ///
    /// The isoelectric point is the pH at which the net charge of the peptidoform is zero. This is determined using a binary
    /// search between pH 0 and 14. The charge at each pH is computed using the Henderson-Hasselbalch equation with pKa values
    /// from the provided `PKaSource`, considering N-terminal, C-terminal, and sidechain ionizable groups.
    ///
    /// # Example
    /// ```rust
    /// # use rustyms::sequence::{Peptidoform, pka::{PKaSource, PKaLide1991}};
    /// // Create a SemiAmbiguous Peptidoform for glutamic acid (E) and Alanine (A)
    /// let peptidoform = Peptidoform::pro_forma(&"EMEVEESPEK", None).unwrap().into_semi_ambiguous().unwrap();
    /// let pi = peptidoform.isoelectic_point::<PKaLide1991>();
    /// // The calculated pI is approximately 3.57 based on Lide 1991 pKa values
    /// assert_eq!(pi.map(|v| (v * 100.0).round() / 100.0), Some(3.57));
    /// ```
    ///
    /// # Shortcomings
    /// - **Naive Approach**: Does not account for interactions between ionizable groups.
    /// - **Modifications Ignored**: Modifications affecting pKa are not considered.
    /// - **Environmental Factors**: Assumes pKa values are independent of sequence and environment.
    ///
    /// Get the calculated pKa value for the given peptidoform, or None if any of the sequence elements do not have a defined pKa.
    #[allow(non_snake_case)]
    pub fn isoelectic_point<Source: PKaSource<AminoAcid>>(&self) -> Option<f64> {
        const EPSILON: f64 = 0.0001;

        let sequence = self.sequence();
        if sequence.is_empty() {
            return None;
        }

        // Collect all ionizable groups with their pKa values
        let mut ionizable = Vec::with_capacity(sequence.len() + 2);

        // Handle N-terminal
        let first = sequence.first()?;
        ionizable.push((
            ChargeClass::Positive,
            Source::pKa(
                first.aminoacid.aminoacid(),
                first.modifications.iter().filter_map(|m| m.simple()),
                Some(self.get_n_term().iter().filter_map(|m| m.simple())),
                (self.len() == 1).then_some(self.get_c_term().iter().filter_map(|m| m.simple())),
            )?
            .n_term(),
        )); // N-terminal is always positive

        // Handle C-terminal
        let last = sequence.last()?;
        ionizable.push((
            ChargeClass::Negative,
            Source::pKa(
                last.aminoacid.aminoacid(),
                last.modifications.iter().filter_map(|m| m.simple()),
                (self.len() == 1).then_some(self.get_n_term().iter().filter_map(|m| m.simple())),
                Some(self.get_c_term().iter().filter_map(|m| m.simple())),
            )?
            .c_term(),
        )); // C-terminal is always negative

        // Handle sidechains
        for (index, aa) in sequence.iter().enumerate() {
            if let Some(sidechain) = Source::pKa(
                aa.aminoacid.aminoacid(),
                aa.modifications.iter().filter_map(|m| m.simple()),
                (index == 0).then_some(self.get_n_term().iter().filter_map(|m| m.simple())),
                (index == self.len() - 1)
                    .then_some(self.get_n_term().iter().filter_map(|m| m.simple())),
            )?
            .sidechain()
            {
                let charge_class = aa.aminoacid.aminoacid().charge_class();
                match charge_class {
                    ChargeClass::Positive | ChargeClass::Negative => {
                        ionizable.push((charge_class, sidechain));
                    }
                    ChargeClass::Unknown => return None,
                    ChargeClass::Uncharged => (),
                }
            }
        }

        // Binary search between pH 0-14 to find isoelectric point
        let mut low = 0.0;
        let mut high = 14.0;
        let mut new_pi = 7.775;

        #[allow(clippy::while_float)]
        while (high - low) > EPSILON {
            new_pi = (low + high) / 2.0;
            let charge = calculate_charge(new_pi, &ionizable);

            if charge > 0.0 {
                low = new_pi;
            } else {
                high = new_pi;
            }
        }

        Some(new_pi)
    }
}

#[allow(non_snake_case)]
fn calculate_charge(pH: f64, ionizable: &[(ChargeClass, f64)]) -> f64 {
    let mut charge = 0.0;

    for (class, pka) in ionizable {
        match class {
            ChargeClass::Positive => charge += 1.0 / (10.0_f64.powf(pH - pka) + 1.0),
            ChargeClass::Negative => charge -= 1.0 / (10.0_f64.powf(pka - pH) + 1.0),
            _ => {}
        }
    }

    charge
}
/// The pKa for a specific Amino Acid
#[derive(Clone, Copy, Debug, Deserialize, PartialEq, PartialOrd, Serialize)]
pub struct AminoAcidPKa {
    n_term: f64,
    sidechain: Option<f64>,
    c_term: f64,
}

impl AminoAcidPKa {
    const fn new(n_term: f64, sidechain: Option<f64>, c_term: f64) -> Self {
        Self {
            n_term,
            sidechain,
            c_term,
        }
    }

    /// Get the pKa value for the n-term of the Amino acid
    pub const fn n_term(self) -> f64 {
        self.n_term
    }

    /// Get the pKa value for the side-chain group of the Amino acid
    pub const fn sidechain(self) -> Option<f64> {
        self.sidechain
    }

    /// Get the pKa value for the c-term of the Amino acid
    pub const fn c_term(self) -> f64 {
        self.c_term
    }
}

/// pKa values from Lide, D. R. (1991). Handbook of Chemistry and Physics: A Ready Reference Book of Chemical and Physical Data.
#[derive(Clone, Copy, Debug)]
pub struct PKaLide1991;

impl PKaSource<AminoAcid> for PKaLide1991 {
    fn pKa(
        amino_acid: AminoAcid,
        mut side_chain_modifications: impl Iterator<Item = impl AsRef<SimpleModificationInner>>,
        n_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
        c_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
    ) -> Option<AminoAcidPKa> {
        if side_chain_modifications.next().is_some()
            || n_terminal_modifications.is_some_and(|mut m| m.next().is_some())
            || c_terminal_modifications.is_some_and(|mut m| m.next().is_some())
        {
            return None;
        }
        match amino_acid {
            AminoAcid::Arginine => Some(AminoAcidPKa::new(9.00, Some(12.10), 2.03)),
            AminoAcid::Histidine => Some(AminoAcidPKa::new(9.09, Some(6.04), 1.70)),
            AminoAcid::Lysine => Some(AminoAcidPKa::new(9.16, Some(10.67), 2.15)),
            AminoAcid::AsparticAcid => Some(AminoAcidPKa::new(9.66, Some(3.71), 1.95)),
            AminoAcid::GlutamicAcid => Some(AminoAcidPKa::new(9.58, Some(4.15), 2.16)),
            AminoAcid::Tyrosine => Some(AminoAcidPKa::new(9.04, Some(10.10), 2.24)),
            AminoAcid::Cysteine => Some(AminoAcidPKa::new(10.28, Some(8.14), 1.91)),
            AminoAcid::Alanine => Some(AminoAcidPKa::new(9.71, None, 2.33)),
            AminoAcid::Glycine => Some(AminoAcidPKa::new(9.58, None, 2.34)),
            AminoAcid::Proline => Some(AminoAcidPKa::new(10.47, None, 1.95)),
            AminoAcid::Serine => Some(AminoAcidPKa::new(9.05, None, 2.13)),
            AminoAcid::Threonine => Some(AminoAcidPKa::new(8.96, None, 2.20)),
            AminoAcid::Methionine => Some(AminoAcidPKa::new(9.08, None, 2.16)),
            AminoAcid::Phenylalanine => Some(AminoAcidPKa::new(9.09, None, 2.18)),
            AminoAcid::Tryptophan => Some(AminoAcidPKa::new(9.34, None, 2.38)),
            AminoAcid::Valine => Some(AminoAcidPKa::new(9.52, None, 2.27)),
            AminoAcid::Isoleucine => Some(AminoAcidPKa::new(9.60, None, 2.26)),
            AminoAcid::Leucine => Some(AminoAcidPKa::new(9.58, None, 2.32)),
            AminoAcid::Glutamine => Some(AminoAcidPKa::new(9.00, None, 2.18)),
            AminoAcid::Asparagine => Some(AminoAcidPKa::new(8.73, None, 2.16)),
            _ => None,
        }
    }
}

/// pKa values from Lehninger, A. L., Nelson, D. L., & Cox, M. M. (2005). Lehninger Principles of Biochemistry. Macmillan.
#[derive(Clone, Copy, Debug)]
pub struct PKaLehninger;

impl PKaSource<AminoAcid> for PKaLehninger {
    fn pKa(
        amino_acid: AminoAcid,
        mut side_chain_modifications: impl Iterator<Item = impl AsRef<SimpleModificationInner>>,
        n_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
        c_terminal_modifications: Option<impl Iterator<Item = impl AsRef<SimpleModificationInner>>>,
    ) -> Option<AminoAcidPKa> {
        if side_chain_modifications.next().is_some()
            || n_terminal_modifications.is_some_and(|mut m| m.next().is_some())
            || c_terminal_modifications.is_some_and(|mut m| m.next().is_some())
        {
            return None;
        }
        match amino_acid {
            AminoAcid::Arginine => Some(AminoAcidPKa::new(9.04, Some(12.48), 2.17)),
            AminoAcid::Histidine => Some(AminoAcidPKa::new(9.17, Some(6.00), 1.82)),
            AminoAcid::Lysine => Some(AminoAcidPKa::new(8.95, Some(10.53), 2.18)),
            AminoAcid::AsparticAcid => Some(AminoAcidPKa::new(9.60, Some(3.65), 1.88)),
            AminoAcid::GlutamicAcid => Some(AminoAcidPKa::new(9.67, Some(4.25), 2.19)),
            AminoAcid::Tyrosine => Some(AminoAcidPKa::new(9.11, Some(10.07), 2.20)),
            AminoAcid::Cysteine => Some(AminoAcidPKa::new(10.28, Some(8.18), 1.96)),
            AminoAcid::Alanine => Some(AminoAcidPKa::new(9.69, None, 2.34)),
            AminoAcid::Glycine => Some(AminoAcidPKa::new(9.60, None, 2.34)),
            AminoAcid::Proline => Some(AminoAcidPKa::new(10.96, None, 1.99)),
            AminoAcid::Serine => Some(AminoAcidPKa::new(9.15, None, 2.21)),
            AminoAcid::Threonine => Some(AminoAcidPKa::new(9.62, None, 2.11)),
            AminoAcid::Methionine => Some(AminoAcidPKa::new(9.21, None, 2.28)),
            AminoAcid::Phenylalanine => Some(AminoAcidPKa::new(9.13, None, 1.83)),
            AminoAcid::Tryptophan => Some(AminoAcidPKa::new(9.39, None, 2.38)),
            AminoAcid::Valine => Some(AminoAcidPKa::new(9.62, None, 2.32)),
            AminoAcid::Isoleucine => Some(AminoAcidPKa::new(9.68, None, 2.36)),
            AminoAcid::Leucine => Some(AminoAcidPKa::new(9.60, None, 2.36)),
            AminoAcid::Glutamine => Some(AminoAcidPKa::new(9.13, None, 2.17)),
            AminoAcid::Asparagine => Some(AminoAcidPKa::new(8.80, None, 2.02)),
            _ => None,
        }
    }
}

#[cfg(test)]
#[expect(clippy::float_cmp, clippy::missing_panics_doc)]
mod tests {
    use super::*;
    use crate::sequence::SimpleModification;

    // Helper to create a Peptidoform from a list of amino acids
    fn create_peptidoform(aas: &str) -> Peptidoform<SemiAmbiguous> {
        Peptidoform::pro_forma(aas, None)
            .unwrap()
            .into_semi_ambiguous()
            .unwrap()
    }

    // Helper function to test pKa values for a given source
    fn test_pka<Source: PKaSource<AminoAcid>>(
        test_cases: &[(AminoAcid, Option<(f64, Option<f64>, f64)>)],
    ) {
        for (aa, maybe_values) in test_cases {
            if let Some((n_term, sidechain, c_term)) = maybe_values {
                let pka = Source::pKa(
                    *aa,
                    std::iter::empty::<SimpleModification>(),
                    None::<std::iter::Empty<SimpleModification>>,
                    None::<std::iter::Empty<SimpleModification>>,
                )
                .unwrap_or_else(|| panic!("Missing pKa for {aa:?}"));
                let round = |v: f64| (v * 100.0).round() / 100.0;

                assert_eq!(round(pka.n_term()), *n_term, "N-term mismatch for {aa:?}");
                assert_eq!(
                    pka.sidechain().map(round),
                    *sidechain,
                    "Sidechain mismatch for {aa:?}"
                );
                assert_eq!(round(pka.c_term()), *c_term, "C-term mismatch for {aa:?}");
            } else {
                assert!(maybe_values.is_none(), "Expected None for {aa:?}");
            }
        }
    }

    // Helper function to test an isoelectric point value given a source
    fn test_isoelectric_point<Source: PKaSource<AminoAcid>>(cases: &[(&str, Option<f64>)]) {
        for &(seq, expected) in cases {
            let peptide = create_peptidoform(seq);
            let round = |v: f64| (v * 100.0).round() / 100.0;
            let iso = peptide.isoelectic_point::<Source>();
            assert_eq!(
                iso.map(round),
                expected,
                "Isoelectric point mismatch for peptide: {seq}"
            );
        }
    }

    #[test]
    fn test_pka_lide1991() {
        let test_cases = [
            (AminoAcid::Arginine, Some((9.00, Some(12.10), 2.03))),
            (AminoAcid::GlutamicAcid, Some((9.58, Some(4.15), 2.16))),
            (AminoAcid::Alanine, Some((9.71, None, 2.33))),
            (AminoAcid::Histidine, Some((9.09, Some(6.04), 1.70))),
            (AminoAcid::Unknown, None),
        ];

        test_pka::<PKaLide1991>(&test_cases);
    }

    #[test]
    fn test_pka_lehninger() {
        let test_cases = [
            (AminoAcid::Cysteine, Some((10.28, Some(8.18), 1.96))),
            (AminoAcid::AsparticAcid, Some((9.60, Some(3.65), 1.88))),
            (AminoAcid::Isoleucine, Some((9.68, None, 2.36))),
            (AminoAcid::Tryptophan, Some((9.39, None, 2.38))),
            (AminoAcid::Selenocysteine, None),
        ];

        test_pka::<PKaLehninger>(&test_cases);
    }

    #[test]
    fn test_isoelectric_point_lide1991() {
        let test_cases = [
            ("E", Some(3.16)),
            ("A", Some(6.02)),
            ("DE", Some(2.85)),
            ("HR", Some(10.6)),
            ("KDEH", Some(5.17)),
            ("AXRT", None),
            ("AXRT[Oxidation]", None),
        ];

        test_isoelectric_point::<PKaLide1991>(&test_cases);
    }

    #[test]
    fn test_isoelectric_point_lehninger() {
        let test_cases = [
            ("G", Some(5.97)),
            ("Y", Some(5.65)),
            ("CQ", Some(6.23)),
            ("KP", Some(9.74)),
            ("FIVS", Some(5.67)),
            ("TKLB", None),
            ("TK[Oxidation]LB", None),
        ];

        test_isoelectric_point::<PKaLehninger>(&test_cases);
    }
}
