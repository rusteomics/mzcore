use serde::{Deserialize, Serialize};

use crate::{
    aminoacid_properties::ChargeClass, aminoacids::IsAminoAcid,
    modification::SimpleModificationInner, AminoAcid, AtMax, Peptidoform, SemiAmbiguous,
};

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
    /// Get the calculated pKa value for the given peptidoform, or None if any of the sequence elements do not have a defined pKa.
    #[allow(non_snake_case)]
    pub fn isoelectic_point<Source: PKaSource<AminoAcid>>(&self) -> Option<f64> {
        let sequence = self.sequence();
        if sequence.is_empty() {
            return None;
        }

        // Collect all ionizable groups with their pKa values
        let mut ionizable = Vec::with_capacity(sequence.len() + 2);

        // N-terminal
        let first = sequence.first().unwrap();
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

        // C-terminal
        let last = sequence.last().unwrap();
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
        const EPSILON: f64 = 0.0001;

        while (high - low) > EPSILON {
            new_pi = (low + high) / 2.0;
            let charge = calculate_charge(new_pi, &ionizable);

            if charge > 0.0 {
                low = new_pi;
            } else {
                high = new_pi;
            }
        }

        Some((new_pi * 100.0).round() / 100.0)
    }
}

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
#[derive(Copy, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
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
