use anyhow::*;
use std::borrow::Cow;
use std::slice::Iter;

use crate::chemistry::api::*;
use crate::chemistry::amino_acid::AminoAcidDefinition;
use crate::chemistry::constants::{WATER_MONO_MASS, WATER_AVERAGE_MASS};
use crate::chemistry::peptide::LinearPeptide;
use crate::chemistry::table::AminoAcidTable;
use crate::ms::MassType;

impl IsAminoAcidSeq for &str {
    fn amino_acids_as_bytes(&self) -> Iter<u8> {
        self.as_bytes().iter()
    }

    fn length(&self) -> usize {
        self.len()
    }
}

impl IsAminoAcidSeq for &Cow<'_, [u8]> {
    fn amino_acids_as_bytes(&self) -> Iter<u8> {
        self.iter()
    }

    fn length(&self) -> usize {
        self.len()
    }
}

pub trait AASeqMassCalc<AA: IsAminoAcid, S: IsAminoAcidSeq>: AminoAcidFactory<AA> {
    fn calc_mass_from_aa_seq(&self, aa_seq: S, mass_type: MassType) -> Result<f64> {
        let mut seq_mass = 0.0;
        for aa in aa_seq.amino_acids_as_bytes() {
            seq_mass += match mass_type {
                MassType::Monoisotopic => self.aa_from_byte(aa)?.mono_mass(),
                MassType::Average => self.aa_from_byte(aa)?.average_mass().ok_or(anyhow!("undefined average mass"))?,
            };
        }

        let water_mass = match mass_type {
            MassType::Monoisotopic => WATER_MONO_MASS,
            MassType::Average => WATER_AVERAGE_MASS,
        };

        Ok(seq_mass + water_mass)
    }
}

impl AASeqMassCalc<AminoAcidDefinition, &Cow<'_, [u8]>> for AminoAcidTable {}
impl AASeqMassCalc<AminoAcidDefinition, LinearPeptide> for AminoAcidTable {}
impl AASeqMassCalc<AminoAcidDefinition, &str> for AminoAcidTable {}

/*
impl AASeqMassCalc<AminoAcidDefinition> for &str {
    fn calc_mass_from_aa_str(&self, mass_type: MassType) -> Result<f64> {
        self.chars().calc_mass_from_aa_chars(aa_table, mass_type)
    }
}

pub trait SeqCharsMassCalc {
    fn calc_mass_from_aa_chars(&self, aa_table: &AminoAcidTable, mass_type: MassType) -> Result<f64>;
}

impl<I> SeqCharsMassCalc for I
where
    I: Iterator<Item = char>,
    I: Clone,
{
    fn calc_mass_from_aa_chars(&self, aa_table: &AminoAcidTable, mass_type: MassType) -> Result<f64> {

        let aa_by_code1 = &aa_table.aa_by_code1;

        let seq_mass = self.clone().try_fold(0.0, |acc, aa_as_char| {
            let aa_opt = aa_by_code1.get(&aa_as_char);
            let aa = aa_opt.ok_or(anyhow!(
                "can't find amino acid '{}' in the provided table",
                aa_as_char
            ))?;

            let mass = match mass_type {
                MassType::Monoisotopic => aa.mono_mass,
                MassType::Average => aa.average_mass,
            };

            Ok(acc + mass)
        })?;

        let water_mass = match mass_type {
            MassType::Monoisotopic => WATER_MONO_MASS,
            MassType::Average => WATER_AVERAGE_MASS,
        };

        Ok(seq_mass + water_mass)
    }
}*/


#[cfg(test)]
mod tests {

    use super::*;
    use crate::chemistry::table::proteinogenic_amino_acid_table;

    const MAX_MASS_DIFF: f64 = 0.0001;
    
    #[test]
    fn cal_pep_seq_mass() -> Result<()> {
        let pep_seq = "INTERSTELLAR";

        let default_aa_table = proteinogenic_amino_acid_table();

        let mono_mass_from_table = default_aa_table.calc_mass_from_aa_seq(
            pep_seq,
            MassType::Monoisotopic
        )?;
        
        let expected_mono_mass = 1401.75759215;
        assert!( (mono_mass_from_table - expected_mono_mass).abs() < MAX_MASS_DIFF, "monoisotopic mass of '{}' should be {}", pep_seq, mono_mass_from_table);

        let mono_mass_from_seq = pep_seq.amino_acids_as_bytes().fold(0.0, |mass_sum, aa_byte| mass_sum + aa_byte.mono_mass()) + WATER_MONO_MASS;
        assert!( (mono_mass_from_table - mono_mass_from_seq).abs() < MAX_MASS_DIFF, "monoisotopic mass differs between '{}' and {}", mono_mass_from_table, mono_mass_from_seq);

        Ok(())
    }
}