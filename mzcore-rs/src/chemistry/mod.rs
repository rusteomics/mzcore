#![allow(dead_code)]

pub mod api;
pub mod amino_acid;
pub mod atom;
pub mod composition;
pub mod constants;
pub mod element;
pub mod glycan;
pub mod isotope;
pub mod ptm;
pub mod peptide;
pub mod table;
pub mod unimod;


#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use anyhow::*;

    use crate::chemistry::api::AminoAcidFactory;
    use crate::chemistry::constants::WATER_MONO_MASS;
    use crate::chemistry::table::proteinogenic_amino_acid_table;

    #[test]
    fn pep_seq_to_amino_acids() -> Result<()> {
        let pep_seq = "INTERSTELLAR";

        let default_aa_table = proteinogenic_amino_acid_table();
        let pep_seq_bytes = &Cow::from(pep_seq.as_bytes());
        let aa_iter = default_aa_table.aa_iter_from_bytes(pep_seq_bytes);

        let mut mono_mass = WATER_MONO_MASS;
        for aa_res in aa_iter {
            let aa = aa_res?;
            mono_mass += aa.mono_mass;
        }

        let expected_mono_mass =  1401.75759215;
        let max_mass_diff = 0.0001;
        assert!( (mono_mass - expected_mono_mass).abs() < max_mass_diff, "monoisotopic mass of '{}' should be {} but got {}", pep_seq, expected_mono_mass, mono_mass);

        Ok(())
    }

}
