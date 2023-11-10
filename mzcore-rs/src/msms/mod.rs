pub mod annotator;
pub mod fragmentation;
pub mod model;

#[cfg(test)]
mod tests {
    use std::borrow::Cow;
    use anyhow::*;
    use crate::chemistry::table::*;
    use crate::msms::fragmentation::FragmentationTableFactory;
    use crate::msms::model::FragmentIonSeries::{b,y};

    #[test]
    fn pep_seq_to_frag_table() -> Result<()> {
        let pep_seq = "INTERSTELLAR";

        let default_aa_table = proteinogenic_amino_acid_table();

        let frag_table = default_aa_table.compute_frag_table_without_mods(
            &Cow::from(pep_seq.as_bytes()),
            &[b,y],
            &vec![1],
        )?;

        let expected_frag_table_len = 2;
        assert_eq!(frag_table.len(), expected_frag_table_len, "unexpected frag table length got {}, expected {}", frag_table.len(), expected_frag_table_len);

        Ok(())
    }


}
