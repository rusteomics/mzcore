//! Fuzz target for parsing mzPAF sequences
use afl::*;
use itertools::Itertools;
use mzannotate::fragment::ToMzPAF;

fn main() {
    fuzz!(|data: &[u8]| {
        if let Ok(data) = std::str::from_utf8(data)
            && let Ok((res, _)) = mzannotate::fragment::Fragment::mz_paf_strict(
                data,
                &mzcore::ontology::STATIC_ONTOLOGIES,
                &[],
            )
        {
            let back = res.iter().map(|a| a.to_mz_paf_string()).join(",");
            let res_back = mzannotate::fragment::Fragment::mz_paf_strict(
                &back,
                &mzcore::ontology::STATIC_ONTOLOGIES,
                &[],
            );
            match res_back {
                Ok((res_back, warnings)) => {
                    let back_back = res_back.iter().map(|a| a.to_mz_paf_string()).join(",");
                    assert_eq!(
                        back, back_back,
                        "{back} != {back_back} (from input: {data})",
                    );
                    if !warnings.is_empty() {
                        println!("{warnings:?}");
                        panic!("Output was not strictly adherent to the standard");
                    }
                }
                Err(err) => {
                    println!("Failed: '{data}' was exported as '{back}'");
                    println!("{err:?}");
                    panic!("Failed test")
                }
            }
        }
    });
}
