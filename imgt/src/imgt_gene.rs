use std::collections::HashMap;
use std::fmt::Display;

use crate::{AnnotatedSequence, Gene};
use itertools::Itertools;
use mzcore::sequence::{AminoAcid, Annotation, CheckedAminoAcid, Region};

use crate::structs::{Location, SequenceRegion, SingleSeq};

#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct IMGTGene {
    pub acc: String,
    pub key: String,
    pub location: Location,
    pub allele: String,
    pub regions: HashMap<String, crate::structs::Region>,
}

impl IMGTGene {
    pub(crate) fn finish(self) -> Result<SingleSeq, String> {
        let (regions, additional_annotations) = self.get_regions()?;

        let sequence: Vec<AminoAcid> = regions.iter().flat_map(|reg| reg.1.0.clone()).collect();
        let dna: String = regions.iter().map(|reg| reg.1.2.clone()).collect();
        let region_lengths = regions
            .iter()
            .map(|reg| (reg.0.clone(), reg.1.0.len()))
            .collect();
        let conserved_map = HashMap::from([
            ("1st-CYS", Annotation::Conserved),
            ("2nd-CYS", Annotation::Conserved),
            ("CONSERVED-TRP", Annotation::Conserved),
            ("J-PHE", Annotation::Conserved),
            ("J-TRP", Annotation::Conserved),
        ]);
        let mut conserved = self
            .regions
            .iter()
            .filter(|(key, _)| {
                ["1st-CYS", "2nd-CYS", "CONSERVED-TRP", "J-PHE", "J-TRP"].contains(&key.as_str())
            })
            .map(|(key, region)| {
                region
                    .location
                    .find_aa_location(&regions)
                    .map(|index| (conserved_map[key.as_str()].clone(), index))
                    .ok_or_else(|| format!("Cannot find location of '{key}' '{region}'"))
            })
            .collect::<Result<Vec<_>, _>>()?;
        conserved.extend(
            find_possible_n_glycan_locations(&sequence)
                .iter()
                .map(|i| (Annotation::NGlycan, *i)),
        );
        conserved.extend(additional_annotations);
        let (name, allele) = Gene::from_imgt_name_with_allele(self.allele.as_str())?;
        Ok(SingleSeq {
            name,
            allele,
            acc: self.acc,
            sequence: AnnotatedSequence::new(
                sequence
                    .iter()
                    .copied()
                    .map(CheckedAminoAcid::new)
                    .map(|p| p.into_unambiguous().unwrap())
                    .collect(),
                region_lengths,
                conserved,
            ),
            dna,
        })
    }

    fn get_regions(&self) -> Result<(Vec<SequenceRegion>, Vec<(Annotation, usize)>), String> {
        let mut additional_annotations = Vec::new();

        let regions = match self.key.as_str() {
            "V-GENE" => {
                vec![
                    self.get_region(&Region::Framework(1), "FR1-IMGT")?,
                    self.get_region(&Region::ComplementarityDetermining(1), "CDR1-IMGT")?,
                    self.get_region(&Region::Framework(2), "FR2-IMGT")?,
                    self.get_region(&Region::ComplementarityDetermining(2), "CDR2-IMGT")?,
                    self.get_region(&Region::Framework(3), "FR3-IMGT")?,
                    self.get_region(&Region::ComplementarityDetermining(3), "CDR3-IMGT")?,
                ]
            }
            "J-GENE" => {
                let j = self.get_region(&Region::Framework(4), "J-REGION")?;
                let motif = j.1.0.iter().tuple_windows().position(|(a, b, _, d)| {
                    (*a == AminoAcid::Tryptophan || *a == AminoAcid::Phenylalanine)
                        && *b == AminoAcid::Glycine
                        && *d == AminoAcid::Glycine
                });
                if let Some(motif_start) = motif {
                    let j = fix_j(j.1, motif_start);
                    additional_annotations.extend(j.1);
                    j.0
                } else {
                    vec![j]
                }
            }
            "D-GENE" => {
                vec![self.get_region(&Region::ComplementarityDetermining(3), "D-REGION")?]
            }
            "C-GENE" => {
                let mut seq = Vec::new();

                // Heavy chain
                seq.extend(self.get_region(&Region::ConstantHeavy(1), "CH1").ok());
                // Try to detect the best H/CH2
                if self.regions.contains_key("H") && self.regions.contains_key("CH2") {
                    seq.extend(self.get_region(&Region::Hinge(None), "H").ok());
                    seq.extend(self.get_region(&Region::ConstantHeavy(2), "CH2").ok());
                } else if self.regions.contains_key("H-CH2") {
                    seq.extend(
                        self.get_region(
                            &Region::Joined(vec![Region::Hinge(None), Region::ConstantHeavy(2)]),
                            "H-CH2",
                        )
                        .ok(),
                    );
                } else {
                    seq.extend(self.get_region(&Region::Hinge(Some(1)), "H1").ok());
                    seq.extend(self.get_region(&Region::Hinge(Some(2)), "H2").ok());
                    seq.extend(self.get_region(&Region::Hinge(Some(3)), "H3").ok());
                    seq.extend(self.get_region(&Region::Hinge(Some(4)), "H4").ok());
                    seq.extend(self.get_region(&Region::ConstantHeavy(2), "CH2").ok());
                }
                let mut secretory_added = false;

                for (region, with_chs) in [
                    (
                        Region::ConstantHeavy(3),
                        Region::Joined(vec![Region::ConstantHeavy(3), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(4),
                        Region::Joined(vec![Region::ConstantHeavy(4), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(5),
                        Region::Joined(vec![Region::ConstantHeavy(5), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(6),
                        Region::Joined(vec![Region::ConstantHeavy(6), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(7),
                        Region::Joined(vec![Region::ConstantHeavy(7), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(8),
                        Region::Joined(vec![Region::ConstantHeavy(8), Region::SecratoryTail]),
                    ),
                    (
                        Region::ConstantHeavy(9),
                        Region::Joined(vec![Region::ConstantHeavy(9), Region::SecratoryTail]),
                    ),
                ] {
                    if self.regions.contains_key(&region.to_string())
                        && self.regions.contains_key("CHS")
                    {
                        // If this region is stored separately and the CHS is available
                        seq.extend(self.get_region(&region, &region.to_string()).ok());
                    } else if self.regions.contains_key(&with_chs.to_string()) {
                        // If this region is stored together with the CHS
                        seq.extend(self.get_region(&with_chs, &with_chs.to_string()).ok());
                        secretory_added = true;
                    } else {
                        // If this region is stored separately but some later region is combined with CHS
                        seq.extend(self.get_region(&region, &region.to_string()).ok());
                    }
                }
                if !secretory_added {
                    seq.extend(self.get_region(&Region::SecratoryTail, "CHS").ok());
                }
                // possibly_add(Region::M, "M")?; // TODO: Potentially find a way to also allow the membrane bound version to be stored
                // possibly_add(Region::M1, "M1")?;
                // possibly_add(Region::M2, "M2")?;

                // Otherwise assume light chain
                if seq.is_empty() {
                    seq.extend(self.get_region(&Region::ConstantLight, "CL").ok());
                }
                if seq.is_empty() {
                    seq.extend(self.get_region(&Region::ConstantLight, "C-REGION").ok());
                }

                if seq.is_empty() {
                    return Err("Empty C sequence".to_string());
                }
                seq
            }
            _ => Vec::new(),
        };

        Ok((regions, additional_annotations))
    }

    fn get_region(&self, region: &Region, key: &str) -> Result<SequenceRegion, String> {
        self.regions
            .get(key)
            .ok_or_else(|| format!("Could not find {key}"))
            .and_then(|region| {
                region
                    .found_seq
                    .as_ref()
                    .map(|seq| {
                        let mut final_seq = region
                            .splice_aa
                            .map(|aa| vec![aa])
                            .filter(|_| region.shift != 2)
                            .unwrap_or_default();
                        final_seq.extend(seq.1.0.clone());
                        (final_seq, region.location.clone(), seq.0.clone())
                    })
                    .map_err(ToOwned::to_owned)
            })
            .map(|res| (region.clone(), res))
    }
}

impl Display for IMGTGene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}\t{}\t{}", self.key, self.location, self.allele)?;
        for region in self.regions.values().sorted_by_key(|reg| &reg.key) {
            writeln!(f, "  R {region}")?;
        }
        Ok(())
    }
}

fn find_possible_n_glycan_locations(sequence: &[AminoAcid]) -> Vec<usize> {
    let mut result = Vec::new();
    for (index, aa) in sequence.windows(3).enumerate() {
        if let (AminoAcid::Asparagine, AminoAcid::Serine | AminoAcid::Threonine) = (aa[0], aa[2])
            && aa[1] != AminoAcid::Proline
        {
            result.push(index);
        }
    }
    result
}

fn fix_j(
    j: (Vec<AminoAcid>, Location, String),
    cdr3_length: usize,
) -> (Vec<SequenceRegion>, Vec<(Annotation, usize)>) {
    let (cdr3_loc, fr4_loc) =
        j.1.splice(cdr3_length)
            .expect("CDR3 should fit in full FR4 of J gene");
    let cdr3 = (
        j.0[..cdr3_length].to_vec(),
        cdr3_loc,
        j.2[..cdr3_length].to_owned(),
    );
    let fr4 = (
        j.0[cdr3_length..].to_vec(),
        fr4_loc,
        j.2[cdr3_length..].to_owned(),
    );

    let mut annotations = Vec::new();
    if fr4.0[0] == AminoAcid::Tryptophan || fr4.0[0] == AminoAcid::Phenylalanine {
        annotations.push((Annotation::Conserved, cdr3_length));
    }
    if fr4.0[1] == AminoAcid::Glycine {
        annotations.push((Annotation::Conserved, cdr3_length + 1));
    }
    if fr4.0[3] == AminoAcid::Glycine {
        annotations.push((Annotation::Conserved, cdr3_length + 3));
    }

    (
        vec![
            (Region::ComplementarityDetermining(3), cdr3),
            (Region::Framework(4), fr4),
        ],
        annotations,
    )
}
