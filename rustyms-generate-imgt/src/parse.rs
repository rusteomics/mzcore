use crate::structs::{AASequence, DataItem, Location, Region};
use crate::{complement, imgt_gene::IMGTGene, translate};
use itertools::Itertools;
use rustyms::imgt::Species;
use rustyms::sequence::AminoAcid;
use std::collections::HashMap;
use std::io::{BufRead, BufReader};

/// Parse the IMGT file
pub fn parse_dat<T: std::io::Read>(
    reader: BufReader<T>,
) -> impl Iterator<Item = Result<DataItem, String>> {
    reader
        .lines()
        .batching(|f| {
            let mut data = PreDataItem::default();
            for line in f.filter_map(|i| i.ok()) {
                if parse_dat_line(&mut data, &line) {
                    return Some(data);
                }
            }
            None
        })
        .filter(|pre| {
            pre.kw.contains(&"immunoglobulin (IG)".to_string())
                && (pre.kw.contains(&"functional".to_string())
                    || pre.kw.contains(&"germline".to_string())
                    || pre.kw.contains(&"productive".to_string()))
                && pre.os.is_some()
        })
        .map(DataItem::new)
}

/// Parse a data item line and return if it is finished or not.
fn parse_dat_line(data: &mut PreDataItem, line: &str) -> bool {
    if line.len() < 2 {
        return false;
    }
    match &line[..2] {
        "//" => return true,
        "ID" => data.id = line.to_string(),
        "KW" => data.kw.extend(
            line[5..]
                .split(';')
                .map(|s| s.trim().to_string())
                .filter(|s| !s.is_empty()),
        ),
        "FH" if line.starts_with("FH   Key") => {
            data.ft_key_width = line.find("Location").expect("Incorrect FH line") - 5
        }
        "FT" => data.ft.push(line[5..].to_string()),
        "OS" if data.os.is_none() => {
            data.os = Species::from_imgt(line[5..].trim()).unwrap_or_else(|()| {
                println!("Not a species name: `{line}`");
                None
            })
        }
        "  " => data.sq.extend(
            line.chars()
                .filter(|c| *c == 'c' || *c == 'a' || *c == 't' || *c == 'g' || *c == 'n'),
        ),
        _ => (),
    }
    false
}

impl DataItem {
    fn new(data: PreDataItem) -> Result<Self, String> {
        // println!("{}", data.id);
        let mut result = Self {
            id: data.id[5..].split(';').next().unwrap().to_string(),
            species: data.os.ok_or("No species found")?,
            sequence: data.sq,
            genes: Vec::new(),
            regions: Vec::new(),
        };
        let mut current: Option<Region> = None;
        let mut is_sequence = false;

        for line in data.ft {
            if !line.starts_with(' ') || current.is_none() {
                if let Some(region) = current.take() {
                    result.add_region(region);
                }
                if let Ok(location) = line[data.ft_key_width..].parse() {
                    let (key, location) = (line[..data.ft_key_width].trim().to_string(), location);
                    current = Some(Region {
                        acc: result.id.clone(),
                        key,
                        location,
                        reported_seq: String::new(),
                        found_seq: Err("Not loaded".to_string()),
                        allele: String::new(),
                        functional: true,
                        partial: false,
                        shift: 0,
                        splice_aa: None,
                    });
                }
                continue;
            }
            if let Some(current) = &mut current {
                DataItem::parse_ft_line(&line, current, &mut is_sequence)?;
            }
        }
        if let Some(region) = current.take() {
            result.add_region(region);
        }

        Ok(result)
    }

    fn parse_ft_line(
        line: &str,
        current: &mut Region,
        is_sequence: &mut bool,
    ) -> Result<(), String> {
        let trimmed = line.trim();
        let split = trimmed
            .split_once('=')
            .map(|(key, tail)| (key.to_ascii_lowercase(), tail));

        match split.as_ref().map(|(key, tail)| (key.as_str(), tail)) {
            Some(("/translation", tail)) => {
                current.reported_seq += tail.trim_matches('\"');
                *is_sequence = !tail.ends_with('\"');
            }
            Some(("/imgt_allele", tail)) => {
                current.allele = tail.trim_matches('\"').to_string();
            }
            Some(("/codon_start", tail)) => {
                current.shift = tail
                    .parse::<usize>()
                    .map_err(|_| format!("Not a valid codon_start: '{tail}'"))?
                    - 1;
            }
            Some(("/splice-expectedcodon", tail)) => {
                if let Some(i) = tail.find(']') {
                    current.splice_aa = AminoAcid::try_from(tail.as_bytes()[i - 1]).ok();
                }
            }
            Some(("/functional", _)) => {
                current.functional = true;
            }
            Some(("/note" | "/imgt_note", s)) if s.to_ascii_lowercase().contains("functional") => {
                current.functional = true;
            }
            Some(("/pseudo", _)) => {
                current.functional = false;
            }
            Some(("/partial", _)) => current.partial = true,
            None if *is_sequence => {
                current.reported_seq += trimmed.trim_end_matches('\"');
                *is_sequence = !trimmed.ends_with('\"');
            }
            _ => (),
        }
        Ok(())
    }

    fn add_region(&mut self, mut region: Region) {
        // println!("AR: {region}");
        // Get the actual sequence
        region.found_seq = self.get_sequence(&region.location, region.shift);

        // Determine if what this region is and if is warrants keeping
        if ["V-GENE", "C-GENE", "J-GENE"].contains(&region.key.as_str()) // , "D-GENE"
            && region.functional
            && !region.partial
            && region.allele.starts_with("IG")
        {
            self.genes.push(IMGTGene {
                acc: region.acc,
                key: region.key,
                location: region.location,
                allele: region.allele,
                regions: HashMap::new(),
            });
        } else if ["V-REGION", "C-REGION", "J-REGION"].contains(&region.key.as_str()) // , "D-GENE"
            && region.functional
            && !region.partial
            && region.allele.starts_with("IG")
        {
            if let Some(existing) = self.genes.iter_mut().find(|g| g.allele == region.allele) {
                existing.regions.insert(region.key.clone(), region);
            } else {
                self.genes.push(IMGTGene {
                    acc: region.acc,
                    key: region.key,
                    location: region.location,
                    allele: region.allele,
                    regions: HashMap::new(),
                });
            }
        } else if [
            "FR1-IMGT",
            "FR2-IMGT",
            "FR3-IMGT",
            "FR4-IMGT",
            "CDR1-IMGT",
            "CDR2-IMGT",
            "CDR3-IMGT",
            "1st-CYS",
            "2nd-CYS",
            "CONSERVED-TRP",
            "J-REGION",
            // "J-TRP",
            // "J-PHE",
            //"J-MOTIF",
            "CH1",
            "CH2",
            "H-CH2",
            "CH3",
            "CH3-CHS",
            "CH4",
            "CH4-CHS",
            "CH5",
            "CH5-CHS",
            "CH6",
            "CH6-CHS",
            "CH7",
            "CH7-CHS",
            "CH8",
            "CH8-CHS",
            "CH9",
            "CH9-CHS",
            "CHS",
            "CL",
            "C-REGION",
            "H", //"D-REGION",
            "H1",
            "H2",
            "H3",
            "H4",
            "M",
            "M1",
            "M2",
        ]
        .contains(&region.key.as_str())
        {
            if region.key == "CDR3-IMGT" {
                if let Some(gene) = self
                    .genes
                    .iter_mut()
                    .find(|g| g.location.overlaps(&region.location))
                // CDR3 does not have to be fully inside a V-REGION
                {
                    gene.regions.insert(region.key.clone(), region);
                } else {
                    self.regions.push(region)
                }
            } else if let Some(gene) = self
                .genes
                .iter_mut()
                .find(|g| g.location.contains(&region.location))
            {
                gene.regions.insert(region.key.clone(), region);
            } else {
                self.regions.push(region)
            }
        }
    }

    fn get_sequence(&self, slice: &Location, shift: usize) -> Result<(String, AASequence), String> {
        let (inner_shift, shift) = if shift == 2 { (1, 0) } else { (0, shift) };

        translate(
            &match slice {
                Location::Normal(range) => {
                    if *range.start() < inner_shift {
                        return Err("Shift outside of range".to_string());
                    }
                    self.sequence
                        .get(range.start() - inner_shift..=*range.end())
                        .ok_or("Normal outside of range")?
                        .to_string()
                }
                Location::SingleNormal(index) => char::from(
                    *self
                        .sequence
                        .as_bytes()
                        .get(*index)
                        .ok_or("Single normal outside of range")?,
                )
                .to_string(),
                Location::Complement(range) => complement(
                    self.sequence
                        .get(*range.start()..=*range.end() + inner_shift)
                        .ok_or("Complement outside of range")?
                        .to_string(),
                ),
                Location::SingleComplement(index) => complement(
                    char::from(
                        *self
                            .sequence
                            .as_bytes()
                            .get(*index)
                            .ok_or("Single complement outside of range")?,
                    )
                    .to_string(),
                ),
            }[shift..],
        )
        .map(|(s, v)| (s.to_owned(), AASequence(v)))
    }
}

#[derive(Default, Debug)]
struct PreDataItem {
    id: String,
    kw: Vec<String>,
    ft_key_width: usize,
    ft: Vec<String>,
    os: Option<Species>,
    sq: String,
}
