use std::fmt::Display;
use std::ops::RangeInclusive;
use std::str::FromStr;

use crate::imgt_gene::IMGTGene;
use crate::shared::{AnnotatedSequence, Gene, Species};
use rustyms::{AminoAcid, IsAminoAcid};

#[derive(Debug)]
pub struct DataItem {
    pub id: String,
    pub genes: Vec<IMGTGene>,
    pub regions: Vec<Region>,
    pub species: Species,
    pub sequence: String,
}

impl Display for DataItem {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "{}\t{}\n{}", self.id, self.species, self.sequence)?;
        for gene in &self.genes {
            writeln!(f, "G {gene}")?;
        }
        for region in &self.regions {
            writeln!(f, "R {region}")?;
        }
        Ok(())
    }
}

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub struct Region {
    pub acc: String,
    pub key: String,
    pub location: Location,
    pub reported_seq: String,
    pub found_seq: Result<(String, AASequence), String>,
    pub allele: String,
    pub functional: bool,
    pub partial: bool,
    pub shift: usize,
    pub splice_aa: Option<AminoAcid>,
}

impl Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}",
            self.key,
            self.location,
            self.allele,
            self.functional,
            self.partial,
            // self.sequence,
            // dna,
            // self.found_seq.0,
            self.found_seq
                .as_ref()
                .map(|seq| seq
                    .1
                     .0
                    .iter()
                    .map(|a| a.pro_forma_definition())
                    .collect::<String>())
                .unwrap_or_else(|e| format!("<NO SEQ!>: {e}")),
        )
    }
}

pub type SequenceRegion = (
    rustyms::peptidoform::Region,
    (Vec<AminoAcid>, Location, String),
);

#[derive(Clone, Debug, Hash, Eq, PartialEq)]
pub enum Location {
    Normal(RangeInclusive<usize>),
    Complement(RangeInclusive<usize>),
    SingleNormal(usize),
    SingleComplement(usize),
}

impl Location {
    /// Check if a location overlaps or is immediately adjacent to this location.
    /// Used to detect if a CDR3 belongs to a certain V-REGION
    pub fn overlaps(&self, other: &Location) -> bool {
        match (self, other) {
            (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                *s.start() <= o.end() + 1 && s.end() + 1 >= *o.start()
            }
            (Self::Complement(s), Self::SingleComplement(o)) => s.contains(o),
            (Self::Normal(s), Self::SingleNormal(o)) => s.contains(o),
            _ => false,
        }
    }

    pub fn contains(&self, other: &Location) -> bool {
        match (self, other) {
            (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                s.start() <= o.start() && s.end() >= o.end()
            }
            (Self::Complement(s), Self::SingleComplement(o)) => s.contains(o),
            (Self::Normal(s), Self::SingleNormal(o)) => s.contains(o),
            _ => false,
        }
    }

    pub fn find_aa_location(&self, sections: &[SequenceRegion]) -> Option<usize> {
        let mut start = 0;
        for section in sections {
            if let Some(index) = section.1 .1.get_aa_loc(self) {
                return Some(start + index.start());
            }
            start += section.1 .0.len();
        }
        None
    }

    fn get_aa_loc(&self, inner: &Self) -> Option<RangeInclusive<usize>> {
        if !self.contains(inner) {
            None
        } else {
            match (self, inner) {
                (Self::Complement(s), Self::Complement(o)) | (Self::Normal(s), Self::Normal(o)) => {
                    Some((o.start() - s.start()) / 3..=(o.end() - s.start()) / 3)
                }
                (Self::Normal(s), Self::SingleNormal(o))
                | (Self::Complement(s), Self::SingleComplement(o)) => {
                    Some((o - s.start()) / 3..=(o - s.start()) / 3)
                }
                _ => None,
            }
        }
    }

    /// Break the location around the given amino acid index in the location. If the position is outside the range or this location is a single it returns None.
    pub fn splice(&self, position: usize) -> Option<(Self, Self)> {
        match self {
            Self::Normal(s) => {
                let mid_point = *s.start() + position * 3;
                if mid_point >= *s.end() {
                    None
                } else {
                    Some((
                        Self::Complement((*s.start())..=mid_point),
                        Self::Complement(mid_point..=*s.end()),
                    ))
                }
            }
            Self::Complement(s) => {
                let mid_point = *s.end() - position * 3;
                if mid_point <= *s.start() {
                    None
                } else {
                    Some((
                        Self::Complement((*s.start())..=mid_point),
                        Self::Complement(mid_point..=*s.end()),
                    ))
                }
            }
            _ => None,
        }
    }
}

impl Display for Location {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Complement(range) => write!(f, "c{}..{}", range.start(), range.end()),
            Self::Normal(range) => write!(f, "{}..{}", range.start(), range.end()),
            Self::SingleComplement(loc) => write!(f, "c{}", loc),
            Self::SingleNormal(loc) => write!(f, "{}", loc),
        }
    }
}

impl FromStr for Location {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let s = s.trim();

        if s.contains("join") || s.contains('^') {
            return Err("Location is complex, joined or it uses ^".to_string());
        }

        if let Some(tail) = s.strip_prefix("complement(") {
            tail.trim_end_matches(')')
                .split_once("..")
                .map(|(start, end)| {
                    Ok(Self::Complement(
                        start
                            .trim_start_matches('<')
                            .parse::<usize>()
                            .map_err(|err| format!("Invalid start number: {err}"))?
                            - 1
                            ..=end
                                .trim_start_matches('>')
                                .parse::<usize>()
                                .map_err(|err| format!("Invalid end number: {err}"))?
                                - 1,
                    ))
                })
                .unwrap_or_else(|| {
                    Ok(Self::SingleComplement(
                        tail.trim_end_matches(')')
                            .parse()
                            .map_err(|err| format!("Invalid single number: {err}"))?,
                    ))
                })
        } else {
            s.split_once("..")
                .map(|(start, end)| {
                    Ok(Self::Normal(
                        start
                            .trim_start_matches('<')
                            .parse::<usize>()
                            .map_err(|err| format!("Invalid start number: {err}"))?
                            - 1
                            ..=end
                                .trim_start_matches('>')
                                .parse::<usize>()
                                .map_err(|err| format!("Invalid end number: {err}"))?
                                - 1,
                    ))
                })
                .unwrap_or_else(|| {
                    Ok(Self::SingleNormal(
                        s.parse()
                            .map_err(|err| format!("Invalid single number: {err}"))?,
                    ))
                })
        }
    }
}

#[derive(Clone, Hash, Eq, PartialEq)]
pub struct AASequence(pub Vec<AminoAcid>);

impl std::fmt::Debug for AASequence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "[{}]",
            self.0
                .iter()
                .map(|a| a.pro_forma_definition())
                .collect::<String>()
        )
    }
}

#[derive(Debug)]
pub struct SingleSeq {
    pub name: Gene,
    pub allele: usize,
    pub acc: String,
    pub sequence: AnnotatedSequence,
    pub dna: String,
}
