use bincode::{Decode, Encode};
use mzcore::sequence::{Annotation, Peptidoform, Region, UnAmbiguous};
use serde::{Deserialize, Serialize};
use std::{fmt::Display, str::FromStr, sync::Arc};

use crate::Allele;

use super::species::Species;

/// A selection of germlines from a single species.
#[derive(Debug, Decode, Deserialize, Encode, Serialize)]
pub struct Germlines {
    /// The species of this selection of germlines
    pub species: Species,
    /// All heavy chain options
    pub h: Chain,
    /// All kappa light chain options
    pub k: Chain,
    /// All lambda light chain options
    pub l: Chain,
    /// All iota chain options
    pub i: Chain,
}

impl Germlines {
    /// Create a new empty collection
    pub fn new(species: Species) -> Self {
        Self {
            species,
            h: Chain::default(),
            k: Chain::default(),
            l: Chain::default(),
            i: Chain::default(),
        }
    }

    /// Insert a germline into this collection
    pub fn insert(&mut self, germline: Germline) {
        match &germline.name.chain {
            ChainType::Heavy => self.h.insert(germline),
            ChainType::LightKappa => self.k.insert(germline),
            ChainType::LightLambda => self.l.insert(germline),
            ChainType::Iota => self.i.insert(germline),
        }
    }
}

impl<'a> IntoIterator for &'a Germlines {
    type IntoIter = std::array::IntoIter<(ChainType, &'a Chain), 4>;
    type Item = (ChainType, &'a Chain);

    fn into_iter(self) -> Self::IntoIter {
        [
            (ChainType::Heavy, &self.h),
            (ChainType::LightKappa, &self.k),
            (ChainType::LightLambda, &self.l),
            (ChainType::Iota, &self.i),
        ]
        .into_iter()
    }
}

impl Germlines {
    /// Iterate over the chains
    pub fn iter(
        &self,
    ) -> impl DoubleEndedIterator<Item = (ChainType, &Chain)> + ExactSizeIterator + '_ {
        [
            (ChainType::Heavy, &self.h),
            (ChainType::LightKappa, &self.k),
            (ChainType::LightLambda, &self.l),
            (ChainType::Iota, &self.i),
        ]
        .into_iter()
    }
}

#[cfg(feature = "rayon")]
use rayon::prelude::*;
#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Germlines {
    type Iter = rayon::array::IntoIter<(ChainType, &'a Chain), 4>;
    type Item = (ChainType, &'a Chain);

    fn into_par_iter(self) -> Self::Iter {
        [
            (ChainType::Heavy, &self.h),
            (ChainType::LightKappa, &self.k),
            (ChainType::LightLambda, &self.l),
            (ChainType::Iota, &self.i),
        ]
        .into_par_iter()
    }
}

/// The intermediate representation for a chain
#[derive(Debug, Decode, Default, Deserialize, Encode, Serialize)]
pub struct Chain {
    /// All V/variable germlines
    pub variable: Vec<Arc<Germline>>,
    /// All J/joining germlines
    pub joining: Vec<Arc<Germline>>,
    /// All C/constant germlines
    pub c: Vec<Arc<Germline>>,
    /// All A constant germlines
    pub a: Vec<Arc<Germline>>,
    /// All D constant germlines
    pub d: Vec<Arc<Germline>>,
    /// All E constant germlines
    pub e: Vec<Arc<Germline>>,
    /// All G constant germlines
    pub g: Vec<Arc<Germline>>,
    /// All M constant germlines
    pub m: Vec<Arc<Germline>>,
    /// All O constant germlines
    pub o: Vec<Arc<Germline>>,
    /// All T constant germlines
    pub t: Vec<Arc<Germline>>,
}

impl Chain {
    /// # Panics
    /// It panics when it inserts an allele it has already placed (has to be filtered and ranked before)
    pub fn insert(&mut self, mut germline: Germline) {
        let db = match &germline.name.kind {
            GeneType::V => &mut self.variable,
            GeneType::J => &mut self.joining,
            GeneType::C(None) => &mut self.c,
            GeneType::C(Some(Constant::A)) => &mut self.a,
            GeneType::C(Some(Constant::D)) => &mut self.d,
            GeneType::C(Some(Constant::E)) => &mut self.e,
            GeneType::C(Some(Constant::G)) => &mut self.g,
            GeneType::C(Some(Constant::M)) => &mut self.m,
            GeneType::C(Some(Constant::O)) => &mut self.o,
            GeneType::C(Some(Constant::T)) => &mut self.t,
        };

        match db.binary_search_by_key(&germline.name, |g| g.name.clone()) {
            // If there are multiple copies of the same region keep the one with the most annotations + regions
            Ok(index) => {
                match db[index]
                    .alleles
                    .binary_search_by_key(&germline.alleles[0].0, |a| a.0)
                {
                    Ok(_allele_index) => {
                        // if germline.alleles[0].1.sequence
                        //     == db[index].alleles[allele_index].1.sequence
                        // {
                        //     db[index].alleles[allele_index].2 = format!(
                        //         "{}|{}",
                        //         db[index].alleles[allele_index].2, germline.alleles[0].2
                        //     );
                        // } else {
                        //     let mut found = false;
                        //     for existing in &mut db[index].duplicates {
                        //         if existing.1.sequence == germline.alleles[0].1.sequence {
                        //             existing.2 =
                        //                 format!("{}|{}", existing.2, germline.alleles[0].2);
                        //             found = true;
                        //         }
                        //     }
                        //     if !found {
                        //         db[index].duplicates.push(germline.alleles.pop().unwrap())
                        //     }
                        // }
                        // if germline.alleles[0].1.conserved.len()
                        //     + germline.alleles[0].1.regions.len()
                        //     > db[index].alleles[allele_index].1.conserved.len()
                        //         + db[index].alleles[allele_index].1.regions.len()
                        // {
                        //     db[index].alleles[allele_index] = germline.alleles.pop().unwrap()
                        // }
                        panic!(
                            "Not allowed to have multiple sequences for one allele in a germline"
                        )
                    }
                    Err(allele_index) => Arc::get_mut(&mut db[index])
                        .map(|g| {
                            g.alleles
                                .insert(allele_index, germline.alleles.pop().unwrap())
                        })
                        .expect("Multiple copies of Arc while building IMGT structure"),
                }
            }
            Err(index) => db.insert(index, Arc::new(germline)),
        }
    }

    /// Create a row for use in the documentation of the numbers of different sequences available
    pub fn doc_row(&self) -> String {
        format!(
            "|{}/{}|{}/{}|{}/{}|",
            self.variable.len(),
            self.variable.iter().map(|g| g.alleles.len()).sum::<usize>(),
            self.joining.len(),
            self.joining.iter().map(|g| g.alleles.len()).sum::<usize>(),
            self.c.len()
                + self.a.len()
                + self.d.len()
                + self.e.len()
                + self.g.len()
                + self.m.len()
                + self.o.len()
                + self.t.len(),
            self.c.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.a.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.d.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.e.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.g.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.m.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.o.iter().map(|g| g.alleles.len()).sum::<usize>()
                + self.t.iter().map(|g| g.alleles.len()).sum::<usize>(),
        )
    }

    /// Get all constant chains germline regardless of the exact type of constant chain
    pub fn constant(&self) -> impl Iterator<Item = &Arc<Germline>> {
        self.c
            .iter()
            .chain(self.a.iter())
            .chain(self.d.iter())
            .chain(self.e.iter())
            .chain(self.g.iter())
            .chain(self.m.iter())
            .chain(self.o.iter())
            .chain(self.t.iter())
    }
}

impl<'a> IntoIterator for &'a Chain {
    type IntoIter = std::array::IntoIter<(GeneType, &'a [Arc<Germline>]), 10>;
    type Item = (GeneType, &'a [Arc<Germline>]);

    fn into_iter(self) -> Self::IntoIter {
        [
            (GeneType::V, self.variable.as_slice()),
            (GeneType::J, self.joining.as_slice()),
            (GeneType::C(None), self.c.as_slice()),
            (GeneType::C(Some(Constant::A)), self.a.as_slice()),
            (GeneType::C(Some(Constant::D)), self.d.as_slice()),
            (GeneType::C(Some(Constant::E)), self.e.as_slice()),
            (GeneType::C(Some(Constant::G)), self.g.as_slice()),
            (GeneType::C(Some(Constant::M)), self.m.as_slice()),
            (GeneType::C(Some(Constant::O)), self.o.as_slice()),
            (GeneType::C(Some(Constant::T)), self.t.as_slice()),
        ]
        .into_iter()
    }
}

impl Chain {
    /// Iterate over the genes
    pub fn iter(
        &self,
    ) -> impl DoubleEndedIterator<Item = (GeneType, &[Arc<Germline>])> + ExactSizeIterator + '_
    {
        [
            (GeneType::V, self.variable.as_slice()),
            (GeneType::J, self.joining.as_slice()),
            (GeneType::C(None), self.c.as_slice()),
            (GeneType::C(Some(Constant::A)), self.a.as_slice()),
            (GeneType::C(Some(Constant::D)), self.d.as_slice()),
            (GeneType::C(Some(Constant::E)), self.e.as_slice()),
            (GeneType::C(Some(Constant::G)), self.g.as_slice()),
            (GeneType::C(Some(Constant::M)), self.m.as_slice()),
            (GeneType::C(Some(Constant::O)), self.o.as_slice()),
            (GeneType::C(Some(Constant::T)), self.t.as_slice()),
        ]
        .into_iter()
    }
}

#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Chain {
    type Iter = rayon::array::IntoIter<(GeneType, &'a [Arc<Germline>]), 10>;
    type Item = (GeneType, &'a [Arc<Germline>]);

    fn into_par_iter(self) -> Self::Iter {
        [
            (GeneType::V, self.variable.as_slice()),
            (GeneType::J, self.joining.as_slice()),
            (GeneType::C(None), self.c.as_slice()),
            (GeneType::C(Some(Constant::A)), self.a.as_slice()),
            (GeneType::C(Some(Constant::D)), self.d.as_slice()),
            (GeneType::C(Some(Constant::E)), self.e.as_slice()),
            (GeneType::C(Some(Constant::G)), self.g.as_slice()),
            (GeneType::C(Some(Constant::M)), self.m.as_slice()),
            (GeneType::C(Some(Constant::O)), self.o.as_slice()),
            (GeneType::C(Some(Constant::T)), self.t.as_slice()),
        ]
        .into_par_iter()
    }
}

/// Intermediate representation for germline
#[derive(Clone, Debug, Decode, Deserialize, Encode, Eq, PartialEq, Serialize)]
pub struct Germline {
    /// The species for this germline
    pub species: Species,
    /// The name for the germline
    pub name: Gene,
    /// All alleles
    pub alleles: Vec<(usize, AnnotatedSequence)>,
}

impl<'a> IntoIterator for &'a Germline {
    type IntoIter = std::slice::Iter<'a, (usize, AnnotatedSequence)>;
    type Item = &'a (usize, AnnotatedSequence);

    fn into_iter(self) -> Self::IntoIter {
        self.alleles.iter()
    }
}

impl Germline {
    /// Iterate over the sequences
    pub fn iter(&self) -> <&Self as IntoIterator>::IntoIter {
        self.into_iter()
    }

    /// Select a specific allele
    pub fn select_allele(&self, allele: Option<usize>) -> Option<Allele<'_>> {
        allele
            .map_or_else(
                || self.alleles.first(),
                |allele| self.alleles.iter().find(|(i, _)| *i == allele),
            )
            .map(|(i, sequence)| Allele {
                species: self.species,
                gene: std::borrow::Cow::Borrowed(&self.name),
                number: *i,
                sequence: &sequence.sequence,
                regions: &sequence.regions,
                annotations: &sequence.annotations,
            })
    }
}

#[cfg(feature = "rayon")]
impl<'a> IntoParallelIterator for &'a Germline {
    type Iter = rayon::slice::Iter<'a, (usize, AnnotatedSequence)>;
    type Item = &'a (usize, AnnotatedSequence);

    fn into_par_iter(self) -> Self::Iter {
        self.alleles.par_iter()
    }
}

/// Intermediate representation for annotated sequence
#[derive(Clone, Debug, Decode, Deserialize, Encode, Eq, PartialEq, Serialize)]
pub struct AnnotatedSequence {
    /// The sequence
    pub sequence: Peptidoform<UnAmbiguous>,
    /// The different regions in the sequence, defined by their name and length
    pub regions: Vec<(Region, usize)>,
    /// 0 based locations of single amino acid annotations, overlapping with the regions defined above
    pub annotations: Vec<(Annotation, usize)>,
}

impl AnnotatedSequence {
    /// Create a new annotated sequence
    pub fn new(
        sequence: Peptidoform<UnAmbiguous>,
        regions: Vec<(Region, usize)>,
        mut conserved: Vec<(Annotation, usize)>,
    ) -> Self {
        conserved.sort_unstable_by_key(|c| c.1);
        Self {
            sequence,
            regions,
            annotations: conserved,
        }
    }
}

/// A germline gene name, broken up in its constituent parts.
#[derive(
    Clone, Debug, Deserialize, Eq, Ord, PartialEq, PartialOrd, Serialize, Hash, Decode, Encode,
)]
pub struct Gene {
    /// The chain of this gene (heavy/kappa etc)
    pub chain: ChainType,
    /// The kind of gene (V/J/C)
    pub kind: GeneType,
    /// If present the additional number _IGHV_ **(I)**
    pub number: Option<usize>,
    /// The family indicators _IGHV_ **-1D**
    pub family: Vec<(Option<usize>, String)>,
}

impl Display for Gene {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        const fn to_roman(n: usize) -> &'static str {
            [
                "0", "I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X",
            ][n]
        }

        write!(
            f,
            "IG{}{}{}{}",
            self.chain,
            self.kind,
            self.number
                .as_ref()
                .map_or_else(String::new, |n| format!("({})", to_roman(*n))),
            if self.number.is_some() && !self.family.is_empty() {
                "-"
            } else {
                ""
            }
        )?;

        let mut first = true;
        let mut last_str = false;
        for element in &self.family {
            if !first && !last_str {
                write!(f, "-")?;
            }
            write!(
                f,
                "{}{}",
                element.0.map(|i| i.to_string()).unwrap_or_default(),
                element.1
            )?;
            last_str = !element.1.is_empty();
            first = false;
        }
        Ok(())
    }
}

impl Gene {
    /// Get an IMGT name with allele, eg IGHV3-23*03
    /// # Errors
    /// If not recognised as a name, returns a description of the error.
    #[expect(clippy::missing_panics_doc)] // Cannot panic
    pub fn from_imgt_name_with_allele(s: &str) -> Result<(Self, usize), String> {
        let s = s.split(" or ").next().unwrap(); // Just ignore double names
        let (gene, tail) = Self::from_imgt_name_internal(s)?;
        if tail.is_empty() {
            return Ok((gene, 1));
        }
        let allele = tail.strip_prefix('*').map_or_else(
            || Err(format!("Invalid allele spec: `{tail}`")),
            |tail| {
                tail.parse()
                    .map_err(|_| format!("Invalid allele spec: `{}`", &tail))
            },
        )?;
        Ok((gene, allele))
    }

    /// Get an IMGT name, eg IGHV3-23
    /// # Errors
    /// If not recognised as a name, returns a description of the error.
    pub fn from_imgt_name(s: &str) -> Result<Self, String> {
        Self::from_imgt_name_internal(s).map(|(gene, _)| gene)
    }

    /// # Errors
    /// If not recognised as a name, returns a description of the error.
    fn from_imgt_name_internal(s: &str) -> Result<(Self, &str), String> {
        #[expect(clippy::missing_panics_doc)] // Cannot panic
        fn parse_name(s: &str) -> (Option<(Option<usize>, String)>, &str) {
            let num = s
                .chars()
                .take_while(char::is_ascii_digit)
                .collect::<String>();
            let tail = s
                .chars()
                .skip(num.len())
                .take_while(char::is_ascii_alphabetic)
                .collect::<String>();
            let rest = &s[num.len() + tail.len()..];
            if num.is_empty() && tail.is_empty() {
                return (None, s);
            }
            let num = if num.is_empty() {
                None
            } else {
                Some(num.parse().unwrap())
            };
            (Some((num, tail)), rest)
        }

        fn from_roman(s: &str) -> Option<usize> {
            match s {
                "Ⅰ" | "I" => Some(1),
                "Ⅱ" | "II" => Some(2),
                "Ⅲ" | "III" => Some(3),
                "Ⅳ" | "IV" => Some(4),
                "Ⅴ" | "V" => Some(5),
                "Ⅵ" | "VI" => Some(6),
                "Ⅶ" | "VII" => Some(7),
                "Ⅷ" | "VIII" => Some(8),
                "Ⅸ" | "IX" => Some(9),
                "Ⅹ" | "X" => Some(10),
                _ => None,
            }
        }

        if s.starts_with("IG") {
            let chain = s[2..3]
                .parse()
                .map_err(|()| format!("Invalid chain: `{}`", &s[2..3]))?;
            let gene = s[3..4]
                .parse()
                .map_err(|()| format!("Invalid gene: `{}`", &s[3..4]))?;
            let mut start = 4;
            let number = if s.len() > 4 && &s[4..5] == "(" {
                let end = s[5..]
                    .find(')')
                    .ok_or_else(|| format!("Invalid gene number `{}` out of `{}`", &s[4..], s))?;
                start += end + 2;
                Some(from_roman(&s[5..5 + end]).ok_or_else(|| {
                    format!("Invalid roman numeral (or too big) `{}`", &s[5..5 + end])
                })?)
            } else {
                None
            };
            let tail = &s[start..];
            let mut tail = tail.trim_start_matches('-');
            let mut family = Vec::new();
            while let (Some(branch), t) = parse_name(tail) {
                family.push(branch);
                tail = t.trim_start_matches('-');
            }

            Ok((
                Self {
                    chain,
                    kind: gene,
                    number,
                    family,
                },
                tail,
            ))
        } else {
            Err("Gene name does not start with IG")?
        }
    }
}

/// Any chain type of germline
#[derive(
    Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Decode, Encode,
)]
pub enum ChainType {
    /// Heavy chain
    Heavy = 0,
    /// Light kappa chain
    LightKappa,
    /// Light lambda chain
    LightLambda,
    /// Fish I kind
    Iota,
}

impl TryFrom<usize> for ChainType {
    type Error = ();
    fn try_from(i: usize) -> Result<Self, Self::Error> {
        match i {
            0 => Ok(Self::Heavy),
            1 => Ok(Self::LightKappa),
            2 => Ok(Self::LightLambda),
            3 => Ok(Self::Iota),
            _ => Err(()),
        }
    }
}

impl FromStr for ChainType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "H" => Ok(Self::Heavy),
            "κ" | "K" => Ok(Self::LightKappa),
            "λ" | "L" => Ok(Self::LightLambda),
            "ι" | "I" => Ok(Self::Iota),
            _ => Err(()),
        }
    }
}

impl Display for ChainType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Heavy => "H",
                Self::LightKappa => "K",
                Self::LightLambda => "L",
                Self::Iota => "I",
            }
        )
    }
}

/// Any gene in a germline, eg variable, joining
#[derive(
    Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Decode, Encode,
)]
pub enum GeneType {
    /// Variable
    V,
    /// Joining
    J,
    /// Constant, potentially with the type of constant given as well
    C(Option<Constant>),
}

/// Any type of constant gene
#[allow(missing_docs)]
#[derive(
    Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Decode, Encode,
)]
pub enum Constant {
    A,
    D,
    E,
    G,
    M,
    O,
    // DD,
    // MD,
    T,
}

impl FromStr for GeneType {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "V" => Ok(Self::V),
            "J" => Ok(Self::J),
            "C" => Ok(Self::C(None)),
            "α" | "A" => Ok(Self::C(Some(Constant::A))),
            "δ" | "D" => Ok(Self::C(Some(Constant::D))),
            "ε" | "E" => Ok(Self::C(Some(Constant::E))),
            "ɣ" | "G" => Ok(Self::C(Some(Constant::G))),
            "μ" | "M" => Ok(Self::C(Some(Constant::M))),
            "ο" | "O" => Ok(Self::C(Some(Constant::O))),
            "τ" | "T" => Ok(Self::C(Some(Constant::T))),
            // "DD" => Ok(Self::C(Some(Constant::DD))),
            // "MD" => Ok(Self::C(Some(Constant::MD))),
            _ => Err(()),
        }
    }
}

impl Display for GeneType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::V => "V",
                Self::J => "J",
                Self::C(None) => "C",
                Self::C(Some(Constant::A)) => "A",
                Self::C(Some(Constant::D)) => "D",
                Self::C(Some(Constant::E)) => "E",
                Self::C(Some(Constant::G)) => "G",
                Self::C(Some(Constant::M)) => "M",
                Self::C(Some(Constant::O)) => "O",
                // Self::C(Some(Constant::DD)) => "DD",
                // Self::C(Some(Constant::MD)) => "MD",
                Self::C(Some(Constant::T)) => "T",
            }
        )
    }
}

#[expect(clippy::missing_panics_doc)]
#[test]
fn imgt_names() {
    assert_eq!(
        Gene::from_imgt_name_with_allele("IGHV3-23*03")
            .map(|(g, a)| (g.to_string(), a))
            .unwrap(),
        ("IGHV3-23".to_string(), 3)
    );
    assert_eq!(
        Gene::from_imgt_name_with_allele("IGKV6-d*01")
            .map(|(g, a)| (g.to_string(), a))
            .unwrap(),
        ("IGKV6-d".to_string(), 1)
    );
}
