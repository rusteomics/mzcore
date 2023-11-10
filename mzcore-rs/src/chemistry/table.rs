#![allow(dead_code)]

use std::collections::HashMap;
use std::sync::OnceLock;

use anyhow::*;
use serde::{Deserialize, Serialize};

use crate::chemistry::amino_acid::*;
use crate::chemistry::api::AminoAcidFactory;
use crate::chemistry::atom::*;
use crate::chemistry::composition::{AtomicComposition, ElementalComposition};
use crate::chemistry::constants::*;
use crate::chemistry::element::Element as El;
use crate::chemistry::isotope::Isotope;

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AtomTable {
    pub atoms: Vec<Atom>,
    pub atom_by_element: HashMap<El, Atom>
}

impl AtomTable {
    pub fn new(
        atoms: Vec<Atom>,
    ) -> Result<AtomTable> {

        if atoms.is_empty() { bail!("atoms is empty") }

        let n_atoms = atoms.len();
        let mut atom_by_element = HashMap::with_capacity(n_atoms);
        for atom in atoms.to_owned() {
            atom_by_element.insert(atom.element, atom);
        }

        if atom_by_element.len() != n_atoms {
            bail!("provided atoms have duplicated entries")
        }

        Ok(AtomTable {
            atoms,
            atom_by_element,
        })
    }

    pub fn elemental_to_atomic_composition(&self, el_comp: ElementalComposition) -> Result<AtomicComposition> {
        let atoms_res: Result<Vec<(AtomIsotopicVariant, f32)>> = el_comp.element_counts.iter().map(|&elc| {
             match self.atom_by_element.get(&elc.element) {
                Some(atom) => {
                    match atom.isotopes.get(elc.isotope_index as usize) {
                        Some(isotope) => {
                            Ok((AtomIsotopicVariant::new(atom.to_owned(), isotope.to_owned()), elc.count))
                        }
                        None => Err(anyhow!("wrong isotope index {}", elc.isotope_index)),
                    }
                }
                None => Err(anyhow!("unknown element {}", elc.element)),
            }
        }).collect();

        Ok(AtomicComposition {
            atoms: atoms_res?,
            additional_mass:el_comp.additional_mass
        })
    }
}

static BIOMOLECULE_ATOM_TABLE: OnceLock<AtomTable> = OnceLock::new();

pub fn biomolecule_atom_table() -> &'static AtomTable {
    BIOMOLECULE_ATOM_TABLE.get_or_init(|| {
        _create_biomolecule_atom_table().unwrap()
    })
}

fn _create_biomolecule_atom_table() -> Result<AtomTable> {
    AtomTable::new(
        vec![
            Atom::new(El::H,"Hydrogen", vec![
                Isotope::new(1, 1.00782503207, 0.999885).unwrap(),
                Isotope::new(2, 2.0141017778, 0.000115).unwrap(),
                // Isotope(3, 3.016049277725, 0.0f)
            ]).unwrap(),
            // TODO: add the proton to the table ???
            /*Atom( symbol = "H+", name = "Proton", atomicNumber = 1, isotopes = Array(
              Isotope(1, 1.007276466812, 1f)
            )),*/
            Atom::new(El::C,"Carbon", vec![
                Isotope::new(12, 12.0000000, 0.9893).unwrap(),
                Isotope::new(13, 13.0033548378, 0.0107).unwrap(),
                // Isotope(14, 14.003241989, 0.0f)
            ]).unwrap(),
            Atom::new(El::N,"Nitrogen", vec![
                Isotope::new(14, 14.0030740048, 0.99636).unwrap(),
                Isotope::new(15, 15.0001088982, 0.00364).unwrap(),
            ]).unwrap(),
            Atom::new(El::O,"Oxygen", vec![
                Isotope::new(16, 15.99491461956, 0.99757).unwrap(),
                Isotope::new(17, 16.99913170, 0.00038).unwrap(),
                Isotope::new(18, 17.9991610, 0.00205).unwrap(),
            ]).unwrap(),
            Atom::new(El::P,"Phosphorus", vec![
                Isotope::new(31, 30.97376163, 1.0000).unwrap(),
            ]).unwrap(),
            Atom::new(El::S,"Sulfur", vec![
                Isotope::new(32, 31.97207100, 0.9499).unwrap(),
                Isotope::new(33, 32.97145876, 0.0075).unwrap(),
                Isotope::new(34, 33.96786690, 0.0425).unwrap(),
                Isotope::new(35, 35.96708076, 0.0001).unwrap(),
            ]).unwrap(),
            Atom::new(El::Se,"Selenium", vec![
                Isotope::new(74, 73.922475934, 0.0089).unwrap(),
                Isotope::new(76, 75.919213704, 0.0937).unwrap(),
                Isotope::new(77, 76.919914154, 0.0763).unwrap(),
                Isotope::new(78, 77.91730928, 0.2377).unwrap(),
                Isotope::new(80, 79.9165218, 0.4961).unwrap(),
                Isotope::new(82, 81.9166995, 0.0873).unwrap(),
            ]).unwrap(),
        ]
    )
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct AminoAcidTable {
    pub amino_acids: Vec<AminoAcidDefinition>,
    pub aa_by_code1: HashMap<u8, AminoAcidDefinition>
}

impl AminoAcidTable {
    pub fn new(
        amino_acids: Vec<AminoAcidDefinition>,
    ) -> Result<AminoAcidTable> {

        if amino_acids.is_empty() { bail!("amino_acids is empty") }

        let n_aas = amino_acids.len();
        let mut aa_by_code1 = HashMap::with_capacity(n_aas);
        for amino_acid in amino_acids.to_owned() {
            aa_by_code1.insert(amino_acid.code1, amino_acid);
        }

        if aa_by_code1.len() != n_aas {
            bail!("provided amino acids have duplicated entries")
        }

        Ok(AminoAcidTable {
            amino_acids: amino_acids,
            aa_by_code1: aa_by_code1,
        })
    }
}

impl AminoAcidFactory<AminoAcidDefinition> for AminoAcidTable {

    fn aa_from_byte(&self, aa: &u8) -> Result<&AminoAcidDefinition> {
        let aa_ref = self.aa_by_code1.get(&aa).ok_or_else(|| anyhow!("amino acid '{aa}' not found"))?;

        Ok(aa_ref)
    }
}


#[allow(unused_macros)]
macro_rules! svec {
    ($($x:expr),*) => {
        [ $($x),* ].iter().map(|s| s.to_string()).collect()
    }
}

/*macro_rules! s {
    ($x:expr) => { $x.to_string() }
}*/



static STANDARD_AMINO_ACID_TABLE: OnceLock<AminoAcidTable> = OnceLock::new();

pub fn standard_amino_acid_table() -> &'static AminoAcidTable {
    STANDARD_AMINO_ACID_TABLE.get_or_init(|| {
        _create_standard_amino_acid_table().unwrap()
    })
}

// Sources:
// - https://proteomicsresource.washington.edu/tools/masses.php
// - http://www.matrixscience.com/help/aa_help.html
fn _create_standard_amino_acid_table() -> Result<AminoAcidTable> {
    AminoAcidTable::new(
        vec![
            AminoAcidDefinition {
                code1: b'A',
                code3: "Ala".to_string(),
                name: "Alanine".to_string(),
                formula: Some("C(3) H(5) O N".to_string()),
                mono_mass: 71.03711381,
                average_mass: 71.0779,
                occurrence: 0.078,
                pka1: 2.35,
                pka2: 9.87,
                pka3: 0.0,
                pi: 6.01,
                codons: svec!["GCA", "GCC", "GCG", "GCU"]
            },
            AminoAcidDefinition {
                code1: b'R',
                code3: "Arg".to_string(),
                name: "Arginine".to_string(),
                formula: Some("C(6) H(12) O N(4)".to_string()),
                mono_mass: 156.1011111,
                average_mass: 156.18568,
                occurrence: 0.051,
                pka1: 1.82,
                pka2: 8.99,
                pka3: 12.48,
                pi: 10.76,
                codons: svec!["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"]
            },
            AminoAcidDefinition {
                code1: b'N',
                code3: "Asn".to_string(),
                name: "Asparagine".to_string(),
                formula: Some("C(4) H(6) O(2) N(2)".to_string()),
                mono_mass: 114.0429275,
                average_mass: 114.10264,
                occurrence: 0.043,
                pka1: 2.14,
                pka2: 8.72,
                pka3: 0.0,
                pi: 5.41,
                codons: svec!["AAC", "AAU"]
            },
            AminoAcidDefinition {
                code1: b'D',
                code3: "Asp".to_string(),
                name: "Aspartic acid".to_string(),
                formula: Some("C(4) H(5) O(3) N".to_string()),
                mono_mass: 115.0269431,
                average_mass: 115.0874,
                occurrence: 0.053,
                pka1: 1.99,
                pka2: 9.9,
                pka3: 3.9,
                pi: 2.85,
                codons: svec!["GAC", "GAU"]
            },
            AminoAcidDefinition {
                code1: b'C',
                code3: "Cys".to_string(),
                name: "Cysteine".to_string(),
                formula: Some("C(3) H(5) O N S".to_string()),
                mono_mass: 103.0091845,
                average_mass: 103.1429,
                occurrence: 0.019,
                pka1: 1.92,
                pka2: 10.7,
                pka3: 8.18,
                pi: 5.05,
                codons: svec!["UGC", "UGU"]
            },
            AminoAcidDefinition {
                code1: b'E',
                code3: "Glu".to_string(),
                name: "Glutamic acid".to_string(),
                formula: Some("C(5) H(7) O(3) N".to_string()),
                mono_mass: 129.0425931,
                average_mass: 129.11398,
                occurrence: 0.063,
                pka1: 2.1,
                pka2: 9.47,
                pka3: 4.07,
                pi: 3.15,
                codons: svec!["GAA", "GAG"]
            },
            AminoAcidDefinition {
                code1: b'Q',
                code3: "Gln".to_string(),
                name: "Glutamine".to_string(),
                formula: Some("C(5) H(8) O(2) N(2)".to_string()),
                mono_mass: 128.0585775,
                average_mass: 128.12922,
                occurrence: 0.042,
                pka1: 2.17,
                pka2: 9.13,
                pka3: 0.0,
                pi: 5.65,
                codons: svec!["CAA", "CAG"]
            },
            AminoAcidDefinition {
                code1: b'G',
                code3: "Gly".to_string(),
                name: "Glycine".to_string(),
                formula: Some("C(2) H(3) O N".to_string()),
                mono_mass: 57.02146374,
                average_mass: 57.05132,
                occurrence: 0.072,
                pka1: 2.35,
                pka2: 9.78,
                pka3: 0.0,
                pi: 6.06,
                codons: svec!["GGA", "GGC", "GGG", "GGU"]
            },
            AminoAcidDefinition {
                code1: b'H',
                code3: "His".to_string(),
                name: "Histidine".to_string(),
                formula: Some("C(6) H(7) O N(3)".to_string()),
                mono_mass: 137.0589119,
                average_mass: 137.13928,
                occurrence: 0.023,
                pka1: 1.8,
                pka2: 9.33,
                pka3: 6.04,
                pi: 7.6,
                codons: svec!["CAC", "CAU"]
            },
            AminoAcidDefinition {
                code1: b'I',
                code3: "Ile".to_string(),
                name: "Isoleucine".to_string(),
                formula: Some("C(6) H(11) O N".to_string()),
                mono_mass: 113.084064,
                average_mass: 113.15764,
                occurrence: 0.053,
                pka1: 2.32,
                pka2: 9.76,
                pka3: 0.0,
                pi: 6.05,
                codons: svec!["AUA", "AUC", "AUU"]
            },
            AminoAcidDefinition {
                code1: b'L',
                code3: "Leu".to_string(),
                name: "Leucine".to_string(),
                formula: Some("C(6) H(11) O N".to_string()),
                mono_mass: 113.084064,
                average_mass: 113.15764,
                occurrence: 0.091,
                pka1: 2.33,
                pka2: 9.74,
                pka3: 0.0,
                pi: 6.01,
                codons: svec!["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"]
            },
            AminoAcidDefinition {
                code1: b'K',
                code3: "Lys".to_string(),
                name: "Lysine".to_string(),
                formula: Some("C(6) H(12) O N(2)".to_string()),
                mono_mass: 128.0949631,
                average_mass: 128.17228,
                occurrence: 0.059,
                pka1: 2.16,
                pka2: 9.06,
                pka3: 10.54,
                pi: 9.6,
                codons: svec!["AAA", "AAG"]
            },
            AminoAcidDefinition {
                code1: b'M',
                code3: "Met".to_string(),
                name: "Methionine".to_string(),
                formula: Some("C(5) H(9) O N S".to_string()),
                mono_mass: 131.0404846,
                average_mass: 131.19606,
                occurrence: 0.023,
                pka1: 2.13,
                pka2: 9.28,
                pka3: 0.0,
                pi: 5.74,
                codons: svec!["AUG"]
            },
            AminoAcidDefinition {
                code1: b'F',
                code3: "Phe".to_string(),
                name: "Phenylalanine".to_string(),
                formula: Some("C(9) H(9) O N".to_string()),
                mono_mass: 147.0684139,
                average_mass: 147.17386,
                occurrence: 0.039,
                pka1: 2.2,
                pka2: 9.31,
                pka3: 0.0,
                pi: 5.49,
                codons: svec!["UUC", "UUU"]
            },
            AminoAcidDefinition {
                code1: b'P',
                code3: "Pro".to_string(),
                name: "Proline".to_string(),
                formula: Some("C(5) H(7) O N".to_string()),
                mono_mass: 97.05276388,
                average_mass: 97.11518,
                occurrence: 0.052,
                pka1: 1.95,
                pka2: 10.64,
                pka3: 0.0,
                pi: 6.3,
                codons: svec!["CCA", "CCC", "CCG", "CCU"]
            },
            AminoAcidDefinition {
                code1: b'S',
                code3: "Ser".to_string(),
                name: "Serine".to_string(),
                formula: Some("C(3) H(5) O(2) N".to_string()),
                mono_mass: 87.03202844,
                average_mass: 87.0773,
                occurrence: 0.068,
                pka1: 2.19,
                pka2: 9.21,
                pka3: 5.68,
                pi: 5.68,
                codons: svec!["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"]
            },
            AminoAcidDefinition {
                code1: b'T',
                code3: "Thr".to_string(),
                name: "Threonine".to_string(),
                formula: Some("C(4) H(7) O(2) N".to_string()),
                mono_mass: 101.0476785,
                average_mass: 101.10388,
                occurrence: 0.059,
                pka1: 2.09,
                pka2: 9.1,
                pka3: 5.53,
                pi: 5.6,
                codons: svec!["ACA", "ACC", "ACG", "ACU"]
            },
            AminoAcidDefinition {
                code1: b'W',
                code3: "Trp".to_string(),
                name: "Tryptophan".to_string(),
                formula: Some("C(11) H(10) O N(2)".to_string()),
                mono_mass: 186.079313,
                average_mass: 186.2099,
                occurrence: 0.014,
                pka1: 2.46,
                pka2: 9.41,
                pka3: 5.885,
                pi: 5.89,
                codons: svec!["UGG"]
            },
            AminoAcidDefinition {
                code1: b'Y',
                code3: "Tyr".to_string(),
                name: "Tyrosine".to_string(),
                formula: Some("C(9) H(9) O(2) N".to_string()),
                mono_mass: 163.0633286,
                average_mass: 163.17326,
                occurrence: 0.032,
                pka1: 2.2,
                pka2: 9.21,
                pka3: 10.46,
                pi: 5.64,
                codons: svec!["UAC", "UAU"]
            },
            AminoAcidDefinition {
                code1: b'V',
                code3: "Val".to_string(),
                name: "Valine".to_string(),
                formula: Some("C(5) H(9) O N".to_string()),
                mono_mass: 99.06841395,
                average_mass: 99.13106,
                occurrence: 0.066,
                pka1: 2.39,
                pka2: 9.74,
                pka3: 0.0,
                pi: 6.0,
                codons: svec!["GUA", "GUC", "GUG", "GUU"]
            }
        ]
    )
}

static PROTEINOGENIC_AMINO_ACID_TABLE: OnceLock<AminoAcidTable> = OnceLock::new();

pub fn proteinogenic_amino_acid_table() -> &'static AminoAcidTable {
    PROTEINOGENIC_AMINO_ACID_TABLE.get_or_init(|| {
        _create_standard_amino_acid_table().unwrap()
    })
}

// Source: http://en.wikipedia.org/wiki/Proteinogenic_amino_acid
// - https://proteomicsresource.washington.edu/tools/masses.php
// - http://www.matrixscience.com/help/aa_help.html
fn _create_proteinogenic_amino_acid_table() -> Result<AminoAcidTable> {
    AminoAcidTable::new({
        let mut extra_amino_acids = vec![
            AminoAcidDefinition {
                code1: b'B',
                code3: "Asx".to_string(),
                name: "Asn or Asp".to_string(),
                formula: None,
                mono_mass: 114.5349353,
                average_mass: 114.59502,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 0.0,
                pi: 0.0,
                codons: vec![]
            },
            AminoAcidDefinition {
                code1: b'J',
                code3: "Xle".to_string(),
                name: "Ile or Leu".to_string(),
                formula: Some("C(6) H(11) O N".to_string()),
                mono_mass: 113.084064,
                average_mass: 113.15764,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 0.0,
                pi: 0.0,
                codons: vec![]
            },
            AminoAcidDefinition {
                code1: b'O',
                code3: "Pyl".to_string(),
                name: "Pyrrolysine".to_string(),
                formula: Some("C(12) H(21) O(3) N(3)".to_string()),
                mono_mass: 237.1477266,
                average_mass: 237.298143,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 0.0,
                pi: 0.0,
                codons: svec!["UAG"]
            },
            AminoAcidDefinition {
                code1: b'U',
                code3: "Sec".to_string(),
                name: "Selenocysteine".to_string(),
                formula: Some("C(3) H(5) N O Se".to_string()),
                mono_mass: 150.9536353,
                average_mass: 150.0379,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 5.73,
                pi: 5.47,
                codons: svec!["UGA"]
            },
            AminoAcidDefinition {
                code1: b'X',
                code3: "Xaa".to_string(),
                name: "Unknown".to_string(),
                formula: None,
                mono_mass: AVERAGE_AA_MASS,
                average_mass: AVERAGE_AA_MASS,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 0.0,
                pi: 0.0,
                codons: vec![]
            },
            AminoAcidDefinition {
                code1: b'Z',
                code3: "Glx".to_string(),
                name: "Glu or Gln".to_string(),
                formula: None,
                mono_mass: 128.5505853,
                average_mass: 128.6216,
                occurrence: 0.0,
                pka1: 0.0,
                pka2: 0.0,
                pka3: 0.0,
                pi: 0.0,
                codons: vec![]
            },
        ];

        let mut v = standard_amino_acid_table().amino_acids.clone();
        v.append(&mut extra_amino_acids);
        v
    })
}

/*
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize, Deserialize)]
#[repr(u8)]
pub enum AA {
    // The 20 standard AAs
    A(AminoAcidDefinition) = b'A',
    C = b'C',
    D = b'D',
    E = b'E',
    F = b'F',
    G = b'G',
    H = b'H',
    I = b'I',
    K = b'K',
    L = b'L',
    M = b'M',
    N = b'N',
    P = b'P',
    Q = b'Q',
    R = b'R',
    S = b'S',
    T = b'T',
    V = b'V',
    W = b'W',
    Y = b'Y',

    // Proteinogenic AAs
    U = b'U',
    O = b'O',

    // Ambiguous AAs
    B = b'B',
    J = b'J',
    X = b'X',
    Z = b'Z',
}*/