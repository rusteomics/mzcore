use std::slice::Iter;
use std::sync::{Arc, OnceLock};

use anyhow::*;
use serde::{Deserialize, Serialize};

use crate::chemistry::api::{HasMass, IsAminoAcidSeq};

// TODO: implement 3 kinds of peptides: LinearPeptide, ComplexPeptide and SerializablePeptide
// The current peptide definition should match the SerializablePeptide type
// We need more information for the LinearPeptide and ComplexPeptide types

// TODO: implement Eq&Hash
// FIXME: how to implement serde::{Deserialize, Serialize} for Arc<[u8]>???
#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct LinearPeptide {
    pub sequence: Arc<[u8]>,
    pub mods: Vec<SimpleModification>, // (ptm_id, seq_position)
    mono_mass: f64,
    average_mass: Option<f64>
}

/*
pub struct LinearPeptide {
    pub sequence: Arc<[u8]>,
    pub modifications: Vec<Modification>,
    mono_mass: f64,
    average_mass: f64,
}
 */

impl LinearPeptide {
    pub fn new(
        sequence: Arc<[u8]>,
        mods: Vec<SimpleModification>,
        mono_mass: f64,
        average_mass: Option<f64>
    ) -> Result<LinearPeptide> {

        if sequence.is_empty() { bail!("sequence is empty") }
        if mono_mass <= 0.0 { bail!("mono_mass must be a strictly positive number") }
        if average_mass.is_some() && average_mass.unwrap() <= 0.0 { bail!("average_mass must be a strictly positive number") }

        Ok(LinearPeptide {
            sequence: sequence,
            mods: mods,
            mono_mass: mono_mass,
            average_mass: average_mass,
        })
    }
}

impl HasMass for LinearPeptide {
    fn mono_mass(&self) -> f64 {
        self.mono_mass
    }

    fn average_mass(&self) -> Option<f64> {
        self.average_mass
    }
}

impl IsAminoAcidSeq for LinearPeptide {
    fn amino_acids_as_bytes(&self) -> Iter<u8> {
        self.sequence.iter()
    }
    fn length(&self) -> usize {
        self.sequence.len()
    }
}

trait IsModification: HasMass {
    fn id(&self) -> i64;
    fn position(&self) -> Option<i32>;
}

// TODO: put somewhere else?
// TODO: implement Eq&Hash
#[derive(Clone, Debug, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct SimpleModification {
    pub id: i64, // can be used to map to set of pre-defined PTMs
    pub mono_mass: f64,
    pub position: Option<i32>, // -1 means N-term
}

impl SimpleModification {
    fn new(id: i64,  mono_mass: f64, position: Option<i32>) -> SimpleModification {
        SimpleModification {
            id,
            mono_mass,
            position,
        }
    }
    fn from_tuple(tuple: (f64, Option<i32>)) -> SimpleModification {
        SimpleModification {
            id: generate_new_ptm_id(),
            mono_mass: tuple.0,
            position: tuple.1,
        }
    }

}

impl HasMass for SimpleModification {
    fn mono_mass(&self) -> f64 {
        self.mono_mass
    }

    fn average_mass(&self) -> Option<f64> {
        None
    }
}

impl IsModification for SimpleModification {
    fn id(&self) -> i64 {
        self.id
    }

    fn position(&self) -> Option<i32> {
        self.position
    }
}


use std::sync::atomic::{AtomicI64, Ordering};

trait InMemoryIdGen {
    fn generate_new_id(&self) -> i64;
}

impl InMemoryIdGen for AtomicI64 {
    fn generate_new_id(&self) -> i64 {
        self.fetch_sub(1, Ordering::Relaxed)
    }
}

static PTM_ID_GEN: OnceLock<AtomicI64> = OnceLock::new();

fn generate_new_ptm_id() -> i64 {
    PTM_ID_GEN.get_or_init(|| {
        AtomicI64::new(0)
    }).generate_new_id()
}

/*
macro_rules! create_id_generator {
    ($name:ident, $func_name:ident) => {
        static $name: OnceLock<AtomicI64> = OnceLock::new();

        fn $func_name() -> i64 {
            $name.get_or_init(|| {
                AtomicI64::new(0)
            })
            .fetch_sub(1, Ordering::Relaxed)
        }
    };
}

create_id_generator!(PTM_ID_GEN, generate_new_ptm_id);
*/

