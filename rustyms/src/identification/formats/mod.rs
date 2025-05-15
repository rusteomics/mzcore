//! Read in the annotations from peptide identification sources

mod basic_csv;
mod deepnovofamily;
mod fasta;
mod fragpipe;
mod instanovo;
mod maxquant;
mod mztab;
mod novob;
mod novor;
mod opair;
mod peaks;
mod pepnet;
mod plgs;
mod plink;
mod powernovo;
mod sage;
mod ssl;

use crate::*;
pub use basic_csv::*;
pub use deepnovofamily::*;
pub use fasta::*;
pub use fragpipe::*;
pub use instanovo::*;
pub use maxquant::*;
pub use mztab::*;
pub use novob::*;
pub use novor::*;
pub use opair::*;
pub use peaks::*;
pub use pepnet::*;
pub use plgs::*;
pub use plink::*;
pub use powernovo::*;
pub use sage::*;
pub use ssl::*;

#[cfg(test)]
mod deepnovofamily_tests;
#[cfg(test)]
mod fragpipe_tests;
#[cfg(test)]
mod instanovo_tests;
#[cfg(test)]
mod maxquant_tests;
#[cfg(test)]
mod mztab_test;
#[cfg(test)]
mod novob_tests;
#[cfg(test)]
mod novor_tests;
#[cfg(test)]
mod opair_tests;
#[cfg(test)]
mod peaks_tests;
#[cfg(test)]
mod pepnet_tests;
#[cfg(test)]
mod plgs_tests;
#[cfg(test)]
mod plink_tests;
#[cfg(test)]
mod powernovo_tests;
#[cfg(test)]
mod sage_tests;
#[cfg(test)]
mod ssl_tests;
