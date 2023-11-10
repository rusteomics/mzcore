#![allow(dead_code)]

// Source: http://pdg.lbl.gov/2012/reviews/rpp2012-rev-phys-constants.pdf
pub const AVERAGE_AA_MASS: f64 = 111.1254; // TODO: marco => why difference with 111.10523866044295 by computation
// TODO: refine this value and put a source reference here (publication ?)
pub const AVERAGE_PEPTIDE_ISOTOPE_MASS_DIFF: f64 = 1.0027;

pub const ELECTRON_MASS: f64 = 0.00054857990946; // Source: NIST 2010 CODATA
pub const PROTON_MASS: f64 = 1.007276466812; // Source: NIST 2010 CODATA

pub const CO_MONO_MASS: f64 = 27.99491461956;
pub const CO2_MONO_MASS: f64 = 0.0;  // FIXME
pub const H2O_MONO_MASS: f64 = 18.010565;
pub const NH3_MONO_MASS: f64 = 17.02654910101;

pub const WATER_MONO_MASS: f64 = 18.010565;
pub const WATER_AVERAGE_MASS: f64 = 18.01525697318;


// --- Sage definition --- //
// FIXME: some letters are sometimes used for undetermined amino acids (should we add them)
// FIXME: is this really needed with the previously defined enum?
// Example: B for D/N, J for I/L, Z for E/Q, X for any
pub const VALID_AA: [u8; 22] = [
    b'A', b'C', b'D', b'E', b'F', b'G', b'H', b'I', b'K', b'L', b'M', b'N', b'P', b'Q', b'R', b'S',
    b'T', b'V', b'W', b'Y', b'U', b'O',
];

// --- Sage definition --- //
pub const MONOISOTOPIC_AA_MASSES: [f64; 26] = [
    71.03711, 0.0, 103.00919, 115.02694, 129.04259, 147.0684, 57.02146, 137.05891, 113.08406, 0.0,
    128.09496, 113.08406, 131.0405, 114.04293, 237.14774, 97.05276, 128.05858, 156.1011, 87.03203,
    101.04768, 150.95363, 99.06841, 186.07932, 0.0, 163.06332, 0.0,
];

pub const fn monoisotopic_aa_mass(aa: u8) -> f64 {
    if aa.is_ascii_uppercase() {
        MONOISOTOPIC_AA_MASSES[(aa - b'A') as usize]
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {

    use super::{monoisotopic_aa_mass, VALID_AA};

    #[test]
    fn valid_aa() {
        for ch in VALID_AA {
            assert!(monoisotopic_aa_mass(ch) > 0.0);
        }
    }

}
