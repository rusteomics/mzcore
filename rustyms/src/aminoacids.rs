use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::{
    formula::MolecularFormula,
    fragment::{Fragment, FragmentType, PeptidePosition},
    model::*,
    molecular_charge::CachedCharge,
    system::Mass,
    MassMode, Multi, MultiChemical, NeutralLoss, SequencePosition,
};

use std::borrow::Cow;

/// A general trait to define amino acids.
pub trait IsAminoAcid {
    /// The full name for this amino acid.
    fn name(&self) -> Cow<'_, str>;
    /// The three letter code for this amino acid. Or None if there is no common three letter
    /// definition for this amino acid.
    fn three_letter_code(&self) -> Option<Cow<'_, str>>;
    /// The one letter code for this amino acid. Or None if there is no common single character
    /// definition for this amino acid.
    #[doc(alias = "code")]
    fn one_letter_code(&self) -> Option<char>;
    /// The ProForma definition for this amino acid. If this is not a simple amino acid it can be
    /// defined as an amino acid with an additional modification. For example `X[H9C2N2]` could be
    /// used if Arginine was not defined as `R` in ProForma.
    fn pro_forma_definition(&self) -> Cow<'_, str>;
    /// The full molecular formula for this amino acid. It allows multiple molecular formulas to
    /// allow ambiguous amino acids such as B and Z.
    fn formulas(&self) -> Cow<'_, Multi<MolecularFormula>>;
    /// The monoisotopic mass of this amino acid. Should be redefined for better performance.
    fn monoisotopic_mass(&self) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(
            self.formulas()
                .iter()
                .map(MolecularFormula::monoisotopic_mass)
                .collect(),
        )
    }
    /// The average weight of this amino acid. Should be redefined for better performance.
    fn average_weight(&self) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(
            self.formulas()
                .iter()
                .map(MolecularFormula::average_weight)
                .collect(),
        )
    }
    /// The mass with a given mass mode for this amino acid. Should be redefined for better performance.
    fn mass(&self, mode: MassMode) -> Cow<'_, Multi<Mass>> {
        Cow::Owned(self.formulas().iter().map(|f| f.mass(mode)).collect())
    }
    /// The molecular formula of the side chain of the amino acid.
    fn side_chain(&self) -> Cow<'_, Multi<MolecularFormula>>;
    /// The molecular formulas that can fragment for satellite ions (d and w). Commonly the fragment
    /// after the second carbon into the side chain. `MolecularFormula::default()` can be returned
    /// if no satellite ions are possible.
    fn satellite_ion_fragments(&self) -> Option<Cow<'_, Multi<MolecularFormula>>>;
    /// Common neutral losses for the immonium ion of this amino acid.
    fn immonium_losses(&self) -> Cow<'_, [NeutralLoss]>;
}

include!("shared/aminoacid.rs");

impl AminoAcid {
    /// All amino acids with a unique mass (no I/L in favour of J, no B, no Z, and no X)
    pub const UNIQUE_MASS_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::AmbiguousLeucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
        Self::Selenocysteine,
        Self::Pyrrolysine,
    ];

    /// All 20 canonical amino acids
    pub const CANONICAL_AMINO_ACIDS: &'static [Self] = &[
        Self::Glycine,
        Self::Alanine,
        Self::Arginine,
        Self::Asparagine,
        Self::AsparticAcid,
        Self::Cysteine,
        Self::Glutamine,
        Self::GlutamicAcid,
        Self::Histidine,
        Self::Leucine,
        Self::Isoleucine,
        Self::Lysine,
        Self::Methionine,
        Self::Phenylalanine,
        Self::Proline,
        Self::Serine,
        Self::Threonine,
        Self::Tryptophan,
        Self::Tyrosine,
        Self::Valine,
    ];

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    /// # Panics
    /// When the sequence index is terminal.
    pub(crate) fn satellite_ion_fragments(
        self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Multi<MolecularFormula> {
        let crate::SequencePosition::Index(sequence_index) = sequence_index else {
            panic!("Not allowed to call satellite ion fragments with a terminal sequence index")
        };
        match self {
            Self::Alanine
            | Self::Glycine
            | Self::Histidine
            | Self::Phenylalanine
            | Self::Proline
            | Self::Tryptophan
            | Self::Tyrosine
            | Self::Unknown => Multi::default(),
            Self::Arginine => molecular_formula!(H 9 C 2 N 2).into(),
            Self::Asparagine => molecular_formula!(H 2 C 1 N 1 O 1).into(),
            Self::AsparticAcid => molecular_formula!(H 1 C 1 O 2).into(),
            Self::AmbiguousAsparagine => vec![
                molecular_formula!(H 2 C 1 N 1 O 1 (crate::AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index})),
                molecular_formula!(H 1 C 1 O 2 (crate::AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Cysteine => molecular_formula!(H 1 S 1).into(),
            Self::Glutamine => molecular_formula!(H 4 C 2 N 1 O 1).into(),
            Self::GlutamicAcid => molecular_formula!(H 3 C 2 O 2).into(),
            Self::AmbiguousGlutamine => vec![
                molecular_formula!(H 4 C 2 N 1 O 1 (crate::AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index})),
                molecular_formula!(H 3 C 2 O 2 (crate::AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Isoleucine => vec![
                molecular_formula!(H 3 C 1),
                molecular_formula!(H 5 C 2),
            ]
            .into(),
            Self::Leucine => molecular_formula!(H 7 C 3).into(),
            Self::AmbiguousLeucine => vec![
                molecular_formula!(H 3 C 1 (crate::AmbiguousLabel::AminoAcid{option: Self::Isoleucine, sequence_index, peptidoform_index})),
                molecular_formula!(H 5 C 2 (crate::AmbiguousLabel::AminoAcid{option: Self::Isoleucine, sequence_index, peptidoform_index})),
                molecular_formula!(H 7 C 3 (crate::AmbiguousLabel::AminoAcid{option: Self::Leucine, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Lysine => molecular_formula!(H 8 C 3 N 1).into(),
            Self::Methionine => molecular_formula!(H 5 C 2 S 1).into(),
            Self::Pyrrolysine => molecular_formula!(H 15 C 9 N 2 O 1).into(),
            Self::Selenocysteine => molecular_formula!(Se 1).into(),
            Self::Serine => molecular_formula!(H 1 O 1).into(),
            Self::Threonine => vec![
                molecular_formula!(H 1 O 1),
                molecular_formula!(H 3 C 1),
            ]
            .into(),
            Self::Valine => molecular_formula!(H 3 C 1).into(), // Technically two options, but both have the same mass
        }
    }

    // TODO: generalise over used storage type, so using molecularformula, monoisotopic mass, or average mass, also make sure that AAs can return these numbers in a const fashion
    #[expect(clippy::too_many_lines, clippy::too_many_arguments)]
    pub(crate) fn fragments(
        self,
        n_term: &Multi<MolecularFormula>,
        c_term: &Multi<MolecularFormula>,
        modifications: &Multi<MolecularFormula>,
        charge_carriers: &mut CachedCharge,
        sequence_index: SequencePosition,
        sequence_length: usize,
        ions: &PossibleIons,
        peptidoform_ion_index: usize,
        peptidoform_index: usize,
        allow_terminal: (bool, bool),
    ) -> Vec<Fragment> {
        let mut base_fragments = Vec::with_capacity(ions.size_upper_bound());
        let n_pos = PeptidePosition::n(sequence_index, sequence_length);
        let c_pos = PeptidePosition::c(sequence_index, sequence_length);

        if allow_terminal.0 {
            if let Some(settings) = &ions.a {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications - molecular_formula!(H 1 C 1 O 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::a(n_pos, 0),
                    n_term,
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.b {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications - molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::b(n_pos, 0),
                    n_term,
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.c {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications + molecular_formula!(H 2 N 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::c(n_pos, 0),
                    n_term,
                    charge_carriers,
                    settings,
                ));
            }
            for (aa, distance) in &ions.d.0 {
                let satellite_fragments =
                    aa.satellite_ion_fragments(sequence_index - *distance, peptidoform_index);
                if !(satellite_fragments.is_empty()
                    || satellite_fragments.len() == 0 && satellite_fragments[0].is_empty())
                {
                    base_fragments.extend(Fragment::generate_series(
                        &(-satellite_fragments
                            * modifications
                            * self.formulas_inner(sequence_index, peptidoform_index)
                            + molecular_formula!(H 1 C 1 O 1)),
                        peptidoform_ion_index,
                        peptidoform_index,
                        &FragmentType::d(n_pos, *aa, *distance, 0),
                        n_term,
                        charge_carriers,
                        &ions.d.1,
                    ));
                }
            }
        }
        if allow_terminal.1 {
            for (aa, distance) in &ions.v.0 {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * -aa.formulas_inner(sequence_index + *distance, peptidoform_index)
                        + molecular_formula!(H 3 C 2 N 1 O 1)),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::v(c_pos, *aa, *distance, 0),
                    c_term,
                    charge_carriers,
                    &ions.v.1,
                ));
            }
            for (aa, distance) in &ions.w.0 {
                let satellite_fragments =
                    aa.satellite_ion_fragments(sequence_index - *distance, peptidoform_index);
                if !(satellite_fragments.is_empty()
                    || satellite_fragments.len() == 0 && satellite_fragments[0].is_empty())
                {
                    base_fragments.extend(Fragment::generate_series(
                        &(-satellite_fragments
                            * modifications
                            * self.formulas_inner(sequence_index, peptidoform_index)
                            + molecular_formula!(H 2 N 1)),
                        peptidoform_ion_index,
                        peptidoform_index,
                        &FragmentType::w(c_pos, *aa, *distance, 0),
                        c_term,
                        charge_carriers,
                        &ions.w.1,
                    ));
                }
            }
            if let Some(settings) = &ions.x {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications + molecular_formula!(C 1 O 1) - molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::x(c_pos, 0),
                    c_term,
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.y {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications + molecular_formula!(H 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::y(c_pos, 0),
                    c_term,
                    charge_carriers,
                    settings,
                ));
            }
            if let Some(settings) = &ions.z {
                base_fragments.extend(Fragment::generate_series(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications - molecular_formula!(H 2 N 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::z(c_pos, 0),
                    c_term,
                    charge_carriers,
                    settings,
                ));
            }
        }

        if allow_terminal.0 && allow_terminal.1 {
            if let Some((charge, losses)) = &ions.immonium {
                base_fragments.extend(Fragment::generate_all(
                    &(self.formulas_inner(sequence_index, peptidoform_index)
                        * (modifications - molecular_formula!(C 1 O 1))),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::Immonium(n_pos, self.into()), // TODO: get the actual sequenceelement here
                    &Multi::default(),
                    &losses
                        .iter()
                        .filter(|(aa, _)| aa.contains(&self))
                        .flat_map(|(_, l)| l.iter())
                        .map(|l| vec![l.clone()])
                        .collect_vec(),
                    charge_carriers,
                    *charge,
                ));
            }
        }
        base_fragments
    }

    /// Get the single letter representation of the amino acid
    pub const fn char(self) -> char {
        match self {
            Self::Alanine => 'A',
            Self::AmbiguousAsparagine => 'B',
            Self::Cysteine => 'C',
            Self::AsparticAcid => 'D',
            Self::GlutamicAcid => 'E',
            Self::Phenylalanine => 'F',
            Self::Glycine => 'G',
            Self::Histidine => 'H',
            Self::Isoleucine => 'I',
            Self::AmbiguousLeucine => 'J',
            Self::Lysine => 'K',
            Self::Leucine => 'L',
            Self::Methionine => 'M',
            Self::Asparagine => 'N',
            Self::Pyrrolysine => 'O',
            Self::Proline => 'P',
            Self::Glutamine => 'Q',
            Self::Arginine => 'R',
            Self::Serine => 'S',
            Self::Threonine => 'T',
            Self::Selenocysteine => 'U',
            Self::Valine => 'V',
            Self::Tryptophan => 'W',
            Self::Unknown => 'X',
            Self::Tyrosine => 'Y',
            Self::AmbiguousGlutamine => 'Z',
        }
    }

    /// Get the 3 letter code for the amino acid
    pub const fn code(self) -> &'static str {
        match self {
            Self::Alanine => "Ala",
            Self::AmbiguousAsparagine => "Asx",
            Self::Cysteine => "Cys",
            Self::AsparticAcid => "Asp",
            Self::GlutamicAcid => "Glu",
            Self::Phenylalanine => "Phe",
            Self::Glycine => "Gly",
            Self::Histidine => "His",
            Self::Isoleucine => "Ile",
            Self::AmbiguousLeucine => "Xle",
            Self::Lysine => "Lys",
            Self::Leucine => "Leu",
            Self::Methionine => "Met",
            Self::Asparagine => "Asn",
            Self::Pyrrolysine => "Pyl",
            Self::Proline => "Pro",
            Self::Glutamine => "Gln",
            Self::Arginine => "Arg",
            Self::Serine => "Ser",
            Self::Threonine => "Thr",
            Self::Selenocysteine => "Sec",
            Self::Valine => "Val",
            Self::Tryptophan => "Trp",
            Self::Unknown => "Xaa",
            Self::Tyrosine => "Tyr",
            Self::AmbiguousGlutamine => "Glx",
        }
    }

    /// Get the full name for the amino acid
    pub const fn name(self) -> &'static str {
        match self {
            Self::Alanine => "Alanine",
            Self::AmbiguousAsparagine => "AmbiguousAsparagine",
            Self::Cysteine => "Cysteine",
            Self::AsparticAcid => "AsparticAcid",
            Self::GlutamicAcid => "GlutamicAcid",
            Self::Phenylalanine => "Phenylalanine",
            Self::Glycine => "Glycine",
            Self::Histidine => "Histidine",
            Self::Isoleucine => "Isoleucine",
            Self::AmbiguousLeucine => "AmbiguousLeucine",
            Self::Lysine => "Lysine",
            Self::Leucine => "Leucine",
            Self::Methionine => "Methionine",
            Self::Asparagine => "Asparagine",
            Self::Pyrrolysine => "Pyrrolysine",
            Self::Proline => "Proline",
            Self::Glutamine => "Glutamine",
            Self::Arginine => "Arginine",
            Self::Serine => "Serine",
            Self::Threonine => "Threonine",
            Self::Selenocysteine => "Selenocysteine",
            Self::Valine => "Valine",
            Self::Tryptophan => "Tryptophan",
            Self::Unknown => "Unknown",
            Self::Tyrosine => "Tyrosine",
            Self::AmbiguousGlutamine => "AmbiguousGlutamine",
        }
    }

    /// Check if two amino acids are considered identical. X is identical to anything, J to IL, B to ND, Z to EQ.
    pub(crate) fn canonical_identical(self, rhs: Self) -> bool {
        match (self, rhs) {
            (a, b) if a == b => true,
            (Self::Unknown, _)
            | (_, Self::Unknown)
            | (Self::AmbiguousLeucine, Self::Leucine | Self::Isoleucine)
            | (Self::Leucine | Self::Isoleucine, Self::AmbiguousLeucine)
            | (Self::AmbiguousAsparagine, Self::Asparagine | Self::AsparticAcid)
            | (Self::Asparagine | Self::AsparticAcid, Self::AmbiguousAsparagine)
            | (Self::AmbiguousGlutamine, Self::Glutamine | Self::GlutamicAcid)
            | (Self::Glutamine | Self::GlutamicAcid, Self::AmbiguousGlutamine) => true,
            _ => false,
        }
    }
}

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.char())
    }
}

#[cfg(test)]
#[expect(clippy::unreadable_literal, clippy::missing_panics_doc)]
mod tests {
    use super::*;

    #[test]
    fn mass() {
        let weight_ala = AminoAcid::Alanine.formulas()[0].average_weight();
        let mass_ala = AminoAcid::Alanine.formulas()[0].monoisotopic_mass();
        assert_ne!(weight_ala, mass_ala);
        assert!((weight_ala.value - 71.07793).abs() < 1e-5);
        assert!((mass_ala.value - 71.037113783).abs() < 1e-5);
    }

    #[test]
    fn mass_lysine() {
        let weight_lys = AminoAcid::Lysine.formulas()[0].average_weight();
        let mass_lys = AminoAcid::Lysine.formulas()[0].monoisotopic_mass();
        assert_ne!(weight_lys, mass_lys);
        assert!((weight_lys.value - 128.17240999999999).abs() < 1e-5);
        assert!((mass_lys.value - 128.094963010536).abs() < 1e-5);
    }

    #[test]
    fn masses() {
        let known = &[
            ('A', 71.03711, 71.08),
            ('R', 156.10111, 156.2),
            ('N', 114.04293, 114.1),
            ('D', 115.02694, 115.1),
            ('C', 103.00919, 103.1),
            ('E', 129.04259, 129.1),
            ('Q', 128.05858, 128.1),
            ('G', 57.02146, 57.05),
            ('H', 137.05891, 137.1),
            ('I', 113.08406, 113.2),
            ('L', 113.08406, 113.2),
            ('K', 128.09496, 128.2),
            ('M', 131.04049, 131.2),
            ('F', 147.06841, 147.2),
            ('P', 97.05276, 97.12),
            ('S', 87.03203, 87.08),
            ('T', 101.04768, 101.1),
            ('W', 186.07931, 186.2),
            ('Y', 163.06333, 163.2),
            ('V', 99.06841, 99.13),
        ];

        for (aa, mono_mass, average_weight) in known {
            let aa = AminoAcid::try_from(*aa).unwrap();
            let (mono, weight) = (
                aa.formulas()[0].monoisotopic_mass().value,
                aa.formulas()[0].average_weight().value,
            );
            println!(
                "{}: {} {} {} {}",
                aa.char(),
                mono,
                mono_mass,
                weight,
                average_weight
            );
            assert!((mono - *mono_mass).abs() < 1e-5);
            assert!((weight - *average_weight).abs() < 1e-1);
        }
    }

    #[test]
    fn read_aa() {
        assert_eq!(
            AminoAcid::try_from('B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(
            AminoAcid::try_from(b'B').unwrap(),
            AminoAcid::AmbiguousAsparagine
        );
        assert_eq!(AminoAcid::try_from('c'), Ok(AminoAcid::Cysteine));
        assert_eq!(AminoAcid::try_from('🦀'), Err(()));
    }
}
