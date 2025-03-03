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
pub trait IsAminoAcid: MultiChemical {
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
    /// The molecular formula of the side chain of the amino acid.  The `sequence_index` and
    /// `peptidoform_index` are used to keep track of ambiguous amino acids.
    fn side_chain(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Cow<'_, Multi<MolecularFormula>>;
    /// The molecular formulas that can fragment for satellite ions (d and w). Commonly the fragment
    /// after the second carbon into the side chain. `MolecularFormula::default()` can be returned
    /// if no satellite ions are possible. The `sequence_index` and `peptidoform_index` are used to
    /// keep track of ambiguous amino acids.
    fn satellite_ion_fragments(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Option<Cow<'_, Multi<MolecularFormula>>>;
    /// Common neutral losses for the immonium ion of this amino acid.
    fn immonium_losses(&self) -> Cow<'_, [NeutralLoss]>;
}

impl std::fmt::Display for dyn IsAminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pro_forma_definition())
    }
}

include!("shared/aminoacid.rs");

impl std::fmt::Display for AminoAcid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.pro_forma_definition())
    }
}

impl IsAminoAcid for AminoAcid {
    /// Get the single letter representation of the amino acid
    fn one_letter_code(&self) -> Option<char> {
        Some(match self {
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
        })
    }

    fn pro_forma_definition(&self) -> Cow<'_, str> {
        Cow::Borrowed(match self {
            Self::Alanine => "A",
            Self::AmbiguousAsparagine => "B",
            Self::Cysteine => "C",
            Self::AsparticAcid => "D",
            Self::GlutamicAcid => "E",
            Self::Phenylalanine => "F",
            Self::Glycine => "G",
            Self::Histidine => "H",
            Self::Isoleucine => "I",
            Self::AmbiguousLeucine => "J",
            Self::Lysine => "K",
            Self::Leucine => "L",
            Self::Methionine => "M",
            Self::Asparagine => "N",
            Self::Pyrrolysine => "O",
            Self::Proline => "P",
            Self::Glutamine => "Q",
            Self::Arginine => "R",
            Self::Serine => "S",
            Self::Threonine => "T",
            Self::Selenocysteine => "U",
            Self::Valine => "V",
            Self::Tryptophan => "W",
            Self::Unknown => "X",
            Self::Tyrosine => "Y",
            Self::AmbiguousGlutamine => "Z",
        })
    }

    /// Get the 3 letter code for the amino acid
    fn three_letter_code(&self) -> Option<Cow<'_, str>> {
        Some(Cow::Borrowed(match self {
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
        }))
    }

    /// Get the full name for the amino acid
    fn name(&self) -> Cow<'_, str> {
        Cow::Borrowed(match self {
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
        })
    }

    fn side_chain(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Cow<'_, Multi<MolecularFormula>> {
        let crate::SequencePosition::Index(sequence_index) = sequence_index else {
            return Cow::Owned(Multi::default());
        };
        Cow::Owned(match self {
            Self::Alanine => molecular_formula!(H 3 C 1).into(),
            Self::Arginine => molecular_formula!(H 10 C 4 N 3).into(), // One of the H's counts as the charge carrier and is added later
            Self::Asparagine => molecular_formula!(H 4 C 2 O 1 N 1).into(),
            Self::AsparticAcid => molecular_formula!(H 3 C 2 O 2).into(),
            Self::AmbiguousAsparagine => vec![
                molecular_formula!(H 4 C 2 O 1 N 1 (crate::AmbiguousLabel::AminoAcid{option: Self::Asparagine, sequence_index, peptidoform_index})),
                molecular_formula!(H 3 C 2 O 2 (crate::AmbiguousLabel::AminoAcid{option: Self::AsparticAcid, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Cysteine => molecular_formula!(H 3 C 1 S 1).into(),
            Self::Glutamine => molecular_formula!(H 6 C 3 O 1 N 1).into(),
            Self::GlutamicAcid => molecular_formula!(H 5 C 3 O 2).into(),
            Self::AmbiguousGlutamine => vec![
                molecular_formula!(H 6 C 3 O 1 N 1 (crate::AmbiguousLabel::AminoAcid{option: Self::Glutamine, sequence_index, peptidoform_index})),
                molecular_formula!(H 5 C 3 O 2 (crate::AmbiguousLabel::AminoAcid{option: Self::GlutamicAcid, sequence_index, peptidoform_index})),
            ]
            .into(),
            Self::Glycine => molecular_formula!(H 1).into(),
            Self::Histidine => molecular_formula!(H 5 C 4 N 2).into(),
            Self::AmbiguousLeucine | Self::Isoleucine | Self::Leucine => {
                molecular_formula!(H 9 C 4).into()
            }
            Self::Lysine => molecular_formula!(H 10 C 4 N 1).into(),
            Self::Methionine => molecular_formula!(H 7 C 3 S 1).into(),
            Self::Phenylalanine => molecular_formula!(H 7 C 7).into(),
            Self::Proline => molecular_formula!(H 5 C 3).into(),
            Self::Pyrrolysine => molecular_formula!(H 17 C 9 O 1 N 2).into(),
            Self::Selenocysteine => molecular_formula!(H 3 C 1 Se 1).into(),
            Self::Serine => molecular_formula!(H 3 C 1 O 1).into(),
            Self::Threonine => molecular_formula!(H 5 C 2 O 1).into(),
            Self::Tryptophan => molecular_formula!(H 8 C 9 N 1).into(),
            Self::Tyrosine => molecular_formula!(H 7 C 7 O 1).into(),
            Self::Valine => molecular_formula!(H 7 C 3).into(),
            Self::Unknown => molecular_formula!().into(),
        })
    }

    // TODO: Take side chain mutations into account (maybe define pyrrolysine as a mutation)
    fn satellite_ion_fragments(
        &self,
        sequence_index: SequencePosition,
        peptidoform_index: usize,
    ) -> Option<Cow<'_, Multi<MolecularFormula>>> {
        let crate::SequencePosition::Index(sequence_index) = sequence_index else {
            return None;
        };
        Some(Cow::Owned(
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
        }))
    }

    /// All losses from the base immonium ions. Compiled from the sources below.
    ///
    /// | AA    | [Wikipedia](https://upload.wikimedia.org/wikipedia/commons/thumb/0/01/Amino_acid_fragment_ions.png/400px-Amino_acid_fragment_ions.png) |  0.1016/1044-0305(93)87006-X  | [ionsource](https://www.ionsource.com/Card/immon/immon.htm) | [10.1002/chin.199624319](http://dx.doi.org/10.1002/chin.199624319) | [Prospector   (MS-Comp)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=mscomp) | [10.1186/1477-5956-9-2](http://dx.doi.org/10.1186/1477-5956-9-2) |  10.1016/j.ymeth.2004.08.013   | 10.1385/1597452750 (table 5)  | 10.1021/ac902712f  | [Prospector   (MS-Product)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msproduct) | [ThermoFisher](https://tools.thermofisher.com/content/sfs/brochures/cms_040030.pdf) | 10.1074/mcp.O113.035915 | 10.1074/mcp.O113.035915 | 10.1021/ac902712f | [Prospector   (MS-Product)](https://prospector.ucsf.edu/prospector/cgi-bin/msform.cgi?form=msproduct) |  | 10.1385/1597452750 (table 5) | Sources | Best mass | Best formula | Loss     | Loss formula | Interpreted loss | Interpreted formula     | Final      |
    /// |-------|---------------------------------------------------------------------------------------------------------------------------|-----------------------------|------------------------------------------------|------------------------------------------|-----------------------------------------------------------------------|-----------------------------------------|-----------------------------|------------------------------|-------------------|--------------------------------------------------------------------------|---------------------------------------------------------------------|-------------------------|-------------------------|-------------------|--------------------------------------------------------------------------|------------------------------|---------|----------:|--------------|----------|--------------|------------------|-------------------------|------------|
    /// | A     | 44                                                                                                                        | 44                          |                                                | 44                                       |                                                                       | 44                                      |                             | 44.05                        |                   |                                                                          | 44.0500                                                             |                         |                         |                   |                                                                          |                              | 6       |   44.0500 |              |          |              |                  |                         |            |
    /// | R     | 129                                                                                                                       | 129                         |                                                | 129                                      |                                                                       | 129                                     |                             | 129.11                       |                   |                                                                          | 129.1140                                                            | 129.1135                | C5H13N4+                |                   |                                                                          |                              | 8       |  129.1138 | C5H13N4+     |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             | 185                                            |                                          | 185                                                                   |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 2       |       185 |              | -55.8862 |              | C-2O-2           |                         | C-2O-2     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 115.09                       |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |    115.09 |              | 14.0238  |              | C1H2             |                         | C1H2       |
    /// |       | 112                                                                                                                       | 112                         | 112                                            | 112                                      | 112                                                                   | 112                                     | 112.09                      | 112.09                       |                   | 112.0869                                                                 | 112.0875                                                            |                         |                         |                   | C5H10N3+                                                                 | C5H10N3+                     | 12      |  112.0872 | C5H10N3+     | 17.0266  | H3N1         |                  |                         | H3N1       |
    /// |       | 100                                                                                                                       | 100                         | 100                                            | 100                                      | 100                                                                   | 100                                     | 100.09                      |                              |                   | 100.0869                                                                 | 100.0875                                                            |                         |                         |                   | C4H10N3+                                                                 |                              | 10      |  100.0872 | C4H10N3+     | 29.0266  | C1H3N1       |                  |                         | C1H3N1     |
    /// |       | 87                                                                                                                        | 87                          | 87                                             | 87                                       | 87                                                                    | 87                                      | 87.09                       | 87.09                        | 87.0922           | 87.0917                                                                  |                                                                     |                         |                         | C4H11N2+          | C4H11N2+                                                                 |                              | 12      |   87.0920 | C4H11N2+     | 42.0218  | C2H2N2       |                  |                         | C2H2N2     |
    /// |       | 73                                                                                                                        | 73                          | 73                                             |                                          | 73                                                                    | 72                                      | 73.00                       |                              | 73.0640           |                                                                          |                                                                     |                         |                         | C2H7N3+           |                                                                          |                              | 8       |   73.0640 | C2H7N3+      | 56.0498  | C3H6N1       |                  |                         | C3H6N1     |
    /// |       | 70                                                                                                                        | 70                          | 70                                             | 70                                       | 70                                                                    | 70                                      | 70.07                       | 70.07                        | 70.0657           | 70.0651                                                                  | 70.0657                                                             |                         |                         | C4H8N1+           | C4H8N1+                                                                  |                              | 13      |   70.0655 | C4H8N1+      | 59.0483  | C1H5N3       |                  |                         | C1H5N3     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 60.06                        |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |     60.06 |              | 69.0538  |              | C3H4N2O-1        | C2H6N1O1+               | C3H4N2O-1  |
    /// |       |                                                                                                                           | 59                          |                                                |                                          |                                                                       | 59                                      |                             |                              | 59.0483           |                                                                          |                                                                     |                         |                         | CH5N3+            |                                                                          |                              | 4       |   59.0483 | CH5N3+       | 70.0655  | C4H8N1       |                  |                         | C4H8N1     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              | 43.0296           |                                                                          |                                                                     |                         |                         | C1H3N2+           |                                                                          |                              | 2       |   43.0296 | C1H3N2+      | 86.0842  | C4H10N2      |                  |                         | C4H10N2    |
    /// |       | 29                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        29 |              | 100.1138 |              |                  | H1N2/C1H1O1/C1H3N1/C2H5 |            |
    /// | N     | 87                                                                                                                        | 87                          | 87                                             | 87                                       | 87                                                                    | 87                                      | 87.09                       | 87.06                        |                   | 87.0553                                                                  | 87.0558                                                             |                         |                         |                   | C3H7N2O1+                                                                |                              | 11      |   87.0556 | C3H7N2O1+    |          |              |                  |                         |            |
    /// |       | 70                                                                                                                        | 70                          | 70                                             | 70                                       |                                                                       | 70                                      |                             | 70.03                        |                   |                                                                          | 70.0293                                                             |                         |                         |                   |                                                                          | C3H4N1O1+                    | 8       |   70.0293 | C3H4N1O1+    | 17.0263  | H3N1         |                  |                         | H3N1       |
    /// | D     | 88                                                                                                                        | 88                          | 88                                             | 88                                       | 88                                                                    | 88                                      | 88.04                       | 88.04                        | 88.0399           | 88.0393                                                                  | 88.0399                                                             |                         |                         | C3H6N1O2+         | C3H6N1O2+                                                                |                              | 13      |   88.0397 | C3H6N1O2+    |          |              |                  |                         |            |
    /// |       | 70                                                                                                                        |                             | 70                                             | 70                                       |                                                                       | 70                                      |                             | 70.03                        |                   |                                                                          | 70.0293                                                             |                         |                         |                   |                                                                          | C3H4N1O1+                    | 7       |   70.0293 | C3H4N1O1+    | 18.0104  | H2O1         |                  |                         | H2O1       |
    /// | C     | 76                                                                                                                        |                             |                                                | 76                                       |                                                                       | 76                                      |                             |                              |                   |                                                                          | 76.0221                                                             |                         |                         |                   |                                                                          |                              | 4       |   76.0221 |              |          |              |                  |                         |            |
    /// | E     | 102                                                                                                                       | 102                         |                                                | 102                                      | 102                                                                   | 102                                     | 102.06                      | 102.05                       | 102.0555          | 102.0550                                                                 | 102.0555                                                            | 102.0550                | C4H8N1O2+               | C4H8N1O2+         | C4H8N1O2+                                                                |                              | 14      |  102.0553 | C4H8N1O2+    |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             |                                                | 91                                       |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        91 |              | 11.0553  |              |                  |                         |            |
    /// |       |                                                                                                                           |                             |                                                | 84                                       |                                                                       |                                         |                             | 84.04                        |                   |                                                                          | 84.0449                                                             |                         |                         |                   |                                                                          | C4H6N1O1+                    | 4       |   84.0449 | C4H6N1O1+    | 18.0104  | H2O1         |                  |                         | H2O1       |
    /// | Q     | 101                                                                                                                       | 101                         | 101                                            | 101                                      | 101                                                                   | 101                                     | 101.11                      | 101.11                       |                   | 101.0709                                                                 | 101.0715                                                            | 101.0709                | C4H9N2O1+               |                   | C4H9N2O1+                                                                |                              | 13      |  101.0711 | C4H9N2O1+    |          |              |                  |                         |            |
    /// |       | 129                                                                                                                       | 129                         | 129                                            | 129                                      | 129                                                                   | 129                                     | 129.1                       | 129.11                       |                   | 129.0659                                                                 | 129.1028                                                            |                         |                         |                   | C5H9N2O2+                                                                |                              | 11      |  129.0844 | C5H9N2O2+    | -28.0133 | C-1O-1       |                  |                         | C-1O-1     |
    /// |       | 84                                                                                                                        | 84                          | 84                                             | 84                                       | 84                                                                    | 84                                      | 84.08                       | 84.04                        | 84.0813           | 84.0444                                                                  | 84.0449                                                             |                         |                         | C5H10N1+          | C4H6N1O1+                                                                | C4H6N1O1+                    | 14      |   84.0569 | C5H10N1+     | 17.0142  | H3N1         |                  |                         | H3N1       |
    /// |       | 56                                                                                                                        |                             |                                                | 56                                       |                                                                       | 56                                      |                             | 56.05                        |                   |                                                                          | 56.0500                                                             |                         |                         |                   |                                                                          |                              | 5       |   56.0500 |              | 45.0211  |              | C1H3N1O1         |                         | C1H3N1O1   |
    /// | G     | 30                                                                                                                        | 30                          |                                                | 30                                       |                                                                       | 30                                      |                             | 30.03                        | 30.0344           |                                                                          | 30.0344                                                             |                         |                         | C1H4N1+           |                                                                          |                              | 8       |   30.0344 | C1H4N1+      |          |              |                  |                         |            |
    /// | H     | 110                                                                                                                       | 110                         | 110                                            | 110                                      | 110                                                                   | 110                                     | 110.07                      | 110.07                       | 110.0718          | 110.0713                                                                 | 110.0718                                                            | 110.0713                | C5H8N3+                 | C5H8N3+           | C5H8N3+                                                                  |                              | 15      |  110.0716 | C5H8N3+      |          |              |                  |                         |            |
    /// |       | 166                                                                                                                       | 166                         |                                                |                                          |                                                                       | 166                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       166 |              | -55.9284 |              | C-2O-2           |                         | C-2O-2     |
    /// |       | 138                                                                                                                       | 138                         |                                                |                                          |                                                                       | 138                                     | 138.07                      |                              |                   | 138.0662                                                                 |                                                                     |                         |                         |                   | C6H8N3O1+                                                                |                              | 6       |  138.0662 | C6H8N3O1+    | -27.9946 | C-1O-1       |                  |                         | C-1O-1     |
    /// |       | 123                                                                                                                       | 123                         |                                                |                                          |                                                                       | 123                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       123 |              | -12.9284 |              | H3O-1            |                         | H3O-1      |
    /// |       | 121                                                                                                                       | 121                         |                                                |                                          |                                                                       | 121                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |       121 |              | -10.9284 |              | H5O-1            |                         | H5O-1      |
    /// |       | 82                                                                                                                        | 82                          |                                                |                                          |                                                                       | 82                                      |                             |                              | 82.0531           |                                                                          |                                                                     |                         |                         | C4H6N2+           |                                                                          |                              | 5       |   82.0531 | C4H6N2+      | 28.0185  | C1H2N1       |                  |                         | C1H2N1     |
    /// | I/L/J | 86                                                                                                                        | 86                          | 86                                             | 86                                       | 86                                                                    | 86                                      | 86.1                        | 86.10                        | 86.0970           | 86.0964                                                                  | 86.0970                                                             |                         |                         | C5H12N+           | C5H12N1+                                                                 |                              | 13      |   86.0968 | C5H12N+      |          |              |                  |                         |            |
    /// |       | 72                                                                                                                        | 72                          |                                                | 72                                       |                                                                       | 72                                      |                             |                              |                   |                                                                          | 72.0449                                                             |                         |                         |                   |                                                                          |                              | 5       |   72.0449 |              | 14.0519  |              | C1H2             |                         | C1H2       |
    /// |       | 44                                                                                                                        |                             |                                                | 44                                       |                                                                       | 44                                      |                             |                              |                   |                                                                          | 44.0500                                                             |                         |                         |                   |                                                                          |                              | 4       |   44.0500 |              | 42.0468  |              | C3H6             |                         | C3H6       |
    /// | K     | 101                                                                                                                       | 101                         | 101                                            | 101                                      | 101                                                                   | 101                                     | 101.11                      | 101.11                       |                   | 101.1073                                                                 | 101.1079                                                            |                         |                         |                   | C5H13N2+                                                                 |                              | 11      |  101.1076 | C5H13N2+     |          |              |                  |                         |            |
    /// |       | 129                                                                                                                       | 129                         | 129                                            | 129                                      | 129                                                                   | 129                                     | 129.1                       | 129.11                       |                   | 129.1022                                                                 |                                                                     | 129.1022                | C6H13N2O1+              |                   | C6H13N2O1+                                                               |                              | 12      |  129.1022 | C6H13N2O1+   | -27.9946 | C-1O-1       |                  |                         | C-1O-1     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   | 126.0913                                                                 |                                                                     |                         |                         |                   | C7H12N1O1+                                                               |                              | 2       |  126.0913 | C7H12N1O1+   | -24.9837 | C-2H1N1O-1   |                  |                         | C-2H1N1O-1 |
    /// |       | 112                                                                                                                       | 112                         |                                                | 112                                      |                                                                       | 112                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 4       |       112 |              | -10.8924 |              | H5O-1            |                         | H5O-1      |
    /// |       | 84                                                                                                                        | 84                          | 84                                             | 84                                       | 84                                                                    | 84                                      | 84.08                       | 84.08                        | 84.0813           | 84.0808                                                                  | 84.0813                                                             |                         |                         | C5H10N1+          | C5H10N1+                                                                 | C5H10N1+                     | 14      |   84.0811 | C5H10N1+     | 17.0265  | H3N1         |                  |                         | H3N1       |
    /// |       | 70                                                                                                                        | 70                          |                                                |                                          |                                                                       | 70                                      |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 3       |        70 |              | 31.1076  |              | C1H5N1           |                         | C1H5N1     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             | 56.05                        |                   |                                                                          | 56.0500                                                             |                         |                         |                   |                                                                          |                              | 2       |   56.0500 |              | 45.0576  |              | C2H7N1           |                         | C2H7N1     |
    /// | M     | 104                                                                                                                       | 104                         | 104                                            | 104                                      | 104                                                                   | 104                                     | 104.05                      | 104.06                       |                   | 104.0528                                                                 | 104.0534                                                            |                         |                         |                   | C4H10N1S1+                                                               |                              | 11      |  104.0531 | C4H10N1S1+   |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       | 70                                      |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              |         |        70 |              | 34.0531  |              | H2S1             |                         | H2S1       |
    /// |       | 61                                                                                                                        | 61                          |                                                |                                          |                                                                       | 61                                      |                             |                              | 61.0112           |                                                                          |                                                                     |                         |                         | C2H5S1+           |                                                                          |                              | 5       |   61.0112 | C2H5S1+      | 43.0419  | C2H3N1       |                  |                         | C2H3N1     |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C3H6N1+                      | 1       |        ?? | C3H6N1+      | ??       | C1H4S1       |                  |                         | C1H4S1     |
    /// | F     | 120                                                                                                                       | 120                         | 120                                            | 120                                      | 120                                                                   | 120                                     | 120.08                      | 120.08                       | 120.0813          | 120.0808                                                                 | 120.0813                                                            | 120.0808                | C8H10N+                 | C8H10N1+          | C8H10N1+                                                                 |                              | 15      |  120.0811 | C8H10N+      |          |              |                  |                         |            |
    /// |       | 91                                                                                                                        | 91                          |                                                | 91                                       |                                                                       | 91                                      |                             |                              |                   |                                                                          | 91.0548                                                             |                         |                         |                   |                                                                          |                              | 5       |   91.0548 |              | 29.0263  |              | C1H3N1           | C7H7+                   | C1H3N1     |
    /// | P     | 70                                                                                                                        | 70                          | 70                                             | 70                                       | 70                                                                    | 70                                      | 70.07                       | 70.07                        | 70.0657           | 70.0651                                                                  | 70.0657                                                             |                         |                         | C4H8N1+           | C4H8N1+                                                                  |                              | 13      |   70.0655 | C4H8N1+      |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             | 126                                            |                                          | 126                                                                   |                                         | 126.06                      |                              |                   | 126.055                                                                  |                                                                     |                         |                         |                   | C6H8N1O2+                                                                |                              | 5       |  126.0550 | C6H8N1O2+    | -55.9895 | C-2O-2       |                  |                         | C-2O-2     |
    /// | S     | 60                                                                                                                        | 60                          | 60                                             | 60                                       | 60                                                                    | 60                                      | 60.04                       | 60.04                        |                   | 60.0444                                                                  | 60.0449                                                             |                         |                         |                   | C2H6N1O1+                                                                |                              | 11      |   60.0447 | C2H6N1O1+    |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C2H4N1+                      | 1       |        ?? | C2H4N1+      | ??       | H2O1         |                  |                         | H2O1       |
    /// | T     | 74                                                                                                                        | 74                          | 74                                             | 74                                       | 74                                                                    | 74                                      |                             | 74.06                        | 74.0606           | 74.0600                                                                  | 74.0606                                                             |                         |                         | C3H8N1O1+         | C3H8N1O1+                                                                |                              | 12      |   74.0604 | C3H8N1O1+    |          |              |                  |                         |            |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          | C3H6N1O1+                    | 1       |        ?? | C3H6N1O1+    | ??       | H2N1         |                  |                         | H2N1       |
    /// | W     | 159                                                                                                                       | 159                         | 159                                            | 159                                      | 159                                                                   | 159                                     | 159.09                      | 159.09                       |                   | 159.0917                                                                 | 159.0922                                                            | 159.0917                | C10H11N2+               |                   | C10H11N2+                                                                |                              | 13      |  159.0919 | C10H11N2+    |          |              |                  |                         |            |
    /// |       |                                                                                                                           | 171                         |                                                |                                          |                                                                       | 171                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 2       |       171 |              | -11.9081 |              | H4O-1            |                         | H4O-1      |
    /// |       | 170                                                                                                                       | 170                         | 170                                            |                                          | 170                                                                   | 170                                     |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 5       |       170 |              | -10.9081 |              | H5O-1            |                         | H5O-1      |
    /// |       | 132                                                                                                                       |                             |                                                | 132                                      |                                                                       | 132                                     |                             | 132.08                       |                   |                                                                          | 132.0813                                                            |                         |                         |                   |                                                                          |                              | 5       |  132.0813 |              | 27.0106  |              | C1H1N1           |                         | C1H1N1     |
    /// |       | 130                                                                                                                       | 130                         | 130                                            | 130                                      | 130                                                                   | 130                                     |                             | 130.07                       |                   |                                                                          | 130.0657                                                            |                         |                         |                   |                                                                          |                              | 8       |  130.0657 |              | 29.0262  |              | C1H3N1           |                         | C1H3N1     |
    /// |       | 117                                                                                                                       | 117                         | 117                                            | 117                                      | 117                                                                   | 117                                     |                             |                              |                   |                                                                          | 117.0578                                                            |                         |                         |                   |                                                                          |                              | 7       |  117.0578 |              | 42.0341  |              | C2H4N1           |                         | C2H4N1     |
    /// |       | 100                                                                                                                       |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |       100 |              | 59.0919  |              | C3H9N1/C2H7N2    |                         |            |
    /// |       |                                                                                                                           |                             |                                                | 77                                       |                                                                       | 77                                      |                             |                              |                   |                                                                          | 77.0391                                                             |                         |                         |                   |                                                                          |                              | 3       |   77.0391 |              | 82.0528  |              | C4H6N2           | C6H5                    | C4H6N2     |
    /// |       | 11                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        11 |              | 148.0919 |              |                  |                         |            |
    /// | Y     | 136                                                                                                                       | 136                         | 136                                            | 136                                      | 136                                                                   | 136                                     | 136.08                      | 136.08                       | 136.0762          | 136.0757                                                                 | 136.0762                                                            | 136.0757                | C8H10N1O1+              | C8H10N1O1+        | C8H10N1O1+                                                               |                              | 15      |  136.0760 | C8H10N1O1+   |          |              |                  |                         |            |
    /// |       | 107                                                                                                                       | 107                         |                                                | 107                                      |                                                                       | 107                                     |                             |                              |                   |                                                                          | 107.0497                                                            |                         |                         |                   |                                                                          |                              | 5       |  107.0497 |              | 29.0263  |              | C1H3N1           |                         | C1H3N1     |
    /// |       | 91                                                                                                                        | 91                          |                                                | 91                                       |                                                                       | 91                                      |                             |                              |                   |                                                                          | 91.0548                                                             |                         |                         |                   |                                                                          |                              | 5       |   91.0548 |              | 45.0212  |              | C1H3N1O1         |                         | C1H3N1O1   |
    /// |       |                                                                                                                           |                             |                                                |                                          |                                                                       |                                         |                             |                              | 55.0184           |                                                                          |                                                                     |                         |                         | C3H3O1+           |                                                                          |                              | 2       |   55.0184 | C3H3O1+      | 81.0576  | C5H7N1       |                  |                         | C5H7N1     |
    /// | V     | 72                                                                                                                        | 72                          | 72                                             | 72                                       | 72                                                                    | 72                                      | 72.08                       | 72.08                        | 72.0813           | 72.0808                                                                  | 72.0813                                                             |                         |                         | C4H10N1+          | C4H10N1+                                                                 |                              | 13      |   72.0811 | C4H10N1+     |          |              |                  |                         |            |
    /// |       | 69                                                                                                                        |                             |                                                | 69                                       |                                                                       | 69                                      |                             |                              |                   |                                                                          | 69.0704                                                             |                         |                         |                   |                                                                          |                              | 4       |   69.0704 |              | 3.0107   |              | C1H1O-1          |                         | C1H1O-1    |
    /// |       | 55                                                                                                                        |                             |                                                | 55                                       |                                                                       | 55                                      |                             |                              |                   |                                                                          | 55.0548                                                             |                         |                         |                   |                                                                          |                              | 4       |   55.0548 |              | 17.0263  |              | H3N1             |                         | H3N1       |
    /// |       | 44                                                                                                                        |                             |                                                |                                          |                                                                       |                                         |                             |                              |                   |                                                                          |                                                                     |                         |                         |                   |                                                                          |                              | 1       |        44 |              | 28.0811  |              | C1H2N1           |                         | C1H2N1     |
    /// |       |                                                                                                                           |                             |                                                | 41                                       |                                                                       | 41                                      |                             |                              |                   |                                                                          | 41.0391                                                             |                         |                         |                   |                                                                          |                              | 3       |   41.0391 |              | 31.0420  |              | C1H5N1           |                         | C1H5N1     |
    fn immonium_losses(&self) -> Cow<'_, [NeutralLoss]> {
        // TODO: For B/Z there are common immonium ions, but the mass is the same (meaning the loss is different), find a way of representing that
        Cow::Owned(match self {
            Self::Arginine => vec![
                NeutralLoss::Gain(molecular_formula!(C 2 O 2)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 2)),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 2 H 2 N 2)),
                NeutralLoss::Loss(molecular_formula!(C 3 H 6 N 2)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 3)),
                NeutralLoss::Loss(molecular_formula!(C 3 H 4 N 2 O -1)),
                NeutralLoss::Loss(molecular_formula!(C 4 H 8 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 4 H 10 N 2)),
            ],
            Self::Asparagine => vec![NeutralLoss::Loss(molecular_formula!(H 3 N 1))],
            Self::AsparticAcid | Self::GlutamicAcid | Self::Serine => {
                vec![NeutralLoss::Loss(molecular_formula!(H 2 O 1))]
            }
            Self::Glutamine => vec![
                NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1 O 1)),
            ],
            Self::Histidine => vec![
                NeutralLoss::Gain(molecular_formula!(C 2 O 2)),
                NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
                NeutralLoss::Loss(molecular_formula!(H 3 O -1)),
                NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 2 N 1)),
            ],
            Self::Leucine | Self::Isoleucine | Self::AmbiguousLeucine => vec![
                NeutralLoss::Loss(molecular_formula!(C 1 H 2)),
                NeutralLoss::Loss(molecular_formula!(C 3 H 6)),
            ],
            Self::Lysine => vec![
                NeutralLoss::Gain(molecular_formula!(C 1 O 1)),
                NeutralLoss::Loss(molecular_formula!(C -2 H 1 N 1 O -1)),
                NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 2 H 7 N 1)),
            ],
            Self::Methionine => vec![
                NeutralLoss::Loss(molecular_formula!(H 2 S 1)),
                NeutralLoss::Loss(molecular_formula!(C 2 H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 4 S 1)),
            ],
            Self::Phenylalanine => vec![NeutralLoss::Gain(molecular_formula!(C 2 O 2))],
            Self::Threonine => vec![NeutralLoss::Loss(molecular_formula!(H 2 N 1))],
            Self::Tryptophan => vec![
                NeutralLoss::Loss(molecular_formula!(H 4 O -1)),
                NeutralLoss::Loss(molecular_formula!(H 5 O -1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 1 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 2 H 4 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 4 H 6 N 2)),
            ],
            Self::Tyrosine => vec![
                NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 3 N 1 O 1)),
                NeutralLoss::Loss(molecular_formula!(C 5 H 7 N 1)),
            ],
            Self::Valine => vec![
                NeutralLoss::Loss(molecular_formula!(C 1 H 1 O -1)),
                NeutralLoss::Loss(molecular_formula!(H 3 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 2 N 1)),
                NeutralLoss::Loss(molecular_formula!(C 1 H 5 N 1)),
            ],
            _ => Vec::new(),
        })
    }
}

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

        if ions.a.0 && allow_terminal.0 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications - molecular_formula!(H 1 C 1 O 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::a(n_pos),
                n_term,
                ions.a.1,
                charge_carriers,
                ions.a.2,
            ));
        }
        if ions.b.0 && allow_terminal.0 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications - molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::b(n_pos),
                n_term,
                ions.b.1,
                charge_carriers,
                ions.b.2,
            ));
        }
        if ions.c.0 && allow_terminal.0 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications + molecular_formula!(H 2 N 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::c(n_pos),
                n_term,
                ions.c.1,
                charge_carriers,
                ions.c.2,
            ));
        }
        if ions.d.0 && allow_terminal.0 {
            if let Some(satellite_ion_fragments) =
                self.satellite_ion_fragments(sequence_index, peptidoform_index)
            {
                base_fragments.extend(Fragment::generate_all(
                    &(-satellite_ion_fragments.as_ref()
                        * modifications
                        * self.formulas_inner(sequence_index, peptidoform_index)
                        + molecular_formula!(H 1 C 1 O 1)),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::d(n_pos),
                    n_term,
                    ions.d.1,
                    charge_carriers,
                    ions.d.2,
                ));
            }
        }
        if ions.v.0 && allow_terminal.1 {
            base_fragments.extend(Fragment::generate_all(
                &molecular_formula!(H 3 C 2 N 1 O 1).into(),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::v(c_pos),
                c_term,
                ions.v.1,
                charge_carriers,
                ions.v.2,
            ));
        }
        if ions.w.0 && allow_terminal.1 {
            if let Some(satellite_ion_fragments) =
                self.satellite_ion_fragments(sequence_index, peptidoform_index)
            {
                base_fragments.extend(Fragment::generate_all(
                    &(-satellite_ion_fragments.as_ref()
                        * modifications
                        * self.formulas_inner(sequence_index, peptidoform_index)
                        + molecular_formula!(H 2 N 1)),
                    peptidoform_ion_index,
                    peptidoform_index,
                    &FragmentType::w(c_pos),
                    c_term,
                    ions.w.1,
                    charge_carriers,
                    ions.w.2,
                ));
            }
        }
        if ions.x.0 && allow_terminal.1 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications + molecular_formula!(C 1 O 1) - molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::x(c_pos),
                c_term,
                ions.x.1,
                charge_carriers,
                ions.x.2,
            ));
        }
        if ions.y.0 && allow_terminal.1 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications + molecular_formula!(H 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::y(c_pos),
                c_term,
                ions.y.1,
                charge_carriers,
                ions.y.2,
            ));
        }
        if ions.z.0 && allow_terminal.1 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications - molecular_formula!(H 2 N 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::z(c_pos),
                c_term,
                ions.z.1,
                charge_carriers,
                ions.z.2,
            ));
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications - molecular_formula!(H 1 N 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::z·(c_pos),
                c_term,
                ions.z.1,
                charge_carriers,
                ions.z.2,
            ));
        }

        if ions.immonium.0 && allow_terminal.0 && allow_terminal.1 {
            base_fragments.extend(Fragment::generate_all(
                &(self.formulas_inner(sequence_index, peptidoform_index)
                    * (modifications - molecular_formula!(C 1 O 1))),
                peptidoform_ion_index,
                peptidoform_index,
                &FragmentType::Immonium(n_pos, self.into()), // TODO: get the actual sequenceelement here
                &Multi::default(),
                self.immonium_losses().as_ref(),
                charge_carriers,
                ions.immonium.1,
            ));
        }
        base_fragments
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
                aa.pro_forma_definition(),
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
