use crate::{
    annotation::model::{ChargeRange, FragmentationModel, get_all_sidechain_losses},
    fragment::{NeutralLoss, PeptidePosition},
    sequence::{AminoAcid, Peptidoform, SequencePosition},
};

/// The possibilities for primary ions, a list of all allowed neutral losses, all charge options, and all allowed variant ions
pub type PossiblePrimaryIons<'a> = (Vec<Vec<NeutralLoss>>, ChargeRange, &'a [i8]);
/// The possibilities for satellite ions, a list of all satellite ions with their amino acid and
/// distance from the parent backbone cleavage, as well as all ion settings as for primary series.
pub type PossibleSatelliteIons<'a> = (Vec<(AminoAcid, u8)>, PossiblePrimaryIons<'a>);

/// A struct to handle all possible fragments that could be generated on a single location
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
#[non_exhaustive]
pub struct PossibleIons<'a> {
    /// a series ions
    pub a: Option<PossiblePrimaryIons<'a>>,
    /// b series ions
    pub b: Option<PossiblePrimaryIons<'a>>,
    /// c series ions
    pub c: Option<PossiblePrimaryIons<'a>>,
    /// d series ions (side chain fragmentation from a)
    pub d: PossibleSatelliteIons<'a>,
    /// v series ions (full side chain broken off from y)
    pub v: PossibleSatelliteIons<'a>,
    /// w series ions (side chain fragmentation from z)
    pub w: PossibleSatelliteIons<'a>,
    /// x series ions
    pub x: Option<PossiblePrimaryIons<'a>>,
    /// y series ions
    pub y: Option<PossiblePrimaryIons<'a>>,
    /// z series ions
    pub z: Option<PossiblePrimaryIons<'a>>,
    /// immonium
    pub immonium: Option<(ChargeRange, &'a [(Vec<AminoAcid>, Vec<NeutralLoss>)])>,
}

impl PossibleIons<'_> {
    /// Give an upper bound for the number of theoretical fragment for these possible ions
    pub fn size_upper_bound(&self) -> usize {
        self.a
            .as_ref()
            .map(|o| (o.0.len() + 1) * o.2.len())
            .unwrap_or_default()
            + self
                .b
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .c
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self.d.0.len() * 2 * (self.d.1.0.len() + 1) * self.d.1.2.len()
            + self.v.0.len() * (self.v.1.0.len() + 1) * self.v.1.2.len()
            + self.w.0.len() * 2 * (self.w.1.0.len() + 1) * self.w.1.2.len()
            + self
                .x
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .y
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + self
                .z
                .as_ref()
                .map(|o| (o.0.len() + 1) * o.2.len())
                .unwrap_or_default()
            + usize::from(self.immonium.is_some())
    }
}

impl FragmentationModel {
    /// Give all possible ions for the given position.
    /// # Panics
    /// If the position is a terminal position.
    pub fn ions<Complexity>(
        &self,
        position: PeptidePosition,
        peptidoform: &Peptidoform<Complexity>,
    ) -> PossibleIons<'_> {
        let SequencePosition::Index(sequence_index) = position.sequence_index else {
            panic!("Not allowed to call possible with a terminal PeptidePosition")
        };

        let get_neutral_losses = |neutral_losses: &Vec<NeutralLoss>,
                                  amino_acid_specific: &Vec<(Vec<AminoAcid>, Vec<NeutralLoss>)>,
                                  amino_acid_side_chains: &(u8, Option<Vec<AminoAcid>>),
                                  c_terminal| {
            let peptide_slice = &peptidoform[if c_terminal {
                sequence_index..peptidoform.len()
            } else {
                0..sequence_index
            }];
            neutral_losses
                .iter()
                .chain(
                    peptide_slice
                        .iter()
                        .flat_map(|seq| {
                            amino_acid_specific.iter().filter_map(|(rule, loss)| {
                                rule.contains(&seq.aminoacid.aminoacid()).then_some(loss)
                            })
                        })
                        .flatten(),
                )
                .map(|l| vec![l.clone()])
                .chain(get_all_sidechain_losses(
                    peptide_slice,
                    amino_acid_side_chains,
                ))
                .collect()
        };

        let c_position = position.flip_terminal();

        PossibleIons {
            a: self.a.location.possible(position).then_some((
                get_neutral_losses(
                    &self.a.neutral_losses,
                    &self.a.amino_acid_neutral_losses,
                    &self.a.amino_acid_side_chain_losses,
                    false,
                ),
                self.a.charge_range,
                self.a.allowed_variants.as_slice(),
            )),
            b: self.b.location.possible(position).then_some((
                get_neutral_losses(
                    &self.b.neutral_losses,
                    &self.b.amino_acid_neutral_losses,
                    &self.b.amino_acid_side_chain_losses,
                    false,
                ),
                self.b.charge_range,
                self.b.allowed_variants.as_slice(),
            )),
            c: self.c.location.possible(position).then_some((
                get_neutral_losses(
                    &self.c.neutral_losses,
                    &self.c.amino_acid_neutral_losses,
                    &self.c.amino_acid_side_chain_losses,
                    false,
                ),
                self.c.charge_range,
                self.c.allowed_variants.as_slice(),
            )),
            d: (
                self.d.location.possible(position, peptidoform, false),
                (
                    get_neutral_losses(
                        &self.d.neutral_losses,
                        &self.d.amino_acid_neutral_losses,
                        &self.d.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.d.charge_range,
                    self.d.allowed_variants.as_slice(),
                ),
            ),
            v: (
                self.v.location.possible(c_position, peptidoform, true),
                (
                    get_neutral_losses(
                        &self.v.neutral_losses,
                        &self.v.amino_acid_neutral_losses,
                        &self.v.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.v.charge_range,
                    self.v.allowed_variants.as_slice(),
                ),
            ),
            w: (
                self.w.location.possible(c_position, peptidoform, true),
                (
                    get_neutral_losses(
                        &self.w.neutral_losses,
                        &self.w.amino_acid_neutral_losses,
                        &self.w.amino_acid_side_chain_losses,
                        false,
                    ),
                    self.w.charge_range,
                    self.w.allowed_variants.as_slice(),
                ),
            ),
            x: self.x.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.x.neutral_losses,
                    &self.x.amino_acid_neutral_losses,
                    &self.x.amino_acid_side_chain_losses,
                    false,
                ),
                self.x.charge_range,
                self.x.allowed_variants.as_slice(),
            )),
            y: self.y.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.y.neutral_losses,
                    &self.y.amino_acid_neutral_losses,
                    &self.y.amino_acid_side_chain_losses,
                    false,
                ),
                self.y.charge_range,
                self.y.allowed_variants.as_slice(),
            )),
            z: self.z.location.possible(c_position).then_some((
                get_neutral_losses(
                    &self.z.neutral_losses,
                    &self.z.amino_acid_neutral_losses,
                    &self.z.amino_acid_side_chain_losses,
                    false,
                ),
                self.z.charge_range,
                self.z.allowed_variants.as_slice(),
            )),
            immonium: self.immonium.as_ref().map(|(c, l)| (*c, l.as_slice())),
        }
    }
}
