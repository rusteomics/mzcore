use std::fmt::Debug;

use crate::{
    align::{
        align_type::*, alignment::Score, diagonal_array::DiagonalArray, piece::*, scoring::*, Alignment
    },
    annotation::model::GlycanModel,
    chemistry::{MassMode, MolecularFormula},
    quantities::{Multi, WithinTolerance},
    sequence::{
        AtMax, HasPeptidoform, Peptidoform, SequenceElement, SequencePosition, SimpleLinear,
    },
    system::{dalton, Mass},
};

// TODO: no way of handling terminal modifications yet
// TODO: potentially allow any gap to match to a list of aminoacids also if the mass difference is exactly a common modification
// eg X[mass(W[oxidation]A)] should match WA

/// Create an alignment of two peptides based on mass and homology.
/// The substitution matrix is in the exact same order as the definition of [`AminoAcid`](crate::sequence::AminoAcid).
/// The [`AlignScoring`] sets the rules and exact scores while scoring.
/// The [`AlignType`] controls the alignment behaviour, global/local or anything in between.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger than [`isize::MAX`].
pub fn align<const STEPS: u16, A: HasPeptidoform<SimpleLinear>, B: HasPeptidoform<SimpleLinear>>(
    seq_a: A,
    seq_b: B,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
) -> Alignment<A, B> {
    let peptidoform_a = seq_a.cast_peptidoform();
    let peptidoform_b = seq_b.cast_peptidoform();
    let masses_a: DiagonalArray<Multi<Mass>, STEPS> =
        calculate_masses::<STEPS>(peptidoform_a, scoring.mass_mode);
    let masses_b: DiagonalArray<Multi<Mass>, STEPS> =
        calculate_masses::<STEPS>(peptidoform_b, scoring.mass_mode);
    align_cached::<STEPS, A, B>(seq_a, &masses_a, seq_b, &masses_b, scoring, align_type)
}

/// Do mass based alignment, but with precomputed masses, see [align] for the normal entry point variant.
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger than [`isize::MAX`].
#[expect(clippy::too_many_lines)]
pub(super) fn align_cached<
    const STEPS: u16,
    A: HasPeptidoform<SimpleLinear>,
    B: HasPeptidoform<SimpleLinear>,
>(
    seq_a: A,
    masses_a: &DiagonalArray<Multi<Mass>, STEPS>,
    seq_b: B,
    masses_b: &DiagonalArray<Multi<Mass>, STEPS>,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
) -> Alignment<A, B> {
    let peptidoform_a = seq_a.cast_peptidoform();
    let peptidoform_b = seq_b.cast_peptidoform();
    assert!(isize::try_from(peptidoform_a.len()).is_ok());
    assert!(isize::try_from(peptidoform_b.len()).is_ok());

    let mut matrix = Matrix::new(peptidoform_a.len(), peptidoform_b.len());
    let mut global_highest = (0, 0, 0);

    if align_type.left.global_a() {
        matrix.global_start(true, scoring);
    }
    if align_type.left.global_b() {
        matrix.global_start(false, scoring);
    }

    // These two arrays serve for optimizational purposes.
    // In `score` routine equality of masses is checked
    // between masses_a[index_a - len_a .. index_a] and masses_b[index_b - len_b .. index_b]
    // ranges which is expensive.
    // However, in most of the cases masses are very different there and the check fails.
    // We precompute two following arrays:
    // - ranges_a[index_a - len_a .. index_a] = (min_mass, max_mass)
    // min_mass is the minimal mass that matches at least one mass in masses_a[index_a - len_a .. index_a].
    // max_mass is the same for maximal mass.
    // - ranges_b[index_b - len_b .. index_b] = (min_mass, max_mass)
    // min_mass in the minimal mass across masses_b[index_b - len_b .. index_b].
    // max_mass is the same for maximal mass.
    // Then, before calling `score` we can perform a quick check whether
    // the ranges overlap for given slices. If they do not overlap, we skip calling `score`
    // as it is guaranteed that there are no two masses that are within tolerance.
    let ranges_a = {
        let mut ranges: DiagonalArray<(Mass, Mass), STEPS> = DiagonalArray::new(peptidoform_a.len());
        for i in 0..peptidoform_a.len() {
            for j in 0..=i.min(STEPS as usize) {
                let (min, max) = unsafe { masses_a.get_unchecked([i, j]) }
                    .iter()
                    .fold((Mass::new::<dalton>(f64::INFINITY), Mass::new::<dalton>(f64::NEG_INFINITY)), |(min, max), &m| {
                        let range = scoring.tolerance.bounds(m);
                        (min.min(range.0), max.max(range.1))
                    });
                ranges[[i, j]] = (min, max);
            }
        }
        ranges
    };

    let ranges_b = {
        let mut ranges: DiagonalArray<(Mass, Mass), STEPS> = DiagonalArray::new(peptidoform_b.len());
        for i in 0..peptidoform_b.len() {
            for j in 0..=i.min(STEPS as usize) {
                let (min, max) = unsafe { masses_b.get_unchecked([i, j]) }
                    .iter()
                    .fold((Mass::new::<dalton>(f64::INFINITY), Mass::new::<dalton>(f64::NEG_INFINITY)), |(min, max), &m| {
                        (min.min(m), max.max(m))
                    });
                ranges[[i, j]] = (min, max);
            }
        }
        ranges
    };

    for index_a in 1..=peptidoform_a.len() {
        for index_b in 1..=peptidoform_b.len() {
            // Returns the score for a gap transition.
            // gap_a controls whether the gap is in a or b.
            let score_gap = |gap_a: bool| {
                let prev = if gap_a {
                    unsafe { matrix.get_unchecked([index_a - 1, index_b]) }
                } else {
                    unsafe { matrix.get_unchecked([index_a, index_b - 1]) }
                };

                let is_first_step = prev.step_a == 0 && prev.step_b == 0;
                let is_previous_gap = prev.step_a == 0 && !gap_a || prev.step_b == 0 && gap_a;
                let is_gap_start = is_first_step || !is_previous_gap;
                // First check the score to be used for affine gaps
                let score = scoring.gap_extend as isize
                    + scoring.gap_start as isize * isize::from(is_gap_start);

                let len_a = if gap_a { 1 } else { 0 };
                let len_b = if gap_a { 0 } else { 1 };
                Piece::new(
                    prev.score + score,
                    score,
                    MatchType::Gap,
                    len_a,
                    len_b,
                )
            };

            // First try all gap possibilities.
            let gap_score_a = score_gap(true);
            let gap_score_b = score_gap(false);

            let mut highest = if gap_score_a.score >= gap_score_b.score {
                gap_score_a
            } else {
                gap_score_b
            };

            for len_a in 1..=index_a.min(STEPS as usize) {
                let range_a = unsafe { ranges_a.get_unchecked([index_a - 1, len_a - 1]) };

                for len_b in 1..=index_b.min(STEPS as usize) {
                    let range_b = unsafe { ranges_b.get_unchecked([index_b - 1, len_b - 1]) };

                    let prev = unsafe { matrix.get_unchecked([index_a - len_a, index_b - len_b]) };
                    let base_score = prev.score;

                    // len_a and b are always <= STEPS
                    let piece = if len_a == 1 && len_b == 1 {
                        Some(score_pair(
                            unsafe {
                                (
                                    peptidoform_a.sequence().get_unchecked(index_a - 1),
                                    masses_a.get_unchecked([index_a - 1, 0]),
                                )
                            },
                            unsafe {
                                (
                                    peptidoform_b.sequence().get_unchecked(index_b - 1),
                                    masses_b.get_unchecked([index_b - 1, 0]),
                                )
                            },
                            scoring,
                            base_score,
                        ))
                    // Ranges do not overlap, skip scoring.
                    } else if range_a.0 > range_b.1 || range_b.0 > range_a.1 {
                        None
                    } else {
                        score(
                            unsafe {
                                (
                                    peptidoform_a
                                        .sequence()
                                        .get_unchecked((index_a - len_a)..index_a),
                                    masses_a.get_unchecked([index_a - 1, len_a - 1])
                                )
                            },
                            unsafe {
                                (
                                    peptidoform_b
                                        .sequence()
                                        .get_unchecked((index_b - len_b)..index_b),
                                    masses_b.get_unchecked([index_b - 1, len_b - 1])
                                )
                            },
                            scoring,
                            base_score,
                        )
                    };
                    if let Some(p) = piece && p.score > highest.score {
                        highest = p;
                    }
                }
            }
            if highest.score >= global_highest.0 {
                global_highest = (highest.score, index_a, index_b);
            }
            if align_type.left.global() || highest.score > 0 {
                unsafe {
                    *matrix.get_unchecked_mut([index_a, index_b]) = highest;
                }
            }
        }
    }
    let (start_a, start_b, path) = matrix.trace_path(align_type, global_highest);
    let score = determine_final_score(
        peptidoform_a,
        peptidoform_b,
        start_a,
        start_b,
        &path,
        scoring,
    );

    Alignment {
        seq_a,
        seq_b,
        score,
        path,
        start_a,
        start_b,
        align_type,
        maximal_step: STEPS,
    }
}

pub(super) fn determine_final_score<A, B>(
    seq_a: &Peptidoform<A>,
    seq_b: &Peptidoform<B>,
    start_a: usize,
    start_b: usize,
    path: &[Piece],
    scoring: AlignScoring<'_>,
) -> Score {
    let maximal_score = isize::midpoint(seq_a.sequence()
        [start_a..start_a + path.iter().map(|p| p.step_a as usize).sum::<usize>()]
        .iter()
        .map(|a| {
            scoring.matrix[a.aminoacid.aminoacid() as usize][a.aminoacid.aminoacid() as usize]
                as isize
        })
        .sum::<isize>(), seq_b.sequence()
            [start_b..start_b + path.iter().map(|p| p.step_b as usize).sum::<usize>()]
            .iter()
            .map(|a| {
                scoring.matrix[a.aminoacid.aminoacid() as usize][a.aminoacid.aminoacid() as usize]
                    as isize
            })
            .sum::<isize>());
    let absolute_score = path.last().map(|p| p.score).unwrap_or_default();
    Score {
        absolute: absolute_score,
        normalised: if maximal_score == 0 {
            ordered_float::OrderedFloat::default()
        } else {
            ordered_float::OrderedFloat((absolute_score as f64 / maximal_score as f64).min(1.0))
        },
        max: maximal_score,
    }
}

/// Score a pair of sequence elements (AA + mods)
pub(super) fn score_pair<A: AtMax<SimpleLinear>, B: AtMax<SimpleLinear>>(
    a: (&SequenceElement<A>, &Multi<Mass>),
    b: (&SequenceElement<B>, &Multi<Mass>),
    scoring: AlignScoring<'_>,
    score: isize,
) -> Piece {
    match (
        a.0.aminoacid.aminoacid() == b.0.aminoacid.aminoacid(),
        scoring.tolerance.within(a.1, b.1),
    ) {
        (true, true) => {
            let local = scoring.matrix[a.0.aminoacid.aminoacid() as usize]
                [b.0.aminoacid.aminoacid() as usize] as isize;
            Piece::new(score + local, local, MatchType::FullIdentity, 1, 1)
        }
        (true, false) => {
            // The assumption is that if the peptide has modifications and the mass do not match
            // this element in the database this is caused by artefacts of some kind. While if
            // there is a modification on the database this is encoded in the genome (or similar)
            // and has to be present in order for the peptide to match properly.
            //
            // If both have modifications this is assumed to be caused by artefacts on the
            // peptides side.
            if (scoring.pair == PairMode::DatabaseToPeptidoform && !b.0.modifications.is_empty())
                || (scoring.pair == PairMode::PeptidoformToDatabase
                    && !a.0.modifications.is_empty())
            {
                let local = scoring.matrix[a.0.aminoacid.aminoacid() as usize]
                    [b.0.aminoacid.aminoacid() as usize] as isize
                    + scoring.mass_mismatch as isize;
                Piece::new(score + local, local, MatchType::IdentityMassMismatch, 1, 1)
            } else {
                Piece::new(
                    score + scoring.mismatch as isize,
                    scoring.mismatch as isize,
                    MatchType::Mismatch,
                    1,
                    1,
                )
            }
        }
        (false, true) => Piece::new(
            score + scoring.mass_base as isize + scoring.isobaric as isize,
            scoring.mass_base as isize + scoring.isobaric as isize,
            MatchType::Isobaric,
            1,
            1,
        ),
        (false, false) => Piece::new(
            score + scoring.mismatch as isize,
            scoring.mismatch as isize,
            MatchType::Mismatch,
            1,
            1,
        ),
    }
}

/// Score two sets of aminoacids (it will only be called when at least one of a and b has len > 1)
/// Returns none if no sensible explanation can be made
fn score<A: AtMax<SimpleLinear>, B: AtMax<SimpleLinear>>(
    a: (&[SequenceElement<A>], &Multi<Mass>),
    b: (&[SequenceElement<B>], &Multi<Mass>),
    scoring: AlignScoring<'_>,
    score: isize,
) -> Option<Piece> {
    if scoring.tolerance.within(a.1, b.1) {
        let rotated = {
            a.0.len() == b.0.len() && {
                let mut b_copy = vec![false; b.0.len()];
                a.0.iter().all(|el| {
                    b_copy
                        .iter()
                        .enumerate()
                        .position(|(index, used)| !used && b.0[index] == *el)
                        .is_some_and(|pos| {
                            b_copy[pos] = true;
                            true
                        })
                })
            }
        };
        #[expect(clippy::cast_possible_wrap)]
        let local = scoring.mass_base as isize
            + if rotated {
                scoring.rotated as isize * a.0.len() as isize
            } else {
                scoring.isobaric as isize * (a.0.len() + b.0.len()) as isize / 2
            };
        Some(Piece::new(
            score + local,
            local,
            if rotated {
                MatchType::Rotation
            } else {
                MatchType::Isobaric
            },
            a.0.len() as u16,
            b.0.len() as u16,
        ))
    } else {
        None
    }
}

/// Get the masses of all sequence elements
pub(super) fn calculate_masses<const STEPS: u16>(
    sequence: &Peptidoform<impl AtMax<SimpleLinear>>,
    mass_mode: MassMode,
) -> DiagonalArray<Multi<Mass>, STEPS> {
    let mut array = DiagonalArray::new(sequence.len());
    for i in 0..sequence.len() {
        for j in 0..=i.min(STEPS as usize) {
            array[[i, j]] = sequence.sequence()[i - j..=i]
                .iter()
                .map(|p| {
                    p.formulas_all(
                        &[],
                        &[],
                        &mut Vec::new(),
                        false,
                        SequencePosition::Index(i),
                        0,
                        0,
                        &GlycanModel::DISALLOW,
                    )
                    .0
                })
                .sum::<Multi<MolecularFormula>>()
                .iter()
                .map(|f| f.mass(mass_mode))
                .collect();
        }
    }
    array
}

struct Matrix {
    value: Vec<Vec<Piece>>,
    a: usize,
    b: usize,
}

impl Debug for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use std::fmt::Write;
        for column in &self.value {
            let mut line_0 = String::new();
            let mut line_1 = String::new();
            for cell in column {
                let top = format!("{}/{} {:2}", cell.step_a, cell.step_b, cell.local_score);
                let bottom = format!(
                    "{} {:3}",
                    match cell.match_type {
                        MatchType::FullIdentity => "FI",
                        MatchType::Gap => "G ",
                        MatchType::IdentityMassMismatch => "IM",
                        MatchType::Isobaric => "I ",
                        MatchType::Rotation => "R ",
                        MatchType::Mismatch => "M ",
                    },
                    cell.score
                );
                write!(&mut line_0, "⎡{top:0$}⎤", top.len().max(bottom.len()))?;
                write!(&mut line_1, "⎣{bottom:0$}⎦", top.len().max(bottom.len()))?;
            }
            writeln!(f, "{line_0}")?;
            writeln!(f, "{line_1}")?;
        }
        Ok(())
    }
}

impl Matrix {
    pub(super) fn new(a: usize, b: usize) -> Self {
        Self {
            value: vec![vec![Piece::default(); b + 1]; a + 1],
            a,
            b,
        }
    }

    #[expect(clippy::cast_possible_wrap)]
    pub(super) fn global_start(&mut self, is_a: bool, scoring: AlignScoring<'_>) {
        let max = if is_a { self.a } else { self.b };
        for index in 0..=max {
            self.value[if is_a { index } else { 0 }][if is_a { 0 } else { index }] = Piece::new(
                match index {
                    0 => 0,
                    _ => {
                        scoring.gap_start as isize + (index as isize) * scoring.gap_extend as isize
                    }
                },
                match index {
                    0 => 0,
                    1 => scoring.gap_start as isize + scoring.gap_extend as isize,
                    _ => scoring.gap_extend as isize,
                },
                MatchType::Gap,
                if is_a { u16::from(index != 0) } else { 0 },
                if is_a { 0 } else { u16::from(index != 0) },
            );
        }
    }

    pub(super) fn trace_path(
        &self,
        ty: AlignType,
        high: (isize, usize, usize),
    ) -> (usize, usize, Vec<Piece>) {
        let mut path = Vec::new();
        let mut high = self.find_end(ty, high);

        // Loop back to left side
        while ty.left.global() || !(high.1 == 0 && high.2 == 0) {
            let value = self.value[high.1][high.2];
            if value.step_a == 0 && value.step_b == 0 || !ty.left.global() && value.score < 0 {
                break;
            }
            high = (
                0,
                high.1 - value.step_a as usize,
                high.2 - value.step_b as usize,
            );
            path.push(value);
        }
        (high.1, high.2, path.into_iter().rev().collect())
    }

    fn find_end(&self, ty: AlignType, high: (isize, usize, usize)) -> (isize, usize, usize) {
        if ty.right.global_a() && ty.right.global_a() {
            (self.value[self.a][self.b].score, self.a, self.b)
        } else if ty.right.global_b() {
            let value = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, value.0, self.b)
        } else if ty.right.global_a() {
            let value = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            (value.1, self.a, value.0)
        } else if ty.right.global() {
            let value_a = (0..=self.a)
                .map(|v| (v, self.value[v][self.b].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            let value_b = (0..=self.b)
                .map(|v| (v, self.value[self.a][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            if value_a.1 >= value_b.1 {
                (value_a.1, value_a.0, self.b)
            } else {
                (value_b.1, self.a, value_b.0)
            }
        } else {
            high
        }
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`Vec::get_unchecked`].
    /// A debug assertion hold up this promise on debug builds.
    pub(super) unsafe fn get_unchecked(&self, index: [usize; 2]) -> &Piece {
        debug_assert!(self.value.len() > index[0]);
        debug_assert!(self.value[index[0]].len() > index[1]);
        unsafe { self.value.get_unchecked(index[0]).get_unchecked(index[1]) }
    }

    /// # Safety
    /// This function assumes the index to be valid. Not upholding this does an out of bounds unsafe [`Vec::get_unchecked_mut`].
    /// A debug assertion hold up this promise on debug builds.
    pub(super) unsafe fn get_unchecked_mut(&mut self, index: [usize; 2]) -> &mut Piece {
        debug_assert!(self.value.len() > index[0]);
        debug_assert!(self.value[index[0]].len() > index[1]);
        unsafe {
            self.value
                .get_unchecked_mut(index[0])
                .get_unchecked_mut(index[1])
        }
    }
}

impl std::ops::Index<[usize; 2]> for Matrix {
    type Output = Piece;
    fn index(&self, index: [usize; 2]) -> &Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &self.value[index[0]][index[1]]
    }
}
impl std::ops::IndexMut<[usize; 2]> for Matrix {
    fn index_mut(&mut self, index: [usize; 2]) -> &mut Self::Output {
        assert!(index[0] <= self.a + 1);
        assert!(index[1] <= self.b + 1);
        &mut self.value[index[0]][index[1]]
    }
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::score;
    use crate::{
        align::scoring::AlignScoring,
        annotation::model::GlycanModel,
        chemistry::MolecularFormula,
        quantities::Multi,
        sequence::{CheckedAminoAcid, SequenceElement, SequencePosition},
    };

    #[test]
    fn pair() {
        let a = [SequenceElement::new(CheckedAminoAcid::N, None)];
        let b = [
            SequenceElement::new(CheckedAminoAcid::G, None),
            SequenceElement::new(CheckedAminoAcid::G, None),
        ];
        let pair = dbg!(score(
            (
                &a,
                &a.iter()
                    .map(|p| p
                        .formulas_all(
                            &[],
                            &[],
                            &mut Vec::new(),
                            false,
                            SequencePosition::default(),
                            0,
                            0,
                            &GlycanModel::DISALLOW,
                        )
                        .0)
                    .sum::<Multi<MolecularFormula>>()[0]
                    .monoisotopic_mass()
                    .into()
            ),
            (
                &b,
                &b.iter()
                    .map(|p| p
                        .formulas_all(
                            &[],
                            &[],
                            &mut Vec::new(),
                            false,
                            SequencePosition::default(),
                            0,
                            0,
                            &GlycanModel::DISALLOW,
                        )
                        .0)
                    .sum::<Multi<MolecularFormula>>()[0]
                    .monoisotopic_mass()
                    .into()
            ),
            AlignScoring::default(),
            0,
        ));
        assert!(pair.is_some());
    }
}
