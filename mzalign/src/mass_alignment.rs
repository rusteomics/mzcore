use crate::{
    Alignment, align_matrix::Matrix, align_type::*, alignment::Score,
    diagonal_array::DiagonalArray, piece::*, scoring::*,
};
use mzcore::{
    chemistry::{MassMode, MolecularFormula},
    prelude::MultiChemical,
    quantities::{Multi, WithinTolerance},
    sequence::{AtMax, HasPeptidoform, Linear, Peptidoform, SequenceElement, SequencePosition},
    system::{Mass, dalton},
};

// TODO: potentially allow any gap to match to a list of aminoacids also if the mass difference is exactly a common modification
// eg X[mass(W[oxidation]A)] should match WA

/// Create an alignment of two peptides based on mass and homology.
/// The substitution matrix is in the exact same order as the definition of [`AminoAcid`](mzcore::sequence::AminoAcid).
/// The [`AlignScoring`] sets the rules and exact scores while scoring.
/// The [`AlignType`] controls the alignment behaviour, global/local or anything in between.
///
/// # Note
/// * Terminal modifications are added to the mass of the terminal amino acid.
/// * This ignores ambiguous aminoacids `(?AA)` because these kinds of errors are why the mass
///   alignment was built in the first place.
/// * Ambiguous modifications are handled by allowing those positions to have the mass with and
///   without the modification. It does not guarantee that the modification is placed exactly once.
/// * Labile modifications and charge carriers are ignored.
///
/// # Panics
/// It panics when the length of `seq_a` or `seq_b` is bigger than [`isize::MAX`].
pub fn align<const STEPS: u16, A: HasPeptidoform<Linear>, B: HasPeptidoform<Linear>>(
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
#[allow(clippy::similar_names)]
pub(super) fn align_cached<
    const STEPS: u16,
    A: HasPeptidoform<Linear>,
    B: HasPeptidoform<Linear>,
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
        let mut ranges: DiagonalArray<(Mass, Mass), STEPS> =
            DiagonalArray::new(peptidoform_a.len());
        for i in 0..peptidoform_a.len() {
            for j in 0..=i.min(STEPS as usize) {
                let (min, max) = unsafe { masses_a.get_unchecked([i, j]) }.iter().fold(
                    (
                        Mass::new::<dalton>(f64::INFINITY),
                        Mass::new::<dalton>(f64::NEG_INFINITY),
                    ),
                    |(min, max), &m| {
                        let range = scoring.tolerance.bounds(m);
                        (min.min(range.0), max.max(range.1))
                    },
                );
                ranges[[i, j]] = (min, max);
            }
        }
        ranges
    };

    let ranges_b = {
        let mut ranges: DiagonalArray<(Mass, Mass), STEPS> =
            DiagonalArray::new(peptidoform_b.len());
        for i in 0..peptidoform_b.len() {
            for j in 0..=i.min(STEPS as usize) {
                let (min, max) = unsafe { masses_b.get_unchecked([i, j]) }.iter().fold(
                    (
                        Mass::new::<dalton>(f64::INFINITY),
                        Mass::new::<dalton>(f64::NEG_INFINITY),
                    ),
                    |(min, max), &m| (min.min(m), max.max(m)),
                );
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

                let len_a = u16::from(gap_a);
                let len_b = u16::from(!gap_a);
                Piece::new(prev.score + score, score, MatchType::Gap, len_a, len_b)
            };

            // First try all gap possibilities.
            let gap_score_a = score_gap(true);
            let gap_score_b = score_gap(false);

            let mut highest = if gap_score_a.score >= gap_score_b.score {
                gap_score_a
            } else {
                gap_score_b
            };

            // Now try matching single aminoacids.
            let prev = unsafe { matrix.get_unchecked([index_a - 1, index_b - 1]) };
            let pair_score = score_pair(
                (
                    unsafe { peptidoform_a.sequence().get_unchecked(index_a - 1) },
                    unsafe { masses_a.get_unchecked([index_a - 1, 0]) },
                ),
                (
                    unsafe { peptidoform_b.sequence().get_unchecked(index_b - 1) },
                    unsafe { masses_b.get_unchecked([index_b - 1, 0]) },
                ),
                scoring,
                prev.score,
            );
            if pair_score.score > highest.score {
                highest = pair_score;
            }

            // Now try matching longer sequences.
            if highest.match_type != MatchType::FullIdentity {
                for len_a in 1..=index_a.min(STEPS as usize) {
                    let range_a = unsafe { ranges_a.get_unchecked([index_a - 1, len_a - 1]) };

                    let min_len_b = if len_a == 1 { 2 } else { 1 };

                    for len_b in min_len_b..=index_b.min(STEPS as usize) {
                        let range_b = unsafe { ranges_b.get_unchecked([index_b - 1, len_b - 1]) };
                        // Note that ranges are already expanded by tolerance, so
                        // exact comparison is fine here.
                        if range_a.0 > range_b.1 || range_b.0 > range_a.1 {
                            continue;
                        }
                        // len_a and b are always <= STEPS
                        let match_score = {
                            let prev =
                                unsafe { matrix.get_unchecked([index_a - len_a, index_b - len_b]) };
                            let base_score = prev.score;

                            score(
                                unsafe {
                                    (
                                        peptidoform_a
                                            .sequence()
                                            .get_unchecked((index_a - len_a)..index_a),
                                        masses_a.get_unchecked([index_a - 1, len_a - 1]),
                                    )
                                },
                                unsafe {
                                    (
                                        peptidoform_b
                                            .sequence()
                                            .get_unchecked((index_b - len_b)..index_b),
                                        masses_b.get_unchecked([index_b - 1, len_b - 1]),
                                    )
                                },
                                scoring,
                                base_score,
                            )
                        };
                        if let Some(p) = match_score
                            && p.score > highest.score
                        {
                            highest = p;
                        }
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
    let maximal_score = isize::midpoint(
        seq_a.sequence()[start_a..start_a + path.iter().map(|p| p.step_a as usize).sum::<usize>()]
            .iter()
            .map(|a| {
                scoring.matrix[a.aminoacid.aminoacid() as usize][a.aminoacid.aminoacid() as usize]
                    as isize
            })
            .sum::<isize>(),
        seq_b.sequence()[start_b..start_b + path.iter().map(|p| p.step_b as usize).sum::<usize>()]
            .iter()
            .map(|a| {
                scoring.matrix[a.aminoacid.aminoacid() as usize][a.aminoacid.aminoacid() as usize]
                    as isize
            })
            .sum::<isize>(),
    );
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
pub(super) fn score_pair<A: AtMax<Linear>, B: AtMax<Linear>>(
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
                let local = scoring.matrix[a.0.aminoacid.aminoacid() as usize]
                    [b.0.aminoacid.aminoacid() as usize] as isize
                    + scoring.mismatch as isize;
                Piece::new(score + local, local, MatchType::Mismatch, 1, 1)
            }
        }
        (false, true) => Piece::new(
            score + scoring.mass_base as isize + scoring.isobaric as isize,
            scoring.mass_base as isize + scoring.isobaric as isize,
            MatchType::Isobaric,
            1,
            1,
        ),
        (false, false) => {
            let local = scoring.matrix[a.0.aminoacid.aminoacid() as usize]
                [b.0.aminoacid.aminoacid() as usize] as isize
                + scoring.mismatch as isize;
            Piece::new(score + local, local, MatchType::Mismatch, 1, 1)
        }
    }
}

/// Score two sets of aminoacids (it will only be called when at least one of a and b has len > 1)
/// Returns none if no sensible explanation can be made
pub(super) fn score<A: AtMax<Linear>, B: AtMax<Linear>>(
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

/// Get the masses of all sequence elements. This also adds the N and C terminal modifications.
/// # Panics
/// If the global isotope modifications are invalid
pub(super) fn calculate_masses<const STEPS: u16>(
    sequence: &Peptidoform<impl AtMax<Linear>>,
    mass_mode: MassMode,
) -> DiagonalArray<Multi<Mass>, STEPS> {
    let mut array = DiagonalArray::new(sequence.len());
    let n = sequence
        .get_n_term()
        .iter()
        .map(mzcore::sequence::Modification::formula)
        .sum::<MolecularFormula>();
    let c = sequence
        .get_c_term()
        .iter()
        .map(mzcore::sequence::Modification::formula)
        .sum::<MolecularFormula>();
    for i in 0..sequence.len() {
        for j in 0..=i.min(STEPS as usize) {
            let mut seq = sequence.sequence()[i - j..=i]
                .iter()
                .map(|p| p.formulas_inner(SequencePosition::Index(i), 0, 0))
                .sum::<Multi<MolecularFormula>>();
            if i - j == 0 {
                seq += n.clone();
            }
            if i == sequence.len() - 1 {
                seq += c.clone();
            }
            array[[i, j]] = seq
                .iter()
                .map(|f| {
                    f.with_global_isotope_modifications(mzcore::sequence::HiddenInternalMethods::get_global(sequence))
                        .expect("Invalid global isotope modifications while calculating masses for alignment")
                        .mass(mass_mode)
                })
                .collect();
        }
    }
    array
}

#[cfg(test)]
#[expect(clippy::missing_panics_doc)]
mod tests {
    use super::score;
    use crate::scoring::AlignScoring;
    use mzcore::{
        chemistry::MolecularFormula,
        prelude::MultiChemical,
        quantities::Multi,
        sequence::{CheckedAminoAcid, SequenceElement},
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
                    .map(MultiChemical::formulas)
                    .sum::<Multi<MolecularFormula>>()[0]
                    .monoisotopic_mass()
                    .into()
            ),
            (
                &b,
                &b.iter()
                    .map(MultiChemical::formulas)
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
