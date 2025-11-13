use std::{any::TypeId, borrow::Cow, fmt::Display, ops::RangeBounds};

use itertools::Itertools;

use crate::{AlignType, Alignment, MatchType, Piece};
use mzcore::{
    chemistry::{Chemical, MolecularFormula, MultiChemical},
    prelude::IsAminoAcid,
    sequence::{AtMax, Linear, Peptidoform, SequenceElement},
    system::Mass,
};
use std::fmt::Write;

/// An alignment of multiple peptides
#[derive(Debug, Clone, PartialEq)]
pub struct MultipleSequenceAlignment<Complexity> {
    /// The sequences
    pub sequences: Vec<MSAPlacement<Complexity>>,
    /// The alignment type
    pub ty: AlignType,
}

/// The placement of a single peptide in a multiple sequence alignment
#[derive(Debug, Clone, PartialEq)]
pub struct MSAPlacement<Complexity> {
    /// The sequence
    pub sequence: Peptidoform<Complexity>,
    /// The start position in the sequence where the alignment starts
    pub start: usize,
    /// The path the alignment follows
    pub path: Vec<MSAPosition>,
    /// The absolute score
    pub score: usize,
    /// The normalised score
    pub normalised_score: f64,
}

impl<Complexity> MSAPlacement<Complexity> {
    pub fn len(&self) -> usize {
        self.path.iter().map(|p| p.len()).sum()
    }
}

/// A single position in a multiple sequence alignment
#[derive(Debug, Clone, Hash, PartialEq, Eq)]
pub enum MSAPosition {
    /// A gap one sequence position spanning one sequence position
    Gap,
    /// A selection of sequence positions with its width, the maximal number of other sequence positions it spans
    Placed(MatchType, usize, usize),
}

impl<Complexity> MultipleSequenceAlignment<Complexity> {
    /// Normalise the alignment by making sure all steps are expanded (no steps that squish any piece)
    /// and that any location that is expanded for every sequence is squished back.
    pub(super) fn normalise(&mut self) {
        fn fix<C>(
            msa: &mut MultipleSequenceAlignment<C>,
            search_index: usize,
            expand: usize,
            except: usize,
        ) {
            for sequence_index in (0..msa.sequences.len()).filter(|i| *i != except) {
                let mut index = msa.sequences[sequence_index].start;
                for path_index in 0..msa.sequences[sequence_index].path.len() {
                    if let MSAPosition::Placed(_, _, b) =
                        &mut msa.sequences[sequence_index].path[path_index]
                    {
                        if index + *b >= search_index {
                            *b += expand;
                            break;
                        }
                        index += *b;
                    } else {
                        if index + 1 >= search_index {
                            for _ in 0..expand {
                                msa.sequences[sequence_index]
                                    .path
                                    .insert(path_index, MSAPosition::Gap);
                            }
                            break;
                        }
                        index += 1;
                    }
                }
            }
        }

        for sequence_index in 0..self.sequences.len() {
            let mut index = self.sequences[sequence_index].start;
            for path_index in 0..self.sequences[sequence_index].path.len() {
                if let MSAPosition::Placed(_, a, b) =
                    self.sequences[sequence_index].path[path_index]
                {
                    if a > b {
                        fix(self, index, a - b, sequence_index);
                    }
                    index += b;
                } else {
                    index += 1;
                }
            }
        }
    }

    /// Generate HTML
    pub fn html(&self) -> String {
        let mut res = String::new();
        write!(
            &mut res,
            "<div class='msa alignment' style='--length:{}'>",
            self.total_length()
        )
        .unwrap();
        for sequence in &self.sequences {
            write!(&mut res, "<div class='peptide' title='path: {} norm: {:.4} abs: {}' data-normalised-score='{1}' style='--start:{};--length:{};--score:{1}'>", 
        sequence
        .path
        .iter()
        .map(ToString::to_string)
        .collect::<String>(),
        sequence.normalised_score,
        sequence.score,
        sequence.start,
        sequence.len(),
    ).unwrap();
            let mut start = sequence.start;
            for piece in &sequence.path {
                match piece {
                    MSAPosition::Gap => write!(&mut res, "<del></del>").unwrap(),
                    MSAPosition::Placed(ty, a, b) => {
                        write!(
                            &mut res,
                            "<span class='{}' style='--wa:{};--wb:{};'><span>{}</span></span>",
                            match ty {
                                MatchType::Rotation => "rotation",
                                MatchType::Mismatch => "mismatch",
                                MatchType::Isobaric => "isobaric",
                                MatchType::Gap => "gap",
                                MatchType::FullIdentity => "identity",
                                MatchType::IdentityMassMismatch => "massmismatch",
                            },
                            a,
                            b,
                            render_seq(
                                &sequence.sequence,
                                start..(start + b).min(sequence.sequence.len()),
                                *ty
                            )
                        )
                        .unwrap();
                        start += a;
                    }
                }
            }
            write!(&mut res, "</div>").unwrap();
        }
        write!(&mut res, "</div>").unwrap();
        res
    }

    // If this multiple sequence alignment consists of only a single pair, create a pairwise mass alignment struct.
    // fn assume_single(&self) -> Option<Alignment> {
    //     todo!();
    // }

    /// Determine the scores for all alignments (should be called after the full alignment is done)
    fn determine_scores(&mut self) {
        // For each location find the highest score it has for all other sets
        todo!();
    }
}

fn render_seq<Complexity>(
    seq: &Peptidoform<Complexity>,
    sel: std::ops::Range<usize>,
    ty: MatchType,
) -> String {
    use std::fmt::Write;
    let mut res = String::new();
    for s in &seq.sequence()[sel] {
        if s.modifications.is_empty() {
            write!(&mut res, "{}", s.aminoacid.pro_forma_definition()).unwrap();
        } else {
            write!(
                &mut res,
                "<span class='modified{}'>{}</span>",
                match ty {
                    MatchType::IdentityMassMismatch => " massmismatch",
                    MatchType::Mismatch => " mismatch",
                    _ => "",
                },
                s.aminoacid.pro_forma_definition()
            )
            .unwrap();
        }
    }
    res
}

impl<Complexity> Display for MultipleSequenceAlignment<Complexity> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for sequence in &self.sequences {
            write!(f, "{}", " ".repeat(sequence.start))?;
            let mut seq = sequence.sequence.sequence().iter().skip(sequence.start);
            for step in &sequence.path {
                match step {
                    MSAPosition::Gap => write!(f, "-")?,
                    MSAPosition::Placed(_, steps, width) => write!(
                        f,
                        "{}{}",
                        (&mut seq)
                            .take(*steps)
                            .map(|s| s.aminoacid.pro_forma_definition())
                            .collect::<String>(),
                        "·".repeat(width.saturating_sub(*steps)) // TODO: handle the cases where steps is too big
                    )?,
                }
            }
            writeln!(
                f,
                " [Score: {}, Normalised score: {:.3}, Path: {}]",
                sequence.score,
                sequence.normalised_score,
                sequence
                    .path
                    .iter()
                    .map(ToString::to_string)
                    .collect::<String>()
            )?;
        }
        Ok(())
    }
}

impl MSAPosition {
    pub const fn len(&self) -> usize {
        match self {
            Self::Gap => 1,
            Self::Placed(_, a, _) => *a,
        }
    }
    pub const fn ty(&self) -> MatchType {
        match self {
            Self::Gap => MatchType::Gap,
            Self::Placed(ty, _, _) => *ty,
        }
    }
}

impl Display for MSAPosition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::Gap => "[-]".to_string(),
                Self::Placed(ty, a, b) => format!(
                    "[{}{a},{b}]",
                    match ty {
                        MatchType::FullIdentity => '=',
                        MatchType::IdentityMassMismatch => '≅',
                        MatchType::Isobaric => 'i',
                        MatchType::Rotation => 'r',
                        MatchType::Mismatch => 'x',
                        MatchType::Gap => 'g',
                    },
                ),
            }
        )
    }
}

impl<Complexity: AtMax<Linear>> MassAlignable<Complexity>
    for MultipleSequenceAlignment<Complexity>
{
    fn total_length(&self) -> usize {
        self.sequences
            .iter()
            .map(|s| s.len() + s.start)
            .max()
            .unwrap_or(0)
    }
    fn number_of_sequences(&self) -> usize {
        self.sequences.len()
    }
    fn index(&self, index: usize, sequence_index: usize) -> &SequenceElement<Complexity> {
        &self.sequences[sequence_index].sequence.sequence()[index]
    }
    // TODO: has to be Option output? Also has to be aligned with the path?
    fn index_slice(
        &self,
        index: impl std::ops::RangeBounds<usize>,
        sequence_index: usize,
    ) -> std::borrow::Cow<[SequenceElement<Complexity>]> {
        std::borrow::Cow::Borrowed(
            &self.sequences[sequence_index].sequence.sequence()
                [(index.start_bound().cloned(), index.end_bound().cloned())],
        )
    }
    fn sequence_bounds(&self) -> Vec<(usize, usize)> {
        self.sequences
            .iter()
            .map(|s| (s.start, s.len() + s.start))
            .collect()
    }
    ///Calculate all masses for all steps beforehand. First dimension is the index in the multiple sequences
    /// (or just a single sequence if available). Second dimension is the index into the sequence (offset by
    /// the start and the length is the length of this particular sequence, does not have to span to the very
    /// end). Third dimension is the masses for that number of steps - 1 (a step of 1 is at index 0). Finally
    /// the result is None when that selection given you an invalid step (like getting the mass halfway in an
    /// aligned aminoacid)
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<[Option<Mass>; STEPS]>> {
        self.sequences
            .iter()
            .enumerate()
            .map(|(sequence_index, sequence)| {
                (0..sequence.len())
                    .map(|index| {
                        std::array::from_fn(|i| {
                            let size = i + 1;
                            if index < size {
                                None
                            } else {
                                Some(
                                    self.index_slice(index - size..index, sequence_index)
                                        .iter()
                                        .map(|s| s.formulas())
                                        .sum::<MolecularFormula>()
                                        .monoisotopic_mass()
                                        .unwrap(),
                                )
                            }
                        })
                    })
                    .collect()
            })
            .collect()
    }
    fn sequences_with_path(
        &self,
        is_a: bool,
        start: usize,
        path: &[super::Piece],
    ) -> Vec<MSAPlacement<Complexity>> {
        todo!();
    }
}

/// The shared behaviour needed to be mass alignable into a MSA
pub trait MassAlignable<Complexity> {
    /// Score for gap extend
    const GAP_EXTEND: i8 = -1;
    /// Score for gap start
    const GAP_START: i8 = -5;
    /// Additional score when there is a mass mismatch but AA identity
    const MASS_MISMATCH: i8 = -1;
    /// Isomass score (per AA)
    const ISOMASS: i8 = 2;
    /// Additional score on top of the matrix score for mismatch
    const MISMATCH: i8 = -1;
    /// Set rotation (swap) score (per AA)
    const SWITCHED: i8 = 3;
    /// Additional score on top of the matrix score for identity
    const IDENTITY: i8 = 0;

    /// Index the underlying structure on one location
    fn index(&self, index: usize, sequence_index: usize) -> &SequenceElement<Complexity>;
    /// Get a slice of the underlying structure
    fn index_slice(
        &self,
        index: impl RangeBounds<usize>,
        sequence_index: usize,
    ) -> Cow<[SequenceElement<Complexity>]>;

    /// Total length of the structure (number of AA for single, length from start to end for multiple)
    fn total_length(&self) -> usize;
    /// Give the total number of sequences stored in this structure
    fn number_of_sequences(&self) -> usize;
    /// For all sequences in this structure give their bounds (start index, end index)
    fn sequence_bounds(&self) -> Vec<(usize, usize)>;
    /// Calculate all masses for all steps beforehand.
    /// First dimension is the index in the multiple sequences (or just a single sequence if available).
    /// Second dimension is the index into the sequence (offset by the start and the length is the length of this particular sequence, does not have to span to the very end).
    /// Third dimension is the masses for that number of steps - 1 (a step of 1 is at index 0).
    /// Finally the result is None when that selection given you an invalid step (like getting the mass halfway in an aligned aminoacid)
    fn calculate_masses<const STEPS: usize>(&self) -> Vec<Vec<[Option<Mass>; STEPS]>>;

    /// Get the [`MSAPlacement`] of the underlying sequences given the path resulting from the alignment
    fn sequences_with_path(
        &self,
        is_a: bool,
        start: usize,
        path: &[Piece],
    ) -> Vec<MSAPlacement<Complexity>>;

    /// Align two sequences
    fn align<const STEPS: usize, A: MassAlignable<Complexity>, B: MassAlignable<Complexity>>(
        seq_a: &A,
        seq_b: &B,
        alphabet: &[&[i8]],
        tolerance: Tolerance,
        ty: AlignType,
    ) -> MultipleSequenceAlignment<Complexity> {
        assert!(isize::try_from(seq_a.total_length()).is_ok());
        assert!(isize::try_from(seq_b.total_length()).is_ok());
        let mut matrix =
            vec![vec![Piece::default(); seq_b.total_length() + 1]; seq_a.total_length() + 1];
        let mut high = (0, 0, 0);
        let masses_a = seq_a.calculate_masses::<STEPS>();
        let masses_b = seq_b.calculate_masses::<STEPS>();
        let bounds_a = seq_a.sequence_bounds();
        let bounds_b = seq_b.sequence_bounds();

        if ty.bind_start_a {
            #[allow(clippy::cast_possible_wrap)]
            // b is always less than seq_b
            for index_b in 0..=seq_b.total_length() {
                matrix[0][index_b] = Piece::new(
                    (index_b as isize).saturating_sub(1) * Self::GAP_EXTEND as isize
                        + Self::GAP_START as isize,
                    if index_b == 0 {
                        Self::GAP_START
                    } else {
                        Self::GAP_EXTEND
                    },
                    MatchType::Gap,
                    0,
                    u8::from(index_b != 0),
                );
            }
        }
        if ty.bind_start_b {
            #[allow(clippy::cast_possible_wrap)]
            // a is always less than seq_a
            for (index_a, row) in matrix.iter_mut().enumerate() {
                row[0] = Piece::new(
                    (index_a as isize).saturating_sub(1) * Self::GAP_EXTEND as isize
                        + Self::GAP_START as isize,
                    if index_a == 0 {
                        Self::GAP_START
                    } else {
                        Self::GAP_EXTEND
                    },
                    MatchType::Gap,
                    u8::from(index_a != 0),
                    0,
                );
            }
        }

        // Main loop
        let mut values = Vec::with_capacity(
            STEPS * STEPS * seq_a.number_of_sequences() * seq_b.number_of_sequences() + 2,
        );
        for index_a in 0..seq_a.total_length() {
            for index_b in 0..seq_b.total_length() {
                values.clear();
                for sequence_index_a in 0..seq_a.number_of_sequences() {
                    // If this sequence does not have info on this location skip it
                    if index_a < bounds_a[sequence_index_a].0
                        || index_a > bounds_a[sequence_index_a].1
                    {
                        continue;
                    }
                    for sequence_index_b in 0..seq_b.number_of_sequences() {
                        // If this sequence does not have info on this location skip it
                        if index_b < bounds_b[sequence_index_b].0
                            || index_b > bounds_b[sequence_index_b].1
                        {
                            continue;
                        }
                        for len_a in 0..=STEPS {
                            for len_b in 0..=STEPS {
                                if len_a == 0 && len_b != 1
                                    || len_a != 1 && len_b == 0
                                    || len_a > index_a + 1
                                    || len_b > index_b + 1
                                {
                                    continue; // Do not allow double gaps, any double gaps will be counted as two gaps after each other
                                }
                                let base_score =
                                    matrix[index_a + 1 - len_a][index_b + 1 - len_b].score;
                                let piece = if len_a == 0 || len_b == 0 {
                                    // First check the score to be used for affine gaps
                                    let prev = &matrix[index_a + 1 - len_a][index_b + 1 - len_b];
                                    let score = if prev.step_a == 0 && len_a == 0
                                        || prev.step_b == 0 && len_b == 0
                                    {
                                        Self::GAP_EXTEND
                                    } else {
                                        Self::GAP_START
                                    };
                                    Some(Piece::new(
                                        base_score + score as isize,
                                        score,
                                        MatchType::Gap,
                                        len_a as u8,
                                        len_b as u8,
                                    ))
                                } else if len_a == 1 && len_b == 1 {
                                    masses_a[sequence_index_a][index_a][0].and_then(|mass_a| {
                                        masses_b[sequence_index_b][index_b][0].map(|mass_b| {
                                            Self::score_pair(
                                                seq_a.index(index_a, sequence_index_a),
                                                mass_a,
                                                seq_b.index(index_b, sequence_index_b),
                                                mass_b,
                                                alphabet,
                                                base_score,
                                                tolerance,
                                            )
                                        })
                                    })
                                } else {
                                    masses_a[sequence_index_a][index_a][len_a - 1].and_then(
                                        |mass_a| {
                                            masses_b[sequence_index_b][index_b][len_b - 1].and_then(
                                                |mass_b| {
                                                    Self::score(
                                                        &seq_a.index_slice(
                                                            index_a + 1 - len_a..=index_a,
                                                            sequence_index_a,
                                                        ),
                                                        mass_a,
                                                        &seq_b.index_slice(
                                                            index_b + 1 - len_b..=index_b,
                                                            sequence_index_b,
                                                        ),
                                                        mass_b,
                                                        base_score,
                                                        tolerance,
                                                    )
                                                },
                                            )
                                        },
                                    )
                                };
                                if let Some(p) = piece {
                                    values.push(p);
                                }
                            }
                        }
                    }
                }
                // Determine the best step
                let value = values
                    .iter()
                    .max_by(|x, y| x.score.cmp(&y.score))
                    .cloned()
                    .expect("No possible steps recorded");
                // Keep track of the highest scoring cell
                if value.score >= high.0 {
                    high = (value.score, index_a + 1, index_b + 1);
                }
                // If local and the score is too low just do not store the result
                if ty.is_global() || value.score > 0 {
                    matrix[index_a + 1][index_b + 1] = value;
                }
            }
        }

        // Find the end cell
        let mut target = high;
        if ty.bind_end_a && ty.bind_end_b {
            target = (
                matrix[seq_a.total_length()][seq_b.total_length()].score,
                seq_a.total_length(),
                seq_b.total_length(),
            );
        } else if ty.bind_end_a {
            let value = (0..=seq_a.total_length())
                .map(|v| (v, matrix[v][seq_b.total_length()].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, value.0, seq_b.total_length());
        } else if ty.bind_end_b {
            let value = (0..=seq_b.total_length())
                .map(|v| (v, matrix[seq_a.total_length()][v].score))
                .max_by(|a, b| a.1.cmp(&b.1))
                .unwrap_or_default();
            target = (value.1, seq_a.total_length(), value.0);
        }

        // Walk the path back
        let mut path = Vec::new();
        while !(target.1 == 0 && target.2 == 0) {
            let value = matrix[target.1][target.2].clone();
            if value.step_a == 0 && value.step_b == 0 {
                break;
            }
            target = (
                0,
                target.1 - value.step_a as usize,
                target.2 - value.step_b as usize,
            );
            path.push(value);
        }

        let path: Vec<Piece> = path.into_iter().rev().collect();
        let mut sequences = seq_a.sequences_with_path(true, target.1, &path);
        sequences.extend(seq_b.sequences_with_path(false, target.2, &path));
        let mut msa = MultipleSequenceAlignment { sequences, ty };
        msa.normalise();
        msa
    }

    /// Determine the score for an alignment of two AAs (with modifications)
    fn score_pair(
        a: &SequenceElement,
        mass_a: Mass,
        b: &SequenceElement,
        mass_b: Mass,
        alphabet: &[&[i8]],
        score: isize,
        tolerance: MassTolerance,
    ) -> Piece {
        match (a == b, tolerance.within(mass_a, mass_b)) {
            (true, true) => {
                let local = Self::IDENTITY + alphabet[a.aminoacid as usize][b.aminoacid as usize];
                Piece::new(score + local as isize, local, MatchType::FullIdentity, 1, 1)
            }
            (true, false) => {
                let local =
                    alphabet[a.aminoacid as usize][b.aminoacid as usize] + Self::MASS_MISMATCH;
                Piece::new(
                    score + local as isize,
                    local,
                    MatchType::IdentityMassMismatch,
                    1,
                    1,
                )
            }
            (false, true) => Piece::new(
                score + Self::ISOMASS as isize,
                Self::ISOMASS,
                MatchType::Isobaric,
                1,
                1,
            ), // TODO: I/L/J is now also scored as isobaric, which is correct but the score is considerably lower then in previous iterations
            (false, false) => Piece::new(
                score + Self::MISMATCH as isize,
                Self::MISMATCH,
                MatchType::Mismatch,
                1,
                1,
            ),
        }
    }

    /// Score two sets of aminoacids (it will only be called when at least one of a and b has len > 1)
    /// Returns none if no sensible explanation can be made
    fn score(
        a: &[SequenceElement],
        mass_a: Mass,
        b: &[SequenceElement],
        mass_b: Mass,
        score: isize,
        tolerance: MassTolerance,
    ) -> Option<Piece> {
        if tolerance.within(mass_a, mass_b) {
            let mut b_copy = b.to_owned();
            let switched = a.len() == b.len()
                && a.iter().all(|el| {
                    b_copy.iter().position(|x| x == el).map_or(false, |pos| {
                        b_copy.remove(pos);
                        true
                    })
                });
            #[allow(clippy::cast_possible_wrap)]
            let local = if switched {
                Self::SWITCHED * a.len() as i8
            } else {
                Self::ISOMASS * (a.len() + b.len()) as i8 / 2
            };
            Some(Piece::new(
                score + local as isize,
                local,
                if switched {
                    MatchType::Switched
                } else {
                    MatchType::Isobaric
                },
                a.len() as u8,
                b.len() as u8,
            ))
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn print_msma() {
        let msma = MultipleSequenceAlignment {
            sequences: vec![
                MSAPlacement {
                    sequence: ComplexPeptide::pro_forma("WGGD").unwrap().assume_linear(),
                    start: 0,
                    path: vec![
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 1),
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                    ],
                    score: 0,
                    normalised_score: 0.0,
                },
                MSAPlacement {
                    sequence: ComplexPeptide::pro_forma("WND").unwrap().assume_linear(),
                    start: 0,
                    path: vec![
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                        MSAPosition::Placed(MatchType::Isobaric, 1, 2),
                        MSAPosition::Placed(MatchType::FullIdentity, 1, 1),
                    ],
                    score: 0,
                    normalised_score: 0.0,
                },
            ],
            ty: AlignmentType {
                bind_start_a: true,
                bind_start_b: true,
                bind_end_a: true,
                bind_end_b: true,
            },
        };
        println!("{msma}");
        assert_eq!(
            msma.to_string(),
            "WGGD [Score: 0, Normalised score: 0.000]\nWN·D [Score: 0, Normalised score: 0.000]\n"
        );
    }
}
