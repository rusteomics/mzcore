#![allow(dead_code)]

use serde::{Deserialize, Serialize};

use crate::{
    AlignIndex, AlignScoring, AlignType, MatchType, Score,
    diagonal_array::DiagonalArray,
    mass_alignment::{align_cached, calculate_masses},
};
use mzcore::{
    prelude::*,
    quantities::Multi,
    sequence::{HasPeptidoform, SimpleLinear},
    system::Mass,
};

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiAlignment<Sequence> {
    lines: Vec<MultiAlignmentLine<Sequence>>,
    score: Score,
    maximal_step: u16,
    align_type: AlignType,
}

impl<Sequence: HasPeptidoform<SimpleLinear>> MultiAlignment<Sequence> {
    fn debug_display(&self, mut w: impl std::fmt::Write) {
        writeln!(
            w,
            "Multi score: {} ({}/{}) max_step: {}",
            self.score.normalised, self.score.absolute, self.score.max, self.maximal_step
        )
        .unwrap();
        for line in &self.lines {
            line.debug_display(&mut w);
        }
    }

    fn new<const STEPS: u16>(
        sequences: Vec<MultiAlignmentLineTemp<'_, Sequence, STEPS>>,
        align_type: AlignType,
    ) -> Self {
        Self {
            lines: sequences
                .into_iter()
                .map(|temp| MultiAlignmentLine {
                    original_index: temp.original_index,
                    sequence: temp.sequence,
                    path: temp.path,
                    start: temp.start,
                })
                .collect(),
            score: Score {
                // TODO: build a scoring function that gets a score from a MultiAlignment
                normalised: 0.0.into(),
                absolute: 0,
                max: 0,
            },
            align_type,
            maximal_step: STEPS,
        }
    }
}

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiAlignmentLine<Sequence> {
    original_index: usize,
    sequence: Sequence,
    path: Vec<MultiPiece>,
    start: usize,
}

#[derive(Debug)]
struct MultiAlignmentLineTemp<'a, Sequence, const STEPS: u16> {
    original_index: usize,
    sequence: Sequence,
    masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    path: Vec<MultiPiece>,
    start: usize,
}

impl<'a, Sequence: HasPeptidoform<SimpleLinear>, const STEPS: u16>
    MultiAlignmentLineTemp<'a, Sequence, STEPS>
{
    fn single(
        original_index: usize,
        sequence: Sequence,
        masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    ) -> Self {
        let len = sequence.cast_peptidoform().len();
        Self {
            original_index,
            sequence,
            masses,
            path: std::iter::repeat_n(
                MultiPiece {
                    match_type: MatchType::FullIdentity,
                    ref_step: 1,
                    seq_step: 1,
                },
                len,
            )
            .collect(),
            start: 0,
        }
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiPiece {
    match_type: MatchType,
    /// Ref seq is required to always be at least `seq_step`
    ref_step: u16,
    seq_step: u16,
}

impl<Sequence: HasPeptidoform<SimpleLinear>> MultiAlignmentLine<Sequence> {
    fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let mut index = self.start;
        for piece in &self.path {
            let subseq = &sequence[index..index + piece.seq_step as usize];
            // Obviously misses mods now
            let display = subseq
                .iter()
                .map(|s| s.aminoacid.one_letter_code().unwrap_or('X'))
                .collect::<String>();
            write!(
                w,
                "{display}{}",
                "·".repeat((piece.ref_step - piece.seq_step) as usize)
            )
            .unwrap();
            index += piece.seq_step as usize;
        }
        writeln!(w).unwrap();
    }
}

impl<const STEPS: u16, Sequence: HasPeptidoform<SimpleLinear> + Clone> AlignIndex<STEPS, Sequence> {
    fn multi_align(
        &self,
        maximal_distance: Option<f64>, // Maximal distance to join clusters, or join all if set to None
        scoring: AlignScoring<'_>,
        align_type: AlignType, // TODO: figure out what to do for align_type
    ) -> Vec<MultiAlignment<Sequence>> {
        assert!(self.sequences.len() > 2);
        // Create an outgroup of 'random' sequence with the length of the average of all used sequences
        let outgroup: Peptidoform<SimpleLinear> = CheckedAminoAcid::CANONICAL_AMINO_ACIDS
            .iter()
            .cycle()
            .take({
                let (n, l) = self.sequences.iter().fold((0, 0), |(n, l), s| {
                    (n + 1, l + s.0.cast_peptidoform().len())
                });
                l / n
            })
            .map(|a| SequenceElement::new(*a, None))
            .collect();
        let outgroup_masses = calculate_masses::<STEPS>(&outgroup, self.mode);

        // Make distance matrix
        let mut matrix = DiagonalArray::<f64, { u16::MAX }>::new(self.sequences.len() + 1);
        let mut nodes = Vec::with_capacity(self.sequences.len() + 1);
        for (i, s) in self.sequences.iter().enumerate() {
            matrix[[i, i]] = 0.0;
            for (o, so) in self.sequences.iter().enumerate().take(i) {
                matrix[[i, o]] = align_cached(
                    s.0.cast_peptidoform(),
                    &s.1,
                    so.0.cast_peptidoform(),
                    &so.1,
                    scoring,
                    align_type,
                )
                .distance();
            }
            matrix[[self.sequences.len(), i]] = align_cached(
                s.0.cast_peptidoform(),
                &s.1,
                &outgroup,
                &outgroup_masses,
                scoring,
                align_type,
            )
            .distance();
            nodes.push(vec![MultiAlignmentLineTemp::single(i, s.0.clone(), &s.1)]);
        }
        matrix[[self.sequences.len(), self.sequences.len()]] = 0.0;
        // Remember that the last node is not made (the outgroup node) so do not try to merge that one

        println!("{}", matrix.to_csv());
        // Create tree
        // It might work out to now rely on the fact that the outgroup is always the furthest distance from all other nodes that we can just immediately start merging the nodes to produce the multi alignment

        // Progressively merge
        while nodes.len() > 2 {
            let (distance, [far_index, close_index]) = matrix.min::<true>();
            if maximal_distance.is_some_and(|m| distance > m)
                || far_index == nodes.len()
                || close_index == nodes.len()
            {
                // Stop if this distance is more than the threshold, or if the outgroup would be merged
                break;
            }
            // Merge nodes (this is where the actual alignment happens)
            let close_node = std::mem::take(&mut nodes[close_index]);
            let far_node = nodes.remove(far_index);
            nodes[close_index] = merge(close_node, far_node, scoring, align_type);

            // Update matric to merge these columns
            matrix = matrix.merge_columns(close_index, far_index, distance);
            println!(
                "Merged {close_index} with {far_index} left nodes {} left matrix {}",
                nodes.len(),
                matrix.len(),
            );
        }

        nodes
            .into_iter()
            .map(|seqs| MultiAlignment::new(seqs, align_type))
            .collect()
    }
}

/// Do the actual alignment
fn merge<'a, Sequence, const STEPS: u16>(
    mut a: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    mut b: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
) -> Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>> {
    a.append(&mut b); // TODO
    a
}

enum MultiAlignTree<'a, Sequence, const STEPS: u16> {
    Branch {
        length: f64,
        left: Box<MultiAlignTree<'a, Sequence, STEPS>>,
        right: Box<MultiAlignTree<'a, Sequence, STEPS>>,
    },
    Leaf {
        length: f64,
        sequences: Vec<MultiAlignmentLineTemp<'a, Sequence, STEPS>>,
    },
    // To easily get nodes out of the list for processing
    Empty,
}

#[cfg(test)]
mod tests {
    use super::*;
    use itertools::Itertools;
    use mzcore::ontology::STATIC_ONTOLOGIES;

    //     #[test]
    //     fn simple() {
    //         // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
    //         let seq = |def: &str| {
    //             Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
    //                 .unwrap()
    //                 .0
    //                 .into_simple_linear()
    //                 .unwrap()
    //         };
    //         let sequences = vec![seq("AGGWHD"), seq("ANWHN[Deamidated]"), seq("AHYDH")];
    //         let alignment = multi_align::<4, Peptidoform<SimpleLinear>>(
    //             sequences,
    //             AlignScoring::default(),
    //             AlignType::GLOBAL,
    //         );
    //         let mut buf = String::new();
    //         alignment.debug_display(&mut buf);
    //         buf = buf.split('\n').skip(1).join("\n");
    //         assert_eq!(
    //             buf,
    //             "AGGWHD
    // AN·WHN
    // AHY·DH
    // "
    //         );
    //     }

    #[test]
    fn many() {
        // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
        let seq = |def: &str| {
            Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_simple_linear()
                .unwrap()
        };
        let sequences = vec![
            seq("WRGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("WRNDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWNDGFYW[U:Oxidation]DYWGQG"),
            seq("RWGGN[U:Deamidated]GFYW[U:Oxidation]DYWGQG"),
        ];
        let index =
            AlignIndex::<4, Peptidoform<SimpleLinear>>::new(sequences, MassMode::Monoisotopic);
        let alignment = index.multi_align(None, AlignScoring::default(), AlignType::GLOBAL);
        let mut buf = String::new();
        alignment[0].debug_display(&mut buf);
        buf = buf.split('\n').skip(1).join("\n");
        assert_eq!(
            buf,
            "AGGWHD
AN·WHN
AHY·DH
"
        );
    }
}
