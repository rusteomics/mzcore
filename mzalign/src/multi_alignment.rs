#![allow(dead_code)]
use std::borrow::Cow;

use serde::{Deserialize, Serialize};

use crate::{
    AlignIndex, AlignScoring, AlignType, MatchType, Score,
    diagonal_array::DiagonalArray,
    mass_alignment::{align_cached, calculate_masses},
};
use mzcore::{
    prelude::*,
    quantities::Multi,
    sequence::{HasPeptidoform, Linear, SimpleLinear, UnAmbiguous},
    system::Mass,
};

#[derive(Clone, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiAlignment<Sequence> {
    lines: Vec<MultiAlignmentLine<Sequence>>,
    score: Score,
    maximal_step: u16,
    align_type: AlignType,
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignment<Sequence> {
    fn debug_display(&self, mut w: impl std::fmt::Write) {
        writeln!(
            w,
            "Multi score: {} ({}/{}) max_step: {}",
            self.score.normalised, self.score.absolute, self.score.max, self.maximal_step
        );
        for line in &self.lines {
            line.debug_display(&mut w);
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
struct MultiAlignmentLineTemp<'a, const STEPS: u16> {
    original_index: usize,
    sequence: &'a Peptidoform<SimpleLinear>,
    masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    path: Vec<MultiPiece>,
    start: usize,
}

impl<'a, const STEPS: u16> MultiAlignmentLineTemp<'a, STEPS> {
    fn single(
        original_index: usize,
        sequence: &'a Peptidoform<SimpleLinear>,
        masses: &'a DiagonalArray<Multi<Mass>, STEPS>,
    ) -> Self {
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
                sequence.len(),
            )
            .collect(),
            start: 0,
        }
    }
}

#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
struct MultiPiece {
    match_type: MatchType,
    /// Ref seq is required to always be at least seq_step
    ref_step: u16,
    seq_step: u16,
}

impl<Sequence: HasPeptidoform<Linear>> MultiAlignmentLine<Sequence> {
    fn debug_display(&self, mut w: impl std::fmt::Write) {
        let sequence = self.sequence.cast_peptidoform().sequence();
        let index = self.start;
        for piece in self.path.iter() {
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
            );
        }
        writeln!(w);
    }
}

impl<const STEPS: u16, Sequence: HasPeptidoform<SimpleLinear>> AlignIndex<STEPS, Sequence> {
    fn multi_align(
        &self,
        scoring: AlignScoring<'_>,
        align_type: AlignType, // TODO: figure out what to do for align_type
    ) -> MultiAlignment<Sequence> {
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
            matrix[[i, i]] = 1.0;
            for (o, so) in self.sequences.iter().enumerate().skip(i + 1) {
                matrix[[o, i]] = align_cached(
                    s.0.cast_peptidoform(),
                    &s.1,
                    so.0.cast_peptidoform(),
                    &so.1,
                    scoring,
                    align_type,
                )
                .normalised_score();
            }
            matrix[[self.sequences.len(), i]] = align_cached(
                s.0.cast_peptidoform(),
                &s.1,
                &outgroup,
                &outgroup_masses,
                scoring,
                align_type,
            )
            .normalised_score();
            nodes.push(MultiAlignTree::Leaf {
                length: 0.0,
                sequences: vec![MultiAlignmentLineTemp::single(
                    i,
                    s.0.cast_peptidoform(),
                    &s.1,
                )],
            });
        }
        matrix[[self.sequences.len(), self.sequences.len()]] = 1.0;
        // Maybe this trickery with the temp alignment is not needed if this outgroup is just never given a node, but then it needs to be special cased in the rest of the code to not use it anywhere
        nodes.push(MultiAlignTree::Leaf {
            length: 0.0,
            sequences: vec![MultiAlignmentLineTemp::single(
                self.sequences.len(),
                &outgroup,
                &outgroup_masses,
            )],
        });

        println!("{}", matrix.to_csv());
        // Create tree
        // It might work out to now rely on the fact that the outgroup is always the furthest distance from all other nodes that we can just immediately start merging the nodes to produce the multi alignment
        // Root the tree

        // Progressively merge

        MultiAlignment {
            lines: Vec::new(),
            score: Score {
                normalised: 0.0.into(),
                absolute: 0,
                max: 0,
            },
            maximal_step: STEPS,
            align_type,
        }
    }
}

enum MultiAlignTree<'a, const STEPS: u16> {
    Branch {
        length: f64,
        left: Box<MultiAlignTree<'a, STEPS>>,
        right: Box<MultiAlignTree<'a, STEPS>>,
    },
    Leaf {
        length: f64,
        sequences: Vec<MultiAlignmentLineTemp<'a, STEPS>>,
    },
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
        let alignment = index.multi_align(AlignScoring::default(), AlignType::GLOBAL);
        let mut buf = String::new();
        alignment.debug_display(&mut buf);
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
