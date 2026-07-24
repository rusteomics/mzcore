#![allow(clippy::missing_panics_doc)]

use std::println;

use itertools::Itertools;
use mzcore::{ontology::STATIC_ONTOLOGIES, prelude::*, quantities::Multi, sequence::Linear};

use super::*;
use crate::{
    AlignIndex, AlignScoring, MatchType, Piece, mass_alignment::calculate_masses,
    multi_alignment::calculate::MultiAlignmentLineTemp,
};

fn seq(def: &str) -> Peptidoform<Linear> {
    Peptidoform::pro_forma(def, &STATIC_ONTOLOGIES)
        .unwrap()
        .0
        .into_linear()
        .unwrap()
}

#[test]
fn simple() {
    // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
    let sequences = vec![seq("AGGWHD"), seq("ANWHN[Deamidated]"), seq("AHYDH")];
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
    let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::GLOBAL);
    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "AGGWHD\nAN·WHN\nAH·YDH\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn gaps() {
    // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
    let sequences = vec![
        seq("WRGGDGFYAMDYWGQG"),
        seq("RWGDGFYAMYDWGQG"),
        seq("RWNDGYAMEDYWGQG"),
        seq("RWGGDGFMDYWGQG"),
    ];
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
    let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::GLOBAL);
    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "WRGGDGFYAM-DYWGQG\nRWGGDGF--M-DYWGQG\nRWG-DGFYAM-YDWGQG\nRWN·DG-YAMEDYWGQG\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn either_global() {
    // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
    let sequences = vec![
        seq("WRGGDGFYAMDYWGQG"),
        seq("RWGDGFYAMYD"),
        seq("RWNDGYAMDYWGQG"),
        seq("GDGFMDYWGQG"),
    ];
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
    let alignment = index.multi_align(None, AlignScoring::default(), MultiAlignType::EITHER_GLOBAL);
    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "WRGGDGFYAMDYWGQG\n---GDGF--MDYWGQG\nRWN·DG-YAMDYWGQG\nRWG-DGFYAMYD----\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn either_global_2() {
    let sequences = vec![
        seq("CSRWRGGDGF"),
        seq("WRNDGFYAM"),
        seq("RWGGDGFYAMDYWG"),
        seq("YW[U:Oxidation]DYWGQG"),
        seq("RWGGNGFYW[U:Oxidation]DYWGQG"),
    ];
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
    let alignment = index.multi_align(
        None,
        AlignScoring {
            tolerance: mzcore::quantities::Tolerance::new_ppm(20.0),
            ..Default::default()
        },
        MultiAlignType::EITHER_GLOBAL,
    );
    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "CSRWRGGDGF---------\n---WRN·DGFYAM------\n---RWGGDGFYAMDYWG--\n----------YW·DYWGQG\n---RWGGNGFYW·DYWGQG\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn crash() {
    let sequences = vec![
        seq("HYTTPPTFGQGT"),
        seq("WGG"),
        seq("VTC[U:Carboxymethyl]QGLSSPKSL"),
    ];
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(sequences, MassMode::Monoisotopic);
    let scoring = AlignScoring::<'_> {
        tolerance: mzcore::quantities::Tolerance::Relative(
            mzcore::system::Ratio::new::<mzcore::system::ratio::ppm>(20.0).into(),
        ),
        ..Default::default()
    };
    let alignment = index.multi_align(None, scoring, MultiAlignType::EITHER_GLOBAL);
    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "----HYTTPPTFGQGT\n-----------WG-G-\nVTCQGLSSPKSL----\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn properties() {
    let sequence = seq("AGGWHD");
    let masses = calculate_masses::<4>(&sequence, MassMode::Monoisotopic);
    let mut line = MultiAlignmentLineTemp::single(0, &sequence, &masses);
    assert_eq!(line.get_sequence_index(4), (4, 0));
    assert_eq!(line.get_sequence_index(5), (5, 0));
    assert_eq!(
        line.get_aa_at(5),
        (
            &SequenceElement::<Linear>::new(AminoAcid::Histidine.into(), None),
            &Multi::from(CheckedAminoAcid::Histidine.formula().monoisotopic_mass())
        )
    );
    assert_eq!(
        line.get_block_at(5, 2),
        Some((
            [
                SequenceElement::<Linear>::new(AminoAcid::Tryptophan.into(), None),
                SequenceElement::<Linear>::new(AminoAcid::Histidine.into(), None),
            ]
            .as_slice(),
            &Multi::from(
                (CheckedAminoAcid::Tryptophan.formula() + CheckedAminoAcid::Histidine.formula())
                    .monoisotopic_mass()
            )
        ))
    );
    let updated = line.clone().update(
        &[
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // A
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
            Piece::new(0, 0, MatchType::Isobaric, 1, 2),     // W
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // H
            Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
        ],
        true,
        0,
    );
    line.path[3].aligned_length = 2; // The W is bigger
    assert_eq!(updated, line);
    assert_eq!(line.get_sequence_index(4), (4, 1));
    assert_eq!(line.get_sequence_index(5), (4, 0));
    assert_eq!(line.get_sequence_index(6), (5, 0));
    assert_eq!(
        line.get_aa_at(5),
        (
            &SequenceElement::<Linear>::new(AminoAcid::Histidine.into(), None),
            &Multi::from(CheckedAminoAcid::Histidine.formula().monoisotopic_mass())
        )
    );
    assert_eq!(
        line.get_block_at(5, 2),
        Some((
            [
                SequenceElement::<Linear>::new(AminoAcid::Tryptophan.into(), None),
                SequenceElement::<Linear>::new(AminoAcid::Histidine.into(), None),
            ]
            .as_slice(),
            &Multi::from(
                (CheckedAminoAcid::Tryptophan.formula() + CheckedAminoAcid::Histidine.formula())
                    .monoisotopic_mass()
            )
        ))
    );
}

#[test]
fn properties_gaps() {
    let sequence = seq("WRGGDGFYAMDYWGQG");
    let sequence2 = seq("RWNDGYAMEDYWGQG");
    let path = &[
        Piece::new(0, 0, MatchType::Rotation, 2, 2),     // WR
        Piece::new(0, 0, MatchType::Isobaric, 2, 1),     // GG N
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
        Piece::new(0, 0, MatchType::Gap, 1, 0),          // F
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Y
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // A
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // M
        Piece::new(0, 0, MatchType::Gap, 0, 1),          // E
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // D
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Y
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // W
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // Q
        Piece::new(0, 0, MatchType::FullIdentity, 1, 1), // G
    ];
    let masses = calculate_masses::<4>(&sequence, MassMode::Monoisotopic);
    let masses2 = calculate_masses::<4>(&sequence2, MassMode::Monoisotopic);
    let line = MultiAlignmentLineTemp::single(0, &sequence, &masses);
    let line2 = MultiAlignmentLineTemp::single(0, &sequence2, &masses2);
    let updated = line.clone().update(path, true, 0);
    let updated2 = line2.update(path, false, 0);
    let mut buf1 = String::new();
    updated.debug_display(&mut buf1);
    let mut buf2 = String::new();
    updated2.debug_display(&mut buf2);
    assert_eq!(buf1, "WRGGDGFYAM-DYWGQG\n");
    assert_eq!(buf2, "RWN·DG-YAMEDYWGQG\n");
}

#[ignore = "Known issue"]
#[test]
fn inherit_aligned_length() {
    // N = GG, N[Deamidated] = D, WGG = HY, HD = DH
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(
        [
            seq("WRGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAM[U:Oxidation]DYWGQG"),
            seq("WRNDGFYAM[U:Oxidation]DYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWGGDGFYAMDYWGQG"),
            seq("RWNDGFYW[U:Oxidation]DYWGQG"),
            seq("RWGGN[U:Deamidated]GFYW[U:Oxidation]DYWGQG"),
        ],
        MassMode::Monoisotopic,
    );
    let alignment = index.multi_align(
        None,
        AlignScoring {
            tolerance: mzcore::quantities::Tolerance::new_ppm(20.0),
            ..Default::default()
        },
        MultiAlignType::GLOBAL,
    );

    for line in &alignment[0] {
        let mut b = String::new();
        line.debug_display(&mut b);
        println!("{} {}", b.trim_ascii_end(), line.short());
    }

    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "WRGGDGFYAMDYWGQG\nWRN·DGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWGGDGFYAMDYWGQG\nRWN·DGFYW·DYWGQG\nRWGGNGFYW·DYWGQG\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn double_isobaric() {
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(
        [seq("RWGGDGFYAMDYWGQG"), seq("RWNDGFYW[U:Oxidation]DYWGQG")],
        MassMode::Monoisotopic,
    );
    let alignment = index.multi_align(
        None,
        AlignScoring {
            tolerance: mzcore::quantities::Tolerance::new_ppm(20.0),
            ..Default::default()
        },
        MultiAlignType::GLOBAL,
    );

    for line in &alignment[0] {
        let mut b = String::new();
        line.debug_display(&mut b);
        println!("{} {}", b.trim_ascii_end(), line.short());
    }

    let mut buf = String::new();
    alignment[0].debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "RWGGDGFYAMDYWGQG\nRWN·DGFYW·DYWGQG\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}

#[test]
fn bad_either_global() {
    let index = AlignIndex::<4, Peptidoform<Linear>>::new(
        [
            seq("HPVYTFG"),
            seq("[U:Carboxymethyl]-AARVTDJDJGVYYAKJGVYYAFQGSHVPYTFGG"),
            seq("MPKTGNRDJGVYYC[U:Carboxymethyl]FQGSHVPYTFGG"),
            seq("DJGVYYC[U:Carboxymethyl]YQGSHVPYTFGG"),
            seq("AGH[U:Oxidation]VQFMJJKJGVYYC[U:Carboxymethyl]FQGSHVPYTFGG"),
        ],
        MassMode::Monoisotopic,
    );
    let mmsa = index
        .multi_align(None, AlignScoring::default(), MultiAlignType {
            left: MultiAlignSide::EitherGlobal,
            right: MultiAlignSide::Global,
        })
        .pop()
        .unwrap();
    let mut buf = String::new();
    mmsa.debug_display(&mut buf);
    println!("{buf}");
    buf = buf.split('\n').skip(1).join("\n");
    let expected = "-------------------------HPVYTF-G\nMPK-----TGNR--DJGVYYCFQGSHVPYTFGG\n--------------DJGVYYCYQGSHVPYTFGG\n---AGHVQFMJJ--KJGVYYCFQGSHVPYTFGG\nAARVTDJDJGVYYAKJGVYYAFQGSHVPYTFGG\n";
    assert_eq!(buf, expected, "Expected:\n{expected}");
}
