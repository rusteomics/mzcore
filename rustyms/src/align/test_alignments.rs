#![allow(clippy::missing_panics_doc)]

use crate::{
    align::{align, scoring::AlignScoring, AlignType, Alignment},
    sequence::{Peptidoform, SimpleLinear},
};

#[test]
fn global_b() {
    test_alignment(
        "AEDTAVYYC[Carboxymethyl]SRWGGDGFYA",
        "AEDTAVYYM[oxidation]QSW[oxidation]SR",
        AlignScoring::default(),
        AlignType::GLOBAL_B,
        "8=1X3I2=",
    );
}

#[test]
fn global_a_1() {
    test_alignment(
        "HJKASFHLKJFHAS",
        "JHSLAFKJHFLKSFJ",
        AlignScoring::default(),
        AlignType::GLOBAL_A,
        "1=2X1=1X2I2r2=1X1=1X2D",
    );
}

#[test]
fn global_a_2() {
    test_alignment(
        "AAASSS",
        "ASSS",
        AlignScoring::default(),
        AlignType::GLOBAL_A,
        "1=2D3=",
    );
}

#[test]
fn either_global_1() {
    test_alignment(
        "HHHHHHAA",
        "AAHHHHHHH",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "6=",
    );
}

#[test]
fn either_global_2() {
    test_alignment(
        "EPELVESGGGLAQPGTSLRLSCEASGFTFDDYAMHWVRQAPGRALEWVSGISWSSDNLAYSDSVEGRFTISRDNAKNSLYLQMNSLRLDDTAFYYCAKDVRPYYDFWAFDSWGRGTPVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        "DLQLVESGGGLVQPGGSLRLSC[carboxymethyl]AASGF",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "3X8=1X3=1X5=2i4=",
    );
}

#[test]
fn either_global_3() {
    test_alignment(
        "EPELVESGGGLAQPGTSLRLSCEASGFTFDDYAMHWVRQAPGRALEWVSGISWSSDNLAYSDSVEGRFTISRDNAKNSLYLQMNSLRLDDTAFYYCAKDVRPYYDFWAFDSWGRGTPVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKRVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        "TSGVHTAPFVLQSSGLY",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "6=3r8=",
    );
}

#[test]
fn either_global_4() {
    test_alignment(
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        "SGFGGLKN[deamidated]TYLHW",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "3=1:2i1i1=1i2=1i2=",
    );
}

#[test]
fn either_global_5() {
    test_alignment(
        "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKVEPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPG",
        "DTLM[oxidation]LSRTPEVTC[carbamidomethyl]VVVN[deamidated]VSHEN[deamidated]PEVK",
        AlignScoring::default(),
        AlignType::EITHER_GLOBAL,
        "3=1m1i7=1m3=1i4=1i4=",
    );
}

#[test]
fn normal_1() {
    test_alignment(
        "ANNA",
        "AGGGGA",
        AlignScoring {
            mismatch: 12,
            ..Default::default()
        },
        AlignType::LOCAL,
        "4X",
    );
}

#[test]
fn normal_2() {
    test_alignment(
        "ADA",
        "AN[deamidated]A",
        AlignScoring::default(),
        AlignType::LOCAL,
        "1=1i1=",
    );
}

#[test]
fn overextended_rotation() {
    test_alignment(
        "IVQEVS",
        "LEVQVES",
        AlignScoring::default(),
        AlignType::GLOBAL,
        "1i1I2=2r1=",
    );
}

#[test]
fn no_detected_rotation() {
    test_alignment(
        "AGGHT",
        "ANTH",
        AlignScoring::default(),
        AlignType::GLOBAL,
        "1=2:1i2r",
    );
}

#[test]
fn no_detected_rotation_2() {
    test_alignment(
        "AGGHTK",
        "ANTHK",
        AlignScoring::default(),
        AlignType::GLOBAL,
        "1=2:1i2r1=",
    );
}

/// Test if the given alignment is as expected and can be recreated
/// # Errors
/// When the alignment is not identical to path and when the alignment cannot be recreated from the path.
fn test_alignment(
    seq_a: &str,
    seq_b: &str,
    scoring: AlignScoring<'_>,
    align_type: AlignType,
    path: &str,
) {
    const MAXIMAL_STEP: u16 = 4;
    let first_peptide = Peptidoform::pro_forma(seq_a, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let second_peptide = Peptidoform::pro_forma(seq_b, None)
        .unwrap()
        .into_simple_linear()
        .unwrap();
    let alignment = align::<MAXIMAL_STEP, SimpleLinear, SimpleLinear>(
        &first_peptide,
        &second_peptide,
        scoring,
        align_type,
    );
    assert_eq!(
        alignment.short(),
        path,
        "Alignment of {seq_a} vs {seq_b} did not produce the expected alignment"
    );

    let recreated = Alignment::create_from_path(
        &first_peptide,
        &second_peptide,
        alignment.start_a,
        alignment.start_b,
        &alignment.short(),
        scoring,
        align_type,
        MAXIMAL_STEP,
    )
    .unwrap();
    // First compare paths because thta gives smaller error messages
    assert_eq!(recreated.path, alignment.path);
    assert_eq!(recreated, alignment);
}
