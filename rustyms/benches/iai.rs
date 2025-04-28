use std::hint::black_box;

use rustyms::SimpleLinear;
use rustyms::align::*;
use rustyms::*;

use iai_callgrind::{
    LibraryBenchmarkConfig, Tool, ValgrindTool, library_benchmark, library_benchmark_group, main,
};

#[inline(never)]
fn setup(a: &str, b: &str) -> (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>) {
    let _force_elements_init = black_box(AminoAcid::Alanine.formulas());
    (
        Peptidoform::pro_forma(a, None)
            .unwrap()
            .into_simple_linear()
            .unwrap(),
        Peptidoform::pro_forma(b, None)
            .unwrap()
            .into_simple_linear()
            .unwrap(),
    )
}

#[inline(never)]
fn setup_simple() -> (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>) {
    setup("ANAGRS", "AGGQRS")
}

#[inline(never)]
fn setup_igha() -> (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>) {
    setup(
        "ASPTSPKVFPLSLDSTPQDGNVVVACLVQGFFPQEPLSVTWSESGQNVTARNFPPSQDASGDLYTTSSQLTLPATQCPDGKSVTCHVKHYTNSSQDVTVPCRVPPPPPCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGATFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAQPWNHGETFTCTAAHPELKTPLTANITKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTYAVTSILRVAAEDWKKGETFSCMVGHEALPLAFTQKTIDRMAGSCCVADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGKREGPQY",
        "ASPTSPKVFPLSLCSTQPDGNVVIACLVQGFFPQEPLSVTWSESGQGVTARNFPPSQDASGDLYTTSSQLTLPATQCLAGKSVTCHVKHYTNPSQDVTVPCPVPSTPPTPSPSTPPTPSPSCCHPRLSLHRPALEDLLLGSEANLTCTLTGLRDASGVTFTWTPSSGKSAVQGPPERDLCGCYSVSSVLPGCAEPWNHGKTFTCTAAYPESKTPLTATLSKSGNTFRPEVHLLPPPSEELALNELVTLTCLARGFSPKDVLVRWLQGSQELPREKYLTWASRQEPSQGTTTFAVTSILRVAAEDWKKGDTFSCMVGHEALPLAFTQKTIDRLADWQMPPPYVVLDLPQETLEEETPGANLWPTTITFLTLFLLSLFYSTALTVTSVRGPSGNREGPQY",
    )
}

#[library_benchmark]
#[bench::simple_1(setup_simple())]
#[bench::igha_1(setup_igha())]
pub fn align_1(setup: (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>)) {
    align::<1, SimpleLinear, SimpleLinear>(
        &setup.0,
        &setup.1,
        AlignScoring::default(),
        AlignType::GLOBAL,
    );
}

#[library_benchmark]
#[bench::simple_4(setup_simple())]
#[bench::igha_4(setup_igha())]
#[bench::ambiguous_not(setup("ANQRS", "ANQRS"))]
#[bench::ambiguous_a(setup("ANZRS", "ANQRS"))]
#[bench::ambiguous_b(setup("ANQRS", "ABQRS"))]
#[bench::ambiguous_ab(setup("ANZRS", "ABQRS"))]
// #[bench::igha_8(setup_igha(Some(8)))]
pub fn align_4(setup: (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>)) {
    align::<4, SimpleLinear, SimpleLinear>(
        &setup.0,
        &setup.1,
        AlignScoring::default(),
        AlignType::GLOBAL,
    );
}

#[library_benchmark]
#[bench::simple_unbounded(setup_simple())]
// #[bench::igha_8(setup_igha(Some(8)))]
pub fn align_unbounded(setup: (Peptidoform<SimpleLinear>, Peptidoform<SimpleLinear>)) {
    align::<{ u16::MAX }, SimpleLinear, SimpleLinear>(
        &setup.0,
        &setup.1,
        AlignScoring::default(),
        AlignType::GLOBAL,
    );
}

library_benchmark_group!(name = alignment; benchmarks = align_1, align_4, align_unbounded);

main!(config = LibraryBenchmarkConfig::default()
.tool(Tool::new(ValgrindTool::DHAT)).tool(Tool::new(ValgrindTool::Massif)); library_benchmark_groups = alignment);
