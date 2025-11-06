<style>
.align-cli {
  display: block;
  white-space: pre-wrap;
  font-family: monospace;
}

.r {color: #d3735d}
.dr {color: #9e432f}
.y {color: #d9d36d}
.dy {color: #9d962a}
.g {color: #00ab92}
.dg {color: #007963}
.lb {color: #84d6ff}
.b {color: #4699da}
.db {color: #0268a5}
.lp {color: #f5b2ff}
.p {color: #b677c0}
.dp {color: #84488d}
.n {color: #7E7E7E}
</style>

# Mass alignment

If you want to know how two peptidoforms/proteoforms are related and want to keep all mass-based mistakes in mind you need to align the peptidoforms with mass-alignment[^1].

## Library features

 - [Align two peptidoforms](crate::align)
 - With [custom scoring](crate::AlignScoring)
 - [Global](crate::AlignType::GLOBAL), [local](crate::AlignType::LOCAL), [either-global](crate::AlignType::EITHER_GLOBAL) and mixtures of these all
 - [Indexed alignment](crate::AlignIndex) to speed up repeated alignments against the same database
 - [Consecutive alignment](crate::consecutive_align) (only when the `imgt` feature is turned on) to create a DomainGapAlign[^2] like alignment with the mass-based alignment

## The alignment itself

A mass-based alignment handles the case in which multiple amino acids are wrong, but the total mass
of this set of amino acids is equal to the mass of a set of different amino acids on the other peptide.
This is quite common in mass spectrometry where mistakes based on mass coincidences are very common.
For example `N` has the same mass as `GG`, so if we want to make a mass spectrometry faithful alignment
of `ANA` with `AGGA` the result should reflect this fact:

<pre class="align-cli">Identity:&nbsp;0.500&nbsp;<span class="n">(2/4)</span>, Similarity:&nbsp;0.750&nbsp;<span class="n">(3/4)</span>, Gaps:&nbsp;0.000&nbsp;<span class="n">(0/4)</span>, Score:&nbsp;0.706&nbsp;<span class="n">(12/17)</span>, Equal mass, Tolerance:&nbsp;10&nbsp;ppm, Alignment:&nbsp;global
Start: A 0 B 0, Path: <span class="n">1=1:2i1=</span>
A<span class="y">N·</span>A <span class="n">A</span>
A<span class="y">GG</span>A <span class="n">B</span>
&nbsp;<span class="y">╶╴</span></pre>

_Generated using this algorithm bound to a cli tool: <https://github.com/snijderlab/align-cli>_
```rust
use rustyms::{prelude::*, sequence::SimpleLinear, align::*};
let a = Peptidoform::pro_forma("ANA", None).unwrap().into_simple_linear().unwrap();
let b = Peptidoform::pro_forma("AGGA", None).unwrap().into_simple_linear().unwrap();
let alignment = align::<4, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(&a, &b, AlignScoring::default(), AlignType::GLOBAL);
assert_eq!(alignment.short(), "1=1:2i1=");
assert_eq!(alignment.ppm().value, 0.0);
```

This extends the Smith-Waterman/Needleman-Wunsch alignment states (Match, Mismatch, Insertion, and Deletion) with the following additional states:
1. IndentityMassMismatch: the same amino acid but the mass of the sequence element is different, this happens when either of the options has a modification that is not present on the other, this could be an error or could be fine depending on the use of the alignment, use [`PairMode`](crate::PairMode) to control this behaviour.
1. Rotation: the same amino acids but in a different order, example: `AHK` on `KAH`.
1. Isobaric: different amino acids (or different modifications) that have the same mass, example `N` on `GG`, or `M[Oxidation]` on `F`.

As is visible in the examples the Rotation errors can be longer than just one amino acids, and Isobaric errors can even be of different lengths for the two peptidoforms. This means that these errors need special care to be visualised properly. Also because these need to handle these different lengths the main loop of the alignment needs to loop over all combinations of lengths from 1 to the maximum chosen length. This means that algorothmically speaking this alignment is slower than SW/NW.

## AlignTypes

Alignments can be done global (as is Needleman-Wunsch) or local (as is Smith-Waterman). Global means that both peptidoforms have to be fully matched, local means that both peptidoforms can have unmatched sequences on both sides of the match. Either global is a variation that enforces that at least one of the two peptidoforms has to be fully matched or said inversely only one peptidoform can have unmatched sequence. This can be used when it is known that sequences are related but neither is known _a priori_ to be a full superset of the other sequence. For example if you are matching a peptide to a database and the database does not contain the full protein sequence but only a part of the sequence (as is common in antibodies, those are built from three/four separate genes spliced together).

<pre class="align-cli">Identity:&nbsp;0.176&nbsp;<span class="n">(3/17)</span>, Mass&nbsp;similarity:&nbsp;0.353&nbsp;<span class="n">(6/17)</span>, Similarity:&nbsp;0.176&nbsp;<span class="n">(3/17)</span>, Gaps:&nbsp;0.118&nbsp;<span class="n">(2/17)</span>, Score:&nbsp;0.100&nbsp;<span class="n">(9/90)</span>, Mass&nbsp;difference:&nbsp;240.101&nbsp;Da&nbsp;120.775&nbsp;‰
Start: A 0 B 0, Path: <span class="n">1i3X2D2=1:2i1=6X</span>
<span class="n">Tolerance: 10 ppm, Alignment: global (==), Maximal isobaric step: 4</span>

<span class="y">I</span><span class="r">PIW</span><span class="y">KT</span>LK<span class="y">N·</span>H<span class="r">GHDDWF</span> <span class="n">A</span>
<span class="y">L</span><span class="r">FSA</span><span class="y">--</span>LK<span class="y">GG</span>H<span class="r">EIYFER</span> <span class="n">B</span>
<span class="y">─</span><span class="r">xxx</span><span class="y">++</span>  <span class="y">╶╴</span> <span class="r">xxxxxx</span></pre>
_Global alignment, the full peptidoforms are matched even if that match is quite poor_

<pre class="align-cli">Identity:&nbsp;0.600&nbsp;<span class="n">(3/5)</span>, Mass&nbsp;similarity:&nbsp;1.000&nbsp;<span class="n">(5/5)</span>, Similarity:&nbsp;0.600&nbsp;<span class="n">(3/5)</span>, Gaps:&nbsp;0.000&nbsp;<span class="n">(0/5)</span>, Score:&nbsp;0.808&nbsp;<span class="n">(21/26)</span>, Equal&nbsp;mass
Start: A 6 B 4, <span class="n">Path: 2=1:2i1=</span>
<span class="n">Tolerance: 10 ppm, Alignment: local (⁐==⁐), Maximal isobaric step: 4</span>

<span class="n">IPIWKT</span>LK<span class="y">N·</span>H<span class="n">GHDDWF A</span>
&nbsp;&nbsp;<span class="n">LFSA</span>LK<span class="y">GG</span>H<span class="n">EIYFER B</span>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<span class="y">╶╴</span></pre>
_Local alignment, the best parts of the peptidoforms are aligned_

<pre class="align-cli">Identity:&nbsp;0.357&nbsp;<span class="n">(5/14)</span>, Mass&nbsp;similarity:&nbsp;0.429&nbsp;<span class="n">(6/14)</span>, Similarity:&nbsp;0.357&nbsp;<span class="n">(5/14)</span>, Gaps:&nbsp;0.071&nbsp;<span class="n">(1/14)</span>, Score:&nbsp;0.250&nbsp;<span class="n">(20/80)</span>, Mass&nbsp;difference:&nbsp;315.108&nbsp;Da&nbsp;177.239&nbsp;‰
Start: A 2 B 0, Path: <span class="n">1i3X2=1X1D2=3X1=</span>
<span class="n">Tolerance: 10 ppm, Alignment: either global (-==-), Maximal isobaric step: 4</span>

<span class="n">IP</span><span class="y">I</span><span class="r">WKT</span>LK<span class="r">N</span><span class="y">H</span>GH<span class="r">DDW</span>F   <span class="n">A</span>
&nbsp;&nbsp;<span class="y">L</span><span class="r">FSA</span>LK<span class="r">G</span><span class="y">-</span>GH<span class="r">EIY</span>F<span class="n">ER B</span>
&nbsp;&nbsp;<span class="y">─</span><span class="r">xxx</span>  <span class="r">x</span><span class="y">+</span>  <span class="r">xxx</span></pre>
_Either global, for both sides (left and right) only one peptidoform can have unmatched sequence_

In this crate [`AlignType`](crate::AlignType) can be used to set these types of alingments. Note that this allows setting one of these types per side and per peptidoform. This allows for encoding as much knowledge into the alignment as possible. For example if one would build a tool to align an antibody sequence to germline V genes you could make it fixed to a global left alignment (emaning both sequences have to start together) and end with an either global right alignment if it is not known if the user will supply full antibody sequences or might even supply truncated V gene sequences as shown below.

<pre class="align-cli">Selected: Mus musculus House mouse IGHV1-18-26*01 / IgV1-18-26*01
Identity:&nbsp;0.847&nbsp;<span class="n">(83/98)</span>, Mass&nbsp;similarity:&nbsp;0.888&nbsp;<span class="n">(87/98)</span>, Similarity:&nbsp;0.847&nbsp;<span class="n">(83/98)</span>, Gaps:&nbsp;0.000&nbsp;<span class="n">(0/98)</span>, Score:&nbsp;0.848&nbsp;<span class="n">(440/519)</span>, Mass difference:&nbsp;127.979&nbsp;Da&nbsp;11.715&nbsp;‰
Start: IGHV1-18-26*01 0 Query 0, Path: <span class="n">26=1X4=3X6=1X4=1X1=1X6=4r3=1X2=1X8=1X2=1X21=</span>
<span class="n">Tolerance: 10 ppm, Alignment: global left, either global right (==-), Maximal isobaric step: 4</span>

FR1<span class="n">      10       20     </span>CDR1    FR2<span class="n">  40        50</span>
EVQLQQSGPELVKPGASVKISCKASG<span class="r">Y</span>TFTD<span class="r">YNM</span>HWVKQS<span class="r">H</span>GKSL<span class="r">E</span>W<span class="r">I</span>GY <span class="n">IGHV1-18-26*01</span>
EVQLQQSGPELVKPGASVKISCKASG<span class="r">F</span>TFTD<span class="r">FSI</span>HWVKQS<span class="r">Q</span>GKSL<span class="r">D</span>W<span class="r">V</span>GY <span class="n">Query</span>
                          <span class="r">x    xxx      x    x x</span>
CDR2    FR3<span class="n">       70        80        90      </span>C3
IYPY<span class="y">NGGT</span>GYN<span class="r">Q</span>KF<span class="r">K</span>SKATLTVD<span class="r">N</span>SS<span class="r">S</span>TAYMELRSLTSEDSAVYYCAR   <span class="n">IGHV1-18-26*01</span>
IYPY<span class="y">TGGN</span>GYN<span class="r">L</span>KF<span class="r">Q</span>SKATLTVD<span class="r">T</span>SS<span class="r">T</span>TAYMELRSLTSEDSAVYYCAR<span class="n">RE Query</span>
    <span class="y">╶──╴</span>   <span class="r">x  x        x  x</span>
<span class="n">GNFVGAMDYWGQGTSVTVSSAKTTAPSVYPLAPVCGDTTGSSVTLGCLVK Query
GYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSI Query
TCNVAHPASSTKVDKKIEPRGPTIKPCPPCKCPAPNLLGGPSVFIFPPKI Query
KDVLMISLSPIVTCVVVDVSEDDPDVQISWFVNNVEVHTAQTQTHREDYN Query
STLRVVSALPIQHQDWMSGKEFKCKVNNKDLPAPIERTISKPKGSVRAPQ Query
VYVLPPPEEEMTKKQVTLTCMVTDFMPEDIYVEWTNNGKTELNYKNTEPV Query
LDSDGSYFMYSKLRVEKKNWVERNSYSCSVVHEGLHNHHTTKSFSRTPGK Query</span></pre>
_Aligning a V gene with a full antibody sequence_

<pre class="align-cli">Selected: Mus musculus House mouse IGHV1-18-26*01 / IgV1-18-26*01
Identity:&nbsp;0.759&nbsp;<span class="n">(63/83)</span>, Mass&nbsp;similarity:&nbsp;0.807&nbsp;<span class="n">(67/83)</span>, Similarity:&nbsp;0.759&nbsp;<span class="n">(63/83)</span>, Gaps:&nbsp;0.060&nbsp;<span class="n">(5/83)</span>, Score:&nbsp;0.770&nbsp;<span class="n">(331/430)</span>, Mass difference:&nbsp;725.292&nbsp;Da&nbsp;78.644&nbsp;‰
Start: IGHV1-18-26*01 0 Query 0, Path: <span class="n">5D21=1X4=3X6=1X4=1X1=1X6=4r3=1X2=1X8=1X2=1X6=</span>
<span class="n">Tolerance: 10 ppm, Alignment: global left, either global right (==-), Maximal isobaric step: 4</span>

FR1<span class="n">      10       20     </span>CDR1    FR2<span class="n">  40        50</span>
<span class="y">EVQLQ</span>QSGPELVKPGASVKISCKASG<span class="r">Y</span>TFTD<span class="r">YNM</span>HWVKQS<span class="r">H</span>GKSL<span class="r">E</span>W<span class="r">I</span>GY <span class="n">IGHV1-18-26*01</span>
<span class="y">-----</span>QSGPELVKPGASVKISCKASG<span class="r">F</span>TFTD<span class="r">FSI</span>HWVKQS<span class="r">Q</span>GKSL<span class="r">D</span>W<span class="r">V</span>GY <span class="n">Query</span>
<span class="y">     </span>                     <span class="r">x    xxx      x    x x</span>
CDR2    FR3<span class="n">       70        80        90      </span>C3
IYPY<span class="y">NGGT</span>GYN<span class="r">Q</span>KF<span class="r">K</span>SKATLTVD<span class="r">N</span>SS<span class="r">S</span>TAYMEL<span class="n">RSLTSEDSAVYYCAR   IGHV1-18-26*01</span>
IYPY<span class="y">TGGN</span>GYN<span class="r">L</span>KF<span class="r">Q</span>SKATLTVD<span class="r">T</span>SS<span class="r">T</span>TAYMEL                  <span class="n">Query</span>
    <span class="y">╶──╴</span>   <span class="r">x  x        x  x</span></pre>
_Aligning a V gene with a beginning and end truncated antibody sequence, showing that the missing start is scored negatively as a leading deletion but the missing end is scored neutral as an allowed truncation_

## Example usage

```rust
# fn main() -> Result<(), context_error::BoxedError<'static, context_error::BasicKind>> {
use mzcore::{prelude::*, sequence::SimpleLinear};
use mzalign::{prelude::*, align};
// Check how this peptide compares to a similar peptide (using the feature `align`)
let first_peptide = Peptidoform::pro_forma("IVQEVT", None)?.into_simple_linear().unwrap();
let second_peptide = Peptidoform::pro_forma("LVQVET", None)?.into_simple_linear().unwrap();
// Align the two peptides using mass based alignment
// IVQEVT A
// LVQVET B
// ─  ╶╴
let alignment = align::<4, &Peptidoform<SimpleLinear>, &Peptidoform<SimpleLinear>>(
                  &first_peptide,
                  &second_peptide,
                  AlignScoring::default(),
                  AlignType::GLOBAL);
# dbg!(&alignment);
// Calculate some more statistics on this alignment
let stats = alignment.stats();
assert_eq!(stats.mass_similar, 6); // 6 out of the 6 positions are mass similar
# Ok(()) }
```

## Compilation features

* `imgt` - enables access to the IMGT database of antibodies germline sequences, with annotations. This also turns on the use of consecutive alignments.
* `rayon` - enables parallel iterators using rayon, this also turns on rayon in `imgt`.

[^1]: Schulte, D.; Snijder, J. A Handle on Mass Coincidence Errors in De Novo Sequencing of Antibodies by Bottom-up Proteomics. <https://doi.org/10.1021/acs.jproteome.4c00188>.
[^2]: Ehrenmann, F.; Kaas, Q.; Lefranc, M.-P. IMGT/3Dstructure-DB and IMGT/DomainGapAlign: A Database and a Tool for Immunoglobulins or Antibodies, T Cell Receptors, MHC, IgSF and MhcSF. <https://doi.org/10.1093/nar/gkp946>.
