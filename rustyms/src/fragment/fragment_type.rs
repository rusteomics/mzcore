use std::{
    borrow::Cow,
    cmp::Ordering,
    fmt::{Debug, Display},
};

use itertools::Itertools;
use serde::{Deserialize, Serialize};

#[cfg(feature = "glycan-render")]
use crate::glycan::GlycanSelection;
use crate::{
    fragment::{DiagnosticPosition, PeptidePosition, SatelliteLabel},
    glycan::{GlycanPosition, MonoSaccharide},
    sequence::{AminoAcid, IsAminoAcid, SemiAmbiguous, SequenceElement, SequencePosition},
};

/// The possible types of fragments
#[derive(Clone, Debug, Default, Deserialize, Eq, Hash, PartialEq, Serialize)]
#[expect(non_camel_case_types)]
pub enum FragmentType {
    /// a
    a(PeptidePosition, i8),
    /// b
    b(PeptidePosition, i8),
    /// c
    c(PeptidePosition, i8),
    /// d a position, originating amino acid, distance from a break, variant,
    d(PeptidePosition, AminoAcid, u8, i8, SatelliteLabel),
    /// v
    v(PeptidePosition, AminoAcid, u8, i8),
    /// w
    w(PeptidePosition, AminoAcid, u8, i8, SatelliteLabel),
    /// x
    x(PeptidePosition, i8),
    /// y
    y(PeptidePosition, i8),
    /// z
    z(PeptidePosition, i8),
    // glycan A fragment (Never generated)
    //A(GlycanPosition),
    /// glycan B fragment
    // B(GlycanPosition),
    // glycan C fragment (Never generated)
    //C(GlycanPosition),
    // glycan X fragment (Never generated)
    //X(GlycanPosition),
    /// glycan Y fragment, generated by one or more branches broken
    Y(Vec<GlycanPosition>),
    // glycan Z fragment (Never generated)
    // Z(GlycanPosition),
    /// B glycan fragment, potentially with additional Y breakages
    B {
        /// The root break
        b: GlycanPosition,
        /// The branch breakages
        y: Vec<GlycanPosition>,
        /// All branches that are not broken
        end: Vec<GlycanPosition>,
    },
    /// A B or internal glycan fragment for a glycan where only the composition is known, also saves the attachment (AA + sequence index)
    BComposition(
        Vec<(MonoSaccharide, isize)>,
        Option<(AminoAcid, SequencePosition)>,
    ),
    /// A B or internal glycan fragment for a glycan where only the composition is known, also saves the attachment (AA + sequence index)
    YComposition(
        Vec<(MonoSaccharide, isize)>,
        Option<(AminoAcid, SequencePosition)>,
    ),
    /// Immonium ion
    Immonium(Option<PeptidePosition>, SequenceElement<SemiAmbiguous>),
    /// Precursor with amino acid side chain loss
    PrecursorSideChainLoss(PeptidePosition, AminoAcid),
    /// Diagnostic ion for a given position
    Diagnostic(DiagnosticPosition),
    /// An internal fragment, potentially with the named bonds that resulted in this fragment
    Internal(
        Option<(BackboneNFragment, BackboneCFragment)>,
        PeptidePosition,
        PeptidePosition,
    ),
    /// An unknown series, with potentially the series number
    Unknown(Option<usize>),
    /// precursor
    #[default]
    Precursor,
}

impl Ord for FragmentType {
    #[allow(clippy::match_same_arms)] // This ordering is helpful for human readers
    fn cmp(&self, other: &Self) -> Ordering {
        // Sort of type first (precursor/abcxyz/dw/v)
        match (self, other) {
            // Peptide
            (Self::Precursor, Self::Precursor) => Ordering::Equal,
            (Self::Precursor, _) => Ordering::Less,
            (_, Self::Precursor) => Ordering::Greater,
            (Self::a(s, sv), Self::a(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::a(_, _), _) => Ordering::Less,
            (_, Self::a(_, _)) => Ordering::Greater,
            (Self::b(s, sv), Self::b(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::b(_, _), _) => Ordering::Less,
            (_, Self::b(_, _)) => Ordering::Greater,
            (Self::c(s, sv), Self::c(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::c(_, _), _) => Ordering::Less,
            (_, Self::c(_, _)) => Ordering::Greater,
            (Self::x(s, sv), Self::x(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::x(_, _), _) => Ordering::Less,
            (_, Self::x(_, _)) => Ordering::Greater,
            (Self::y(s, sv), Self::y(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::y(_, _), _) => Ordering::Less,
            (_, Self::y(_, _)) => Ordering::Greater,
            (Self::z(s, sv), Self::z(o, ov)) => s.cmp(o).then(sv.cmp(ov)),
            (Self::z(_, _), _) => Ordering::Less,
            (_, Self::z(_, _)) => Ordering::Greater,
            (Self::d(s, _, sd, sv, sl), Self::d(o, _, od, ov, ol)) => {
                s.cmp(o).then(sd.cmp(od)).then(sv.cmp(ov)).then(sl.cmp(ol))
            }
            (Self::d(_, _, _, _, _), _) => Ordering::Less,
            (_, Self::d(_, _, _, _, _)) => Ordering::Greater,
            (Self::w(s, _, sd, sv, sl), Self::w(o, _, od, ov, ol)) => {
                s.cmp(o).then(sd.cmp(od)).then(sv.cmp(ov)).then(sl.cmp(ol))
            }
            (Self::w(_, _, _, _, _), _) => Ordering::Less,
            (_, Self::w(_, _, _, _, _)) => Ordering::Greater,
            (Self::v(s, _, sd, sv), Self::v(o, _, od, ov)) => {
                s.cmp(o).then(sd.cmp(od)).then(sv.cmp(ov))
            }
            (Self::v(_, _, _, _), _) => Ordering::Less,
            (_, Self::v(_, _, _, _)) => Ordering::Greater,
            (Self::Immonium(s, _), Self::Immonium(o, _)) => s.cmp(o),
            (Self::Immonium(_, _), _) => Ordering::Less,
            (_, Self::Immonium(_, _)) => Ordering::Greater,
            (Self::PrecursorSideChainLoss(s, _), Self::PrecursorSideChainLoss(o, _)) => s.cmp(o),
            (Self::PrecursorSideChainLoss(_, _), _) => Ordering::Less,
            (_, Self::PrecursorSideChainLoss(_, _)) => Ordering::Greater,
            (Self::Internal(st, sa, sb), Self::Internal(ot, oa, ob)) => {
                sa.cmp(oa).then(sb.cmp(ob)).then(st.cmp(ot))
            }
            (Self::Internal(_, _, _), _) => Ordering::Less,
            (_, Self::Internal(_, _, _)) => Ordering::Greater,
            // Glycans
            (Self::B { b: sb, y: sy, .. }, Self::B { b: ob, y: oy, .. }) => {
                sy.len().cmp(&oy.len()).then(sb.cmp(ob))
            }
            (Self::Y(s), Self::Y(o)) => s.len().cmp(&o.len()),
            (Self::B { y: sy, .. }, Self::Y(o)) => {
                (sy.len() + 1).cmp(&o.len()).then(Ordering::Greater)
            }
            (Self::Y(s), Self::B { y: oy, .. }) => {
                s.len().cmp(&(oy.len() + 1)).then(Ordering::Less)
            }
            (Self::B { .. }, _) => Ordering::Less,
            (_, Self::B { .. }) => Ordering::Greater,
            (Self::Y(_), _) => Ordering::Less,
            (_, Self::Y(_)) => Ordering::Greater,
            (Self::BComposition(s, sl), Self::BComposition(o, ol))
            | (Self::YComposition(s, sl), Self::YComposition(o, ol)) => {
                s.len().cmp(&o.len()).then(sl.cmp(ol))
            }
            (Self::BComposition(s, sl), Self::YComposition(o, ol)) => s
                .len()
                .cmp(&o.len())
                .then(sl.cmp(ol))
                .then(Ordering::Greater),
            (Self::YComposition(s, sl), Self::BComposition(o, ol)) => {
                s.len().cmp(&o.len()).then(sl.cmp(ol)).then(Ordering::Less)
            }
            (Self::BComposition(_, _), _) => Ordering::Less,
            (_, Self::BComposition(_, _)) => Ordering::Greater,
            (Self::YComposition(_, _), _) => Ordering::Less,
            (_, Self::YComposition(_, _)) => Ordering::Greater,
            // Other
            (Self::Diagnostic(s), Self::Diagnostic(o)) => s.cmp(o),
            (Self::Diagnostic(_), _) => Ordering::Less,
            (_, Self::Diagnostic(_)) => Ordering::Greater,
            (Self::Unknown(s), Self::Unknown(o)) => s.cmp(o),
        }
    }
}

impl PartialOrd for FragmentType {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl FragmentType {
    /// Get a main ion series fragment with the specified variant, or pass the fragment type through unchanged
    #[must_use]
    pub fn with_variant(&self, variant: i8) -> Self {
        match self {
            Self::a(p, _) => Self::a(*p, variant),
            Self::b(p, _) => Self::b(*p, variant),
            Self::c(p, _) => Self::c(*p, variant),
            Self::d(p, a, d, _, l) => Self::d(*p, *a, *d, variant, *l),
            Self::v(p, a, d, _) => Self::v(*p, *a, *d, variant),
            Self::w(p, a, d, _, l) => Self::w(*p, *a, *d, variant, *l),
            Self::x(p, _) => Self::x(*p, variant),
            Self::y(p, _) => Self::y(*p, variant),
            Self::z(p, _) => Self::z(*p, variant),
            other => other.clone(),
        }
    }

    /// Get the position of this ion (or None if it is not known)
    pub const fn position(&self) -> Option<&PeptidePosition> {
        match self {
            Self::a(n, _)
            | Self::b(n, _)
            | Self::c(n, _)
            | Self::d(n, _, _, _, _)
            | Self::v(n, _, _, _)
            | Self::w(n, _, _, _, _)
            | Self::x(n, _)
            | Self::y(n, _)
            | Self::z(n, _)
            | Self::Diagnostic(DiagnosticPosition::Peptide(n, _))
            | Self::PrecursorSideChainLoss(n, _) => Some(n),
            Self::Immonium(n, _) => n.as_ref(),
            _ => None,
        }
    }

    /// Get the root glycan position of this ion (or None if not applicable), Y is not defined as it does not have a root break
    pub const fn glycan_position(&self) -> Option<&GlycanPosition> {
        match self {
            Self::Diagnostic(DiagnosticPosition::Glycan(b, _)) | Self::B { b, .. } => Some(b),
            _ => None,
        }
    }

    /// Get the glycan break positions of this ion (or None if not applicable), gives the sequence index, the root break, and the branch breaks.
    /// Only available with feature 'glycan-render'.
    #[cfg(feature = "glycan-render")]
    pub fn glycan_break_positions(
        &self,
    ) -> Option<(Option<SequencePosition>, GlycanSelection<'_>)> {
        match self {
            Self::Diagnostic(DiagnosticPosition::Glycan(n, _)) => Some((
                n.attachment.map(|(_, p)| p),
                GlycanSelection::SingleSugar(n),
            )),
            Self::Y(breaks) => Some((
                breaks.first().and_then(|p| p.attachment.map(|(_, p)| p)),
                GlycanSelection::Subtree(None, breaks),
            )),
            Self::B { b, y, .. } => Some((
                b.attachment.map(|(_, p)| p),
                GlycanSelection::Subtree(Some(b), y),
            )),
            _ => None,
        }
    }

    /// Get the position label, unless it is a precursor ion
    pub fn position_label(&self) -> Option<String> {
        match self {
            Self::a(n, _)
            | Self::b(n, _)
            | Self::c(n, _)
            | Self::d(n, _, _, _, _)
            | Self::v(n, _, _, _)
            | Self::w(n, _, _, _, _)
            | Self::x(n, _)
            | Self::y(n, _)
            | Self::z(n, _)
            | Self::Diagnostic(DiagnosticPosition::Peptide(n, _))
            | Self::PrecursorSideChainLoss(n, _) => Some(n.series_number.to_string()),
            Self::Immonium(n, _) => n.map(|n| n.series_number.to_string()),
            Self::Diagnostic(DiagnosticPosition::Glycan(n, _)) => Some(n.label()),
            Self::Y(bonds) => Some(bonds.iter().map(GlycanPosition::label).join("Y")),
            Self::B { b, y, end } => Some(
                b.label()
                    + "Y"
                    + &y.iter()
                        .chain(end.iter())
                        .map(GlycanPosition::label)
                        .join("Y"),
            ),
            Self::YComposition(sugars, _) | Self::BComposition(sugars, _) => Some(
                sugars
                    .iter()
                    .map(|(sugar, amount)| format!("{sugar}{amount}"))
                    .join(""),
            ),
            Self::Internal(_, pos1, pos2) => {
                Some(format!("{}:{}", pos1.sequence_index, pos2.sequence_index,))
            }
            Self::Precursor
            | Self::Unknown(_)
            | Self::Diagnostic(
                DiagnosticPosition::Labile(_)
                | DiagnosticPosition::GlycanCompositional(_, _)
                | DiagnosticPosition::Reporter,
            ) => None,
        }
    }

    /// Get the label for this fragment type, the first argument is the optional superscript prefix, the second is the main label
    pub fn label(&self) -> (Option<String>, Cow<str>) {
        let get_label = |ion: &'static str, v: i8| {
            if v == 0 {
                Cow::Borrowed(ion)
            } else {
                Cow::Owned(format!(
                    "{ion}{}",
                    if v < 0 {
                        "\'".repeat((-v) as usize)
                    } else {
                        "·".repeat(v as usize)
                    }
                ))
            }
        };

        match self {
            Self::a(_, v) => (None, get_label("a", *v)),
            Self::b(_, v) => (None, get_label("b", *v)),
            Self::c(_, v) => (None, get_label("c", *v)),
            Self::d(_, _, n, v, l) => (
                (*n != 0).then_some(n.to_string()),
                Cow::Owned(format!(
                    "d{l}{}",
                    if *v < 0 {
                        "\'".repeat((-v) as usize)
                    } else {
                        "·".repeat(*v as usize)
                    }
                )),
            ),
            Self::v(_, _, n, v) => ((*n != 0).then_some(n.to_string()), get_label("v", *v)),
            Self::w(_, _, n, v, l) => (
                (*n != 0).then_some(n.to_string()),
                Cow::Owned(format!(
                    "w{l}{}",
                    if *v < 0 {
                        "\'".repeat((-v) as usize)
                    } else {
                        "·".repeat(*v as usize)
                    }
                )),
            ),
            Self::x(_, v) => (None, get_label("x", *v)),
            Self::y(_, v) => (None, get_label("y", *v)),
            Self::z(_, v) => (None, get_label("z", *v)),
            Self::B { .. } | Self::BComposition(_, _) => (None, Cow::Borrowed("B")),
            Self::Y(_) | Self::YComposition(_, _) => (None, Cow::Borrowed("Y")),
            Self::Diagnostic(DiagnosticPosition::Peptide(_, aa)) => (
                None,
                Cow::Owned(
                    aa.one_letter_code()
                        .map(|c| format!("d{c}"))
                        .or_else(|| aa.three_letter_code().map(|c| format!("d{c}")))
                        .unwrap_or_else(|| format!("d{}", aa.name())),
                ),
            ),
            Self::Diagnostic(DiagnosticPosition::Reporter) => (None, Cow::Borrowed("r")),
            Self::Diagnostic(DiagnosticPosition::Labile(m)) => (None, Cow::Owned(format!("d{m}"))),
            Self::Diagnostic(
                DiagnosticPosition::Glycan(_, sug)
                | DiagnosticPosition::GlycanCompositional(sug, _),
            ) => (None, Cow::Owned(format!("d{sug}"))),
            Self::Immonium(_, aa) => (
                None,
                Cow::Owned(
                    aa.aminoacid
                        .one_letter_code()
                        .map(|c| format!("i{c}"))
                        .or_else(|| aa.aminoacid.three_letter_code().map(|c| format!("i{c}")))
                        .unwrap_or_else(|| format!("i{}", aa.aminoacid.name())),
                ),
            ),
            Self::PrecursorSideChainLoss(_, aa) => (
                None,
                Cow::Owned(
                    aa.one_letter_code()
                        .map(|c| format!("p-s{c}"))
                        .or_else(|| aa.three_letter_code().map(|c| format!("p-s{c}")))
                        .unwrap_or_else(|| format!("p-s{}", aa.name())),
                ),
            ),
            Self::Precursor => (None, Cow::Borrowed("p")),
            Self::Internal(fragmentation, _, _) => (
                None,
                Cow::Owned(format!(
                    "m{}",
                    fragmentation.map_or(String::new(), |(n, c)| format!("{n}:{c}")),
                )),
            ),
            Self::Unknown(series) => (
                None,
                Cow::Owned(format!(
                    "?{}",
                    series.map_or(String::new(), |s| s.to_string()),
                )),
            ),
        }
    }

    /// Get the kind of fragment, easier to match against
    pub const fn kind(&self) -> FragmentKind {
        match self {
            Self::a(_, _) => FragmentKind::a,
            Self::b(_, _) => FragmentKind::b,
            Self::c(_, _) => FragmentKind::c,
            Self::d(_, _, _, _, _) => FragmentKind::d,
            Self::v(_, _, _, _) => FragmentKind::v,
            Self::w(_, _, _, _, _) => FragmentKind::w,
            Self::x(_, _) => FragmentKind::x,
            Self::y(_, _) => FragmentKind::y,
            Self::z(_, _) => FragmentKind::z,
            Self::Y(_) | Self::YComposition(_, _) => FragmentKind::Y,
            Self::Diagnostic(
                DiagnosticPosition::Glycan(_, _) | DiagnosticPosition::GlycanCompositional(_, _),
            )
            | Self::B { .. }
            | Self::BComposition(_, _) => FragmentKind::B,
            Self::Diagnostic(_) => FragmentKind::diagnostic,
            Self::Immonium(_, _) => FragmentKind::immonium,
            Self::PrecursorSideChainLoss(_, _) => FragmentKind::precursor_side_chain_loss,
            Self::Precursor => FragmentKind::precursor,
            Self::Internal(_, _, _) => FragmentKind::internal,
            Self::Unknown(_) => FragmentKind::unknown,
        }
    }
}

impl Display for FragmentType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let (sup, label) = self.label();
        write!(
            f,
            "{}{}{}",
            sup.unwrap_or_default(),
            label,
            self.position_label().unwrap_or_default()
        )
    }
}

/// The possible kinds of N terminal backbone fragments.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[expect(non_camel_case_types)]
pub enum BackboneNFragment {
    /// a
    a,
    /// b
    b,
    /// c
    c,
}

impl Display for BackboneNFragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a => "a",
                Self::b => "b",
                Self::c => "c",
            }
        )
    }
}

/// The possible kinds of C terminal backbone fragments.
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[expect(non_camel_case_types)]
pub enum BackboneCFragment {
    /// x
    x,
    /// y
    y,
    /// z and z·
    z,
}

impl Display for BackboneCFragment {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::x => "x",
                Self::y => "y",
                Self::z => "z",
            }
        )
    }
}

/// The possible kinds of fragments, same options as [`FragmentType`] but without any additional data
#[derive(Clone, Copy, Debug, Deserialize, Eq, Hash, Ord, PartialEq, PartialOrd, Serialize)]
#[expect(non_camel_case_types)]
pub enum FragmentKind {
    /// a
    a,
    /// b
    b,
    /// c
    c,
    /// d
    d,
    /// v
    v,
    /// w
    w,
    /// x
    x,
    /// y
    y,
    /// z and z·
    z,
    /// glycan Y fragment, generated by one or more branches broken
    Y,
    /// B or glycan diagnostic ion or Internal glycan fragment, meaning both a B and Y breakages (and potentially multiple of both), resulting in a set of monosaccharides
    B,
    /// Immonium ion
    immonium,
    /// Precursor with amino acid side chain loss
    precursor_side_chain_loss,
    /// Diagnostic ion for a given position
    diagnostic,
    /// Internal ion
    internal,
    /// precursor
    precursor,
    /// unknown fragment
    unknown,
}

impl Display for FragmentKind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Self::a => "a",
                Self::b => "b",
                Self::c => "c",
                Self::d => "d",
                Self::x => "x",
                Self::y => "y",
                Self::v => "v",
                Self::w => "w",
                Self::z => "z",
                Self::Y => "Y",
                Self::B => "oxonium",
                Self::immonium => "immonium",
                Self::precursor_side_chain_loss => "precursor side chain loss",
                Self::diagnostic => "diagnostic",
                Self::internal => "m",
                Self::precursor => "precursor",
                Self::unknown => "unknown",
            }
        )
    }
}
