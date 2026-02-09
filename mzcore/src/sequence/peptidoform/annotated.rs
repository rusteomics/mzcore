use std::ops::Bound;

use crate::{
    prelude::Peptidoform,
    sequence::{AtMax, HasPeptidoformImpl, Linear},
    space::UsedSpace,
};

/// An annotated peptidoform
pub trait AnnotatedPeptidoform: HasPeptidoformImpl {
    /// Get the regions, as a list of the regions in order with the length of each region, these are
    /// required to be as long as the full peptidoform.
    fn regions(&self) -> &[(Region, usize)];
    /// Get the annotations, specified as the annotation and the into the peptidoform where it is
    /// located.
    fn annotations(&self) -> &[(Annotation, usize)];

    /// Get the region for a specific index into the sequence, None if outside range,
    /// the additional bool indicates if this is the starting position for the region.
    fn get_region(&self, index: usize) -> Option<(&Region, bool)> {
        let regions = self.regions();
        let mut left = index;
        let mut regions_index = 0;
        let mut next = &regions[regions_index];
        while left > next.1 {
            left -= next.1;
            regions_index += 1;
            if regions_index == regions.len() {
                return None;
            }
            next = &regions[regions_index];
        }
        Some((&next.0, left == 1))
    }

    /// Get all annotations for this position
    fn get_annotations(&self, index: usize) -> impl Iterator<Item = &Annotation> + '_ {
        self.annotations()
            .iter()
            .filter(move |a| a.1 == index)
            .map(|a| &a.0)
    }

    /// Get a sub peptidoform from this annotated peptidoform. It returns None if the range is outside of bounds of the sequence.
    fn sub_peptidoform(
        &self,
        range: impl std::ops::RangeBounds<usize> + Clone,
    ) -> Option<(
        Peptidoform<Self::Complexity>,
        Vec<(Region, usize)>,
        Vec<(Annotation, usize)>,
    )>
    where
        Self::Complexity: AtMax<Linear>,
    {
        self.peptidoform()
            .sub_peptidoform(range.clone())
            .map(|peptidoform| {
                let start = match range.start_bound() {
                    Bound::Excluded(e) => e + 1,
                    Bound::Included(i) => *i,
                    Bound::Unbounded => 0,
                };
                let end = match range.end_bound() {
                    Bound::Excluded(e) => e - 1,
                    Bound::Included(i) => *i,
                    Bound::Unbounded => 0,
                };
                let mut index = 0;
                let mut regions = Vec::new();

                for (r, len) in self.regions() {
                    if index + len < start {
                        // Do nothing
                    } else if index < start && index + len > start {
                        if index + len > end {
                            regions.push((r.clone(), end - start));
                            break;
                        } else {
                            regions.push((r.clone(), index + len - start));
                        }
                    } else if index + len <= end {
                        regions.push((r.clone(), *len));
                    } else {
                        regions.push((r.clone(), len - (index + len - start)));
                        break;
                    }
                    index += len;
                }

                let res: (
                    Peptidoform<<Self as HasPeptidoformImpl>::Complexity>,
                    Vec<(Region, usize)>,
                    Vec<(Annotation, usize)>,
                ) = (
                    peptidoform,
                    regions,
                    self.annotations()
                        .iter()
                        .filter(|(_, i)| range.contains(i))
                        .cloned()
                        .collect(),
                );

                debug_assert_eq!(res.0.len(), res.1.iter().map(|(_, l)| l).sum::<usize>());
                res
            })
    }
}

use bincode::{Decode, Encode};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

/// A region on an antibody
#[expect(missing_docs)]
#[derive(Clone, Debug, Decode, Deserialize, Encode, Eq, Hash, PartialEq, Serialize)]
pub enum Region {
    Framework(usize),
    ComplementarityDetermining(usize),
    Hinge(Option<usize>),
    ConstantHeavy(usize),
    ConstantLight,
    SecratoryTail,
    MembraneTail(Option<usize>),
    Other(String),
    /// When multiple regions are joined, or when the exact boundary is not known
    Joined(Vec<Region>),
    None,
}

impl crate::space::Space for Region {
    fn space(&self) -> UsedSpace {
        match self {
            Self::Framework(d) => d.space(),
            Self::ComplementarityDetermining(d) => d.space(),
            Self::Hinge(d) => d.space(),
            Self::ConstantHeavy(d) => d.space(),
            Self::ConstantLight | Self::None | Self::SecratoryTail => UsedSpace::default(),
            Self::MembraneTail(d) => d.space(),
            Self::Other(d) => d.space(),
            Self::Joined(d) => d.space(),
        }
        .set_total::<Self>()
    }
}

/// A sequence annotation
#[derive(Clone, Debug, Decode, Deserialize, Encode, Eq, Hash, PartialEq, Serialize)]
pub enum Annotation {
    /// A conserved residue
    Conserved,
    /// A potential N linked glycan position
    NGlycan,
    /// Any other annotation
    Other(String),
}

impl crate::space::Space for Annotation {
    fn space(&self) -> UsedSpace {
        match self {
            Self::Conserved | Self::NGlycan => UsedSpace::default(),
            Self::Other(d) => d.space(),
        }
        .set_total::<Self>()
    }
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Framework(n) => write!(f, "FR{n}"),
            Self::ComplementarityDetermining(n) => write!(f, "CDR{n}"),
            Self::Hinge(n) => write!(f, "H{}", n.map_or(String::new(), |n| n.to_string())),
            Self::ConstantHeavy(n) => write!(f, "CH{n}"),
            Self::ConstantLight => write!(f, "CL"),
            Self::SecratoryTail => write!(f, "CHS"),
            Self::MembraneTail(n) => write!(f, "M{}", n.map_or(String::new(), |n| n.to_string())),
            Self::Other(o) => write!(f, "{o}"),
            Self::Joined(o) => write!(f, "{}", o.iter().join("-")),
            Self::None => Ok(()),
        }
    }
}

impl std::str::FromStr for Region {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "" => Self::None,
            "CL" => Self::ConstantLight,
            "CHS" => Self::SecratoryTail,
            "H" => Self::Hinge(None),
            "M" => Self::MembraneTail(None),
            cdr if cdr.starts_with("CDR") => cdr[3..].parse::<usize>().map_or_else(
                |_| Self::Other(cdr.to_string()),
                Self::ComplementarityDetermining,
            ),
            fr if fr.starts_with("FR") => fr[2..]
                .parse::<usize>()
                .map_or_else(|_| Self::Other(fr.to_string()), Self::Framework),
            ch if ch.starts_with("CH") => ch[2..]
                .parse::<usize>()
                .map_or_else(|_| Self::Other(ch.to_string()), Self::ConstantHeavy),
            h if h.starts_with('H') => h[1..]
                .parse::<usize>()
                .map_or_else(|_| Self::Other(h.to_string()), |c| Self::Hinge(Some(c))),
            m if m.starts_with('M') => m[1..].parse::<usize>().map_or_else(
                |_| Self::Other(m.to_string()),
                |c| Self::MembraneTail(Some(c)),
            ),
            o => Self::Other(o.to_string()),
        })
    }
}

impl std::str::FromStr for Annotation {
    type Err = ();
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Ok(match s {
            "C" | "Conserved" => Self::Conserved,
            "N" | "NGlycan" => Self::NGlycan,
            o => Self::Other(o.to_string()),
        })
    }
}

impl std::fmt::Display for Annotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Conserved => write!(f, "Conserved"),
            Self::NGlycan => write!(f, "NGlycan"),
            Self::Other(o) => write!(f, "{o}"),
        }
    }
}
