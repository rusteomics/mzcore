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
                    Bound::Excluded(e) => e.saturating_add(1),
                    Bound::Included(i) => *i,
                    Bound::Unbounded => 0,
                };
                let end = match range.end_bound() {
                    Bound::Excluded(e) => *e,
                    Bound::Included(i) => i.saturating_add(1),
                    Bound::Unbounded => peptidoform.len(),
                };
                let mut index = 0;
                let mut regions = Vec::new();

                for (r, len) in self.regions() {
                    if index + len <= start {
                        // Before selection, so ignore
                    } else if index <= start && index + len > start {
                        // If the selections starts inside this region
                        if index + len >= end {
                            regions.push((r.clone(), end.saturating_sub(start)));
                            break;
                        }
                        regions.push((r.clone(), len - (start - index)));
                    } else if index + len < end {
                        regions.push((r.clone(), *len));
                    } else {
                        // If this extends beyond the selection
                        regions.push((r.clone(), end - index));
                        break;
                    }
                    index += len;
                }

                (
                    peptidoform,
                    regions,
                    self.annotations()
                        .iter()
                        .filter(|&(_, i)| range.contains(i))
                        .map(|(a, i)| (a.clone(), i - start))
                        .collect(),
                )
            })
    }
}

impl<T: AnnotatedPeptidoform + HasPeptidoformImpl> AnnotatedPeptidoform for &T {
    fn regions(&self) -> &[(Region, usize)] {
        (*self).regions()
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        (*self).annotations()
    }
}

impl<T: AnnotatedPeptidoform + HasPeptidoformImpl> AnnotatedPeptidoform for std::sync::Arc<T> {
    fn regions(&self) -> &[(Region, usize)] {
        self.as_ref().regions()
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        self.as_ref().annotations()
    }
}

impl<T: AnnotatedPeptidoform + HasPeptidoformImpl> AnnotatedPeptidoform for std::rc::Rc<T> {
    fn regions(&self) -> &[(Region, usize)] {
        self.as_ref().regions()
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        self.as_ref().annotations()
    }
}

impl<T: AnnotatedPeptidoform + HasPeptidoformImpl> AnnotatedPeptidoform for Box<T> {
    fn regions(&self) -> &[(Region, usize)] {
        self.as_ref().regions()
    }
    fn annotations(&self) -> &[(Annotation, usize)] {
        self.as_ref().annotations()
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
    Joined(Vec<Self>),
    None,
}

impl crate::space::Space for Region {
    fn space(&self) -> UsedSpace {
        match self {
            Self::Framework(d) | Self::ComplementarityDetermining(d) | Self::ConstantHeavy(d) => {
                d.space()
            }
            Self::ConstantLight | Self::None | Self::SecratoryTail => UsedSpace::default(),
            Self::Hinge(d) | Self::MembraneTail(d) => d.space(),
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

#[cfg(test)]
mod tests {
    use crate::{
        ontology::STATIC_ONTOLOGIES,
        prelude::{HasPeptidoformImpl, Peptidoform},
        sequence::{AnnotatedPeptidoform, Annotation, Linear, Region},
    };

    struct AP {
        peptidoform: Peptidoform<Linear>,
        regions: Vec<(Region, usize)>,
        annotations: Vec<(Annotation, usize)>,
    }
    impl AnnotatedPeptidoform for AP {
        fn regions(&self) -> &[(Region, usize)] {
            &self.regions
        }
        fn annotations(&self) -> &[(Annotation, usize)] {
            &self.annotations
        }
    }

    impl HasPeptidoformImpl for AP {
        type Complexity = Linear;
        fn peptidoform(&self) -> &Peptidoform<Self::Complexity> {
            &self.peptidoform
        }
    }

    #[test]
    fn annotated_sub_peptidoform() {
        let seq = AP {
            peptidoform: Peptidoform::pro_forma("PEPTIDE", &STATIC_ONTOLOGIES)
                .unwrap()
                .0
                .into_linear()
                .unwrap(),
            regions: vec![
                (Region::Framework(1), 4),
                (Region::ComplementarityDetermining(1), 1),
                (Region::Framework(2), 2),
            ],
            annotations: vec![
                (Annotation::Conserved, 4),
                (Annotation::NGlycan, 4),
                (Annotation::Conserved, 6),
            ],
        };

        let sub = seq.sub_peptidoform(..4).unwrap();
        dbg!(&sub);
        assert_eq!(sub.0.len(), 4);
        assert_eq!(sub.1.len(), 1);
        assert_eq!(sub.1.iter().map(|(_, l)| l).sum::<usize>(), 4);
        assert_eq!(sub.2.len(), 0);

        let sub = seq.sub_peptidoform(3..6).unwrap();
        dbg!(&sub);
        assert_eq!(sub.0.len(), 3);
        assert_eq!(sub.1.len(), 3);
        assert_eq!(sub.1.iter().map(|(_, l)| l).sum::<usize>(), 3);
        assert_eq!(sub.2.len(), 2);

        let sub = seq.sub_peptidoform(5..7).unwrap();
        dbg!(&sub);
        assert_eq!(sub.0.len(), 2);
        assert_eq!(sub.1.len(), 1);
        assert_eq!(sub.1.iter().map(|(_, l)| l).sum::<usize>(), 2);
        assert_eq!(sub.2.len(), 1);

        let sub = seq.sub_peptidoform(1..3).unwrap();
        dbg!(&sub);
        assert_eq!(sub.0.len(), 2);
        assert_eq!(sub.1.len(), 1);
        assert_eq!(sub.1.iter().map(|(_, l)| l).sum::<usize>(), 2);
        assert_eq!(sub.2.len(), 0);

        let sub = seq.sub_peptidoform(4..5).unwrap();
        dbg!(&sub);
        assert_eq!(sub.0.len(), 1);
        assert_eq!(sub.1.len(), 1);
        assert_eq!(sub.1.iter().map(|(_, l)| l).sum::<usize>(), 1);
        assert_eq!(sub.2.len(), 2);
    }
}
