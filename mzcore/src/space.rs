use std::{num::NonZero, ops::{Add, AddAssign}};
use thin_vec::ThinVec;
use crate::{ glycan::{Configuration, GlycanSubstituent, HeptoseIsomer, HexoseIsomer, NonoseIsomer, PentoseIsomer, TetroseIsomer}, ontology::Ontology, prelude::{AminoAcid, Element}, sequence::{GnoSubsumption, Position, SimpleModificationInner}, system::OrderedMass};

pub(crate) trait Space {
    fn space(&self) -> UsedSpace;
}

#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct UsedSpace {
    pub(crate) stack: usize,
    pub(crate) padding: usize,
    pub(crate) heap_used: usize,
    pub(crate) heap_padding: usize,
    pub(crate) heap_unused: usize,
}

impl std::fmt::Display for UsedSpace {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result { 
        writeln!(f, "Stack: {}", display_bytes(self.stack))?;
        writeln!(f, "Padding: {}", display_bytes(self.padding))?;
        writeln!(f, "Heap")?;
        writeln!(f, "- Used: {}", display_bytes(self.heap_used))?;
        writeln!(f, "- Padding: {}", display_bytes(self.heap_padding))?;
        writeln!(f, "- Unused: {}", display_bytes(self.heap_unused))?;
        writeln!(f, "Total: {}", display_bytes(self.stack + self.padding + self.heap_used + self.heap_padding + self.heap_unused))
    }
}

fn display_bytes(bytes: usize) -> String {
    if bytes < 1024 {
        format!("{bytes} B")
    } else {
        const GB: usize = 1024 * 1024 * 1024;
        const MB: usize = 1024 * 1024;
        const KB: usize = 1024;
        let human = if bytes > GB {
            format!("{} GiB", bytes / GB)
        } else if bytes > MB {
            format!("{} MiB", bytes / MB)
        } else {
            format!("{} KiB", bytes / KB)
        };
        format!("{human} ({bytes} B)")
    }
}

impl UsedSpace {
    pub(crate) fn stack(stack: usize) -> Self {
        Self {
            stack,
            ..Default::default()
        }
    }

    pub(crate) fn set_total<Total>(self) -> Self {
        Self {padding: size_of::<Total>() - self.stack, ..self}
    }

    pub(crate) fn add_child(self, child: UsedSpace) -> Self {
        Self {
            stack: self.stack,
            padding: self.padding,
            heap_used: self.heap_used + child.heap_used + child.stack,
            heap_padding: self.heap_padding + child.heap_padding + child.padding,
            heap_unused: self.heap_unused + child.heap_unused
        }
    }
}

impl AddAssign for UsedSpace {
    fn add_assign(&mut self, rhs: Self) {
        self.stack += rhs.stack;
        self.padding += rhs.padding;
        self.heap_used += rhs.heap_used;
        self.heap_padding += rhs.heap_padding;
        self.heap_unused += rhs.heap_unused;
    }
}

impl Add for UsedSpace {
    type Output = Self;
    fn add(self, rhs: Self) -> Self{
        Self {
            stack: self.stack + rhs.stack,
            padding: self.padding + rhs.padding,
            heap_used: self.heap_used + rhs.heap_used,
            heap_padding: self.heap_padding + rhs.heap_padding,
            heap_unused: self.heap_unused + rhs.heap_unused
        }
    }
}

impl<T: Space> Space for Vec<T> {
    fn space(&self) -> UsedSpace {
        let mut total = UsedSpace::default();
        for e in self {
            total += e.space();
        }
        UsedSpace { 
            stack: 24, 
            padding: 0, 
            heap_used: total.stack + total.heap_used, 
            heap_padding: total.padding + total.heap_padding, 
            heap_unused: total.heap_unused + (self.capacity() - self.len()) * size_of::<T>() }
    }
}

impl<T: Space> Space for ThinVec<T> {
    fn space(&self) -> UsedSpace {
        let mut total = UsedSpace::default();
        for e in self {
            total += e.space();
        }
        UsedSpace { 
            stack: 8, 
            padding: 0, 
            heap_used: total.stack + total.heap_used, 
            heap_padding: total.padding + total.heap_padding, 
            heap_unused: total.heap_unused + (self.capacity() - self.len()) * size_of::<T>() }
    }
}

impl Space for String{
    fn space(&self) -> UsedSpace {
        UsedSpace { 
            stack: 24, 
            padding: 0, 
            heap_used: self.len(), 
            heap_padding: 0, 
            heap_unused: (self.capacity() - self.len()) }
    }
}

impl<T: Space> Space for Option<T> {
    fn space(&self) -> UsedSpace {
        if let Some(t) = self {
            t.space().set_total::<Self>()
        } else {
            UsedSpace::default().set_total::<Self>()
        }
    }
}

impl<T: Space> Space for std::sync::Arc<T> {
    fn space(&self) -> UsedSpace {
        UsedSpace::stack(size_of::<Self>()).add_child(self.as_ref().space())
    }
}

impl Space for Box<str> {
    fn space(&self) -> UsedSpace {
         UsedSpace {stack: 8, padding: 0, heap_used: self.len(), heap_padding: 0, heap_unused: 0}
    }
}

impl Space for SimpleModificationInner {
    fn space(&self) -> UsedSpace {
        match self {
            SimpleModificationInner::Mass(m) => m.space(),
            SimpleModificationInner::Formula(f) => f.space(),
            SimpleModificationInner::Glycan(g) => g.space(),
            SimpleModificationInner::GlycanStructure(g) => g.space(),
            SimpleModificationInner::Gno { composition, id, structure_score, subsumption_level, motif, taxonomy, glycomeatlas } => composition.space() + id.space() + structure_score.space() + subsumption_level.space() + motif.space() + taxonomy.space() + glycomeatlas.space(),
            SimpleModificationInner::Database { specificities, formula, id } => specificities.space() + formula.space() + id.space(),
            SimpleModificationInner::Linker { specificities, formula, id, length } => specificities.space() + formula.space() + id.space() + length.space(),
        }.set_total::<Self>()
    }
}

impl<A: Space> Space for (A,) {
    fn space(&self) -> UsedSpace {
        self.0.space()
    }
}

impl<A: Space, B: Space> Space for (A, B) {
    fn space(&self) -> UsedSpace {
        (self.0.space() + self.1.space()).set_total::<Self>()
    }
}

impl<A: Space, B: Space, C: Space> Space for (A, B, C) {
    fn space(&self) -> UsedSpace {
        (self.0.space() + self.1.space() + self.2.space()).set_total::<Self>()
    }
}

impl<A: Space, B: Space, C: Space, D: Space> Space for (A, B, C, D) {
    fn space(&self) -> UsedSpace {
        (self.0.space() + self.1.space() + self.2.space() + self.3.space()).set_total::<Self>()
    }
}

macro_rules! simple_space {
    ($ty:ty) => {
        impl Space for $ty {
            fn space(&self) -> UsedSpace {
                UsedSpace::stack(size_of::<Self>())
            }
        }
    };
}

simple_space!(u8);
simple_space!(u16);
simple_space!(u32);
simple_space!(u64);
simple_space!(u128);
simple_space!(usize);
simple_space!(i8);
simple_space!(i16);
simple_space!(i32);
simple_space!(i64);
simple_space!(i128);
simple_space!(isize);
simple_space!(NonZero<u8>);
simple_space!(NonZero<u16>);
simple_space!(NonZero<u32>);
simple_space!(NonZero<u64>);
simple_space!(NonZero<u128>);
simple_space!(NonZero<usize>);
simple_space!(NonZero<i8>);
simple_space!(NonZero<i16>);
simple_space!(NonZero<i32>);
simple_space!(NonZero<i64>);
simple_space!(NonZero<i128>);
simple_space!(NonZero<isize>);
simple_space!(f32);
simple_space!(f64);
simple_space!(ordered_float::OrderedFloat<f32>);
simple_space!(ordered_float::OrderedFloat<f64>);
simple_space!(bool);
simple_space!(Element);
simple_space!(AminoAcid);
simple_space!(GlycanSubstituent);
simple_space!(Configuration);
simple_space!(TetroseIsomer);
simple_space!(PentoseIsomer);
simple_space!(HexoseIsomer);
simple_space!(HeptoseIsomer);
simple_space!(NonoseIsomer);
simple_space!(Position);
simple_space!(Ontology);
simple_space!(GnoSubsumption);
simple_space!(OrderedMass);
simple_space!(mzcv::SynonymScope);

#[test]
fn gnome_space() {
    let space = crate::ontology::STATIC_ONTOLOGIES.gnome().data().fold(UsedSpace::default(), |acc, e| acc + e.space());
    let space_description = crate::ontology::STATIC_ONTOLOGIES.gnome().data().fold(UsedSpace::default(), |acc, e| if let Some(d) = e.description() {d.space() + acc} else {acc});
    let space_structures = crate::ontology::STATIC_ONTOLOGIES.gnome().data().fold(UsedSpace::default(), |acc, e| if let SimpleModificationInner::Gno { composition: crate::sequence::GnoComposition::Topology(t), .. } = &*e {t.space() + acc} else {acc} );
    println!("Total GNOme modifications: {}", crate::ontology::STATIC_ONTOLOGIES.gnome().data().len());
    println!("All used space");
    print!("{}", space);
    println!("All used space for ModificationIds");
    print!("{}", space_description);
    println!("All used space for GlycanStructures");
    print!("{}", space_structures);
    todo!();
}