use std::num::NonZeroU16;

use rand::distr::{Distribution, StandardUniform};

use crate::{
    chemistry::{Element, MolecularFormula},
    glycan::{BaseSugar, GlycanStructure, GlycanSubstituent, MonoSaccharide},
    sequence::SimpleModificationInner,
    system::{Mass, OrderedMass, dalton},
};

impl Distribution<SimpleModificationInner> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> SimpleModificationInner {
        match rng.random_range(0..=3) {
            0 => SimpleModificationInner::Mass(rng.random()),
            1 => SimpleModificationInner::Formula(rng.random()),
            2 => {
                let mut glycans: Vec<(MonoSaccharide, i8)> = Vec::new();
                for _ in 0..rng.random_range(0..32) {
                    glycans.push(rng.random());
                }
                SimpleModificationInner::Glycan(
                    glycans
                        .into_iter()
                        .map(|(ms, number)| (ms, number as isize))
                        .collect(),
                )
            }
            3 => SimpleModificationInner::GlycanStructure(rng.random()),
            _ => todo!(),
        }
    }
}

impl Distribution<GlycanStructure> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> GlycanStructure {
        let mut branches = Vec::new();
        for _ in 0..rng.random_range(0..=2) {
            branches.push(rng.random());
        }
        GlycanStructure::new(rng.random(), branches)
    }
}

impl Distribution<MonoSaccharide> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> MonoSaccharide {
        let mut substituents = Vec::new();
        for _ in 0..rng.random_range(0..5) {
            substituents.push(rng.random());
        }
        MonoSaccharide::new(rng.random(), &substituents)
    }
}

impl Distribution<MolecularFormula> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> MolecularFormula {
        let mut formula = MolecularFormula::default();
        for _ in 0..rng.random_range(0..32) {
            let element: Element = rng.random();
            let isotope: u16 = rng.random();
            let isotope = if rng.random::<f32>() > 0.9 {
                None
            } else {
                NonZeroU16::new(isotope)
            };

            if element.is_valid(isotope) {
                let _ = formula.add((element, isotope, rng.random()));
            }
        }
        formula
    }
}

impl Distribution<GlycanSubstituent> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> GlycanSubstituent {
        match rng.random_range(1..=36) {
            1 => GlycanSubstituent::Acetimidoyl,
            2 => GlycanSubstituent::Acetyl,
            3 => GlycanSubstituent::Acid,
            4 => GlycanSubstituent::Alanyl,
            5 => GlycanSubstituent::Alcohol,
            6 => GlycanSubstituent::Amino,
            7 => GlycanSubstituent::Aric,
            8 => GlycanSubstituent::CargoxyEthylidene,
            9 => GlycanSubstituent::Deoxy,
            10 => GlycanSubstituent::Didehydro,
            11 => GlycanSubstituent::DiMethyl,
            12 => GlycanSubstituent::Element(rng.random()),
            13 => GlycanSubstituent::Ethanolamine,
            14 => GlycanSubstituent::EtOH,
            15 => GlycanSubstituent::Formyl,
            16 => GlycanSubstituent::Glyceryl,
            17 => GlycanSubstituent::Glycolyl,
            18 => GlycanSubstituent::Glycyl,
            19 => GlycanSubstituent::HydroxyButyryl,
            20 => GlycanSubstituent::Lac,
            21 => GlycanSubstituent::Lactyl,
            22 => GlycanSubstituent::Methyl,
            23 => GlycanSubstituent::NAcetyl,
            24 => GlycanSubstituent::NDiMe,
            25 => GlycanSubstituent::NFo,
            26 => GlycanSubstituent::NGlycolyl,
            27 => GlycanSubstituent::OCarboxyEthyl,
            28 => GlycanSubstituent::PCholine,
            29 => GlycanSubstituent::Phosphate,
            30 => GlycanSubstituent::Pyruvyl,
            31 => GlycanSubstituent::Suc,
            32 => GlycanSubstituent::Sulfate,
            33 => GlycanSubstituent::Tauryl,
            34 => GlycanSubstituent::Ulo,
            35 => GlycanSubstituent::Ulof,
            _ => GlycanSubstituent::Water,
        }
    }
}

impl Distribution<BaseSugar> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> BaseSugar {
        match rng.random_range(0..=9) {
            0 => BaseSugar::None,
            1 => BaseSugar::Sugar,
            2 => BaseSugar::Triose,
            3 => BaseSugar::Tetrose(None),
            4 => BaseSugar::Pentose(None),
            5 => BaseSugar::Hexose(None),
            6 => BaseSugar::Heptose(None),
            7 => BaseSugar::Octose,
            8 => BaseSugar::Nonose(None),
            _ => BaseSugar::Decose,
        }
    }
}

impl Distribution<Element> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> Element {
        Element::try_from(rng.random_range(Element::Electron as usize..Element::Og as usize))
            .unwrap()
    }
}

impl Distribution<OrderedMass> for StandardUniform {
    fn sample<R: rand::prelude::Rng + ?Sized>(&self, rng: &mut R) -> OrderedMass {
        Mass::new::<dalton>(rng.random_range(f64::MIN..f64::MAX)).into()
    }
}
