use std::num::NonZeroU16;

use rand::distr::{Distribution, StandardUniform};

use crate::{
    glycan::{BaseSugar, GlycanStructure, GlycanSubstituent, MonoSaccharide},
    modification::SimpleModificationInner,
    system::{dalton, Mass, OrderedMass},
    Element, MolecularFormula,
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
        match rng.random_range(1..=44) {
            1 => GlycanSubstituent::Acetimidoyl,
            2 => GlycanSubstituent::Acetyl,
            3 => GlycanSubstituent::AcetylAlanyl,
            4 => GlycanSubstituent::AcetylGlutaminyl,
            5 => GlycanSubstituent::Acid,
            6 => GlycanSubstituent::Alanyl,
            7 => GlycanSubstituent::Alcohol,
            8 => GlycanSubstituent::Amino,
            9 => GlycanSubstituent::Aric,
            10 => GlycanSubstituent::CargoxyEthylidene,
            11 => GlycanSubstituent::Deoxy,
            12 => GlycanSubstituent::Didehydro,
            13 => GlycanSubstituent::DiHydroxyButyryl,
            14 => GlycanSubstituent::DiMethyl,
            15 => GlycanSubstituent::DiMethylAcetimidoyl,
            16 => GlycanSubstituent::DiMethylGlyceryl,
            17 => GlycanSubstituent::Element(rng.random()),
            18 => GlycanSubstituent::Ethanolamine,
            19 => GlycanSubstituent::EtOH,
            20 => GlycanSubstituent::Formyl,
            21 => GlycanSubstituent::Glyceryl,
            22 => GlycanSubstituent::Glycolyl,
            23 => GlycanSubstituent::Glycyl,
            24 => GlycanSubstituent::HydroxyButyryl,
            25 => GlycanSubstituent::HydroxyMethyl,
            26 => GlycanSubstituent::Lac,
            27 => GlycanSubstituent::Lactyl,
            28 => GlycanSubstituent::Methyl,
            29 => GlycanSubstituent::MethylAcetimidoyl,
            30 => GlycanSubstituent::MethylGlutamyl,
            31 => GlycanSubstituent::NAcetyl,
            32 => GlycanSubstituent::NDiMe,
            33 => GlycanSubstituent::NFo,
            34 => GlycanSubstituent::NGlycolyl,
            35 => GlycanSubstituent::OCarboxyEthyl,
            36 => GlycanSubstituent::PCholine,
            37 => GlycanSubstituent::Phosphate,
            38 => GlycanSubstituent::Pyruvyl,
            39 => GlycanSubstituent::Suc,
            40 => GlycanSubstituent::Sulfate,
            41 => GlycanSubstituent::Tauryl,
            42 => GlycanSubstituent::Ulo,
            43 => GlycanSubstituent::Ulof,
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
