#![allow(dead_code)]

///
/// Some parts of this file originates from [rustyms](https://github.com/snijderlab/rustyms)
/// Copyright (c) 2023 Douwe Schulte and contributors
/// MIT License
///

use std::str::FromStr;
use anyhow::*;
use crate::chemistry::composition::{ElementalComposition, ElementCount};
use crate::chemistry::element::Element;
use crate::chemistry::glycan::{GlycanComposition, MonoSaccharide};

enum Brick {
    Element(Element),
    Formula(ElementalComposition),
    MonoSaccharide(MonoSaccharide),
}

pub fn parse_unimod_composition(composition: &str) -> Result<(ElementalComposition, GlycanComposition)> {
    let mut elem_comp = ElementalComposition::new(&[]);
    let mut monosaccharides: GlycanComposition = Vec::new();

    let mut last_name = String::new();
    let mut last_number = String::new();
    let mut minus = 1;
    for c in composition.bytes() {
        match c {
            b'-' => minus = -1,
            b'(' => (),
            b')' => {
                let parsed_number = last_number.parse::<i16>()?;
                let num = parsed_number * minus;
                match parse_unimod_composition_brick(last_name.as_str())? {
                    Brick::Formula(f) => elem_comp += &(f * num),
                    Brick::Element(e) => elem_comp.add(ElementCount::new_monoisotope(e, num as f32)),
                    Brick::MonoSaccharide(m) => monosaccharides.push((m, num)),
                }
                last_name.clear();
                last_number.clear();
                minus = 1;
            }
            b' ' => {
                if !last_name.is_empty() {
                    match parse_unimod_composition_brick(last_name.as_str())? {
                        Brick::Formula(f) => elem_comp += &f,
                        Brick::Element(e) => elem_comp.add(ElementCount::new_monoisotope(e, 1f32)),
                        Brick::MonoSaccharide(m) => monosaccharides.push((m, 1)),
                    }
                    last_name.clear();
                }
            }
            n if n.is_ascii_digit() => last_number.push(n as char),
            n if n.is_ascii_alphabetic() => last_name.push(n as char),
            _ => panic!("Weird formula composition: {composition}"),
        }
    }
    if !last_name.is_empty() {
        match parse_unimod_composition_brick(last_name.as_str())? {
            Brick::Formula(f) => elem_comp += &f,
            Brick::Element(e) => elem_comp.add(ElementCount::new_monoisotope(e, 1f32)),
            Brick::MonoSaccharide(m) => monosaccharides.push((m, 1)),
        }
    }
    Ok((elem_comp, monosaccharides))
}


fn parse_unimod_composition_brick(name: &str) -> Result<Brick> {
    match name.to_lowercase().as_str() {
        "ac" => Ok(Brick::Formula(ElementalComposition::from_monoisotope_tuples(&[
            (Element::O, 1),
            (Element::C, 2),
            (Element::H, 2),
        ]))),
        "me" => Ok(Brick::Formula(ElementalComposition::from_monoisotope_tuples(&[
            (Element::C, 1),
            (Element::H, 2),
        ]))),
        "kdn" => Ok(Brick::Formula(ElementalComposition::from_monoisotope_tuples(&[
            (Element::C, 9),
            (Element::H, 14),
            (Element::O, 8),
        ]))),
        "kdo" => Ok(Brick::Formula(ElementalComposition::from_monoisotope_tuples(&[
            (Element::C, 8),
            (Element::H, 12),
            (Element::O, 7),
        ]))),
        "sulf" => Ok(Brick::Formula(ElementalComposition::from_monoisotope_tuples(&[(Element::S, 1)]))),
        _ => {
            // Try some fallbacks
            Element::from_str(name)
                .map(|el| Brick::Element(el))
                .or_else(|_e| MonoSaccharide::from_short_iupac(name, 0, 0).map(|(ms, _)| Brick::MonoSaccharide(ms)))
                .map_err(|e| Error::new(e).context(format!("Could not parse unimod brick: `{}`", name)))
        }
    }
}
