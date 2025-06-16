use std::{io::Write, path::Path};

use bincode::config::Configuration;
use rustyms::{
    chemistry::{Element, ElementalData},
    identification::csv::parse_csv,
    system::da,
};

pub(crate) fn build_atomic_masses(out_dir: &Path) {
    let mut atomic_weights = vec![None; 118];
    let mut isotopic_abundances = vec![Vec::new(); 118];
    let mut atomic_masses = vec![Vec::new(); 118];

    let table = parse_csv(
        "rustyms-generate-databases/data/IUPAC-atomic-masses.csv",
        b',',
        None,
    )
    .unwrap();
    for line in table {
        let line = line.unwrap();
        let (nuclide, mass, _uncertainty, year) = (&line[0], &line[1], &line[2], &line[3]);
        if nuclide.starts_with("AME")
            || nuclide.is_empty()
            || nuclide == "nuclide"
            || !year.ends_with("2020</a>")
        {
            continue;
        }
        let isotope = nuclide
            .trim_end_matches(|c: char| c.is_alphabetic())
            .parse::<usize>()
            .map_err(|e| e.to_string())
            .unwrap();
        let element = Element::try_from(nuclide.trim_start_matches(|c: char| c.is_ascii_digit()))
            .map_err(|()| {
                format!("Not a valid isotope+element, could not recognise element: {nuclide}")
            })
            .unwrap();
        let mass = mass
            .parse::<f64>()
            .map_err(|e| format!("{e}@{mass}"))
            .unwrap();
        atomic_masses[element as usize - 1].push((isotope, mass));
    }

    let mut last_element = 0;
    let table = parse_csv(
        "rustyms-generate-databases/data/CIAAW-isotopic-abundances.csv",
        b',',
        Some(vec![
            "z".to_string(),
            "symbol".to_string(),
            "name".to_string(),
            "isotope".to_string(),
            "abundance".to_string(),
            "note".to_string(),
        ]),
    )
    .unwrap();
    for line in table {
        let line = line.unwrap();
        let (element, _element, _name, isotope, abundance, _note) =
            (&line[0], &line[1], &line[2], &line[3], &line[4], &line[5]);
        let mut abundance = abundance.to_owned();
        abundance.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || isotope.starts_with('[') || abundance == "-" {
            continue;
        }

        if !element.is_empty() {
            last_element = element
                .parse::<usize>()
                .map_err(|_| format!("Not a valid number for element Z: {element}"))
                .unwrap();
        }
        let element = last_element;

        let isotope = isotope.parse::<usize>().unwrap();

        isotopic_abundances[element - 1].push((isotope, get_ciaaw_number(&abundance).unwrap()));
    }

    let table = parse_csv(
        "rustyms-generate-databases/data/CIAAW-atomic-weights.csv",
        b',',
        Some(vec![
            "z".to_string(),
            "symbol".to_string(),
            "name".to_string(),
            "weight".to_string(),
            "note".to_string(),
        ]),
    )
    .unwrap();
    for line in table {
        let line = line.unwrap();
        let (element, weight) = (&line[0], &line[3]);
        let mut weight = weight.to_owned();
        weight.retain(|c| !c.is_whitespace()); // Remove any whitespace, including any sneaking non breaking spaces
        if element == "Z" || weight == "—" {
            continue;
        }
        let element = element
            .parse::<usize>()
            .map_err(|_| format!("Not valid element number (not a number): {element}"))
            .unwrap();
        atomic_weights[element - 1] = Some(get_ciaaw_number(&weight).unwrap());
    }

    // Monoisotopic, average weight, isotopes: (num, mass, abundance)
    let combined_data = isotopic_abundances
        .into_iter()
        .zip(atomic_masses)
        .zip(atomic_weights)
        .map(|((isotopic_abundance, isotopic_mass), atomic_weight)| {
            let isotopes = isotopic_mass
                .iter()
                .map(|(i, m)| {
                    (
                        *i,
                        *m,
                        isotopic_abundance
                            .iter()
                            .find(|(ai, _)| ai == i)
                            .map_or(0.0, |(_, a)| *a),
                    )
                })
                .collect::<Vec<_>>();
            (
                isotopes
                    .iter()
                    .fold(
                        (None, 0.0),
                        |acc, (_, m, a)| if *a > acc.1 { (Some(*m), *a) } else { acc },
                    )
                    .0,
                atomic_weight,
                isotopes,
            )
        });

    // Write out the data
    let dest_path = Path::new(&out_dir).join("elements.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let elements = combined_data
        .into_iter()
        .map(|(m, a, i)| {
            (
                m.map(da),
                a.map(da),
                i.into_iter()
                    .map(|(n, m, i)| (n as u16, da(m), i))
                    .collect(),
            )
        })
        .collect();
    file.write_all(
        &bincode::serde::encode_to_vec::<ElementalData, Configuration>(
            elements,
            Configuration::default(),
        )
        .unwrap(),
    )
    .unwrap();
}

fn get_ciaaw_number(text: &str) -> Result<f64, String> {
    let parse = |t: &str| {
        t.parse::<f64>()
            .map_err(|_| format!("Not a valid number: {t}"))
    };
    if text.starts_with('[') {
        let (low, high) = text[1..text.len() - 1]
            .split_once(',')
            .ok_or_else(|| format!("Not a valid range: {text}"))?;
        Ok((parse(low)? + parse(high)?) / 2.0)
    } else if text.ends_with(')') {
        Ok(parse(
            text.split_once('(')
                .ok_or_else(|| format!("Not valid error indication: {text}"))?
                .0,
        )?)
    } else {
        Ok(parse(text)?)
    }
}
