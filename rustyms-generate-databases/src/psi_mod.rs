use std::{io::Write, path::Path};

use bincode::config::Configuration;
use rustyms::{
    chemistry::MolecularFormula,
    ontology::{Ontology, OntologyModificationList},
    sequence::{PlacementRule, Position},
};

use super::{obo::OboOntology, ontology_modification::OntologyModification, ModData};

pub(crate) fn build_psi_mod_ontology(out_dir: &Path) {
    let mods = parse_psi_mod();

    let dest_path = Path::new(&out_dir).join("psimod.dat");
    let mut file = std::fs::File::create(dest_path).unwrap();
    let final_mods = mods.into_iter().map(|m| m.into_mod()).collect::<Vec<_>>();
    println!("Found {} PSI-MOD modifications", final_mods.len());
    file.write_all(
        &bincode::serde::encode_to_vec::<OntologyModificationList, Configuration>(
            final_mods,
            Configuration::default(),
        )
        .unwrap(),
    )
    .unwrap();
}

fn parse_psi_mod() -> Vec<OntologyModification> {
    let obo = OboOntology::from_file("rustyms-generate-databases/data/PSI-MOD-newstyle.obo")
        .expect("Not a valid obo file");
    let mut mods = Vec::new();

    for obj in obo.objects {
        if obj.name != "Term" {
            continue;
        }
        let mut modification = OntologyModification {
            id: obj.lines["id"][0]
                .split_once(':')
                .expect("Incorrect psi mod id, should contain a colon")
                .1
                .parse()
                .expect("Incorrect psi mod id, should be numerical"),
            name: obj.lines["name"][0].to_string(),
            ontology: Ontology::Psimod,
            ..OntologyModification::default()
        };
        if let Some(values) = obj.lines.get("def") {
            assert!(values.len() == 1);
            let line = values[0][1..].split_once('\"').unwrap();
            modification.description = line.0.to_string();
            let ids = line.1.trim();
            modification.cross_ids = ids[1..ids.len() - 1]
                .split(',')
                .filter(|s| !s.is_empty())
                .map(|id| id.trim().split_once(':').unwrap())
                .map(|(r, i)| (r.to_string(), i.to_string()))
                .collect();
        }
        if let Some(values) = obj.lines.get("synonym") {
            for line in values {
                let line = line[1..].split_once('\"').unwrap();
                modification.synonyms.push(line.0.to_string());
            }
        }

        let mut rules = Vec::new();
        let mut origins = Vec::new();
        let mut term = None;
        for (id, value) in obj.property_values {
            if id == "DiffFormula" {
                modification.formula =
                    MolecularFormula::from_psi_mod(&value[0].to_string(), ..).unwrap();
            } else if id == "Origin" {
                origins = value[0]
                    .to_string()
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .collect();
            } else if id == "TermSpec" {
                if value[0].to_string() == "N-term" {
                    term = Some(Position::AnyNTerm);
                } else if value[0].to_string() == "C-term" {
                    term = Some(Position::AnyCTerm);
                } else {
                    panic!("Invalid TermSpec: {}", value[0])
                }
            }
        }
        // If the list of possible origins contains "X" than the mod can be placed on any aminoacid
        // But if there is a TermSpec definition that should still be accounted for
        let all_aminoacids = origins.contains(&"X".to_string());
        if !all_aminoacids {
            for origin in &origins {
                if origin.len() == 1 {
                    rules.push((
                        vec![PlacementRule::AminoAcid(
                            vec![origin.try_into().unwrap()],
                            term.unwrap_or(Position::Anywhere),
                        )],
                        Vec::new(),
                        Vec::new(),
                    ));
                } else {
                    rules.push((
                        vec![PlacementRule::PsiModification(
                            origin
                                .split_once(':')
                                .expect("Incorrect psi mod id, should contain a colon")
                                .1
                                .parse()
                                .expect("Incorrect psi mod id, should be numerical"),
                            term.unwrap_or(Position::Anywhere),
                        )],
                        Vec::new(),
                        Vec::new(),
                    ));
                }
            }
        }
        if origins.is_empty() || all_aminoacids {
            if let Some(term) = term {
                rules.push((vec![PlacementRule::Terminal(term)], Vec::new(), Vec::new()));
            }
        }
        modification.data = ModData::Mod {
            specificities: rules,
        };
        mods.push(modification);
    }

    mods
}

#[cfg(test)]
mod tests {
    use rustyms::{chemistry::MolecularFormula, molecular_formula};

    #[test]
    fn parse_molecular_formula() {
        assert_eq!(
            MolecularFormula::from_psi_mod("(12)C -5 (13)C 5 H 0 N 0 O 0 S 0", ..).unwrap(),
            molecular_formula!([12 C -5] [13 C 5] H 0 N 0 O 0 S 0)
        );
        assert_eq!(
            MolecularFormula::from_psi_mod("(12)C -9 (13)C 9", ..).unwrap(),
            molecular_formula!([12 C -9] [13 C 9])
        );
    }
}
