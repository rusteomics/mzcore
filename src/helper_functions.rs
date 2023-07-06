use crate::{element, formula::MolecularFormula, ELEMENT_PARSE_LIST};
pub trait ResultExtensions<T, E> {
    fn flat_err(self) -> Result<T, E>;
}

impl<T, E> ResultExtensions<T, E> for Result<T, Result<T, E>> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => Ok(o),
            Err(r) => r,
        }
    }
}

impl<T, E> ResultExtensions<T, E> for Result<Result<T, E>, E> {
    fn flat_err(self) -> Result<T, E> {
        match self {
            Ok(o) => o,
            Err(r) => Err(r),
        }
    }
}

pub fn parse_named_counter<T: Copy>(
    value: &str,
    names: &[(&str, T)],
    allow_negative: bool,
) -> Result<Vec<(T, isize)>, String> {
    let mut index = 0;
    let mut output = Vec::new();
    while index < value.len() {
        if value[index..].starts_with(' ') {
            index += 1;
        } else {
            let mut found = false;
            for name in names {
                if value[index..].starts_with(name.0) {
                    index += name.0.len();
                    let num = &value[index..]
                        .chars()
                        .skip_while(char::is_ascii_whitespace)
                        .take_while(|c| c.is_ascii_digit() || (allow_negative && *c == '-'))
                        .collect::<String>()
                        .trim()
                        .to_string();
                    if num.is_empty() {
                        output.push((name.1, 1));
                    } else {
                        output.push((name.1, num.parse().unwrap()));
                        index += num.len()
                            + value[index..]
                                .chars()
                                .take_while(char::is_ascii_whitespace)
                                .count();
                    }
                    found = true;
                    break; // Names loop
                }
            }
            if !found {
                return Err(format!("Name not recognised {}", &value[index..]));
            }
        }
    }
    Ok(output)
}

// PSI-MOD: (12)C -5 (13)C 5 H 0 N 0 O 0 S 0
pub fn parse_molecular_formula_psi_mod(value: &str) -> Result<MolecularFormula, String> {
    let mut index = 0;
    let mut isotope = 0;
    let mut element = None;
    let bytes = value.as_bytes();
    let mut result = MolecularFormula::default();
    while index < value.len() {
        match bytes[index] {
            b'(' if isotope == 0 => {
                let len = bytes
                    .iter()
                    .skip(index)
                    .position(|c| *c == b')')
                    .ok_or(format!(
                        "No closing round bracket for round bracket at index: {index}"
                    ))?;
                isotope = value[index + 1..index + len]
                    .parse::<u16>()
                    .map_err(|e| e.to_string())?;
                index += len + 1;
            }
            b'-' | b'0'..=b'9' if element.is_some() => {
                let (num, len) = std::str::from_utf8(
                    &bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit() || **c == b'-')
                        .copied()
                        .collect::<Vec<_>>(),
                )
                .map(|v| (v.parse::<i16>().map_err(|e| e.to_string()), v.len()))
                .map_err(|e| e.to_string())?;
                let num = num?;
                if num != 0 {
                    result.add((element.unwrap(), isotope, num));
                }
                element = None;
                isotope = 0;
                index += len;
            }
            b' ' => index += 1,
            _ => {
                let mut found = false;
                for possible in ELEMENT_PARSE_LIST {
                    if value[index..].starts_with(possible.0) {
                        element = Some(possible.1);
                        index += possible.0.len();
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(format!(
                        "Could not parse PSI-MOD elemental formula, broke down at index: {index}"
                    ));
                }
            }
        }
    }
    if isotope != 0 || element.is_some() {
        Err("Last element missed a count".to_string())
    } else {
        Ok(result)
    }
}

// ProForma: [13C2][12C-2]H2N
pub fn parse_molecular_formula_pro_forma(value: &str) -> Result<MolecularFormula, String> {
    let mut index = 0;
    let mut element = None;
    let bytes = value.as_bytes();
    let mut result = MolecularFormula::default();
    while index < value.len() {
        match bytes[index] {
            b'[' => {
                index += 1; // Skip the open square bracket
                let len = bytes
                    .iter()
                    .skip(index)
                    .position(|c| *c == b']')
                    .ok_or(format!(
                        "No closing round bracket for round bracket at index: {index}"
                    ))?;
                let isotope = bytes
                    .iter()
                    .skip(index)
                    .take_while(|c| c.is_ascii_digit())
                    .count();
                let ele = bytes
                    .iter()
                    .skip(index + isotope)
                    .take_while(|c| c.is_ascii_alphabetic())
                    .count();

                for possible in ELEMENT_PARSE_LIST {
                    if value[index + isotope..index + isotope + ele] == *possible.0 {
                        element = Some(possible.1);
                        break;
                    }
                }
                let num = value[index + isotope + ele..index + len]
                    .parse::<i16>()
                    .map_err(|e| e.to_string())?;
                let isotope = value[index..index + isotope]
                    .parse::<u16>()
                    .map_err(|e| e.to_string())?;

                result.add((element.unwrap(), isotope, num));
                element = None;
                index += len + 1;
            }
            b'-' | b'0'..=b'9' if element.is_some() => {
                let (num, len) = std::str::from_utf8(
                    &bytes
                        .iter()
                        .skip(index)
                        .take_while(|c| c.is_ascii_digit() || **c == b'-')
                        .copied()
                        .collect::<Vec<_>>(),
                )
                .map(|v| (v.parse::<i16>().map_err(|e| e.to_string()), v.len()))
                .map_err(|e| e.to_string())?;
                let num = num?;
                if num != 0 {
                    result.add((element.unwrap(), 0, num));
                }
                element = None;
                index += len;
            }
            b' ' => index += 1,
            _ => {
                let mut found = false;
                for possible in ELEMENT_PARSE_LIST {
                    if value[index..].starts_with(possible.0) {
                        element = Some(possible.1);
                        index += possible.0.len();
                        found = true;
                        break;
                    }
                }
                if !found {
                    return Err(format!(
                        "Could not parse ProForma elemental formula, broke down at index: {index}"
                    ));
                }
            }
        }
    }
    if let Some(element) = element {
        result.add((element, 0, 1));
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Element;
    use crate::MolecularFormula;

    #[test]
    fn simple_formulae() {
        assert_eq!(
            parse_molecular_formula_psi_mod("(12)C -5 (13)C 5 H 0 N 0 O 0 S 0").unwrap(),
            molecular_formula!((12)C -5 (13)C 5 H 0 N 0 O 0 S 0)
        );
        assert_eq!(
            parse_molecular_formula_pro_forma("[13C2][12C-2]H2N").unwrap(),
            molecular_formula!((12)C -2 (13)C 2 H 2 N 1)
        );
    }
}
