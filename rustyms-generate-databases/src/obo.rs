use std::{collections::HashMap, fs::File, io::BufReader, path::Path};

use flate2::bufread::GzDecoder;

#[derive(Clone, Debug, Default)]
pub(crate) struct OboOntology {
    pub headers: Vec<(String, String)>,
    pub objects: Vec<OboObject>,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct OboObject {
    pub name: String,
    pub lines: HashMap<String, Vec<String>>,
    pub property_values: HashMap<String, Vec<OboValue>>,
}

#[derive(Clone, Debug)]
pub(crate) enum OboValue {
    String(String),
    Float(f64),
    Integer(isize),
    Boolean(bool),
}

impl std::fmt::Display for OboValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::String(s) => write!(f, "{s}"),
            Self::Float(s) => write!(f, "{s}"),
            Self::Integer(s) => write!(f, "{s}"),
            Self::Boolean(s) => write!(f, "{s}"),
        }
    }
}

impl OboOntology {
    pub(crate) fn from_file(path: impl AsRef<Path>) -> Result<OboOntology, String> {
        let file = File::open(path.as_ref()).map_err(|e| e.to_string())?;
        if path
            .as_ref()
            .extension()
            .is_some_and(|ext| ext.eq_ignore_ascii_case("gz"))
        {
            Self::from_raw(BufReader::new(GzDecoder::new(BufReader::new(file))))
        } else {
            Self::from_raw(BufReader::new(file))
        }
    }

    fn from_raw<T: std::io::BufRead>(reader: T) -> Result<Self, String> {
        let mut obo = Self::default();
        let mut recent_obj = None;

        for line in reader.lines() {
            let line = line.map_err(|e| e.to_string())?.trim_end().to_string();
            if line.is_empty() {
                continue;
            }
            if line.starts_with('[') && line.ends_with(']') {
                if let Some(obj) = recent_obj {
                    obo.objects.push(obj);
                }
                recent_obj = Some(OboObject::new(&line[1..=line.len() - 2]));
            } else if let Some((id, value_line)) = line.split_once(':') {
                if let Some(obj) = &mut recent_obj {
                    if id == "property_value" {
                        let value_line = value_line.trim();
                        let first_space = value_line
                            .char_indices()
                            .find_map(|(i, c)| (c == ' ').then_some(i))
                            .unwrap();
                        let last_space = value_line
                            .char_indices()
                            .rfind(|(_, c)| *c == ' ')
                            .map(|(i, _)| i)
                            .unwrap();
                        if first_space == last_space {
                            let name = value_line[..first_space].trim();
                            let value = value_line[first_space..].trim().trim_matches('"');
                            obj.property_values
                                .entry(name.to_string())
                                .or_insert(Vec::new())
                                .push(OboValue::String(value.to_string()));
                        } else {
                            let name = value_line[..first_space].trim().trim_end_matches(':');
                            let value =
                                value_line[first_space..last_space].trim().trim_matches('"');
                            let unit = value_line[last_space..].trim();
                            obj.property_values
                                .entry(name.to_string())
                                .or_insert(Vec::new())
                                .push(match unit {
                                    "xsd:string" => OboValue::String(value.to_string()),
                                    "xsd:double" | "xsd:float" => {
                                        if !value.starts_with('-') && value.contains('-') {
                                            // Some ontologies use a range
                                            OboValue::String(value.to_string())
                                        } else {
                                            OboValue::Float(value.parse().unwrap_or_else(|err| {
                                                panic!(
                                            "Invalid float: '{value}' ({err}) for line: '{line}'"
                                        )
                                            }))
                                        }
                                    }
                                    "xsd:boolean" => {
                                        OboValue::Boolean(value == "true" || value == "1")
                                    }
                                    "xsd:integer"
                                    | "xsd:nonNegativeInteger"
                                    | "xsd:positiveInteger" => {
                                        OboValue::Integer(value.parse().unwrap_or_else(|err| {
                                            panic!(
                                            "Invalid integer: '{value}' ({err}) for line: '{line}'"
                                        )
                                        }))
                                    }
                                    dt => {
                                        unreachable!("Undefined datatype: '{dt}' in line: '{line}'")
                                    }
                                });
                        }
                    } else {
                        obj.lines
                            .entry(id.trim().to_string())
                            .or_insert(Vec::new())
                            .push(value_line.trim().to_string());
                    }
                } else {
                    obo.headers.push((id.to_string(), value_line.to_string()));
                }
            } else {
                return Err(format!("Invalid line in obo file: {line}"));
            }
        }
        if let Some(obj) = recent_obj {
            obo.objects.push(obj);
        }
        Ok(obo)
    }
}

impl OboObject {
    pub(crate) fn new(name: &str) -> Self {
        Self {
            name: name.to_string(),
            ..Self::default()
        }
    }
}
