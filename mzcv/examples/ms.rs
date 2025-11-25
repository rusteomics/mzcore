//! An example to show how to work with the MS ontology
use std::{borrow::Cow, io::Write, sync::Arc};

use bincode::{Encode, Decode};
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVData, CVError, CVFile, CVIndex, CVSource, CVVersion, ControlledVocabulary,
    HashBufReader, OboOntology, OboStanzaType, SynonymScope,
};

fn main() {
    let (mut ms, _) = CVIndex::<MS>::init();
    if let Err(e) = ms.update_from_path(None, true) {
        eprintln!("Failed to load from cache {e}. Attempting to download from the internet");
        ms.update_from_url(&[]).unwrap();
    }
    println!("Stored at: {}", MS::default_stem().display());
    println!("Found {} MS items", ms.data().len());

    println!(
        "Date: {} Version: {} Hash: {}",
        ms.version().last_updated.map_or_else(
            || "?".to_string(),
            |(y, m, d, h, min)| format!("{y}-{m}-{d} {h:02}:{min:02}")
        ),
        ms.version()
            .version
            .as_ref()
            .map_or_else(|| "?".to_string(), ToString::to_string),
        ms.version().hash_hex(),
    );

    let item = ms.get_by_index(&1_000_016).unwrap();
    assert_eq!(item.name, "scan start time".into());

    println!("Search for terms in the MS ontology (prefix with '==' to search for exact names)");

    loop {
        print!("?:");
        std::io::stdout().flush().unwrap();
        let mut term = String::new();
        std::io::stdin().read_line(&mut term).unwrap();
        let term = term.trim();
        let answers = term.strip_prefix("==").map_or_else(
            || {
                ms.search(term, 10, 6)
                    .into_iter()
                    .map(|(a, _, _)| a)
                    .collect()
            },
            |term| ms.get_by_name(term).map(|v| vec![v]).unwrap_or_default(),
        );
        if answers.is_empty() {
            println!("No matches found");
        } else {
            for a in answers {
                println!(
                    " >{} '{}': {}",
                    a.curie()
                        .map_or_else(|| "-".to_string(), |i| i.to_string()),
                    a.name,
                    a.definition
                );
                for i in a.parents() {
                    let parent = ms.get_by_index(&i).unwrap();
                    println!(
                        "\t{}|{}", parent.curie().unwrap(), parent.name().unwrap()
                    )
                }
            }
        }
    }
}

struct MS {}

#[derive(Clone, Debug, Decode, Default, Encode)]
struct MSData {
    index: Option<usize>,
    name: Box<str>,
    definition: Box<str>,
    synonyms: Vec<(SynonymScope, Box<str>)>,
    cross_ids: Vec<(Option<Box<str>>, Box<str>)>,
    is_a: Vec<usize>,
}

impl CVData for MSData {
    type Index = usize;
    fn index(&self) -> Option<usize> {
        self.index
    }
    fn curie(&self) -> Option<mzcv::Curie> {
        self.index
            .map(|v| ControlledVocabulary::MS.curie(mzcv::AccessionCode::Numeric(v as u32)))
    }
    fn name(&self) -> Option<Cow<'_, str>> {
        Some(Cow::Borrowed(&self.name))
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.synonyms
            .iter()
            .filter_map(|(s, n)| (*s == SynonymScope::Exact).then_some(n.as_ref()))
    }

    fn parents(&self) -> impl Iterator<Item = &Self::Index> {
        self.is_a.iter()
    }
}

impl CVSource for MS {
    type Data = MSData;
    type Structure = Vec<Arc<MSData>>;
    fn cv_name() -> &'static str {
        "MS"
    }
    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "MS",
            extension: "obo",
            url: Some("http://purl.obolibrary.org/obo/ms.obo"),
            compression: mzcv::CVCompression::None,
        }]
    }
    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        None
    }
    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, Self::Structure), Vec<BoxedError<'static, CVError>>> {
        let reader = reader.next().unwrap();
        OboOntology::from_raw(reader)
            .map_err(|e| {
                vec![
                    BoxedError::small(
                        CVError::FileCouldNotBeParsed,
                        e.get_short_description(),
                        e.get_long_description(),
                    )
                    .add_contexts(e.get_contexts().iter().cloned()),
                ]
            })
            .map(|obo| {
                (
                    obo.version(),
                    obo.objects
                        .into_iter()
                        .filter(|o| o.stanza_type == OboStanzaType::Term)
                        .map(|obj| {
                            let mut data = MSData {
                                index: obj.id.1.parse().ok(),
                                name: obj.lines["name"][0].0.clone(),
                                synonyms: obj
                                    .synonyms
                                    .iter()
                                    .map(|s| (s.scope, s.synonym.clone()))
                                    .collect(),
                                ..Default::default()
                            };
                            if let Some((def, ids, _, _)) = obj.definition {
                                data.definition = def;
                                data.cross_ids = ids;
                            }
                            for parent in obj.is_a.iter() {
                                if let Ok(parent_id) = parent.1.parse() {
                                    data.is_a.push(parent_id);
                                }
                            }
                            Arc::new(data)
                        })
                        .collect(),
                )
            })
    }
}
