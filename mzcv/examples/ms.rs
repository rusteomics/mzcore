//! An example to show how to work with the MS ontology
use std::io::{BufReader, Write};

use bincode::{Decode, Encode};
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVData, CVError, CVIndex, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType,
};

fn main() {
    let (mut ms, _) = CVIndex::<MS>::init();
    // ms.update_from_url(None).unwrap();
    ms.update_from_path(None).unwrap();
    println!(
        "Stored at: {}",
        MS::default_stem()
            .with_extension(MS::cv_extension())
            .display()
    );
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

    let item = ms.get_by_index(1_000_016).unwrap();
    assert_eq!(item.name, "scan start time");

    println!("Search for terms in the MS ontology (prefix with '==' to search for exact names)");

    loop {
        print!("?:");
        std::io::stdout().flush().unwrap();
        let mut term = String::new();
        std::io::stdin().read_line(&mut term).unwrap();
        let term = term.trim();
        let answers = term.strip_prefix("==").map_or_else(
            || ms.search(term, 10, 6),
            |term| ms.get_by_name(term).map(|v| vec![v]).unwrap_or_default(),
        );
        if answers.is_empty() {
            println!("No matches found");
        } else {
            for a in answers {
                println!(
                    " > MS:{} '{}': {}",
                    a.index
                        .map_or_else(|| "-".to_string(), |i| format!("{i:07}")),
                    a.name,
                    a.definition
                );
            }
        }
    }
}

struct MS {}

#[derive(Clone, Decode, Encode, Debug, Default)]
struct MSData {
    index: Option<usize>,
    name: String,
    definition: String,
    synonyms: Vec<String>,
    cross_ids: Vec<(Option<String>, String)>,
}

impl CVData for MSData {
    type Index = usize;
    fn index(&self) -> Option<usize> {
        self.index
    }
    fn name(&self) -> Option<&str> {
        Some(&self.name)
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        self.synonyms.iter().map(AsRef::as_ref)
    }
}

impl CVSource for MS {
    type Data = MSData;
    fn cv_name() -> &'static str {
        "MS"
    }
    fn cv_extension() -> &'static str {
        "obo"
    }
    fn cv_url() -> Option<&'static str> {
        Some("http://purl.obolibrary.org/obo/ms.obo")
    }
    fn cv_compression() -> mzcv::CVCompression {
        mzcv::CVCompression::None
    }
    fn static_data() -> Option<(CVVersion, &'static [std::sync::Arc<Self::Data>])> {
        None
    }
    fn parse(
        reader: HashBufReader<impl std::io::Read, impl sha2::Digest>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
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
                                name: obj.lines["name"][0].to_string(),
                                synonyms: obj.synonyms.iter().map(|s| s.synonym.clone()).collect(),
                                ..Default::default()
                            };
                            if let Some((def, ids, _, _)) = obj.definition {
                                data.definition = def;
                                data.cross_ids = ids;
                            }
                            data
                        }),
                )
            })
    }
}
