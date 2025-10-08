//! An example to show how to work with the IMGT ontology which uses server side compression
use std::io::Write;

use bincode::{Decode, Encode};
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVData, CVError, CVIndex, CVSource, CVVersion, HashBufReader, OboOntology, OboStanzaType,
};
use sha2::Digest;

fn main() {
    let (mut imgt, _) = CVIndex::<IMGT>::init();
    imgt.update_from_url(None).unwrap();
    println!(
        "Stored at: {}",
        IMGT::default_stem()
            .with_extension(IMGT::cv_extension())
            .display()
    );
    println!("Found {} IMGT items", imgt.data().len());
    println!(
        "Date: {} Version: {} Hash: {}",
        imgt.version().last_updated.map_or_else(
            || "?".to_string(),
            |(y, m, d, h, min)| format!("{y}-{m}-{d} {h:02}:{min:02}")
        ),
        imgt.version()
            .version
            .as_ref()
            .map_or_else(|| "?".to_string(), ToString::to_string),
        imgt.version().hash_hex(),
    );

    println!("Search for terms in the IMGT ontology (prefix with '==' to search for exact names)");

    loop {
        print!("?:");
        std::io::stdout().flush().unwrap();
        let mut term = String::new();
        std::io::stdin().read_line(&mut term).unwrap();
        let term = term.trim();
        let answers = term.strip_prefix("==").map_or_else(
            || imgt.search(term, 10, 6),
            |term| imgt.get_by_name(term).map(|v| vec![v]).unwrap_or_default(),
        );
        if answers.is_empty() {
            println!("No matches found");
        } else {
            for a in answers {
                println!(
                    " > IMGT:{} '{}': {}",
                    a.index
                        .map_or_else(|| "-".to_string(), |i| format!("{i:07}")),
                    a.name,
                    a.definition
                );
            }
        }
    }
}

struct IMGT {}

#[derive(Clone, Decode, Encode, Debug, Default)]
struct IMGTData {
    index: Option<usize>,
    name: String,
    definition: String,
    synonyms: Vec<String>,
    cross_ids: Vec<(Option<String>, String)>,
}

impl CVData for IMGTData {
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

impl CVSource for IMGT {
    type Data = IMGTData;
    fn cv_name() -> &'static str {
        "IMGT"
    }
    fn cv_extension() -> &'static str {
        "dat"
    }
    fn cv_url() -> Option<&'static str> {
        Some("https://www.imgt.org/download/LIGM-DB/imgt.dat.Z")
    }
    fn cv_compression() -> mzcv::CVCompression {
        mzcv::CVCompression::LZW
    }
    fn static_data() -> Option<(CVVersion, &'static [std::sync::Arc<Self::Data>])> {
        None
    }
    fn parse(
        reader: HashBufReader<impl std::io::Read, impl Digest>,
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
                            let mut data = IMGTData {
                                index: obj.lines["id"][0]
                                    .split_once(':')
                                    .expect("Incorrect psi ms id, should contain a colon")
                                    .1
                                    .parse()
                                    .ok(),
                                name: obj.lines["name"][0].to_string(),
                                synonyms: obj.lines.get("synonym").cloned().unwrap_or_default(),
                                ..Default::default()
                            };
                            let (def, ids) = obj.definition();
                            data.definition = def.to_string();
                            data.cross_ids = ids
                                .into_iter()
                                .map(|(r, i)| (r.map(ToString::to_string), i.to_string()))
                                .collect();
                            data
                        }),
                )
            })
    }
}
