//! An example to show how to work with the IMGT ontology which uses server side compression
use std::io::Write;

use bincode::{Decode, Encode};
use context_error::{BoxedError, CreateError, FullErrorContent, StaticErrorContent};

use mzcv::{
    CVData, CVError, CVFile, CVIndex, CVSource, CVVersion, HashBufReader, OboOntology,
    OboStanzaType,
};
use sha2::Digest;

fn main() {
    let (mut imgt, _) = CVIndex::<IMGT>::init();
    imgt.update_from_url(&[]).unwrap();
    println!(
        "Stored at: {}",
        IMGT::default_stem()
            .with_extension(IMGT::files()[0].extension)
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

#[derive(Clone, Debug, Decode, Default, Encode)]
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
    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "IMGT",
            extension: "dat",
            url: Some("https://www.imgt.org/download/LIGM-DB/imgt.dat.Z"),
            compression: mzcv::CVCompression::LZW,
        }]
    }
    fn static_data() -> Option<(CVVersion, &'static [std::sync::Arc<Self::Data>])> {
        None
    }
    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl Digest>>,
    ) -> Result<(CVVersion, impl Iterator<Item = Self::Data>), Vec<BoxedError<'static, CVError>>>
    {
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
                            let mut data = IMGTData {
                                index: obj.id.1.parse().ok(),
                                name: obj.lines["name"][0].0.to_string(),
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
