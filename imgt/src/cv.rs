//! Code to handle the PSI-MOD ontology
use std::{
    borrow::Cow,
    collections::HashMap,
    sync::{Arc, LazyLock},
};

use context_error::BoxedError;

use mzcv::{CVData, CVError, CVFile, CVIndex, CVSource, CVStructure, CVVersion, HashBufReader};

use crate::{Gene, Germline, Germlines, Species, parse::parse_dat};

/// A single shared static access to the static data in the ontologies for cases where no runtime resolution is needed (like tests).
pub static STATIC_IMGT: LazyLock<CVIndex<IMGT>> = LazyLock::new(CVIndex::init_static);

/// IMGT antibody germlines
#[allow(missing_copy_implementations, missing_debug_implementations)]
pub struct IMGT {}

impl CVData for Germline {
    type Index = (Species, Gene);
    fn index(&self) -> Option<Self::Index> {
        Some((self.species, self.name.clone()))
    }
    fn name(&self) -> Option<Cow<'_, str>> {
        Some(Cow::Owned(self.name.to_string()))
    }
    fn synonyms(&self) -> impl Iterator<Item = &str> {
        std::iter::empty()
    }
}

impl CVSource for IMGT {
    type Data = Germline;
    type Structure = HashMap<Species, Germlines>;
    fn cv_name() -> &'static str {
        "IMGT"
    }

    fn files() -> &'static [CVFile] {
        &[CVFile {
            name: "IMGT",
            extension: "dat",
            url: None, // TODO: Uses .Z compression and so cannot automatically be processed
            compression: mzcv::CVCompression::None,
        }]
    }

    fn static_data() -> Option<(CVVersion, Self::Structure)> {
        #[cfg(not(feature = "internal-no-data"))]
        {
            use bincode::config::Configuration;
            let cache = bincode::decode_from_slice::<(CVVersion, Self::Structure), Configuration>(
                include_bytes!("IMGT.dat"),
                Configuration::default(),
            )
            .unwrap()
            .0;
            Some(cache)
        }
        #[cfg(feature = "internal-no-data")]
        None
    }

    fn parse(
        mut reader: impl Iterator<Item = HashBufReader<Box<dyn std::io::Read>, impl sha2::Digest>>,
    ) -> Result<(CVVersion, Self::Structure), Vec<BoxedError<'static, CVError>>> {
        let mut reader = reader.next().unwrap();
        let data = parse_dat(&mut reader);

        let (grouped, _errors) = crate::combine::combine(data);

        let version = CVVersion {
            hash: reader.hash(),
            ..Default::default()
        };

        Ok((version, grouped))
    }
}

#[expect(
    clippy::implicit_hasher,
    reason = "Gave some issues with default and lifetimes, likely easily fixed but just could not be bothered"
)]
impl CVStructure<Germline> for HashMap<Species, Germlines> {
    fn is_empty(&self) -> bool {
        self.is_empty()
    }
    fn len(&self) -> usize {
        self.values().fold(0, |acc, s| {
            acc + s.iter().fold(0, |acc, (_, g)| {
                acc + g.iter().fold(0, |acc, (_, g)| acc + g.len())
            })
        })
    }
    fn clear(&mut self) {
        self.clear();
    }
    fn add(&mut self, data: Arc<Germline>) {
        self.entry(data.species)
            .or_insert_with(|| Germlines::new(data.species))
            .insert(Arc::unwrap_or_clone(data));
    }
    type Index = (Species, Gene);
    type IterIndexed<'a> = Box<dyn Iterator<Item = (Self::Index, Arc<Germline>)> + 'a>;
    fn iter_indexed(&self) -> Self::IterIndexed<'_> {
        Box::new(self.iter().flat_map(|(species, germlines)| {
            germlines.iter().flat_map(move |(_, germlines)| {
                germlines.iter().flat_map(move |(_, germlines)| {
                    germlines
                        .iter()
                        .map(move |germline| ((*species, germline.name.clone()), germline.clone()))
                })
            })
        }))
    }
    type IterData<'a> = Box<dyn Iterator<Item = Arc<Germline>> + 'a>;
    fn iter_data(&self) -> Self::IterData<'_> {
        Box::new(self.iter().flat_map(|(_, germlines)| {
            germlines.iter().flat_map(|(_, germlines)| {
                germlines
                    .iter()
                    .flat_map(|(_, germlines)| germlines.iter().map(Clone::clone))
            })
        }))
    }
    fn index(&self, index: Self::Index) -> Option<Arc<Germline>> {
        self.get(&index.0)
            .and_then(|germlines| germlines.find_germline(index.1))
    }
    fn remove(&mut self, index: Self::Index) {
        if let Some(germlines) = self.get_mut(&index.0) {
            germlines.remove_germline(index.1);
        }
    }
}
