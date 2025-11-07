# mzcv

Handle ontologies/controlled vocabularies with  getting data from many sources:

1. Download from a URL (only with feature `http`) [`CVIndex::update_from_url`]
2. Parse from a given file [`CVIndex::update_from_path`]
3. Open binary cache from the standardised cache location 
4. Update with values directly in memory [`CVIndex::update`]
5. Use statically included data [`CVSource::static_data`]

When downloading a file it downloads it to the standardised location and compresses it with gzip compression to not take up too much space. When opening a file (such as a downloaded file) it first parses the file. If it succeeds it places the file at the standardised location and stores the parsed data in the binary cache. If it fails to parse the file it will report the errors both to the caller of the method and leave the errors next to the standard file location for end user convenience.

## Lookup

There are three major ways of looking up data: [index](CVIndex::get_by_index), [name](CVIndex::get_by_name), and [fuzzy match](CVIndex::search) search. The first two use HashMaps to do constant time lookups, the second uses a trigram index (when `search-index` is turned on, or loops over all data if not) and loops over all matches using Levenshtein distance to find good enough matches.

## Compilation features

- `http` allow downloading ontologies from the internet
- `serde` allow using serde (de)serialise to store data in the cache
- `search-index` builds a trigram index to speed up fuzzy matching
