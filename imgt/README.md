# IMGT

Handle the IMGT database of antibody germlines easily.

## Library features

 - Access to all germlines from all species in IMGT
 - Access to all isotopes of the germlines
 - Access to the regions (CDRs etc) and annotations (conserved etc)
 - Single threaded and multi threaded access

## Example usage

```rust
use imgt::*;
let selection = Selection::default()
    .species([Species::HomoSapiens])
    .chain([ChainType::Heavy])
    .gene([GeneType::V]);
let first = selection.germlines(&STATIC_IMGT).next().unwrap();
assert_eq!(first.name(), "IGHV1-2*01");
```

## Compilation features

* `rayon` - enables parallel iterators using rayon

## Changelog
### 0.2.0 

- Added more constant chain options
- Added Lzw parsing (in mzcv) to allow automatic download and update
- Fixed some issues with parsing
- Updated static data to version 202631-5 date 2026-07-31