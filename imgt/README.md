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
let first = selection.germlines().next().unwrap();
assert_eq!(first.name(), "IGHV1-2*01");
```

## Compilation features

* `rayon` - enables parallel iterators using rayon
