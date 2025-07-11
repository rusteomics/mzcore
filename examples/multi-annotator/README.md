# Multi annotator

Usage:
```
cargo run -p multi-annotator --release -- --in-path examples/multi-annotator/test.csv --out-path examples/multi-annotator/test-out2.csv
```

This takes a CSV file as input that contains a peptidoform and the rawfile (MGF or mzML) it originated from, it then annotates the spectrum with the theoretical fragmentation from rustyms and delivers some statistics on the annotation in a resulting CSV file. This can be used to get a global impression over a whole dataset, so for example see if a certain fragmentation energy increases or decreases the coverage of a particular ion series (peptidoform or glycan).