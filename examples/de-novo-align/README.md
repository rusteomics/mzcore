# De Novo align

Usage:
```
cargo run --release --bin de-novo-align -- --peptides rustyms\data\200305_HER_test_04_DENOVO.csv.gz --database examples\de-novo-align\database.fasta --out-path out.csv
```

This aligns all peptides from a given PSM file, see rustyms for a list of all supported files, to a list of known proteins. It returns a CSV file with the best alignment for each _de novo_ peptide. This can be used to look into how good the _de novo_ predictions actually are.
