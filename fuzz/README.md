# Fuzzing

Uses [cargo-afl](https://crates.io/crates/cargo-afl) for [fuzz testing](https://en.wikipedia.org/wiki/Fuzzing). Note that this only works on Linux.

## Manually

From root directory
```
cargo-afl afl system-config
cargo afl build --release -p rustyms-fuzz
cargo afl fuzz -i fuzz/in_pro_forma -o out_pro_forma target/release/pro_forma
cargo afl fuzz -i fuzz/in_pro_forma -o out_sloppy_pro_forma target/release/sloppy_pro_forma
cargo afl fuzz -i fuzz/in_peaks -o out_peaks target/release/peaks
```
Several fuzz targets are defined: `pro_forma`, `sloppy_pro_forma`, `peaks`. The all targets have an `in_<target>` directory with input examples.

After running the fuzzer the following commands can be used to easily save all crashes into a single file.
```
open out_pro_forma/default/crashes/* | save crashes.txt -f (nushell)
cat out_pro_forma/default/crashes/* >> crashes.txt (bash)
```

## Scripted

From root directory
```
./fuzz/fuzz.rs <project> <cores>
```

This uses a Rust script so need nightly Rust to be installed. If the script immediately gets through without building the project running the following:

```
cargo afl config --build --force
```

After fuzzing running:

```
./fuzz/fuzz.rs <project> <cores> --minify
```

Will create a `crashes_<target>.rs` file with all found crashes.
