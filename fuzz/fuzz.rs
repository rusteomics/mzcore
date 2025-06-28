#!/usr/bin/env -S cargo +nightly -Zscript

---
[package]
edition = "2024"
[dependencies]
clap = { version = "4.2", features = ["derive"] }
---

use clap::Parser;
use std::io::{BufRead, BufWriter, Write};
use std::iter::Sum;
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::time::Duration;

use std::thread;

#[derive(Parser, Debug)]
#[clap(version)]
struct Args {
    #[clap(help = "Name of the target")]
    target: String,
    #[clap(help = "Number of jobs to launch")]
    jobs: u8,
    #[clap(long, help = "Skip the set up part")]
    skip_setup: bool,
    #[clap(long, help = "Skip the set up part")]
    minify: bool,
}

fn main() {
    let args = Args::parse();
    let target = args.target;
    let input = format!("fuzz/in_{target}");
    let output = format!("out_{target}");
    let bin = format!("target/release/{target}");

    if args.minify {
        minify(&target, args.jobs);
    } else {
        // setup
        if !args.skip_setup {
            println!("Build config");
            Command::new("cargo")
                .args(["afl", "config", "--build"])
                .output()
                .expect("Failed to build config");
            println!("System config");
            Command::new("cargo")
                .args(["afl", "system-config"])
                .output()
                .expect("Failed to config");
            println!("Build targets");
            Command::new("cargo")
                .args(["afl", "build", "--release", "-p", "rustyms-fuzz"])
                .output()
                .expect("Failed to build targets");
        }

        // run
        println!("Run fuzzer");
        let id = AtomicUsize::new(0);
        thread::scope(|s| {
            for i in 0..args.jobs {
                if i == 0 {
                    s.spawn(|| {
                        Command::new("cargo")
                            .args(["afl", "fuzz", "-i", &input, "-o", &output, "-M", "fM", &bin])
                            .output()
                            .expect("Failed to launch main");
                    });
                } else {
                    s.spawn(|| {
                        let i = id.fetch_add(1, Ordering::SeqCst);
                        let name = format!("f{i}");
                        Command::new("cargo")
                            .args([
                                "afl", "fuzz", "-i", &input, "-o", &output, "-S", &name, &bin,
                            ])
                            .output()
                            .expect("Failed to launch child");
                    });
                }
            }

            s.spawn(|| {
                thread::sleep(Duration::from_millis(2000));
                loop {
                    thread::sleep(Duration::from_millis(1000));
                    Stats::get(&target, args.jobs);
                }
            });
        });
    }
}

#[derive(Default)]
struct Stats {
    fuzz_time: usize,
    cycles_done: usize,
    execs_done: usize,
    saved_crashes: usize,
    saved_hangs: usize,
    jobs: u8,
}

impl Stats {
    fn get(target: &str, jobs: u8) -> Option<()> {
        let stats = (-1..(jobs - 1) as i16)
            .map(|i| {
                let name = if i == -1 {
                    "fM".to_string()
                } else {
                    format!("f{i}")
                };
                let mut res = Self::default();
                res.jobs = jobs;
                for line in std::io::BufReader::new(
                    std::fs::File::open(&format!("out_{target}/{name}/fuzzer_stats"))
                        .expect(&format!("out_{target}/{name}/fuzzer_stats")),
                )
                .lines()
                {
                    let line = line.unwrap();
                    let (tag, value) = line.split_once(':').unwrap();
                    match tag.trim() {
                        "fuzz_time" => {
                            res.fuzz_time = value.trim().parse().unwrap();
                        }
                        "cycles_done" => {
                            res.cycles_done = value.trim().parse().unwrap();
                        }
                        "execs_done" => {
                            res.execs_done = value.trim().parse().unwrap();
                        }
                        "saved_crashes" => {
                            res.saved_crashes = value.trim().parse().unwrap();
                        }
                        "saved_hangs" => {
                            res.saved_hangs = value.trim().parse().unwrap();
                        }
                        _ => (),
                    }
                }
                res
            })
            .sum::<Self>();
        stats.print();
        Some(())
    }

    fn print(&self) {
        let duration = Duration::from_secs_f64(self.fuzz_time as f64 / self.jobs as f64);
        print!(
            "\rtime: {}, cycles: {}, execs: {}, crashes: {} hangs: {}",
            print_time(duration),
            self.cycles_done,
            self.execs_done,
            self.saved_crashes,
            self.saved_hangs,
        )
    }

    fn add(&mut self, other: &Self) {
        self.fuzz_time += other.fuzz_time;
        self.cycles_done += other.cycles_done;
        self.execs_done += other.execs_done;
        self.saved_crashes += other.saved_crashes;
        self.saved_hangs += other.saved_hangs;
        self.jobs = other.jobs;
    }
}

impl Sum for Stats {
    fn sum<I: Iterator<Item = Stats>>(iter: I) -> Self {
        let mut res = Self::default();
        for item in iter {
            res.add(&item);
        }
        res
    }
}

fn minify(target: &str, jobs: u8) {
    let mut options = Vec::new();
    for job in -1..(jobs - 1) as i16 {
        let name = if job == -1 {
            "fM".to_string()
        } else {
            format!("f{job}")
        };
        for item in std::fs::read_dir(&format!("out_{target}/{name}/crashes")).unwrap() {
            if let Ok(item) = item {
                if !item.path().extension().is_some_and(|e| e == "txt") {
                    options.push(std::fs::read_to_string(item.path()).unwrap())
                }
            }
        }
    }
    let file = std::fs::File::create(&format!("crashes_{target}.rs")).unwrap();
    let mut writer = BufWriter::new(file);
    let mut count = 0;
    for option in options {
        let mut child = Command::new(format!("target/release/{target}"))
            .stdin(Stdio::piped())
            .stdout(Stdio::null())
            .stderr(Stdio::null())
            .spawn()
            .expect("Failed to launch child");
        write!(child.stdin.take().unwrap(), "{}", &option).unwrap();
        if !child.wait().unwrap().success() {
            writeln!(writer, "test!(r\"{option}\", test_{count});",).unwrap();
            count += 1;
        }
    }
}

fn print_time(time: Duration) -> String {
    const BIG_SUFFIXES: &[(f64, &str)] = &[
        (1.0, "s"),
        (60.0, "m"),
        (60.0 * 60.0, "h"),
        (60.0 * 60.0 * 24.0, "d"),
        (60.0 * 60.0 * 24.0 * 7.0, "w"),
    ];
    const SMALL_SUFFIXES: &[(f64, &str)] = &[
        (1.0, "s"),
        (1e-3, "ms"),
        (1e-6, "μs"),
        (1e-9, "ns"),
        (1e-12, "ps"),
        (1e-15, "fs"),
    ];
    const PRECISION: usize = 2;
    let value = time.as_secs_f64();
    let base = if value == 0.0 {
        (1.0, "s")
    } else if value <= 1.0 {
        SMALL_SUFFIXES[(-((value.abs().log10() / 3.0).floor() as isize)
            .clamp(-(SMALL_SUFFIXES.len() as isize - 1), 0)) as usize]
    } else {
        BIG_SUFFIXES
            .iter()
            .filter(|(v, _)| value > *v)
            .last()
            .unwrap()
            .clone()
    };
    format!("{:.PRECISION$}{}", value / base.0, base.1)
}
