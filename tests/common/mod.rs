//! Shared helpers for the simulation integration tests.
//!
//! These drive the real binary (`CARGO_BIN_EXE_baitbench`) rather than the library,
//! because the simulation pipeline lives in the bin crate: `src/lib.rs` exposes only
//! `sdust` and `thermodynamics`.

#![allow(dead_code)]

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

pub const BIN: &str = env!("CARGO_BIN_EXE_baitbench");

/// Absolute path to `examples/minimal/<name>`.
pub fn example(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("examples/minimal")
        .join(name)
}

fn run(args: &[&str]) {
    let out = Command::new(BIN)
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to spawn {BIN}: {e}"));
    assert!(
        out.status.success(),
        "baitbench {:?} failed ({})\nstdout:\n{}\nstderr:\n{}",
        args,
        out.status,
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr),
    );
}

/// Run `prepare` on the minimal example at a given distractor fraction.
///
/// Returns (combined_reference.fa, weights.txt, distractors.txt).
pub fn prepare(outdir: &Path, distractor_fraction: f64) -> (PathBuf, PathBuf, PathBuf) {
    run(&[
        "prepare",
        "--targets",
        example("targets.fa").to_str().unwrap(),
        "--distractors",
        example("distractors.fa").to_str().unwrap(),
        "--distractor-fraction",
        &distractor_fraction.to_string(),
        "--outdir",
        outdir.to_str().unwrap(),
    ]);
    (
        outdir.join("combined_reference.fa"),
        outdir.join("weights.txt"),
        outdir.join("distractors.txt"),
    )
}

/// Run `simulate` and return the path to the fragments FASTA.
pub fn simulate(
    outdir: &Path,
    reference: &Path,
    weights: &Path,
    num_fragments: usize,
    capture_fraction: f64,
    seed: u64,
) -> PathBuf {
    let output = outdir.join("fragments.fa");
    run(&[
        "simulate",
        "--reference",
        reference.to_str().unwrap(),
        "--weights",
        weights.to_str().unwrap(),
        "--probes",
        example("probes.fa").to_str().unwrap(),
        "--num-fragments",
        &num_fragments.to_string(),
        "--capture-fraction",
        &capture_fraction.to_string(),
        "--seed",
        &seed.to_string(),
        "--output",
        output.to_str().unwrap(),
    ]);
    output
}

/// Count fragments per source sequence, keyed by the `{seq_id}_fragment_{n}` contract.
pub fn fragment_sources(fragments: &Path) -> HashMap<String, usize> {
    let text = std::fs::read_to_string(fragments)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", fragments.display()));
    let mut counts = HashMap::new();
    for line in text.lines().filter(|l| l.starts_with('>')) {
        let id = line[1..].split_whitespace().next().unwrap_or("");
        let src = match id.rfind("_fragment_") {
            Some(p) => &id[..p],
            None => panic!("fragment header does not match the expected contract: {line}"),
        };
        *counts.entry(src.to_string()).or_insert(0) += 1;
    }
    counts
}

/// IDs listed in a one-per-line ID file (e.g. `distractors.txt`).
pub fn id_set(path: &Path) -> Vec<String> {
    std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()))
        .lines()
        .map(|l| l.trim().to_string())
        .filter(|l| !l.is_empty())
        .collect()
}

/// Realized fraction of fragments originating from distractor sequences.
pub fn realized_distractor_fraction(fragments: &Path, distractors: &Path) -> f64 {
    let counts = fragment_sources(fragments);
    let distractor_ids = id_set(distractors);
    let total: usize = counts.values().sum();
    assert!(total > 0, "no fragments were generated");
    let d: usize = counts
        .iter()
        .filter(|(id, _)| distractor_ids.iter().any(|x| x == *id))
        .map(|(_, n)| *n)
        .sum();
    d as f64 / total as f64
}
