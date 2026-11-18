//! Integration tests for the fragment simulation.
//!
//! Two of these are expected to FAIL against the pre-rework model; they are the
//! regressions the simulation rework must turn green (see the plan's Stage 0):
//!
//! * `same_seed_is_reproducible` — defect D6. HashMap iteration order leaks into
//!   sampling, so `--seed` is not reproducible between processes.
//! * `distractor_fraction_is_honoured` — reviewer item 6. The realized distractor
//!   fraction saturates at `1 - capture_fraction`.

mod common;

use common::*;
use tempfile::TempDir;

/// The same `--seed` must produce byte-identical fragments.
///
/// Fixed in Stage 1 (defect D6): sampler key vectors are now sorted, so HashMap
/// iteration order no longer leaks into a seeded run.
#[test]
fn same_seed_is_reproducible() {
    let a = TempDir::new().unwrap();
    let b = TempDir::new().unwrap();

    let (ref_a, w_a, _) = prepare(a.path(), 0.5);
    let frags_a = simulate(a.path(), &ref_a, &w_a, 2000, 0.5, 42);

    let (ref_b, w_b, _) = prepare(b.path(), 0.5);
    let frags_b = simulate(b.path(), &ref_b, &w_b, 2000, 0.5, 42);

    let ta = std::fs::read_to_string(&frags_a).unwrap();
    let tb = std::fs::read_to_string(&frags_b).unwrap();

    if ta != tb {
        // Report the first divergence rather than dumping two whole FASTA files.
        let first_diff = ta
            .lines()
            .zip(tb.lines())
            .enumerate()
            .find(|(_, (x, y))| x != y)
            .map(|(i, (x, y))| format!("line {}:\n    run A: {x}\n    run B: {y}", i + 1))
            .unwrap_or_else(|| {
                format!("same prefix, different length ({} vs {} lines)",
                        ta.lines().count(), tb.lines().count())
            });
        panic!(
            "identical --seed produced different fragments across processes \
             (sampling order depends on HashMap iteration order)\n  {first_diff}"
        );
    }
}

/// Requesting more distractor must actually deliver more distractor, across the whole
/// range — not saturate once `distractor_fraction` exceeds `1 - capture_fraction`.
///
/// KNOWN FAILING (reviewer item 6). Un-ignore in Stage 3.
/// Run with: `cargo test --test simulation -- --ignored`
#[test]
#[ignore = "reviewer item 6: realized distractor fraction saturates at 1 - capture_fraction; fixed in Stage 3"]
fn distractor_fraction_is_honoured() {
    let fractions = [0.1, 0.5, 0.9, 0.99];
    let mut realized = Vec::new();

    for &df in &fractions {
        let dir = TempDir::new().unwrap();
        let (reference, weights, distractors) = prepare(dir.path(), df);
        let frags = simulate(dir.path(), &reference, &weights, 4000, 0.5, 42);
        realized.push(realized_distractor_fraction(&frags, &distractors));
    }

    let report: Vec<String> = fractions
        .iter()
        .zip(&realized)
        .map(|(req, got)| format!("requested {req:.3} -> realized {got:.3}"))
        .collect();

    // Strictly increasing: each step up in request must raise the realized fraction.
    for i in 1..realized.len() {
        assert!(
            realized[i] > realized[i - 1],
            "realized distractor fraction is not monotone in the request\n  {}",
            report.join("\n  ")
        );
    }

    // And the top of the range must actually get close to what was asked for.
    let (req_top, got_top) = (fractions[fractions.len() - 1], realized[realized.len() - 1]);
    assert!(
        (got_top - req_top).abs() < 0.1,
        "at --distractor-fraction {req_top:.3} the realized fraction was {got_top:.3}\n  {}",
        report.join("\n  ")
    );
}

/// Every requested fragment should be produced; shortfalls must not be silent.
/// Covers defect D9 (background sampler drops draws with no retry and no warning).
#[test]
fn requested_fragment_count_is_delivered() {
    let dir = TempDir::new().unwrap();
    let (reference, weights, _) = prepare(dir.path(), 0.9);
    let requested = 3000;
    let frags = simulate(dir.path(), &reference, &weights, requested, 0.5, 7);

    let produced: usize = fragment_sources(&frags).values().sum();
    assert_eq!(
        produced, requested,
        "simulate produced {produced} of {requested} requested fragments"
    );
}

/// The fragment header contract that everything downstream depends on
/// (`io_utils::extract_source_id`, metrics classification, read renaming).
#[test]
fn fragment_headers_follow_the_contract() {
    let dir = TempDir::new().unwrap();
    let (reference, weights, _) = prepare(dir.path(), 0.5);
    let frags = simulate(dir.path(), &reference, &weights, 500, 0.5, 1);

    let text = std::fs::read_to_string(&frags).unwrap();
    let mut n = 0;
    for line in text.lines().filter(|l| l.starts_with('>')) {
        n += 1;
        let id = line[1..].split_whitespace().next().unwrap();
        assert!(
            id.contains("_fragment_"),
            "header missing _fragment_ delimiter: {line}"
        );
        assert!(
            line.contains("start=") && line.contains("length="),
            "header missing start=/length= fields: {line}"
        );
    }
    assert!(n > 0, "no fragments written");
}
