//! Integration tests for the fragment simulation.
//!
//! Two of these were written against the pre-rework model and failed by design; both
//! are now green (see the plan's Stages 0-3):
//!
//! * `same_seed_is_reproducible` — defect D6. HashMap iteration order leaked into
//!   sampling, so `--seed` was not reproducible between processes. Fixed in Stage 1.
//! * `distractor_fraction_is_honoured` — reviewer item 6. The realized distractor
//!   fraction saturated at `1 - capture_fraction`. Fixed in Stage 3.

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
/// This is reviewer item 6. `--distractor-fraction` sets the *input* library
/// composition, so the realized output fraction is not expected to equal it at moderate
/// values — capture enriches the targets, which is the point. What was wrong is that
/// the realized fraction could not exceed `1 - capture_fraction` no matter how extreme
/// the input, so the parameter went silently inert over half its range.
///
/// Measured under the retention model (capture_fraction 0.5, background_retention
/// 1e-3, so roughly 500x enrichment):
///
/// ```text
///     requested 0.1000 -> realized 0.00000
///     requested 0.5000 -> realized 0.00125
///     requested 0.9000 -> realized 0.01700
///     requested 0.9900 -> realized 0.16700
///     requested 0.9990 -> realized 0.66775
///     requested 0.9999 -> realized 0.95500
/// ```
///
/// Against the quota model the same sweep gave 0.043 / 0.242 / 0.450 / 0.495 / 0.495 /
/// 0.495 — compressed into half the range and pinned at the ceiling.
#[test]
fn distractor_fraction_is_honoured() {
    let fractions = [0.1, 0.5, 0.9, 0.99, 0.999, 0.9999];
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
        .map(|(req, got)| format!("requested {req:.4} -> realized {got:.5}"))
        .collect();

    eprintln!("distractor titration:\n  {}", report.join("\n  "));

    // Strictly increasing: each step up in request must raise the realized fraction.
    for i in 1..realized.len() {
        assert!(
            realized[i] > realized[i - 1],
            "realized distractor fraction is not monotone in the request\n  {}",
            report.join("\n  ")
        );
    }

    // The ceiling is the actual regression: the quota model could not put more than
    // `1 - capture_fraction` of the output into distractors whatever the input was.
    let capture_fraction = 0.5;
    let ceiling = 1.0 - capture_fraction;
    let got_top = realized[realized.len() - 1];
    assert!(
        got_top > ceiling,
        "realized distractor fraction never exceeded the old 1 - capture_fraction \
         ceiling of {ceiling:.3} (best {got_top:.3}), so the parameter is still inert \
         over the top half of its range\n  {}",
        report.join("\n  ")
    );

    // And at an input that is almost entirely distractor, the output must be too:
    // enrichment is finite, so it cannot rescue a target present at 1 in 10,000.
    let (req_top, got_top) = (fractions[fractions.len() - 1], realized[realized.len() - 1]);
    assert!(
        (got_top - req_top).abs() < 0.1,
        "at --distractor-fraction {req_top:.4} the realized fraction was {got_top:.3}\n  {}",
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
