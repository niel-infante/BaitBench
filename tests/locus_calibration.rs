//! Calibration harness for per-locus capture efficiency.
//!
//! The Borrelia LOD data gives eight loci with measured relative capture efficiency
//! spanning 14x (nifS 0.070 ... rplB 1.000). Reproducing that spread is the open
//! Stage 3 defect: the shipped model predicts ~0.9 for nifS because equilibrium
//! occupancy saturates and the independent-capture product then flattens what is left.
//!
//! Running the whole pipeline to test a change to the capture law takes hours. This
//! harness instead replays **cached probe placements** against the known gene intervals
//! and scores candidate laws directly, so a functional form can be evaluated in
//! seconds. It is `#[ignore]`d because it needs the external LOD data.
//!
//! ```text
//! cargo test --test locus_calibration -- --ignored --nocapture
//! ```
//!
//! Set `LOD_PROBE_CACHE` to a probe-hit cache built against the BB chromosome
//! (`tests/validation/config.sh` writes one to `$PROBE_CACHE_DIR`).

use std::collections::HashMap;

const FRAG_LEN: usize = 225; // the LOD Illumina prep
const MIN_OVERLAP: f64 = 0.5;

/// Measured relative efficiency, from expected/lod_locus_efficiency.tsv.
const REAL: &[(&str, usize, usize, f64)] = &[
    ("nifS", 79582, 80146, 0.0698),
    ("clpA", 377259, 378043, 0.0900),
    ("clpX", 637979, 640472, 0.5632),
    ("uvrA", 889721, 893774, 0.6158),
    ("pyrG", 587461, 590268, 0.7080),
    ("recG", 595930, 599191, 0.7858),
    ("16S", 443987, 446717, 0.9809),
    ("rplB", 497610, 499647, 1.0000),
];

/// A candidate capture law: label, per-probe ln K at unit concentration, and whether
/// overlapping probes combine by summing rates rather than by the independence product.
type Law = (String, Box<dyn Fn(&Placement) -> f64>, bool);

struct Placement {
    start: usize,
    end: usize,
    pairs: Vec<(u8, u8)>,
}

fn load(path: &str) -> Vec<Placement> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {path}: {e}"));
    let mut out = Vec::new();
    for line in text.lines() {
        if line.starts_with('#') || line.starts_with("probe_name") {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 || !f[1].starts_with("BB:") {
            continue;
        }
        out.push(Placement {
            start: f[2].parse().unwrap(),
            end: f[3].parse().unwrap(),
            pairs: f[5].bytes().zip(f[6].bytes()).collect(),
        });
    }
    out
}

/// Fragment start positions from which a probe hit can be captured.
fn window(p: &Placement) -> (usize, usize) {
    let hit_len = p.end - p.start;
    let m = ((hit_len as f64) * MIN_OVERLAP).ceil().max(1.0) as usize;
    ((p.start + m).saturating_sub(FRAG_LEN), p.end.saturating_sub(m))
}

/// Mean retention over the fragment starts spanning a locus, under a given per-probe
/// propensity and a given way of combining overlapping probes.
fn locus_retention(
    pl: &[Placement],
    idx: &[usize],
    theta: &HashMap<usize, f64>,
    st: usize,
    en: usize,
    sum_rule: bool,
) -> f64 {
    let lo = st.saturating_sub(FRAG_LEN);
    let mut total = 0.0;
    let mut n = 0usize;
    for s in lo..en {
        let (mut prod, mut sum) = (1.0f64, 0.0f64);
        // Only the placements that can reach this locus at all; iterating all 792 for
        // every position made the sweep take minutes per law.
        for &i in idx {
            let (a, b) = window(&pl[i]);
            if a <= s && s <= b {
                let t = theta[&i];
                prod *= 1.0 - t;
                sum += t;
            }
        }
        // Product: independent equilibrium per probe, 1 - prod(1 - theta).
        // Poisson: each probe contributes an independent capture *rate* over a finite
        // hybridisation window, so the count is Poisson(sum theta) and p(>=1) is
        // 1 - exp(-sum). In the low-capture regime this is linear in sum theta, so
        // probe density and affinity multiply instead of both being flattened.
        total += if sum_rule { 1.0 - (-sum).exp() } else { 1.0 - prod };
        n += 1;
    }
    if n > 0 { total / n as f64 } else { 0.0 }
}

/// Mean `p_specific` over probe-covered positions, mirroring production's
/// `capture_efficiency_at`.
///
/// The denominator matters: production averages over **covered positions only**, while
/// averaging over every start spanning a locus (as an earlier version of this harness
/// did) includes uncovered ground and so calibrates to a different concentration. That
/// mismatch is why the harness predicted a per-locus win the pipeline did not show.
///
/// Swept by events, as production does, so it is O(n log n) rather than per-position.
fn covered_mean_p_specific(pl: &[Placement], theta: &HashMap<usize, f64>, sum_rule: bool) -> f64 {
    let mut ev: Vec<(usize, usize, bool)> = Vec::new(); // (pos, idx, is_start)
    for (i, p) in pl.iter().enumerate() {
        let (a, b) = window(p);
        if b >= a {
            ev.push((a, i, true));
            ev.push((b + 1, i, false));
        }
    }
    ev.sort_unstable();
    let (mut sum_log, mut sum_theta, mut n_sat) = (0.0f64, 0.0f64, 0i32);
    let (mut total, mut covered) = (0.0f64, 0.0f64);
    let mut pos = 0usize;
    let mut k = 0usize;
    while k < ev.len() {
        let at = ev[k].0;
        if at > pos && (n_sat > 0 || sum_log < 0.0 || sum_theta > 0.0) {
            let span = (at - pos) as f64;
            let p = if sum_rule {
                1.0 - (-sum_theta).exp()
            } else if n_sat > 0 {
                1.0
            } else {
                1.0 - sum_log.exp()
            };
            if p > 0.0 {
                total += span * p;
                covered += span;
            }
        }
        while k < ev.len() && ev[k].0 == at {
            let t = theta[&ev[k].1].clamp(0.0, 1.0);
            let sign = if ev[k].2 { 1.0 } else { -1.0 };
            sum_theta += sign * t;
            if t >= 1.0 {
                n_sat += if ev[k].2 { 1 } else { -1 };
            } else if t > 0.0 {
                sum_log += sign * (1.0 - t).ln();
            }
            k += 1;
        }
        pos = at;
    }
    if covered > 0.0 { total / covered } else { 0.0 }
}

/// Placements that can contribute to any fragment start spanning a locus.
fn locus_index(pl: &[Placement], st: usize, en: usize) -> Vec<usize> {
    let lo = st.saturating_sub(FRAG_LEN);
    (0..pl.len())
        .filter(|&i| {
            let (a, b) = window(&pl[i]);
            b >= lo && a < en
        })
        .collect()
}

fn spearman(a: &[f64], b: &[f64]) -> f64 {
    let rank = |v: &[f64]| {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&i, &j| v[i].partial_cmp(&v[j]).unwrap());
        let mut r = vec![0.0; v.len()];
        for (k, &i) in idx.iter().enumerate() {
            r[i] = k as f64;
        }
        r
    };
    let (ra, rb) = (rank(a), rank(b));
    let n = ra.len() as f64;
    let (ma, mb) = (ra.iter().sum::<f64>() / n, rb.iter().sum::<f64>() / n);
    let num: f64 = ra.iter().zip(&rb).map(|(x, y)| (x - ma) * (y - mb)).sum();
    let da: f64 = ra.iter().map(|x| (x - ma).powi(2)).sum::<f64>().sqrt();
    let db: f64 = rb.iter().map(|y| (y - mb).powi(2)).sum::<f64>().sqrt();
    if da > 0.0 && db > 0.0 { num / (da * db) } else { f64::NAN }
}

#[test]
#[ignore = "needs the external Borrelia LOD probe cache; see the module docs"]
fn per_locus_capture_laws() {
    use baitbench::thermodynamics::{delta_g, min_window_delta_g, occupancy_from_ln_k, ThermoModel};

    let path = std::env::var("LOD_PROBE_CACHE").unwrap_or_else(|_| {
        format!(
            "{}/baitbench_probe_cache/lod_mtec_with_iscap.tsv",
            std::env::var("TMPDIR").unwrap_or_else(|_| "/tmp".into()).trim_end_matches('/')
        )
    });
    let pl = load(&path);
    assert!(!pl.is_empty(), "no BB placements in {path}");
    println!("\n{} placements on the BB chromosome", pl.len());

    let model = ThermoModel::new(58.0, 1.0);
    let real: Vec<f64> = REAL.iter().map(|r| r.3).collect();
    println!("{:26}  {:>10}  {:>13}   real: nifS={:.3} clpA={:.3} clpX={:.3} uvrA={:.3} recG={:.3} 16S={:.3}",
        "law", "", "", real[0], real[1], real[2], real[3], real[5], real[6]);

    // Cross the two changes: how overlapping probes combine, and whether per-probe
    // affinity comes from the whole duplex or its least stable window.
    const RT: f64 = 0.0019872 * 331.15;
    let mut laws: Vec<Law> = Vec::new();
    for (rule_name, sum_rule) in [("product", false), ("poisson", true)] {
        laws.push((
            format!("{rule_name:7} / global dG"),
            Box::new(move |p: &Placement| -delta_g(&p.pairs, &model) / RT),
            sum_rule,
        ));
        for w in [10usize, 12, 15, 20] {
            laws.push((
                format!("{rule_name:7} / weakest {w}bp"),
                Box::new(move |p: &Placement| -min_window_delta_g(&p.pairs, w, &model) / RT),
                sum_rule,
            ));
        }
    }

    let loci_idx: Vec<Vec<usize>> = REAL
        .iter()
        .map(|&(_, st, en, _)| locus_index(&pl, st, en))
        .collect();

    let targets: Vec<f64> = std::env::var("CAL_TARGETS")
        .unwrap_or_else(|_| "0.5".into())
        .split(',')
        .filter_map(|t| t.trim().parse().ok())
        .collect();

    for target in &targets {
    for (label, ln_k0, sum_rule) in &laws {
        // Calibrate the concentration so mean retention over covered positions is 0.5,
        // matching what the solver does, then report the per-locus ratios.
        let k0: Vec<f64> = pl.iter().map(ln_k0).collect();
        let (mut lo, mut hi) = (-400.0f64, 400.0f64);
        let mut best = Vec::new();
        for _ in 0..50 {
            let mid = 0.5 * (lo + hi);
            let theta: HashMap<usize, f64> = k0
                .iter()
                .enumerate()
                .map(|(i, &k)| (i, occupancy_from_ln_k(mid + k)))
                .collect();
            let eff = covered_mean_p_specific(&pl, &theta, *sum_rule);
            if eff < *target { lo = mid } else { hi = mid }
            best = REAL
                .iter()
                .zip(&loci_idx)
                .map(|(&(_, st, en, _), idx)| {
                    locus_retention(&pl, idx, &theta, st, en, *sum_rule)
                })
                .collect();
        }
        let base = *best.last().unwrap();
        let pred: Vec<f64> = best.iter().map(|v| if base > 0.0 { v / base } else { 0.0 }).collect();
        let err: f64 = pred
            .iter()
            .zip(&real)
            .map(|(p, r)| (p.max(1e-9).ln() - r.ln()).powi(2))
            .sum::<f64>()
            .sqrt();
        println!(
            "t={target:<5} {label:26}  rho={:+.3}  logRMSE={:6.2}   nifS={:.3} clpA={:.3} clpX={:.3} uvrA={:.3} recG={:.3} 16S={:.3}",
            spearman(&pred, &real), err, pred[0], pred[1], pred[2], pred[3], pred[5], pred[6]
        );
    }
    }
}
