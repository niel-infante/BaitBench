//! sDUST algorithm for identifying low-complexity regions in DNA sequences.
//!
//! Implementation based on:
//! Morgulis A, Gertz EM, Schäffer AA, Agarwala R. "A Fast and Symmetric DUST
//! Implementation to Mask Low-Complexity DNA Sequences." J Comput Biol.
//! 2006;13(5):1028-1040. doi:10.1089/cmb.2006.13.1028

use std::collections::VecDeque;

const NUM_TRIPLETS: usize = 64;

/// Default score threshold (T parameter from the paper).
#[allow(dead_code)]
pub const DEFAULT_THRESHOLD: f64 = 2.0;

/// Default window size in bases (W parameter from the paper).
#[allow(dead_code)]
pub const DEFAULT_WINDOW: usize = 64;

#[derive(Clone, Debug)]
struct PerfectInterval {
    start: usize,  // start position in sequence (inclusive)
    finish: usize, // end position in sequence (inclusive)
    score: f64,
}

/// Encode a DNA base as 0–3. Returns None for non-ACGT bases.
#[inline]
fn encode_base(b: u8) -> Option<u8> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Algorithm S1: Update triplet counts when a sequence changes by one base.
// ---------------------------------------------------------------------------

/// Add a triplet: running count increases by the old count, then count increments.
#[inline]
fn add_triplet_info(r: &mut u32, c: &mut [u32; NUM_TRIPLETS], t: u8) {
    let idx = t as usize;
    *r += c[idx];
    c[idx] += 1;
}

/// Remove a triplet: count decrements, then running count decreases by the new count.
#[inline]
fn rem_triplet_info(r: &mut u32, c: &mut [u32; NUM_TRIPLETS], t: u8) {
    let idx = t as usize;
    c[idx] -= 1;
    *r -= c[idx];
}

// ---------------------------------------------------------------------------
// Algorithm S2: Shift the sliding window by one triplet.
// ---------------------------------------------------------------------------

fn shift_window(
    t: u8,
    w: &mut VecDeque<u8>,
    threshold: f64,
    window_size: usize,
    l: &mut usize,
    rw: &mut u32,
    rv: &mut u32,
    cw: &mut [u32; NUM_TRIPLETS],
    cv: &mut [u32; NUM_TRIPLETS],
) {
    let max_triplets = window_size - 2;

    // If the window is full, pop the leftmost triplet.
    if w.len() >= max_triplets {
        let s = w.pop_front().unwrap();
        rem_triplet_info(rw, cw, s);
        if *l > w.len() {
            *l -= 1;
            rem_triplet_info(rv, cv, s);
        }
    }

    // Append the new triplet to the right.
    w.push_back(t);
    *l += 1;
    add_triplet_info(rw, cw, t);
    add_triplet_info(rv, cv, t);

    // If the triplet count in suffix v exceeds 2T, shrink v from the left
    // until the offending triplet is removed (Proposition 1).
    let two_t = (2.0 * threshold) as u32;
    if cv[t as usize] > two_t {
        loop {
            let s = w[w.len() - *l];
            rem_triplet_info(rv, cv, s);
            *l -= 1;
            if s == t {
                break;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Algorithm S3: Save masked regions for intervals leaving the window.
// ---------------------------------------------------------------------------

/// Remove perfect intervals whose start falls before `wstart` and record
/// their regions in `res`, merging adjacent/overlapping entries.
///
/// `perfect` is sorted by descending start, so items at the back have the
/// smallest start and are the first to fall off the window's left edge.
fn save_masked_regions(
    res: &mut Vec<(usize, usize)>,
    perfect: &mut Vec<PerfectInterval>,
    wstart: usize,
) {
    while perfect.last().is_some_and(|p| p.start < wstart) {
        let p = perfect.pop().unwrap();
        if let Some(prev) = res.last_mut() {
            if p.start <= prev.1 + 1 {
                // Adjacent or overlapping — extend.
                prev.1 = prev.1.max(p.finish);
            } else {
                res.push((p.start, p.finish));
            }
        } else {
            res.push((p.start, p.finish));
        }
    }
}

// ---------------------------------------------------------------------------
// Algorithm S4: Find perfect intervals that are suffixes of the window.
// ---------------------------------------------------------------------------

fn find_perfect(
    perfect: &mut Vec<PerfectInterval>,
    w: &VecDeque<u8>,
    threshold: f64,
    wstart: usize,
    l: usize,
    rv: u32,
    cv: &[u32; NUM_TRIPLETS],
) {
    let wlen = w.len();
    if wlen <= l {
        return;
    }

    let mut c = *cv; // work on a copy of the suffix-v counts
    let mut r = rv;
    let mut perf_idx: usize = 0;
    let mut max_score: f64 = 0.0;

    // Iterate suffixes from shortest (just beyond v) to longest (entire window).
    for i in (0..=(wlen - l - 1)).rev() {
        let t = w[i];
        add_triplet_info(&mut r, &mut c, t);

        let denom = wlen - i - 1; // ℓ − 1
        if denom == 0 {
            continue; // single-triplet suffix, score is 0 by definition
        }
        let new_score = r as f64 / denom as f64;

        if new_score > threshold {
            // Advance past existing intervals with start ≥ i + wstart.
            while perf_idx < perfect.len() && perfect[perf_idx].start >= i + wstart {
                if perfect[perf_idx].score > max_score {
                    max_score = perfect[perf_idx].score;
                }
                perf_idx += 1;
            }

            if new_score >= max_score {
                max_score = new_score;
                let new_perf = PerfectInterval {
                    start: i + wstart,
                    finish: wlen + 1 + wstart,
                    score: new_score,
                };
                perfect.insert(perf_idx, new_perf);
                // The element formerly at perf_idx shifted to perf_idx + 1.
                perf_idx += 1;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Flush helper: move all remaining perfect intervals into results.
// ---------------------------------------------------------------------------

fn flush_remaining(res: &mut Vec<(usize, usize)>, perfect: &mut Vec<PerfectInterval>) {
    while let Some(p) = perfect.pop() {
        if let Some(prev) = res.last_mut() {
            if p.start <= prev.1 + 1 {
                prev.1 = prev.1.max(p.finish);
            } else {
                res.push((p.start, p.finish));
            }
        } else {
            res.push((p.start, p.finish));
        }
    }
}

// ---------------------------------------------------------------------------
// Algorithm S5: Top-level sDUST function.
// ---------------------------------------------------------------------------

/// Run the sDUST algorithm on a DNA sequence.
///
/// Returns low-complexity regions as `(start, end)` pairs (0-indexed, start
/// inclusive, end exclusive).
///
/// # Parameters
/// - `seq`: DNA sequence (ACGT, case-insensitive). Non-ACGT bases reset the
///   window (treated as sequence breaks).
/// - `threshold`: Score threshold *T* (default [`DEFAULT_THRESHOLD`] = 2.0).
/// - `window_size`: Window size *W* in bases (default [`DEFAULT_WINDOW`] = 64).
pub fn sdust(seq: &[u8], threshold: f64, window_size: usize) -> Vec<(usize, usize)> {
    let n = seq.len();
    if n < 3 {
        return Vec::new();
    }

    let mut res: Vec<(usize, usize)> = Vec::new(); // inclusive on both ends internally
    let mut perfect: Vec<PerfectInterval> = Vec::new();
    let mut rv: u32 = 0;
    let mut rw: u32 = 0;
    let mut cv = [0u32; NUM_TRIPLETS];
    let mut cw = [0u32; NUM_TRIPLETS];
    let mut w: VecDeque<u8> = VecDeque::new();
    let mut l: usize = 0;

    for wfinish in 2..n {
        // Encode the three bases forming the new triplet.
        let b0 = encode_base(seq[wfinish - 2]);
        let b1 = encode_base(seq[wfinish - 1]);
        let b2 = encode_base(seq[wfinish]);

        if let (Some(b0), Some(b1), Some(b2)) = (b0, b1, b2) {
            let wstart = (wfinish + 1).saturating_sub(window_size);

            save_masked_regions(&mut res, &mut perfect, wstart);

            let t = (b0 << 4) | (b1 << 2) | b2;
            shift_window(
                t,
                &mut w,
                threshold,
                window_size,
                &mut l,
                &mut rw,
                &mut rv,
                &mut cw,
                &mut cv,
            );

            // Proposition 2: skip suffix search when window score is low.
            if rw as f64 > l as f64 * threshold {
                find_perfect(&mut perfect, &w, threshold, wstart, l, rv, &cv);
            }
        } else {
            // Non-ACGT base encountered — flush state and restart.
            flush_remaining(&mut res, &mut perfect);
            w.clear();
            cw = [0u32; NUM_TRIPLETS];
            cv = [0u32; NUM_TRIPLETS];
            rw = 0;
            rv = 0;
            l = 0;
        }
    }

    // Flush any remaining perfect intervals.
    flush_remaining(&mut res, &mut perfect);

    // Convert from inclusive-finish to exclusive-end (Rust convention).
    res.iter().map(|&(s, f)| (s, f + 1)).collect()
}

/// Fraction of bases in `seq` that sDUST identifies as low-complexity.
///
/// Returns a value in \[0.0, 1.0\].
pub fn masked_fraction(seq: &[u8], threshold: f64, window_size: usize) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let regions = sdust(seq, threshold, window_size);
    let masked_bases: usize = regions.iter().map(|&(s, e)| e - s).sum();
    masked_bases as f64 / seq.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty() {
        assert!(sdust(b"", DEFAULT_THRESHOLD, DEFAULT_WINDOW).is_empty());
    }

    #[test]
    fn test_too_short() {
        assert!(sdust(b"AC", DEFAULT_THRESHOLD, DEFAULT_WINDOW).is_empty());
    }

    #[test]
    fn test_poly_a_fully_masked() {
        let seq = b"AAAAAAAAAAAAAAAA"; // 16 A's
        let regions = sdust(seq, DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        let masked: usize = regions.iter().map(|&(s, e)| e - s).sum();
        assert_eq!(masked, 16);
    }

    #[test]
    fn test_complex_not_masked() {
        // ACGTACGTACGTACGT: repeating 4-mer, but DUST score < 2.
        //   14 triplets; ACG×4, CGT×4, GTA×3, TAC×3
        //   r = 6+6+3+3 = 18, S = 18/13 ≈ 1.38
        let seq = b"ACGTACGTACGTACGT";
        let regions = sdust(seq, DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        let masked: usize = regions.iter().map(|&(s, e)| e - s).sum();
        assert_eq!(masked, 0);
    }

    #[test]
    fn test_dinucleotide_repeat_masked() {
        // ATAT… has only two distinct triplets (ATA, TAT) each appearing many
        // times, producing a high DUST score.
        let seq = b"ATATATATATATATATATAT"; // 20 bases
        let regions = sdust(seq, DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        assert!(!regions.is_empty());
    }

    #[test]
    fn test_masked_fraction_poly_a() {
        let frac = masked_fraction(b"AAAAAAAAAA", DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        assert!(
            (frac - 1.0).abs() < 1e-6,
            "poly-A should be ~100% masked, got {}",
            frac
        );
    }

    #[test]
    fn test_masked_fraction_complex() {
        let frac = masked_fraction(b"ACGTACGTACGTACGT", DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        assert!(frac < 0.01, "complex seq should be ~0% masked, got {}", frac);
    }

    #[test]
    fn test_mixed_sequence() {
        // Complex prefix + poly-A suffix.
        let seq = b"ACGTACGTACGTACGTAAAAAAAAAAAAAAAA";
        let regions = sdust(seq, DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        // Only the poly-A portion should be masked.
        assert!(!regions.is_empty());
        let masked: usize = regions.iter().map(|&(s, e)| e - s).sum();
        assert!(masked < seq.len());
        assert!(masked > 0);
    }

    #[test]
    fn test_n_breaks_window() {
        // N in the middle should not cause a panic.
        let seq = b"AAAAAAAANACGTACGT";
        let regions = sdust(seq, DEFAULT_THRESHOLD, DEFAULT_WINDOW);
        // The poly-A portion before N should be masked.
        let _ = regions; // just ensure no panic
    }
}
