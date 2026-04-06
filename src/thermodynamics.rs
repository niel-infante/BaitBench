/// Nearest-neighbor thermodynamic model for DNA hybridization.
///
/// Parameters from SantaLucia (1998) "A unified view of polymer, dumbbell, and
/// oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95(4):1460–1465.
///
/// Nucleotide index convention: A=0, C=1, G=2, T=3

/// Gas constant in kcal/(mol·K)
const R_KCAL: f64 = 1.987_204_258_64e-3;

#[derive(Clone, Copy)]
struct NnParams {
    dh: f64, // kcal/mol
    ds: f64, // cal/(K·mol)
}

impl NnParams {
    const fn new(dh: f64, ds: f64) -> Self {
        NnParams { dh, ds }
    }
}

/// Nearest-neighbor table indexed by [top_base][next_top_base], both Watson-Crick paired.
/// Row = 5'→3' top strand nucleotide; column = the next 5'→3' top strand nucleotide.
/// Values represent the stacking contribution of the dinucleotide step:
///   5'-XY-3' / 3'-X'Y'-5'
/// where X' and Y' are the Watson-Crick complements of X and Y.
///
/// Index: A=0, C=1, G=2, T=3
/// Ordering matches SantaLucia (1998) Table 2.
///
/// Missing entries (non-WC contexts) are never accessed because the
/// delta_g() function only accumulates steps where both the current
/// and previous pair are Watson-Crick.
// SantaLucia (1998) Table 2, with reverse-complement symmetry applied to fill
// all 16 entries.  Index: A=0, C=1, G=2, T=3.
// Complement map: A(0)↔T(3), C(1)↔G(2)  (i.e. comp = idx ^ 0b11).
// RC symmetry: ΔG(XY) = ΔG(comp(Y), comp(X)).
// Explicit entries from Table 2 (10 independent parameters):
//   AA(-7.9,-22.2)  AT(-7.2,-20.4)  TA(-7.2,-21.3)  CA(-8.5,-22.7)
//   GT(-8.4,-22.4)  CT(-7.8,-21.0)  GA(-8.2,-22.2)  CG(-10.6,-27.2)
//   GC(-9.8,-24.4)  GG(-8.0,-19.9)
// Derived via RC symmetry (6 remaining):
//   AC = RC(GT) = (-8.4,-22.4)    AG = RC(CT) = (-7.8,-21.0)
//   CC = RC(GG) = (-8.0,-19.9)    TC = RC(GA) = (-8.2,-22.2)
//   TG = RC(CA) = (-8.5,-22.7)    TT = RC(AA) = (-7.9,-22.2)
static NN_TABLE: [[NnParams; 4]; 4] = [
    //                A                     C                     G                     T
    /* prev A */ [NnParams::new(-7.9, -22.2), NnParams::new(-8.4, -22.4), NnParams::new(-7.8, -21.0), NnParams::new(-7.2, -20.4)],
    /* prev C */ [NnParams::new(-8.5, -22.7), NnParams::new(-8.0, -19.9), NnParams::new(-10.6, -27.2), NnParams::new(-7.8, -21.0)],
    /* prev G */ [NnParams::new(-8.2, -22.2), NnParams::new(-9.8, -24.4), NnParams::new(-8.0, -19.9), NnParams::new(-8.4, -22.4)],
    /* prev T */ [NnParams::new(-7.2, -21.3), NnParams::new(-8.2, -22.2), NnParams::new(-8.5, -22.7), NnParams::new(-7.9, -22.2)],
];

/// Convert an ASCII nucleotide byte to an index (A=0, C=1, G=2, T=3).
/// Returns 4 for any non-ACGT character.
#[inline]
fn nt_to_idx(b: u8) -> u8 {
    match b | 0x20 { // to lowercase
        b'a' => 0,
        b'c' => 1,
        b'g' => 2,
        b't' => 3,
        _ => 4,
    }
}

/// Returns the Watson-Crick complement index for a base index (0..=3).
/// A(0)↔T(3), C(1)↔G(2).
#[inline]
fn wc_complement(idx: u8) -> u8 {
    idx ^ 0b11
}

/// Returns true if (top, bot) are a Watson-Crick base pair (A-T or C-G),
/// expressed as nucleotide indices.
#[inline]
fn is_wc(top: u8, bot: u8) -> bool {
    top < 4 && bot < 4 && bot == wc_complement(top)
}

/// Compute the Gibbs free energy (ΔG, kcal/mol) for an aligned probe-reference pair.
///
/// `aligned_pairs` is a slice of `(probe_base, ref_base)` ASCII bytes, in 5'→3' order
/// of the probe strand. Both bases should be ACGT (upper or lower case).
///
/// The SkipStacking strategy is used: stacking energy accumulates only for consecutive
/// Watson-Crick steps. A mismatch or non-standard base breaks the stacking chain.
///
/// Returns ΔG in kcal/mol. More negative = more stable (higher affinity).
pub fn delta_g(aligned_pairs: &[(u8, u8)], temp_c: f64) -> f64 {
    let t_k = temp_c + 273.15;
    let mut sum_dh = 0.0f64;
    let mut sum_ds = 0.0f64;

    let mut prev_top: Option<u8> = None;
    let mut prev_bot: Option<u8> = None;

    for &(probe_b, ref_b) in aligned_pairs {
        let top = nt_to_idx(probe_b);
        let bot = nt_to_idx(ref_b);

        if let (Some(pt), Some(pb)) = (prev_top, prev_bot) {
            let prev_wc = is_wc(pt, pb);
            let curr_wc = is_wc(top, bot);

            if prev_wc && curr_wc {
                // Both ends of this step are WC pairs — accumulate stacking energy
                let params = NN_TABLE[pt as usize][top as usize];
                sum_dh += params.dh;
                sum_ds += params.ds;
            }
            // If either is a mismatch, skip this step (SkipStacking)
        }

        prev_top = Some(top);
        prev_bot = Some(bot);
    }

    // ΔG = ΔH - T × (ΔS/1000)   [convert cal → kcal for ΔS]
    sum_dh - t_k * (sum_ds / 1000.0)
}

/// Compute the Boltzmann binding score for a probe-reference interaction.
///
/// `dg` is the free energy in kcal/mol (from `delta_g()`).
/// Returns `exp(-ΔG / (R × T))`. Higher score = stronger binding.
pub fn boltzmann_score(dg: f64, temp_c: f64) -> f64 {
    let t_k = temp_c + 273.15;
    (-dg / (R_KCAL * t_k)).exp()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn perfect_match_is_negative_dg() {
        // A fully matched 10-mer probe should yield negative ΔG
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg = delta_g(&pairs, 70.0);
        assert!(dg < 0.0, "Perfect match ΔG should be negative, got {}", dg);
    }

    #[test]
    fn all_mismatches_gives_zero_energy() {
        // All mismatches → no stacking → ΔG = 0
        let pairs: Vec<(u8, u8)> = b"AAAA".iter().zip(b"AAAA".iter()).map(|(&p, &r)| (p, r)).collect();
        let dg = delta_g(&pairs, 70.0);
        assert_eq!(dg, 0.0);
    }

    #[test]
    fn boltzmann_perfect_match_gt_one() {
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg = delta_g(&pairs, 70.0);
        let score = boltzmann_score(dg, 70.0);
        assert!(score > 1.0, "Boltzmann score for negative ΔG should be > 1, got {}", score);
    }

    #[test]
    fn single_pair_no_stacking() {
        // A single base pair has no dinucleotide step → ΔG = 0
        let pairs = vec![(b'A', b'T')];
        assert_eq!(delta_g(&pairs, 70.0), 0.0);
    }
}
