/// Nearest-neighbor thermodynamic model for DNA hybridization.
///
/// Parameters from SantaLucia (1998) "A unified view of polymer, dumbbell, and
/// oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95(4):1460–1465.
///
/// Salt correction from Owczarzy et al. (1997):
///   ΔS([Na+]) = ΔS(1 M) + 0.368 × (n_bp - 1) × ln([Na+])
///
/// Nucleotide index convention: A=0, C=1, G=2, T=3
///
/// # Design notes
///
/// `boltzmann_score` is used as a relative sampling weight — only ratios between
/// scores matter.  Score = 1 (ΔG = 0) is the neutral baseline; > 1 means
/// favorable binding; < 1 is never produced by this model (see initiation note).
///
/// ## Initiation terms
///
/// SantaLucia (1998) initiation parameters capture the cost of nucleating a new
/// duplex from two free strands: loss of translational/rotational entropy on
/// association, plus the instability of terminal base pairs (only partially
/// stacked, fewer constraints than interior pairs).  AT termini are more costly
/// than GC because A-T pairs have two hydrogen bonds vs three for G-C.
///
/// Initiation is applied **only when at least one NN stacking step exists**
/// (≥ 2 consecutive WC pairs).  A single isolated WC pair flanked by mismatches
/// cannot sustain a stable duplex; applying the nucleation penalty to it would
/// give ΔG > 0 and score < 1 — implying the probe is actively repelled, which
/// is unphysical.  When no stacking exists, ΔG = 0 and score = 1 (neutral).
///
/// ## Comparison with RAmpSim
///
/// RAmpSim (Rooney et al. 2025, bioRxiv 2025.12.05.692407) uses the same
/// SantaLucia stacking table and SkipStacking strategy but omits initiation and
/// salt correction.  For relative scoring, the constant part of initiation
/// cancels across probes.  BaitBench preserves the AT vs GC terminal distinction
/// (~2 kcal/mol effect) because it is physically real and affects relative
/// capture probabilities.  The salt correction has a larger impact at typical
/// hybridization conditions (50 mM Na+) and meaningfully re-ranks long vs short
/// probe hits.

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

/// SantaLucia (1998) Table 2: initiation parameters.
/// Applied once for each terminal Watson-Crick pair of the duplex.
/// A duplex has two terminals (5' end and 3' end), so two initiation terms are added.
const INIT_GC: NnParams = NnParams::new(0.1, -2.8);  // terminal G-C or C-G pair
const INIT_AT: NnParams = NnParams::new(2.3, 4.1);   // terminal A-T or T-A pair

/// All thermodynamic model parameters for ΔG calculation.
///
/// Construct with `ThermoModel::new(temp_c, na_conc_m)`.
/// Both `delta_g()` and `boltzmann_score()` operate via this struct.
#[derive(Clone, Copy, Debug)]
pub struct ThermoModel {
    /// Hybridization temperature in °C
    pub temp_c: f64,
    /// Na+ concentration in molar (e.g. 0.05 for 50 mM).
    /// SantaLucia (1998) tables assume 1.0 M; values other than 1.0 apply
    /// the Owczarzy et al. (1997) entropy correction.
    pub na_conc_m: f64,
}

impl ThermoModel {
    pub fn new(temp_c: f64, na_conc_m: f64) -> Self {
        ThermoModel { temp_c, na_conc_m }
    }
}

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

/// Returns the initiation `NnParams` for a terminal base pair, or `None` if
/// the pair is not Watson-Crick.
///
/// A-T and T-A pairs use `INIT_AT`; G-C and C-G pairs use `INIT_GC`.
#[inline]
fn initiation_params(top: u8, bot: u8) -> Option<NnParams> {
    if !is_wc(top, bot) {
        return None;
    }
    // A=0, T=3 → AT pair; C=1, G=2 → GC pair
    if top == 0 || top == 3 {
        Some(INIT_AT)
    } else {
        Some(INIT_GC)
    }
}

/// Compute the Gibbs free energy (ΔG, kcal/mol) for an aligned probe-reference pair.
///
/// `aligned_pairs` is a slice of `(probe_base, ref_base)` ASCII bytes, in 5'→3' order
/// of the probe strand. Both bases should be ACGT (upper or lower case).
///
/// Three contributions are summed:
///
/// 1. **NN stacking** (SantaLucia 1998 Table 2): accumulated for each consecutive
///    Watson-Crick dinucleotide step. A mismatch or non-standard base breaks the
///    stacking chain (SkipStacking strategy).
///
/// 2. **Initiation terms** (SantaLucia 1998 Table 2): added once per terminal
///    Watson-Crick pair. The first and last WC pairs of the aligned duplex are
///    identified; each contributes an AT or GC initiation penalty to ΔH and ΔS.
///
/// 3. **Salt correction** (Owczarzy et al. 1997): adjusts ΔS for the actual Na+
///    concentration relative to the 1 M reference used to derive the NN table:
///      ΔS_corrected = ΔS_1M + 0.368 × (n_wc - 1) × ln([Na+])
///    where n_wc is the count of WC-paired positions and [Na+] is in molar.
///    At 1 M Na+, ln(1) = 0 so no correction is applied.
///
/// Returns ΔG in kcal/mol. More negative = more stable (higher affinity).
pub fn delta_g(aligned_pairs: &[(u8, u8)], model: &ThermoModel) -> f64 {
    let t_k = model.temp_c + 273.15;
    let mut sum_dh = 0.0f64;
    let mut sum_ds = 0.0f64;

    // --- 1. NN stacking (SkipStacking) ---
    let mut prev_top: Option<u8> = None;
    let mut prev_bot: Option<u8> = None;
    let mut has_stacking = false;

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
                has_stacking = true;
            }
            // If either is a mismatch, skip this step (SkipStacking)
        }

        prev_top = Some(top);
        prev_bot = Some(bot);
    }

    // --- 2. Initiation terms: first and last WC pair ---
    // Only applied when at least one stacking step exists (≥ 2 consecutive WC pairs).
    // An isolated single WC pair flanked by mismatches cannot form a stable duplex,
    // so paying the initiation cost for it would produce ΔG > 0 (score < 1) — implying
    // the probe is actively repelled, which is unphysical. RAmpSim omits initiation
    // entirely; we apply it only when there is genuine duplex stacking.
    if has_stacking {
        // Forward scan: first WC pair
        for &(probe_b, ref_b) in aligned_pairs.iter() {
            let top = nt_to_idx(probe_b);
            let bot = nt_to_idx(ref_b);
            if let Some(p) = initiation_params(top, bot) {
                sum_dh += p.dh;
                sum_ds += p.ds;
                break;
            }
        }

        // Reverse scan: last WC pair
        for &(probe_b, ref_b) in aligned_pairs.iter().rev() {
            let top = nt_to_idx(probe_b);
            let bot = nt_to_idx(ref_b);
            if let Some(p) = initiation_params(top, bot) {
                sum_dh += p.dh;
                sum_ds += p.ds;
                break;
            }
        }
    }

    // --- 3. Salt correction on ΔS (Owczarzy et al. 1997) ---
    // ΔS([Na+]) = ΔS(1M) + 0.368 × (n_wc - 1) × ln([Na+])
    // 0.368 is in cal/(mol·K) per base pair — same units as sum_ds
    // Only applied when a duplex actually exists (same guard as initiation).
    // Without this guard, isolated WC pairs at low salt would make sum_ds negative
    // with sum_dh = 0, producing ΔG > 0 and score < 1 — unphysical.
    if has_stacking {
        let n_wc = aligned_pairs
            .iter()
            .filter(|&&(p, r)| is_wc(nt_to_idx(p), nt_to_idx(r)))
            .count();

        if n_wc >= 1 {
            let salt_correction = 0.368 * (n_wc as f64 - 1.0) * model.na_conc_m.ln();
            sum_ds += salt_correction;
        }
    }

    // ΔG = ΔH - T × (ΔS/1000)   [convert cal → kcal for ΔS]
    sum_dh - t_k * (sum_ds / 1000.0)
}

/// Compute the Boltzmann binding score for a probe-reference interaction.
///
/// `dg` is the free energy in kcal/mol (from `delta_g()`).
/// Returns `exp(-ΔG / (R × T))`, clamped to a minimum of 1.0.
///
/// The clamp enforces that probes can only *enrich* fragments relative to
/// background — never deplete them. Score = 1 is the neutral baseline (no
/// preferential binding); score > 1 means favorable hybridization. A raw
/// ΔG > 0 indicates the duplex is thermodynamically too weak to be stable at
/// the given temperature and salt conditions (e.g. a single stacking step at
/// high temperature with costly AT initiation), and is treated the same as
/// no binding.
pub fn boltzmann_score(dg: f64, model: &ThermoModel) -> f64 {
    let t_k = model.temp_c + 273.15;
    (-dg / (R_KCAL * t_k)).exp().max(1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper: model at 70°C / 1 M Na+ (no salt correction; matches original behavior).
    fn model_1m() -> ThermoModel {
        ThermoModel::new(70.0, 1.0)
    }

    /// Helper: model at 70°C / 50 mM Na+ (typical hybridization buffer).
    fn model_50mm() -> ThermoModel {
        ThermoModel::new(70.0, 0.05)
    }

    #[test]
    fn perfect_match_is_negative_dg() {
        // A fully matched 10-mer probe should yield negative ΔG (at 1 M Na+)
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg = delta_g(&pairs, &model_1m());
        assert!(dg < 0.0, "Perfect match ΔG should be negative, got {}", dg);
    }

    #[test]
    fn all_mismatches_gives_zero_energy() {
        // All mismatches → no WC pairs → no stacking, no initiation, no salt correction → ΔG = 0
        let pairs: Vec<(u8, u8)> = b"AAAA".iter().zip(b"AAAA".iter()).map(|(&p, &r)| (p, r)).collect();
        let dg = delta_g(&pairs, &model_1m());
        assert_eq!(dg, 0.0);
    }

    #[test]
    fn boltzmann_perfect_match_gt_one() {
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg = delta_g(&pairs, &model_1m());
        let score = boltzmann_score(dg, &model_1m());
        assert!(score > 1.0, "Boltzmann score for negative ΔG should be > 1, got {}", score);
    }

    #[test]
    fn single_at_pair_no_stacking_no_initiation() {
        // A single A-T pair has no stacking (needs ≥ 2 consecutive WC pairs)
        // and therefore no initiation terms either. ΔG = 0, score = 1.
        let pairs = vec![(b'A', b'T')];
        let dg = delta_g(&pairs, &model_1m());
        assert_eq!(dg, 0.0, "Single AT pair with no stacking should give ΔG = 0, got {}", dg);
    }

    #[test]
    fn single_gc_pair_no_stacking_no_initiation() {
        // A single G-C pair: same reasoning — no stacking, no initiation → ΔG = 0.
        let pairs = vec![(b'G', b'C')];
        let dg = delta_g(&pairs, &model_1m());
        assert_eq!(dg, 0.0, "Single GC pair with no stacking should give ΔG = 0, got {}", dg);
    }

    #[test]
    fn boltzmann_score_never_below_one() {
        // A very weak duplex (single AT stacking step, AT terminals, 50 mM Na+) can
        // have ΔG > 0 because initiation + salt correction outweigh stacking energy.
        // boltzmann_score must clamp to 1.0 — probes can only enrich, not deplete.
        let pairs: Vec<(u8, u8)> = vec![(b'A', b'T'), (b'A', b'T')]; // single AT stack
        let model = ThermoModel::new(70.0, 0.05);
        let dg = delta_g(&pairs, &model);
        assert!(dg > 0.0, "Expected ΔG > 0 for weak duplex, got {}", dg);
        let score = boltzmann_score(dg, &model);
        assert!((score - 1.0).abs() < 1e-9, "Score should be clamped to 1.0, got {}", score);
    }

    #[test]
    fn isolated_wc_pairs_no_stacking_score_neutral_at_low_salt() {
        // Two isolated WC pairs with a mismatch between them — no consecutive stacking.
        // At 50 mM Na+, salt correction must NOT be applied (has_stacking = false),
        // otherwise sum_ds would go negative with sum_dh = 0, giving ΔG > 0 and score < 1.
        let pairs = vec![(b'A', b'T'), (b'A', b'A'), (b'G', b'C')]; // AT, mismatch, GC
        let dg = delta_g(&pairs, &model_50mm());
        assert_eq!(dg, 0.0, "Isolated WC pairs at 50 mM should give ΔG = 0 (no stacking), got {}", dg);
        let score = boltzmann_score(dg, &model_50mm());
        assert!((score - 1.0).abs() < 1e-9, "Score should be 1.0, got {}", score);
    }

    #[test]
    fn salt_correction_1m_no_change() {
        // At 1 M Na+, ln(1.0) = 0 → salt correction = 0
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg = delta_g(&pairs, &model_1m());
        assert!(dg.is_finite());
        assert!(dg < 0.0, "Perfect match ΔG at 1M Na+ should be negative");
    }

    #[test]
    fn salt_50mm_weaker_than_1m() {
        // At lower salt, hybridization is less favorable → ΔG more positive
        let pairs: Vec<(u8, u8)> = b"ACGTACGTAC"
            .iter()
            .zip(b"TGCATGCATG".iter())
            .map(|(&p, &r)| (p, r))
            .collect();
        let dg_1m = delta_g(&pairs, &model_1m());
        let dg_50mm = delta_g(&pairs, &model_50mm());
        assert!(
            dg_50mm > dg_1m,
            "ΔG at 50 mM ({:.4}) should be more positive than at 1 M ({:.4})",
            dg_50mm, dg_1m
        );
    }

    #[test]
    fn initiation_gc_terminal_two_base_duplex() {
        // GC + GC duplex: two GC initiation terms + one GC-GC stacking step
        // stacking: NN_TABLE[G=2][G=2] = (-8.0, -19.9)
        // initiation: 2 × INIT_GC = (0.2, -5.6)
        // salt at 1M: 0.368 * (2-1) * ln(1) = 0
        let pairs = vec![(b'G', b'C'), (b'G', b'C')];
        let dg = delta_g(&pairs, &model_1m());
        let t_k = 70.0 + 273.15;
        let expected = (-8.0 + 0.2) - t_k * ((-19.9 + (-5.6)) / 1000.0);
        assert!(
            (dg - expected).abs() < 1e-9,
            "GC+GC duplex ΔG should be {:.6}, got {:.6}", expected, dg
        );
    }
}
