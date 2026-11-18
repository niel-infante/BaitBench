//! Nearest-neighbor thermodynamic model for DNA hybridization.
//!
//! Parameters from SantaLucia (1998) "A unified view of polymer, dumbbell, and
//! oligonucleotide DNA nearest-neighbor thermodynamics", PNAS 95(4):1460–1465.
//!
//! Salt correction from Owczarzy et al. (1997):
//!   ΔS([Na+]) = ΔS(1 M) + 0.368 × (n_bp - 1) × ln([Na+])
//!
//! Nucleotide index convention: A=0, C=1, G=2, T=3
//!
//! # Design notes
//!
//! `boltzmann_score` is used as a relative sampling weight — only ratios between
//! scores matter.  Score = 1 (ΔG = 0) is the neutral baseline; > 1 means
//! favorable binding; < 1 is never produced by this model (see initiation note).
//!
//! ## Initiation terms
//!
//! SantaLucia (1998) initiation parameters capture the cost of nucleating a new
//! duplex from two free strands: loss of translational/rotational entropy on
//! association, plus the instability of terminal base pairs (only partially
//! stacked, fewer constraints than interior pairs).  AT termini are more costly
//! than GC because A-T pairs have two hydrogen bonds vs three for G-C.
//!
//! Initiation is applied **only when at least one NN stacking step exists**
//! (≥ 2 consecutive WC pairs).  A single isolated WC pair flanked by mismatches
//! cannot sustain a stable duplex; applying the nucleation penalty to it would
//! give ΔG > 0 and score < 1 — implying the probe is actively repelled, which
//! is unphysical.  When no stacking exists, ΔG = 0 and score = 1 (neutral).
//!
//! ## Comparison with RAmpSim
//!
//! RAmpSim (Rooney et al. 2025, bioRxiv 2025.12.05.692407) uses the same
//! SantaLucia stacking table and SkipStacking strategy but omits initiation and
//! salt correction.  For relative scoring, the constant part of initiation
//! cancels across probes.  BaitBench preserves the AT vs GC terminal distinction
//! (~2 kcal/mol effect) because it is physically real and affects relative
//! capture probabilities.  The salt correction has a larger impact at typical
//! hybridization conditions (50 mM Na+) and meaningfully re-ranks long vs short
//! probe hits.

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

/// Internal-mismatch nearest-neighbour parameters.
///
/// Sources (as compiled in SantaLucia & Hicks 2004, *Annu Rev Biophys Biomol Struct*
/// 33:415-440, and transcribed in biopython's `MeltingTemp.DNA_IMM1`):
///   Allawi & SantaLucia (1997) *Biochemistry* 36:10581-10594  (G.T)
///   Allawi & SantaLucia (1998) *Biochemistry* 37:9435-9444    (G.A)
///   Allawi & SantaLucia (1998) *Biochemistry* 37:2170-2179    (C.T)
///   Allawi & SantaLucia (1998) *Nucleic Acids Res* 26:2694-2701 (A.C)
///   Peyret et al. (1999) *Biochemistry* 38:3468-3477          (A.A, C.C, G.G, T.T)
///
/// 51 steps: 48 where exactly one end of the dinucleotide is a mismatch, plus 3
/// published tandem (both-ends-mismatched) G.T steps. All 12 single-mismatch types are
/// covered (A/A A/C A/G C/A C/C C/T G/A G/G G/T T/C T/G T/T).
///
/// The key packs the dinucleotide step as four 2-bit base codes,
/// `[prev_top, prev_bot, curr_top, curr_bot]`, where "top" is the probe strand and
/// "bottom" the strand it hybridizes to — matching the `XY/WZ` convention of the
/// published tables, in which the bottom strand runs 3'->5'.
///
/// Before this table existed the code used "SkipStacking": any step touching a
/// mismatch contributed exactly zero. That treats a weakened stack as a broken one and
/// substantially overestimates the penalty for scattered mismatches — which matters
/// most for probes designed against relatives of the sampled organism, the case the
/// capture model is being rebuilt to handle.
#[inline]
fn mismatch_params(prev_top: u8, prev_bot: u8, curr_top: u8, curr_bot: u8) -> Option<NnParams> {
    if prev_top > 3 || prev_bot > 3 || curr_top > 3 || curr_bot > 3 {
        return None; // N or any other ambiguous base: no published parameters
    }
    let key = (prev_top << 6) | (prev_bot << 4) | (curr_top << 2) | curr_bot;
    match key {
        0b00110000 => Some(NnParams::new(1.2, 1.7)), // AA/TA
        0b00110001 => Some(NnParams::new(2.3, 4.6)), // AA/TC
        0b00110010 => Some(NnParams::new(-0.6, -2.3)), // AA/TG
        0b00110100 => Some(NnParams::new(5.3, 14.6)), // AC/TA
        0b00110101 => Some(NnParams::new(0.0, -4.4)), // AC/TC
        0b00110111 => Some(NnParams::new(0.7, 0.2)), // AC/TT
        0b00111000 => Some(NnParams::new(-0.7, -2.3)), // AG/TA
        0b00111010 => Some(NnParams::new(-3.1, -9.5)), // AG/TG
        0b00111011 => Some(NnParams::new(1.0, 0.9)), // AG/TT
        0b00111101 => Some(NnParams::new(-1.2, -6.2)), // AT/TC
        0b00111110 => Some(NnParams::new(-2.5, -8.3)), // AT/TG
        0b00111111 => Some(NnParams::new(-2.7, -10.8)), // AT/TT
        0b01100000 => Some(NnParams::new(-0.9, -4.2)), // CA/GA
        0b01100001 => Some(NnParams::new(1.9, 3.7)), // CA/GC
        0b01100010 => Some(NnParams::new(-0.7, -2.3)), // CA/GG
        0b01100100 => Some(NnParams::new(0.6, -0.6)), // CC/GA
        0b01100101 => Some(NnParams::new(-1.5, -7.2)), // CC/GC
        0b01100111 => Some(NnParams::new(-0.8, -4.5)), // CC/GT
        0b01101000 => Some(NnParams::new(-4.0, -13.2)), // CG/GA
        0b01101010 => Some(NnParams::new(-4.9, -15.3)), // CG/GG
        0b01101011 => Some(NnParams::new(-4.1, -11.7)), // CG/GT
        0b01101101 => Some(NnParams::new(-1.5, -6.1)), // CT/GC
        0b01101110 => Some(NnParams::new(-2.8, -8.0)), // CT/GG
        0b01101111 => Some(NnParams::new(-5.0, -15.8)), // CT/GT
        0b10010000 => Some(NnParams::new(-2.9, -9.8)), // GA/CA
        0b10010001 => Some(NnParams::new(5.2, 14.2)), // GA/CC
        0b10010010 => Some(NnParams::new(-0.6, -1.0)), // GA/CG
        0b10010100 => Some(NnParams::new(-0.7, -3.8)), // GC/CA
        0b10010101 => Some(NnParams::new(3.6, 8.9)), // GC/CC
        0b10010111 => Some(NnParams::new(2.3, 5.4)), // GC/CT
        0b10011000 => Some(NnParams::new(0.5, 3.2)), // GG/CA
        0b10011010 => Some(NnParams::new(-6.0, -15.8)), // GG/CG
        0b10011011 => Some(NnParams::new(3.3, 10.4)), // GG/CT
        0b10011101 => Some(NnParams::new(5.2, 13.5)), // GT/CC
        0b10011110 => Some(NnParams::new(-4.4, -12.3)), // GT/CG
        0b10011111 => Some(NnParams::new(-2.2, -8.4)), // GT/CT
        0b10111011 => Some(NnParams::new(5.8, 16.3)), // GG/TT  // tandem
        0b10111110 => Some(NnParams::new(4.1, 9.5)), // GT/TG  // tandem
        0b11000000 => Some(NnParams::new(4.7, 12.9)), // TA/AA
        0b11000001 => Some(NnParams::new(3.4, 8.0)), // TA/AC
        0b11000010 => Some(NnParams::new(0.7, 0.7)), // TA/AG
        0b11000100 => Some(NnParams::new(7.6, 20.2)), // TC/AA
        0b11000101 => Some(NnParams::new(6.1, 16.4)), // TC/AC
        0b11000111 => Some(NnParams::new(1.2, 0.7)), // TC/AT
        0b11001000 => Some(NnParams::new(3.0, 7.4)), // TG/AA
        0b11001010 => Some(NnParams::new(1.6, 3.6)), // TG/AG
        0b11001011 => Some(NnParams::new(-0.1, -1.7)), // TG/AT
        0b11001101 => Some(NnParams::new(1.0, 0.7)), // TT/AC
        0b11001110 => Some(NnParams::new(-1.3, -5.3)), // TT/AG
        0b11001111 => Some(NnParams::new(0.2, -1.5)), // TT/AT
        0b11101011 => Some(NnParams::new(-1.4, -6.2)), // TG/GT  // tandem
        _ => None,
    }
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
///    ΔS_corrected = ΔS_1M + 0.368 × (n_wc - 1) × ln([Na+])
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
            } else if let Some(params) = mismatch_params(pt, pb, top, bot) {
                // One (or, for the three published tandem G.T steps, both) ends of the
                // step is a mismatch. A mismatch weakens the stack rather than breaking
                // it, so charge the published internal-mismatch parameters.
                //
                // `has_stacking` is deliberately NOT set here: mismatch stacking alone
                // cannot nucleate a duplex, so it must not unlock the initiation and
                // salt terms. Those still require at least one WC|WC step.
                sum_dh += params.dh;
                sum_ds += params.ds;
            }
            // Steps with no published parameters (ambiguous bases, unlisted tandem
            // mismatches) contribute nothing, as before.
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
        // Two isolated WC pairs with a mismatch between them. The mismatch steps now
        // carry published internal-mismatch parameters, so ΔG is no longer exactly
        // zero — but the invariant that matters is unchanged: without a WC|WC step
        // there is no nucleated duplex, so initiation and the salt correction must stay
        // locked, and the probe must not look like a binder.
        let pairs = vec![(b'A', b'T'), (b'A', b'A'), (b'G', b'C')]; // AT, mismatch, GC

        let dg_1m = delta_g(&pairs, &model_1m());
        let dg_50mm = delta_g(&pairs, &model_50mm());
        assert_eq!(
            dg_1m, dg_50mm,
            "salt correction must not apply without a WC|WC stacking step"
        );

        assert!(
            dg_1m > 0.0,
            "isolated WC pairs flanked by a mismatch must not be a stable duplex, got {dg_1m}"
        );
        let score = boltzmann_score(dg_50mm, &model_50mm());
        assert!(
            (score - 1.0).abs() < 1e-9,
            "an unstable duplex must clamp to the neutral score, got {score}"
        );
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
    /// The Stage 2.5 gate: a probe with one internal mismatch must sit strictly between
    /// a perfect match and a fully mismatched probe. Under SkipStacking the mismatch
    /// deleted both adjacent stacking steps and charged nothing for them, overstating
    /// the penalty.
    fn one_internal_mismatch_scores_between_perfect_and_none() {
        let m = model_50mm();
        let probe = b"ACGTACGTACGTACGTACGT";

        let perfect: Vec<(u8, u8)> = probe
            .iter()
            .map(|&b| (b, super::wc_complement(super::nt_to_idx(b))))
            .map(|(b, c)| (b, b"ACGT"[c as usize]))
            .collect();

        let mut one_mm = perfect.clone();
        one_mm[10].1 = one_mm[10].0; // self-pair: a mismatch at a central position

        let all_mm: Vec<(u8, u8)> = probe.iter().map(|&b| (b, b)).collect();

        let dg_perfect = delta_g(&perfect, &m);
        let dg_one = delta_g(&one_mm, &m);
        let dg_none = delta_g(&all_mm, &m);

        assert!(dg_perfect < dg_one, "one mismatch must be less stable than none");
        assert!(dg_one < dg_none, "one mismatch must be more stable than all mismatches");
    }

    /// Mismatch steps must contribute energy rather than being skipped, but must not
    /// unlock initiation or the salt correction on their own.
    #[test]
    fn mismatch_steps_contribute_without_nucleating() {
        let m = model_1m();
        // A single WC|mismatch step: no WC|WC step anywhere, so has_stacking stays false.
        let pairs = vec![(b'A', b'T'), (b'G', b'T')]; // AG/TT is in the table
        let dg = delta_g(&pairs, &m);
        assert_ne!(dg, 0.0, "a published mismatch step must contribute energy");
        assert_eq!(
            dg,
            delta_g(&pairs, &model_50mm()),
            "salt must stay locked without a WC|WC step"
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
