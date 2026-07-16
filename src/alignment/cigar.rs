/// Shared CIGAR parsing utilities.
///
/// Two representations are provided:
/// - `expand_cigar`: per-position expansion (one `CigarOp` per base) — used by
///   thermodynamic scoring which needs to inspect each position individually.
/// - `parse_cigar`: run-length encoded `(len, CigarOp)` pairs — used by coverage
///   computation which advances reference position by `len` at a time.
/// - `cigar_ref_len`: total reference bases consumed — used when only the span is needed.

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CigarOp {
    Match,
    Ins,
    Del,
    Skip,
    SoftClip,
    HardClip,
    Pad,
    Equal,
    Mismatch,
}

impl CigarOp {
    pub fn from_char(c: char) -> Option<Self> {
        match c {
            'M' => Some(Self::Match),
            'I' => Some(Self::Ins),
            'D' => Some(Self::Del),
            'N' => Some(Self::Skip),
            'S' => Some(Self::SoftClip),
            'H' => Some(Self::HardClip),
            'P' => Some(Self::Pad),
            '=' => Some(Self::Equal),
            'X' => Some(Self::Mismatch),
            _ => None,
        }
    }

    /// True for operations that consume reference bases.
    pub fn consumes_ref(self) -> bool {
        matches!(self, Self::Match | Self::Del | Self::Skip | Self::Equal | Self::Mismatch)
    }
}

/// Parse a CIGAR string into run-length encoded `(length, op)` pairs.
pub fn parse_cigar(cigar: &str) -> Vec<(u32, CigarOp)> {
    let mut ops = Vec::new();
    let mut num_start = 0;
    for (i, c) in cigar.char_indices() {
        if !c.is_ascii_digit() {
            if let (Ok(len), Some(op)) = (cigar[num_start..i].parse::<u32>(), CigarOp::from_char(c)) {
                ops.push((len, op));
            }
            num_start = i + c.len_utf8();
        }
    }
    ops
}

/// Expand a CIGAR string into one `CigarOp` per position.
pub fn expand_cigar(cigar: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num = String::new();
    for ch in cigar.chars() {
        if ch.is_ascii_digit() {
            num.push(ch);
        } else if let Some(op) = CigarOp::from_char(ch) {
            let n: usize = num.parse().unwrap_or(1);
            num.clear();
            for _ in 0..n {
                ops.push(op);
            }
        } else {
            num.clear();
        }
    }
    ops
}

/// Total reference bases consumed by a CIGAR string.
pub fn cigar_ref_len(cigar: &str) -> usize {
    parse_cigar(cigar)
        .iter()
        .filter(|(_, op)| op.consumes_ref())
        .map(|&(len, _)| len as usize)
        .sum()
}