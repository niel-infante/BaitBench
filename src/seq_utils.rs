/// Shared DNA sequence utilities used by probe design algorithms.

/// Complement of a single DNA base (case-insensitive). Non-ACGT → N.
#[inline]
pub fn complement(b: u8) -> u8 {
    match b {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Reverse complement of a DNA byte sequence.
pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| complement(b)).collect()
}

/// Hamming distance where N on either side always counts as a mismatch.
pub fn hamming_n_mismatch(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b.iter())
        .filter(|&(&x, &y)| x == b'N' || y == b'N' || x != y)
        .count()
}