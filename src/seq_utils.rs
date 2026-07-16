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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn complement_all_bases() {
        assert_eq!(complement(b'A'), b'T');
        assert_eq!(complement(b'T'), b'A');
        assert_eq!(complement(b'C'), b'G');
        assert_eq!(complement(b'G'), b'C');
        assert_eq!(complement(b'a'), b'T');
        assert_eq!(complement(b't'), b'A');
        assert_eq!(complement(b'c'), b'G');
        assert_eq!(complement(b'g'), b'C');
    }

    #[test]
    fn complement_non_acgt_returns_n() {
        assert_eq!(complement(b'N'), b'N');
        assert_eq!(complement(b'X'), b'N');
        assert_eq!(complement(b'R'), b'N');
    }

    #[test]
    fn reverse_complement_simple() {
        assert_eq!(reverse_complement(b"ATCG"), b"CGAT");
    }

    #[test]
    fn reverse_complement_single_base() {
        assert_eq!(reverse_complement(b"A"), b"T");
    }

    #[test]
    fn reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    #[test]
    fn reverse_complement_palindrome() {
        assert_eq!(reverse_complement(b"AATT"), b"AATT");
    }

    #[test]
    fn hamming_match() {
        assert_eq!(hamming_n_mismatch(b"ATCG", b"ATCG"), 0);
    }

    #[test]
    fn hamming_one_mismatch() {
        assert_eq!(hamming_n_mismatch(b"ATG", b"ACG"), 1);
    }

    #[test]
    fn hamming_all_mismatch() {
        assert_eq!(hamming_n_mismatch(b"AAA", b"TTT"), 3);
    }

    #[test]
    fn hamming_n_always_mismatch() {
        assert_eq!(hamming_n_mismatch(b"NNN", b"ATG"), 3);
        assert_eq!(hamming_n_mismatch(b"ATG", b"NNN"), 3);
        assert_eq!(hamming_n_mismatch(b"N", b"N"), 1);
    }
}