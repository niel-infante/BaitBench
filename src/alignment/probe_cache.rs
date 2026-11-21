//! Persistent cache of probe placements.
//!
//! Placing a probe panel on a reference is by far the most expensive part of a
//! simulation — on the MTEC panel against a 2.1 GB tick genome it runs for well over an
//! hour — and it depends only on the probes and the reference. Everything downstream
//! (abundance, temperature, salt, capture efficiency, fragment length) changes without
//! changing where the probes land. A titration series, a replicate set or a parameter
//! sweep therefore recomputes an identical alignment once per run.
//!
//! What is cached is the **geometry plus the aligned base pairs**, not the free energy.
//! ΔG is recomputed on load, which costs almost nothing, and that is what keeps a cache
//! valid across `--hybridization-temperature`, `--salt-concentration` and
//! `--simulate-mode`. Caching ΔG instead would have been smaller and simpler, but it
//! would silently narrow the cache to one thermodynamic model and defeat exactly the
//! temperature sweeps `coverage-curve` exists to run.
//!
//! Staleness is treated as a correctness problem, not a convenience one: a cache that
//! silently belongs to a different probe set would produce plausible, wrong output.
//! Every file is stamped with the identity of the inputs it came from and is rejected,
//! loudly and by name, when any of them differ.
//!
//! That identity is the **content** of the parsed sequences, not the file's path, size
//! or modification time. `prepare` rebuilds `combined_reference.fa` inside every run
//! directory, so a path-and-mtime stamp would miss on every run of the pipeline and the
//! cache would never once be used — which is the whole motivating case. Hashing looks
//! expensive against a single run and is not: the sequences are already in memory by
//! this point, so it costs no I/O, and a few seconds of hashing stands against an
//! alignment that can run for over an hour.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;

/// Bump when the on-disk format changes in a way older readers cannot handle.
const CACHE_VERSION: u32 = 1;

/// One probe placement: where it landed and what it paired with.
#[derive(Debug, Clone)]
pub struct Placement {
    pub probe_name: String,
    pub seq_id: String,
    /// 0-based, inclusive start in the reference.
    pub start: usize,
    /// 0-based, exclusive end in the reference.
    pub end: usize,
    /// True when the probe matched as its reverse complement.
    pub reverse: bool,
    /// (probe base, reference base) over the alignment, in reference orientation.
    /// Gap columns are already excluded, so this is exactly what `delta_g` consumes
    /// once the reference side is complemented.
    pub pairs: Vec<(u8, u8)>,
}

/// Identity of the inputs a cache was built from.
#[derive(Debug, Clone, PartialEq)]
pub struct CacheKey {
    pub probes: String,
    pub reference: String,
    pub aligner: String,
}

/// FNV-1a over the parsed sequences, in sorted id order.
///
/// Deliberately not `DefaultHasher`: its output is explicitly not stable across Rust
/// releases, so a toolchain upgrade would silently invalidate every cache on disk. FNV
/// is fixed forever, needs no dependency, and is far faster than the alignment it
/// guards. Ids and lengths are folded in as well as bases, so two panels that differ
/// only in naming or in how sequences are split do not collide.
fn content_hash(seqs: &std::collections::HashMap<String, String>) -> String {
    content_hash_where(seqs, |_| true)
}

/// As `content_hash`, over only the sequences `keep` accepts.
fn content_hash_where(
    seqs: &std::collections::HashMap<String, String>,
    keep: impl Fn(&str) -> bool,
) -> String {
    const OFFSET: u64 = 0xcbf2_9ce4_8422_2325;
    const PRIME: u64 = 0x1000_0000_01b3;
    let mut h = OFFSET;
    let feed = |bytes: &[u8], h: &mut u64| {
        for &b in bytes {
            *h ^= b as u64;
            *h = h.wrapping_mul(PRIME);
        }
    };
    let mut ids: Vec<&String> = seqs.keys().filter(|id| keep(id)).collect();
    ids.sort_unstable();
    feed(&(ids.len() as u64).to_le_bytes(), &mut h);
    let n_ids = ids.len();
    for id in ids {
        let seq = &seqs[id];
        feed(id.as_bytes(), &mut h);
        feed(&[0u8], &mut h);
        feed(&(seq.len() as u64).to_le_bytes(), &mut h);
        feed(seq.as_bytes(), &mut h);
    }
    format!("{:016x}:{}", h, n_ids)
}

impl CacheKey {
    /// Build a key from the already-parsed probes and reference. Both are in memory by
    /// the time a caller needs this, so no file is re-read.
    ///
    /// Only sequences with positive weight enter the reference hash, because only those
    /// are aligned against: the aligners skip zero-weight references entirely. Hashing
    /// the whole file instead would give two runs that differ in `--sample` the same
    /// key while legitimately having different placements, and the second would
    /// silently reuse a set of hits missing every sequence the first had excluded.
    pub fn new(
        probes: &std::collections::HashMap<String, String>,
        reference: &std::collections::HashMap<String, String>,
        weights: &std::collections::HashMap<String, f64>,
        aligner: &str,
    ) -> Self {
        CacheKey {
            probes: content_hash(probes),
            reference: content_hash_where(reference, |id| {
                weights.get(id).copied().unwrap_or(0.0) > 0.0
            }),
            aligner: aligner.to_string(),
        }
    }
}

/// Write placements to `path`, stamped with the inputs they came from.
pub fn write(path: &Path, key: &CacheKey, placements: &[Placement]) -> Result<()> {
    if let Some(dir) = path.parent() {
        if !dir.as_os_str().is_empty() {
            std::fs::create_dir_all(dir).ok();
        }
    }
    let file = File::create(path)
        .with_context(|| format!("cannot create probe-hit cache: {}", path.display()))?;
    let mut w = BufWriter::new(file);
    writeln!(w, "#baitbench-probe-hits\tv{}", CACHE_VERSION)?;
    writeln!(w, "#probes\t{}", key.probes)?;
    writeln!(w, "#reference\t{}", key.reference)?;
    writeln!(w, "#aligner\t{}", key.aligner)?;
    writeln!(w, "#placements\t{}", placements.len())?;
    writeln!(w, "probe_name\tseq_id\tstart\tend\tstrand\tprobe_bases\tref_bases")?;
    for p in placements {
        let probe_bases: String = p.pairs.iter().map(|&(a, _)| a as char).collect();
        let ref_bases: String = p.pairs.iter().map(|&(_, b)| b as char).collect();
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            p.probe_name,
            p.seq_id,
            p.start,
            p.end,
            if p.reverse { '-' } else { '+' },
            probe_bases,
            ref_bases
        )?;
    }
    w.flush()?;
    Ok(())
}

/// Why a cache could not be used. Carried so the caller can say which input changed
/// instead of just recomputing without explanation.
#[derive(Debug)]
pub enum Miss {
    Absent,
    Unreadable(String),
    Stale { field: &'static str, want: String, got: String },
}

impl std::fmt::Display for Miss {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Miss::Absent => write!(f, "no cache file yet"),
            Miss::Unreadable(e) => write!(f, "cache unreadable ({e})"),
            Miss::Stale { field, want, got } => {
                write!(f, "{field} changed since the cache was written (cached {got}, now {want})")
            }
        }
    }
}

/// Read placements from `path`, or explain why they cannot be used.
pub fn read(path: &Path, key: &CacheKey) -> std::result::Result<Vec<Placement>, Miss> {
    if !path.exists() {
        return Err(Miss::Absent);
    }
    let file = File::open(path).map_err(|e| Miss::Unreadable(e.to_string()))?;
    let reader = BufReader::new(file);

    let mut placements = Vec::new();
    let mut saw_header = false;
    for line in reader.lines() {
        let line = line.map_err(|e| Miss::Unreadable(e.to_string()))?;
        if let Some(rest) = line.strip_prefix('#') {
            let mut it = rest.splitn(2, '\t');
            let (k, v) = (it.next().unwrap_or(""), it.next().unwrap_or(""));
            match k {
                "baitbench-probe-hits" => {
                    saw_header = true;
                    if v != format!("v{}", CACHE_VERSION) {
                        return Err(Miss::Stale {
                            field: "cache format version",
                            want: format!("v{}", CACHE_VERSION),
                            got: v.to_string(),
                        });
                    }
                }
                "probes" if v != key.probes => {
                    return Err(Miss::Stale { field: "probes file", want: key.probes.clone(), got: v.to_string() })
                }
                "reference" if v != key.reference => {
                    return Err(Miss::Stale { field: "reference file", want: key.reference.clone(), got: v.to_string() })
                }
                "aligner" if v != key.aligner => {
                    return Err(Miss::Stale { field: "aligner", want: key.aligner.clone(), got: v.to_string() })
                }
                _ => {}
            }
            continue;
        }
        if line.starts_with("probe_name\t") || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            return Err(Miss::Unreadable(format!("short row: {line}")));
        }
        let probe_bases = f[5].as_bytes();
        let ref_bases = f[6].as_bytes();
        if probe_bases.len() != ref_bases.len() {
            return Err(Miss::Unreadable(format!(
                "probe/reference base strings differ in length for {} on {}",
                f[0], f[1]
            )));
        }
        placements.push(Placement {
            probe_name: f[0].to_string(),
            seq_id: f[1].to_string(),
            start: f[2].parse().map_err(|_| Miss::Unreadable(format!("bad start: {}", f[2])))?,
            end: f[3].parse().map_err(|_| Miss::Unreadable(format!("bad end: {}", f[3])))?,
            reverse: f[4] == "-",
            pairs: probe_bases.iter().copied().zip(ref_bases.iter().copied()).collect(),
        });
    }
    if !saw_header {
        return Err(Miss::Unreadable("missing #baitbench-probe-hits header".into()));
    }
    Ok(placements)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;

    fn seqs(pairs: &[(&str, &str)]) -> HashMap<String, String> {
        pairs.iter().map(|(k, v)| (k.to_string(), v.to_string())).collect()
    }

    fn tmp(name: &str) -> std::path::PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("bb_probe_cache_{}_{}", std::process::id(), name));
        p
    }

    fn wts(ids: &[(&str, f64)]) -> HashMap<String, f64> {
        ids.iter().map(|(k, v)| (k.to_string(), *v)).collect()
    }

    fn key(aligner: &str) -> CacheKey {
        CacheKey::new(
            &seqs(&[("p1", "ACGT")]),
            &seqs(&[("ref1", "ACGTACGTACGT")]),
            &wts(&[("ref1", 1.0)]),
            aligner,
        )
    }

    fn sample() -> Vec<Placement> {
        vec![
            Placement {
                probe_name: "p1".into(),
                seq_id: "ref1".into(),
                start: 10,
                end: 14,
                reverse: false,
                pairs: vec![(b'A', b'T'), (b'C', b'G'), (b'G', b'C'), (b'T', b'A')],
            },
            Placement {
                probe_name: "p2".into(),
                seq_id: "ref2".into(),
                start: 0,
                end: 2,
                reverse: true,
                pairs: vec![(b'G', b'G'), (b'T', b'A')],
            },
        ]
    }

    #[test]
    fn round_trips_placements_exactly() {
        let k = key("seed-extend");
        let path = tmp("rt.tsv");
        write(&path, &k, &sample()).unwrap();

        let got = read(&path, &k).expect("cache should be readable");
        assert_eq!(got.len(), 2);
        assert_eq!(got[0].probe_name, "p1");
        assert_eq!(got[0].seq_id, "ref1");
        assert_eq!(got[0].start, 10);
        assert_eq!(got[0].end, 14);
        assert!(!got[0].reverse);
        assert_eq!(got[0].pairs, sample()[0].pairs);
        assert!(got[1].reverse, "strand must survive the round trip");
        assert_eq!(got[1].pairs, sample()[1].pairs);
        let _ = std::fs::remove_file(&path);
    }

    /// The motivating case: `prepare` rewrites combined_reference.fa for every run, so
    /// identical content arriving as a brand-new file must still hit the cache.
    #[test]
    fn identical_content_hits_regardless_of_file_identity() {
        let a = CacheKey::new(
            &seqs(&[("p1", "ACGT")]),
            &seqs(&[("ref1", "ACGTACGTACGT")]),
            &wts(&[("ref1", 1.0)]),
            "seed-extend",
        );
        let b = CacheKey::new(
            &seqs(&[("p1", "ACGT")]),
            &seqs(&[("ref1", "ACGTACGTACGT")]),
            &wts(&[("ref1", 1.0)]),
            "seed-extend",
        );
        assert_eq!(a, b);
        let path = tmp("same_content.tsv");
        write(&path, &a, &sample()).unwrap();
        assert_eq!(read(&path, &b).expect("should hit").len(), 2);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn absent_cache_is_a_miss_not_an_error() {
        match read(Path::new("/nonexistent/probe_cache.tsv"), &key("seed-extend")) {
            Err(Miss::Absent) => {}
            other => panic!("expected Absent, got {other:?}"),
        }
    }

    /// A cache written against different probes must be rejected by name, not used.
    #[test]
    fn changed_probes_invalidate_the_cache() {
        let a = key("seed-extend");
        let b = CacheKey::new(
            &seqs(&[("p1", "ACGTACGTACGT")]), // same name, different sequence
            &seqs(&[("ref1", "ACGTACGTACGT")]),
            &wts(&[("ref1", 1.0)]),
            "seed-extend",
        );
        let path = tmp("stale_probes.tsv");
        write(&path, &a, &sample()).unwrap();
        match read(&path, &b) {
            Err(Miss::Stale { field, .. }) => assert_eq!(field, "probes file"),
            other => panic!("expected a stale-probes miss, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    /// Editing the reference must invalidate even when its length is unchanged.
    #[test]
    fn changed_reference_invalidates_the_cache() {
        let a = key("seed-extend");
        let b = CacheKey::new(
            &seqs(&[("p1", "ACGT")]),
            &seqs(&[("ref1", "ACGTACGTACGA")]), // one base differs, same length
            &wts(&[("ref1", 1.0)]),
            "seed-extend",
        );
        let path = tmp("stale_ref.tsv");
        write(&path, &a, &sample()).unwrap();
        match read(&path, &b) {
            Err(Miss::Stale { field, .. }) => assert_eq!(field, "reference file"),
            other => panic!("expected a stale-reference miss, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    /// Renaming a sequence changes which probe lands where in the output, so the id
    /// must be part of the identity even when the bases are the same.
    #[test]
    fn renaming_a_sequence_invalidates_the_cache() {
        let a = key("seed-extend");
        let b = CacheKey::new(
            &seqs(&[("p1", "ACGT")]),
            &seqs(&[("renamed", "ACGTACGTACGT")]),
            &wts(&[("renamed", 1.0)]),
            "seed-extend",
        );
        assert_ne!(a.reference, b.reference);
    }

    /// Only weighted sequences are aligned, so zeroing one changes the placements and
    /// must change the key. Hashing the whole reference instead let two runs that
    /// differ only in `--sample` share a cache and silently reuse the wrong hits.
    #[test]
    fn zero_weight_sequences_are_outside_the_key() {
        let probes = seqs(&[("p1", "ACGT")]);
        let reference = seqs(&[("a", "ACGTACGT"), ("b", "TTTTGGGG")]);
        let both = CacheKey::new(&probes, &reference, &wts(&[("a", 1.0), ("b", 1.0)]), "seed-extend");
        let only_a = CacheKey::new(&probes, &reference, &wts(&[("a", 1.0), ("b", 0.0)]), "seed-extend");
        assert_ne!(both.reference, only_a.reference, "dropping b must change the key");

        // And a key built over only the weighted subset matches one built from a
        // reference that never contained the extra sequence at all.
        let just_a = CacheKey::new(&probes, &seqs(&[("a", "ACGTACGT")]), &wts(&[("a", 1.0)]), "seed-extend");
        assert_eq!(only_a.reference, just_a.reference);
    }

    /// Switching aligner changes which placements are found, so it must invalidate.
    #[test]
    fn changed_aligner_invalidates_the_cache() {
        let path = tmp("aligner.tsv");
        write(&path, &key("seed-extend"), &sample()).unwrap();
        match read(&path, &key("blast")) {
            Err(Miss::Stale { field, .. }) => assert_eq!(field, "aligner"),
            other => panic!("expected a stale-aligner miss, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    /// A truncated or foreign file must not be silently treated as an empty panel.
    #[test]
    fn file_without_the_magic_header_is_rejected() {
        let path = tmp("foreign.tsv");
        std::fs::write(&path, "probe_name\tseq_id\tstart\tend\tstrand\tprobe_bases\tref_bases\n").unwrap();
        match read(&path, &key("seed-extend")) {
            Err(Miss::Unreadable(_)) => {}
            other => panic!("expected Unreadable, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }
}
