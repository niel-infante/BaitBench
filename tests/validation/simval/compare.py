#!/usr/bin/env python3
"""Compare BaitBench-simulated per-organism read composition with the real
D6331 mock-community PacBio CCS data.

  usage: compare.py <TAG>            e.g. compare.py MOXTAA_res_cf0.3

Real counts come from real/<SAMPLE>_counts.tsv (minimap2 primary alignments to the
21 reference genomes). Simulated counts come from the matching BaitBench run.
Both are collapsed contig -> species with inputs/groups.tsv, so the two sides are
treated identically.

Reported: Spearman rho (rank agreement), Pearson r on log10 proportions, and Lin's
concordance correlation coefficient (agreement, not just correlation). Also a variant
with the five near-identical E. coli strains collapsed, since read assignment among
them is close to arbitrary for both real and simulated data.
"""
import os, sys, math, statistics, glob, csv

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# External data is referenced in place, never copied (see ../config.sh).
SIMVAL_DIR = os.environ.get("SIMVAL_DIR", "")
SIMVAL_REAL = os.environ.get("SIMVAL_REAL", os.path.join(SIMVAL_DIR, "real"))
SIMVAL_GROUPS = os.environ.get("SIMVAL_GROUPS", os.path.join(SIMVAL_DIR, "inputs", "groups.tsv"))
RESULTS_DIR = os.environ.get("VALIDATION_OUT", os.path.join(HERE, "results"))

def load_groups():
    g = {}
    with open(SIMVAL_GROUPS) as fh:
        for line in fh:
            c, s = line.rstrip("\n").split("\t")
            g[c] = s
    return g

def load_real(sample, groups):
    counts = {}
    with open(f"{SIMVAL_REAL}/{sample}_counts.tsv") as fh:
        for line in fh:
            ref, n = line.rstrip("\n").split("\t")
            if ref == "__unmapped__":
                continue
            counts[groups.get(ref, ref)] = counts.get(groups.get(ref, ref), 0) + int(n)
    return counts

def load_sim(tag, groups):
    """Read per-organism read counts from the run's group_detail.tsv.

    BaitBench writes that table when --groups is supplied: one row per group with a
    total_reads column, which is already collapsed to species here.
    """
    hits = sorted(glob.glob(f"{RESULTS_DIR}/{tag}/**/group_detail.tsv", recursive=True))
    if not hits:
        return None, None
    path = hits[-1]
    counts = {}
    with open(path) as fh:
        rdr = csv.DictReader(fh, delimiter="\t")
        for row in rdr:
            name = row["group_name"]
            counts[name] = counts.get(name, 0) + int(row["total_reads"])
    return counts, path

def props(d, keys):
    tot = sum(d.get(k, 0) for k in keys)
    return [d.get(k, 0) / tot if tot else 0.0 for k in keys]

def spearman(a, b):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0]*len(v); i = 0
        while i < len(v):
            j = i
            while j+1 < len(v) and v[order[j+1]] == v[order[i]]: j += 1
            avg = (i+j)/2 + 1
            for k in range(i, j+1): r[order[k]] = avg
            i = j+1
        return r
    ra, rb = rank(a), rank(b)
    return pearson(ra, rb)

def pearson(a, b):
    n = len(a); ma, mb = sum(a)/n, sum(b)/n
    num = sum((x-ma)*(y-mb) for x, y in zip(a, b))
    da = math.sqrt(sum((x-ma)**2 for x in a)); db = math.sqrt(sum((y-mb)**2 for y in b))
    return num/(da*db) if da and db else float("nan")

def ccc(a, b):
    n = len(a); ma, mb = sum(a)/n, sum(b)/n
    va = sum((x-ma)**2 for x in a)/n; vb = sum((y-mb)**2 for y in b)/n
    cov = sum((x-ma)*(y-mb) for x, y in zip(a, b))/n
    d = va + vb + (ma-mb)**2
    return 2*cov/d if d else float("nan")

def collapse_ecoli(d):
    out = {}
    for k, v in d.items():
        kk = "Escherichia_coli" if k.startswith("Escherichia_coli") else k
        out[kk] = out.get(kk, 0) + v
    return out

def report(real, sim, label, rows=None):
    keys = sorted(set(real) | set(sim))
    rp, sp = props(real, keys), props(sim, keys)
    both = [(k, r, s) for k, r, s in zip(keys, rp, sp) if r > 0 or s > 0]
    keys2 = [k for k, _, _ in both]; rp2 = [r for _, r, _ in both]; sp2 = [s for _, _, s in both]

    # Floor a zero proportion at that arm's DETECTION LIMIT — the proportion one read
    # would represent — not at an arbitrary epsilon.
    #
    # This was `eps = 1e-9`. Real proportions here span ~1e-5 to 0.8, so a floored zero
    # landed at log10 = -9, four orders of magnitude below the smallest real
    # observation, and behaved as an enormous leverage point in Pearson r and Lin's CCC.
    # Measured on MOXTAC, which has two organisms with zero real reads: r and CCC read
    # 0.75 against 0.96-0.97 for the replicates that happen to have none, making a
    # counting accident look like a worse prediction. A zero means "fewer than one
    # read", so the detection limit is what it is actually bounded by.
    #
    # Spearman is computed on the raw proportions and is unaffected either way; it is
    # the log-scale statistics that the floor distorts.
    real_lim = 1.0 / max(sum(real.values()), 1)
    sim_lim = 1.0 / max(sum(sim.values()), 1)
    lr = [math.log10(max(x, real_lim)) for x in rp2]
    ls = [math.log10(max(x, sim_lim)) for x in sp2]
    censored = sum(1 for r, s in zip(rp2, sp2) if r <= 0 or s <= 0)
    stats = (spearman(rp2, sp2), pearson(lr, ls), ccc(lr, ls))
    print(f"\n=== {label}  (n={len(keys2)}) ===")
    print(f"  Spearman rho          : {stats[0]:.4f}")
    print(f"  Pearson r (log10)     : {stats[1]:.4f}")
    print(f"  Lin's CCC (log10)     : {stats[2]:.4f}")
    if censored:
        print(f"  censored at detection limit: {censored} "
              f"(real {real_lim:.2e}, sim {sim_lim:.2e})")
    print(f"\n  {'organism':<34}{'real':>10}{'sim':>10}{'log2 s/r':>10}")
    for k, r, s in sorted(both, key=lambda t: -t[1]):
        lf = math.log2(max(s, sim_lim)/max(r, real_lim)) if r > 0 and s > 0 else float('nan')
        print(f"  {k:<34}{r:>10.5f}{s:>10.5f}{lf:>10.2f}")
        if rows is not None:
            rows.append({
                "set": label, "organism": k,
                "real_prop": f"{r:.10f}", "sim_prop": f"{s:.10f}",
                "real_reads": real.get(k, 0), "sim_reads": sim.get(k, 0),
                "spearman_rho": f"{stats[0]:.4f}",
                "pearson_r_log10": f"{stats[1]:.4f}",
                "ccc_log10": f"{stats[2]:.4f}",
            })

def main():
    tag = sys.argv[1]
    # Optional second arg overrides which real run to compare against (used for
    # mechanics checks where the results dir name does not encode a sample).
    sample = sys.argv[2] if len(sys.argv) > 2 else tag.split("_")[0]
    groups = load_groups()
    real = load_real(sample, groups)
    sim, path = load_sim(tag, groups)
    if sim is None:
        print(f"could not find a counts table under {RESULTS_DIR}/{tag}", file=sys.stderr)
        return 1
    print(f"real: {SIMVAL_REAL}/{sample}_counts.tsv  ({sum(real.values()):,} reads)")
    print(f"sim : {path}  ({sum(sim.values()):,} reads)")
    rows = []
    report(real, sim, "per strain", rows)
    report(collapse_ecoli(real), collapse_ecoli(sim), "E. coli strains collapsed", rows)

    # Emit the comparison as a table so it can be plotted and reused, rather than only
    # printed. SIMVAL_TSV names the destination; unset means print-only as before.
    out_tsv = os.environ.get("SIMVAL_TSV")
    if out_tsv and rows:
        cols = ["set", "organism", "real_prop", "sim_prop", "real_reads", "sim_reads",
                "spearman_rho", "pearson_r_log10", "ccc_log10"]
        with open(out_tsv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print(f"\nwrote {out_tsv}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
