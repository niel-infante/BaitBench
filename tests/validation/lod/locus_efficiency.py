#!/usr/bin/env python3
"""Per-locus capture efficiency from a BaitBench LOD run, relative to rplB.

  usage: locus_efficiency.py <condition_outdir> [<condition_outdir> ...]

The real reference values in ../expected/lod_locus_efficiency.tsv were measured with
mosdepth over gene intervals on the B. burgdorferi chromosome. This script measures the
simulated analogue differently, because in genome mode reads are mapped to the per-gene
*target* sequences rather than to the chromosome: for each gene it sums mapped reads
across every Borrelia species variant of that gene and divides by the summed target
length, giving mean depth per gene, then normalises to rplB.

The two bases agree on what they measure — mean depth over a locus, relative to rplB —
but not on the exact denominator, so compare the *ranking* and the order of magnitude of
the spread, not the third decimal place. The gate the LOD data actually constrains is
that the weak loci (nifS, clpA) stay roughly an order of magnitude below rplB and that
the ratios hold across the dilution series.
"""
import os, sys, glob, csv, re, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
EXPECTED = os.path.join(HERE, "..", "expected", "lod_locus_efficiency.tsv")


def load_expected():
    out = {}
    with open(EXPECTED) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4 or f[0] == "gene":
                continue
            # "16S/16s" is one locus spelled two ways in the panel.
            for name in f[0].split("/"):
                out[name.lower()] = float(f[3])
    return out


def load_sim(outdir):
    """Mean depth per gene, summed over all Borrelia variants of that gene."""
    detail = None
    for pat in ("detected_detail.tsv", "group_detail.tsv"):
        hits = glob.glob(os.path.join(outdir, "**", pat), recursive=True)
        if hits:
            detail = hits[0]
            break
    if not detail:
        return None

    reads, length, n_variants = {}, {}, {}
    with open(detail) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            ref = row.get("reference_id") or row.get("group") or ""
            if "|" not in ref or not ref.lower().startswith("borrelia"):
                continue
            gene = ref.split("|", 1)[1].lower()
            try:
                n = float(row.get("reads_assigned") or 0)
                L = float(row.get("ref_length") or 0)
            except ValueError:
                continue
            reads[gene] = reads.get(gene, 0.0) + n
            length[gene] = length.get(gene, 0.0) + L
            n_variants[gene] = n_variants.get(gene, 0) + 1
    # Divide by the MEAN variant length, not the summed length.
    #
    # A gene is represented in the panel by several near-identical species variants
    # (rplB 13, clpX 3), and reads from one genomic locus scatter across whichever
    # variant they match best. Dividing the pooled reads by the summed length of every
    # variant therefore penalises a gene in proportion to how many variants the panel
    # happens to carry: it made clpX look 4.2x deeper than rplB purely because
    # 13/3 = 4.3. The locus exists once in the genome, so the denominator is the length
    # of one copy.
    return {
        g: reads[g] / (length[g] / max(n_variants[g], 1))
        for g in reads
        if length.get(g, 0) > 0
    }


def main(dirs):
    exp = load_expected()
    rows = []
    for d in dirs:
        cond = os.path.basename(d.rstrip("/"))
        depth = load_sim(d)
        if not depth:
            print(f"{cond:10s} NO OUTPUT")
            continue
        base = depth.get("rplb")
        if not base:
            print(f"{cond:10s} no rplB signal — cannot normalise")
            continue
        rows.append((cond, {g: v / base for g, v in depth.items()}))

    if not rows:
        return 1

    genes = [g for g in ("nifs", "clpa", "clpx", "uvra", "pyrg", "recg", "16s", "rplb")
             if any(g in rel for _, rel in rows)]

    hdr = f"{'gene':8s} {'real':>7s}" + "".join(f"{c:>10s}" for c, _ in rows)
    print(hdr)
    print("-" * len(hdr))
    for g in genes:
        line = f"{g:8s} {exp.get(g, float('nan')):7.4f}"
        for _, rel in rows:
            v = rel.get(g)
            line += f"{v:10.4f}" if v is not None else f"{'-':>10s}"
        print(line)

    # Flatness across conditions is the second half of the gate.
    print()
    print(f"{'gene':8s} {'spread (max/min across conditions)':>36s}")
    for g in genes:
        vals = [rel[g] for _, rel in rows if rel.get(g, 0) > 0]
        if len(vals) > 1:
            print(f"{g:8s} {max(vals) / min(vals):36.2f}x")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]) if len(sys.argv) > 1 else 1)
