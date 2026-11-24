#!/usr/bin/env python3
"""Derive BaitBench input weights from the non-enriched arm of the mock community.

Why this exists: the rejected manuscript's validation set `--sample` weights from the
very reads it then predicted (reviewer item 3, circularity). PRJNA1081843 includes six
non-enrichment runs of the same synthetic mock community, so the input composition can
be *measured* on a different library from the one being predicted. See
`improvements/parameterisation.md` A3.

Two outputs:

  sample_nonenriched.tsv  per-sequence weights for `--sample`
  input_share.json        per-organism input share, the denominator for enrichment

The weight conversion is the part that is easy to get wrong. Read share in a shotgun
library is proportional to DNA mass, so to abundance x genome length; BaitBench samples
a reference proportional to `weight x retention_mass`, and retention mass is
proportional to length. So `weight = read_share / group_length`, and every sequence in a
group carries the same weight. Weights come out around 1e-8, which is why the weights
file is written in exponential notation — a fixed six decimals rounds them all to zero.

Usage:
  input_weights.py --real DIR --groups FILE --genomes FILE --outdir DIR [--pattern NEG]
"""
import argparse
import json
import os
import sys
from collections import defaultdict


def read_counts(path):
    counts = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                counts[parts[0]] = float(parts[1])
    return counts


def read_groups(path):
    groups = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 2:
                groups[parts[0]] = parts[1]
    return groups


def seq_lengths(path):
    lengths = {}
    name, n = None, 0
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--real', required=True, help='directory of *_counts.tsv files')
    ap.add_argument('--groups', required=True)
    ap.add_argument('--genomes', required=True)
    ap.add_argument('--outdir', required=True)
    ap.add_argument('--pattern', default='NEG',
                    help='substring selecting the non-enriched runs (default NEG)')
    args = ap.parse_args()

    files = sorted(f for f in os.listdir(args.real)
                   if f.endswith('_counts.tsv') and args.pattern in f)
    if not files:
        sys.exit(f"no count files matching {args.pattern!r} in {args.real}")

    groups = read_groups(args.groups)
    lengths = seq_lengths(args.genomes)
    org = lambda s: groups.get(s, s)

    # Pool the replicates. Pooling rather than averaging shares weights each run by its
    # depth, which is what we want: they are replicate measurements of one composition,
    # and A3 shows they agree at r(log10) >= 0.989 both within and across chemistry.
    pooled = defaultdict(float)
    unmapped = 0.0
    for f in files:
        for seq, c in read_counts(os.path.join(args.real, f)).items():
            # the count tables carry an `__unmapped__` row (~2% of reads). It is not an
            # organism, and leaving it in the denominator would deflate every share.
            if seq.startswith('__'):
                unmapped += c
                continue
            pooled[org(seq)] += c
    total = sum(pooled.values())
    if total <= 0:
        sys.exit("non-enriched runs contain no reads")
    share = {o: c / total for o, c in pooled.items() if c > 0}

    # group length, over the sequences actually present in the reference
    glen = defaultdict(float)
    for seq, l in lengths.items():
        glen[org(seq)] += l

    missing = [o for o in share if glen.get(o, 0) <= 0]
    if missing:
        sys.exit(f"organisms counted but absent from the reference: {sorted(missing)}")

    # weight = share / length, then normalised so the largest is 1.0 — the scale is
    # arbitrary (only ratios drive composition) but keeping it O(1) keeps the weights
    # file readable and away from underflow.
    raw = {o: share[o] / glen[o] for o in share}
    scale = max(raw.values())
    per_org = {o: v / scale for o, v in raw.items()}

    os.makedirs(args.outdir, exist_ok=True)
    sample_path = os.path.join(args.outdir, 'sample_nonenriched.tsv')
    n = 0
    with open(sample_path, 'w') as fh:
        for seq in sorted(lengths):
            w = per_org.get(org(seq))
            if w is None:
                continue
            fh.write(f"{seq}\t{w:.10g}\n")
            n += 1
    share_path = os.path.join(args.outdir, 'input_share.json')
    with open(share_path, 'w') as fh:
        json.dump(share, fh, sort_keys=True)

    print(f"pooled {len(files)} non-enriched runs ({total:.0f} mapped reads, "
          f"{len(share)} organisms; {unmapped:.0f} unmapped excluded, "
          f"{100 * unmapped / (total + unmapped):.1f}%)")
    print(f"  {sample_path}: {n} sequences, weights "
          f"{min(per_org.values()):.3g} to {max(per_org.values()):.3g}")
    print(f"  {share_path}: input shares "
          f"{min(share.values()):.3g} to {max(share.values()):.3g}")


if __name__ == '__main__':
    main()
