#!/usr/bin/env python3
"""Score predicted against measured per-organism enrichment factor.

This is the non-circular validation target (`improvements/parameterisation.md` B0):
the input composition comes from the non-enriched arm, the prediction target is the
enriched arm, and they are different libraries.

    measured_o  = enriched_share_o / nonenriched_share_o
    predicted_o = (weight_o x retention_mass_o) / sum_j(...) / nonenriched_share_o

The predicted side needs no reads. The sampler draws a reference proportional to
`weight x retention_mass`, so an organism's expected output share is its pooled
`weight x mass` — which `fragments.retention.tsv` already contains. A `simulate` run
with a warm `--probe-hits` cache therefore scores in ~90 seconds instead of ~25
minutes, and carries no sampling noise. Pass --reads to score the sampled reads instead
and compare the two.

Note that predicted enrichment is independent of the input composition except through
one scalar common to all organisms (parameterisation.md, under B0b), so the ranking and
all centred statistics test the capture model alone.

Scoring is on log10 enrichment, and reports **centred** RMSE alongside raw: enrichment
factors are compositional, so one badly predicted dominant organism shifts every other
value, and centring separates shape error from that closure effect.

Usage:
  enrichment.py --run DIR --share input_share.json --enriched MOXTAA_counts.tsv \
                --groups groups.tsv [--reads] [--min-fragments N] [--tsv OUT]
"""
import argparse
import csv
import json
import math
import os
import sys
from collections import defaultdict

FLOOR = 1e-4  # keeps log10 finite for an organism the simulation never samples


def read_groups(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if line and not line.startswith('#'):
            p = line.split('\t')
            if len(p) >= 2:
                out[p[0]] = p[1]
    return out


def read_counts(path):
    out = {}
    for line in open(path):
        line = line.strip()
        if line and not line.startswith('#'):
            p = line.split('\t')
            if len(p) >= 2 and not p[0].startswith('__'):
                out[p[0]] = float(p[1])
    return out


def pearson(a, b):
    n = len(a)
    ma, mb = sum(a) / n, sum(b) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(a, b))
    den = math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return num / den if den else float('nan')


def ranks(v):
    n = len(v)
    order = sorted(range(n), key=lambda i: v[i])
    rk = [0.0] * n
    i = 0
    while i < n:
        j = i
        while j + 1 < n and v[order[j + 1]] == v[order[i]]:
            j += 1
        for k in range(i, j + 1):
            rk[order[k]] = (i + j) / 2 + 1
        i = j + 1
    return rk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--run', required=True, help='BaitBench run directory')
    ap.add_argument('--share', required=True, help='input_share.json from input_weights.py')
    ap.add_argument('--enriched', required=True, help='enriched arm *_counts.tsv')
    ap.add_argument('--groups', required=True)
    ap.add_argument('--reads', action='store_true',
                    help='score the sampled reads (group_detail.tsv) instead of the '
                         'analytic retention mass')
    ap.add_argument('--num-fragments', type=float, default=None,
                    help='if given, exclude organisms whose expected fragment count '
                         'falls below --min-fragments')
    ap.add_argument('--min-fragments', type=float, default=1.0)
    ap.add_argument('--min-real-reads', type=float, default=10.0,
                    help='exclude organisms with fewer than this many reads in the '
                         'enriched arm; their measured enrichment is Poisson noise')
    ap.add_argument('--pool-prefix', action='append', default=['Escherichia_coli'],
                    help='collapse organisms sharing this name prefix (repeatable); '
                         'strains of one species are counted together')
    ap.add_argument('--exclude', action='append', default=[],
                    help='organism to drop, with the reason stated in the report')
    ap.add_argument('--tsv', help='write the per-organism table here')
    args = ap.parse_args()

    groups = read_groups(args.groups)

    def pool(name):
        for pre in args.pool_prefix:
            if name.startswith(pre):
                return pre
        return name

    org = lambda s: pool(groups.get(s, s))

    share = defaultdict(float)
    for k, v in json.load(open(args.share)).items():
        share[pool(k)] += v

    enr = defaultdict(float)
    for seq, c in read_counts(args.enriched).items():
        enr[org(seq)] += c
    enr_total = sum(enr.values())
    if enr_total <= 0:
        sys.exit(f"no reads in {args.enriched}")

    pred_mass = defaultdict(float)
    if args.reads:
        path = os.path.join(args.run, 'group_detail.tsv')
        rows = list(csv.DictReader(open(path), delimiter='\t'))
        for r in rows:
            pred_mass[pool(r['group_name'])] += float(r['total_reads'])
    else:
        path = os.path.join(args.run, 'fragments.retention.tsv')
        weights = {}
        for line in open(os.path.join(args.run, 'weights.txt')):
            line = line.strip()
            if line and not line.startswith('#'):
                p = line.split('\t')
                weights[p[0]] = float(p[1])
        for r in csv.DictReader(open(path), delimiter='\t'):
            pred_mass[org(r['seq_id'])] += (weights.get(r['seq_id'], 0.0)
                                            * float(r['retention_mass']))
    pred_total = sum(pred_mass.values())
    if pred_total <= 0:
        sys.exit(f"simulation produced nothing in {path}")

    excluded = {}
    for o in args.exclude:
        excluded[o] = 'excluded on the command line'
    if args.num_fragments:
        for o, s in share.items():
            if s * args.num_fragments < args.min_fragments and o not in excluded:
                excluded[o] = (f"{s * args.num_fragments:.2f} expected fragments at "
                               f"{args.num_fragments:.0f} — unsampleable at this depth")
    # A floor on the measured side as well as the simulated one. At 1-4 reads the
    # measured enrichment is Poisson noise with ~50-100% relative error, and because
    # such an organism can be present in one replicate and absent in the next, it also
    # changes *which* organisms are scored between replicates — which shows up as
    # replicate disagreement that is really counting noise. Measured on the ARG panel:
    # E. faecalis has 4 / 1 / 0 reads across MOXTAA/B/C and C. perfringens 0 / 0 / 1.
    for o in share:
        n_reads = enr.get(o, 0.0)
        if 0 < n_reads < args.min_real_reads and o not in excluded:
            excluded[o] = (f"{n_reads:.0f} reads in the enriched arm — below "
                           f"--min-real-reads {args.min_real_reads:.0f}, so the "
                           f"measured enrichment is counting noise")

    scored, skipped = [], []
    for o in share:
        if share[o] <= 0 or enr.get(o, 0) <= 0:
            continue
        m = enr[o] / enr_total / share[o]
        p = pred_mass.get(o, 0.0) / pred_total / share[o]
        (skipped if o in excluded else scored).append((o, share[o], p, m))
    scored.sort(key=lambda r: -r[3])

    w = max([len(r[0]) for r in scored + skipped] + [8])
    print(f"{'organism':<{w}}  {'input':>9}  {'pred':>8}  {'meas':>8}  {'ratio':>7}")
    for o, s, p, m in scored:
        print(f"{o:<{w}}  {s:9.5f}  {p:8.3f}  {m:8.3f}  {p / m:7.2f}")
    for o, s, p, m in sorted(skipped, key=lambda r: -r[3]):
        print(f"{o:<{w}}  {s:9.5f}  {p:8.3f}  {m:8.3f}  {'skip':>7}")

    if args.tsv:
        with open(args.tsv, 'w') as fh:
            fh.write("organism\tinput_share\tpredicted\tmeasured\tscored\n")
            for o, s, p, m in scored:
                fh.write(f"{o}\t{s:.8g}\t{p:.6g}\t{m:.6g}\ttrue\n")
            for o, s, p, m in skipped:
                fh.write(f"{o}\t{s:.8g}\t{p:.6g}\t{m:.6g}\tfalse\n")

    if skipped:
        print("\nexcluded:")
        for o, _, _, _ in skipped:
            print(f"  {o}: {excluded[o]}")

    n = len(scored)
    if n < 3:
        sys.exit("\nfewer than 3 scorable organisms")
    lp = [math.log10(max(r[2], FLOOR)) for r in scored]
    lm = [math.log10(r[3]) for r in scored]
    mp, mm = sum(lp) / n, sum(lm) / n
    rmse = math.sqrt(sum((a - b) ** 2 for a, b in zip(lp, lm)) / n)
    crmse = math.sqrt(sum(((a - mp) - (b - mm)) ** 2 for a, b in zip(lp, lm)) / n)
    lo, hi = min(r[3] for r in scored), max(r[3] for r in scored)
    print()
    print(f"n                          : {n}  ({'reads' if args.reads else 'analytic'})")
    print(f"Spearman rho               : {pearson(ranks(lp), ranks(lm)):+.4f}")
    print(f"Pearson r (log10)          : {pearson(lp, lm):+.4f}")
    print(f"log10 RMSE                 : {rmse:.3f}  ({10 ** rmse:.2f}x)")
    print(f"log10 RMSE, centred        : {crmse:.3f}  ({10 ** crmse:.2f}x)")
    print(f"log10 bias (pred - meas)   : {mp - mm:+.3f}")
    print(f"measured dynamic range     : {lo:.3f} - {hi:.2f}  ({hi / lo:.0f}x)")


if __name__ == '__main__':
    main()
