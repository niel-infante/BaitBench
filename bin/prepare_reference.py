#!/usr/bin/env python3
"""
Prepare Reference for BaitBench

Combines target and distractor FASTA files into a single reference file
and generates a weights file for read simulation.
"""

import argparse
import sys


def parse_fasta_ids(fasta_file):
    """Parse FASTA file and return list of sequence IDs."""
    ids = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line[1:].strip().split()[0]
                ids.append(seq_id)
    return ids


def combine_fastas(target_fasta, distractor_fasta, output_fasta):
    """Combine two FASTA files into one."""
    with open(output_fasta, 'w') as out:
        # Write targets
        with open(target_fasta, 'r') as f:
            out.write(f.read())

        # Ensure newline separation
        out.write('\n')

        # Write distractors
        with open(distractor_fasta, 'r') as f:
            out.write(f.read())


def generate_weights(target_ids, distractor_ids, distractor_fraction, output_weights):
    """
    Generate weights file.

    Targets get weight 1.0 each.
    Distractors are weighted to produce distractor_fraction of total reads.

    Example: If distractor_fraction=0.9 and there are 2 targets and 10 distractors,
    targets collectively get 10% of reads, distractors get 90%.
    Each target: 0.1/2 = 0.05
    Each distractor: 0.9/10 = 0.09
    Normalized so targets = 1.0:
    Each target: 1.0
    Each distractor: 0.09/0.05 = 1.8
    """
    n_targets = len(target_ids)
    n_distractors = len(distractor_ids)

    if n_targets == 0:
        print("Error: No target sequences found!", file=sys.stderr)
        sys.exit(1)

    # Calculate weights
    # We want: sum(target_weights) / sum(all_weights) = (1 - distractor_fraction)
    # And: sum(distractor_weights) / sum(all_weights) = distractor_fraction

    # Set target weight to 1.0 each
    target_weight = 1.0

    # Calculate distractor weight to achieve desired fraction
    if n_distractors > 0 and distractor_fraction > 0:
        # total_target_weight = n_targets * target_weight
        # total_distractor_weight = n_distractors * distractor_weight
        # distractor_fraction = total_distractor_weight / (total_target_weight + total_distractor_weight)
        # Solving for distractor_weight:
        # distractor_fraction * (n_targets * target_weight + n_distractors * distractor_weight) = n_distractors * distractor_weight
        # distractor_fraction * n_targets * target_weight = n_distractors * distractor_weight * (1 - distractor_fraction)
        # distractor_weight = (distractor_fraction * n_targets * target_weight) / (n_distractors * (1 - distractor_fraction))

        if distractor_fraction >= 1.0:
            print("Error: distractor_fraction must be less than 1.0", file=sys.stderr)
            sys.exit(1)

        distractor_weight = (distractor_fraction * n_targets * target_weight) / (n_distractors * (1 - distractor_fraction))
    else:
        distractor_weight = 0.0

    # Write weights file
    with open(output_weights, 'w') as f:
        f.write("# BaitBench weights file\n")
        f.write(f"# Targets: {n_targets}, Distractors: {n_distractors}\n")
        f.write(f"# Distractor fraction: {distractor_fraction}\n")
        f.write(f"# Target weight: {target_weight}, Distractor weight: {distractor_weight:.6f}\n")
        f.write("#\n")

        for seq_id in target_ids:
            f.write(f"{seq_id}\t{target_weight}\n")

        for seq_id in distractor_ids:
            f.write(f"{seq_id}\t{distractor_weight:.6f}\n")

    return target_weight, distractor_weight


def write_id_lists(target_ids, distractor_ids, targets_list, distractors_list):
    """Write separate ID list files for metrics calculation."""
    with open(targets_list, 'w') as f:
        for seq_id in target_ids:
            f.write(f"{seq_id}\n")

    with open(distractors_list, 'w') as f:
        for seq_id in distractor_ids:
            f.write(f"{seq_id}\n")


def main():
    parser = argparse.ArgumentParser(
        description='Prepare combined reference and weights for BaitBench',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s --targets targets.fa --distractors background.fa --output-fasta combined.fa --output-weights weights.txt
  %(prog)s -t viral.fa -d host.fa -f 0.95 -o combined.fa -w weights.txt

The distractor_fraction parameter controls what fraction of simulated reads
come from distractor sequences (default 0.9 = 90%% distractors, 10%% targets).
        '''
    )

    parser.add_argument('-t', '--targets', required=True,
                        help='Input FASTA file with target sequences')
    parser.add_argument('-d', '--distractors', required=True,
                        help='Input FASTA file with distractor sequences')
    parser.add_argument('-o', '--output-fasta', required=True,
                        help='Output combined FASTA file')
    parser.add_argument('-w', '--output-weights', required=True,
                        help='Output weights file')
    parser.add_argument('-f', '--distractor-fraction', type=float, default=0.9,
                        help='Fraction of reads from distractors (default: 0.9)')
    parser.add_argument('--targets-list',
                        help='Output file listing target IDs (for metrics)')
    parser.add_argument('--distractors-list',
                        help='Output file listing distractor IDs (for metrics)')

    args = parser.parse_args()

    # Parse input FASTAs
    print(f"Reading targets from {args.targets}...", file=sys.stderr)
    target_ids = parse_fasta_ids(args.targets)
    print(f"  Found {len(target_ids)} target sequences", file=sys.stderr)

    print(f"Reading distractors from {args.distractors}...", file=sys.stderr)
    distractor_ids = parse_fasta_ids(args.distractors)
    print(f"  Found {len(distractor_ids)} distractor sequences", file=sys.stderr)

    # Check for overlapping IDs
    overlap = set(target_ids) & set(distractor_ids)
    if overlap:
        print(f"Warning: {len(overlap)} sequence IDs appear in both targets and distractors:", file=sys.stderr)
        for seq_id in list(overlap)[:5]:
            print(f"  {seq_id}", file=sys.stderr)
        if len(overlap) > 5:
            print(f"  ... and {len(overlap) - 5} more", file=sys.stderr)

    # Combine FASTAs
    print(f"Combining FASTAs to {args.output_fasta}...", file=sys.stderr)
    combine_fastas(args.targets, args.distractors, args.output_fasta)

    # Generate weights
    print(f"Generating weights to {args.output_weights}...", file=sys.stderr)
    target_w, distractor_w = generate_weights(
        target_ids, distractor_ids,
        args.distractor_fraction,
        args.output_weights
    )
    print(f"  Target weight: {target_w}", file=sys.stderr)
    print(f"  Distractor weight: {distractor_w:.6f}", file=sys.stderr)

    # Write ID lists if requested
    if args.targets_list or args.distractors_list:
        if args.targets_list:
            print(f"Writing targets list to {args.targets_list}...", file=sys.stderr)
        if args.distractors_list:
            print(f"Writing distractors list to {args.distractors_list}...", file=sys.stderr)
        write_id_lists(
            target_ids, distractor_ids,
            args.targets_list or '/dev/null',
            args.distractors_list or '/dev/null'
        )

    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
