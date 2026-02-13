#!/usr/bin/env python3
"""
FASTA Fragment Sampler for BaitBench

Generates random fragments from multi-sequence FASTA with weighted sampling.
Sequences are weighted by both their specified weight AND sequence length.
"""

import argparse
import random
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)


def parse_fasta(fasta_file):
    """Parse FASTA file and return dictionary of sequences."""
    sequences = {}
    current_id = None
    current_seq = []

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Get first word after '>'
                current_seq = []
            elif line:
                current_seq.append(line.upper())

        # Don't forget the last sequence
        if current_id is not None:
            sequences[current_id] = ''.join(current_seq)

    return sequences


def parse_weights(weights_file):
    """Parse weights file and return dictionary of weights."""
    weights = {}

    with open(weights_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 2:
                    seq_id = parts[0]
                    weight = float(parts[1])
                    weights[seq_id] = weight

    return weights


def weighted_choice(seq_ids, weights):
    """Select sequence ID based on weights."""
    return random.choices(seq_ids, weights=weights, k=1)[0]


def generate_fragments(sequences, weights, num_fragments, output_file,
                       read_length_mean=175, read_length_min=150, read_length_max=200):
    """Generate random fragments with weighted sampling."""

    # Filter sequences that have positive weights
    weighted_seqs = {sid: sequences[sid] for sid in sequences if weights.get(sid, 0) > 0}

    if not weighted_seqs:
        logger.error("No sequences with positive weights found!")
        sys.exit(1)

    # Prepare for weighted sampling (weight by both specified weight AND length)
    seq_ids = list(weighted_seqs.keys())
    seq_weights = [weights[sid] * len(weighted_seqs[sid]) for sid in seq_ids]

    logger.info(f"Loaded {len(sequences)} sequences from FASTA")
    logger.info(f"Found weights for {len(weights)} sequences")
    logger.info(f"Generating {num_fragments} fragments from {len(weighted_seqs)} weighted sequences")
    logger.info(f"Read length: {read_length_min}-{read_length_max} (mean {read_length_mean})")

    # Fragment length parameters (normal distribution)
    sd_length = (read_length_max - read_length_min) / 6  # ~99.7% within 3 standard deviations

    fragments_generated = 0
    source_counts = {}

    with open(output_file, 'w') as out:
        for i in range(num_fragments):
            # Select sequence based on weights
            seq_id = weighted_choice(seq_ids, seq_weights)
            seq = weighted_seqs[seq_id]

            # Track source counts
            source_counts[seq_id] = source_counts.get(seq_id, 0) + 1

            # Generate fragment length (normal distribution, clamped)
            frag_len = int(random.gauss(read_length_mean, sd_length))
            frag_len = max(read_length_min, min(read_length_max, frag_len))

            # Skip if sequence is too short
            if len(seq) < frag_len:
                continue

            # Random start position
            max_start = len(seq) - frag_len
            start = random.randint(0, max_start)

            # Extract fragment
            fragment = seq[start:start + frag_len]

            # Write to output
            header = f">{seq_id}_fragment_{i + 1} start={start} length={frag_len}"
            out.write(f"{header}\n{fragment}\n")

            fragments_generated += 1

            # Progress indicator
            if (i + 1) % 10000 == 0:
                logger.info(f"Generated {i + 1} fragments...")

    logger.info(f"Successfully generated {fragments_generated} fragments to {output_file}")

    # Log source distribution
    logger.info("Source distribution:")
    for sid in sorted(source_counts.keys()):
        logger.info(f"  {sid}: {source_counts[sid]} reads")

    return fragments_generated


def main():
    parser = argparse.ArgumentParser(
        description='Generate random fragments from FASTA with weighted sampling',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  %(prog)s -f sequences.fasta -w weights.txt -n 1000 -o fragments.fasta
  %(prog)s --fasta input.fa --weights abundances.txt --num 5000 --output out.fa

Weights file format (tab or space separated):
  seq1    10.5
  seq2    5.0
  seq3    25.3
        '''
    )

    parser.add_argument('-f', '--fasta', required=True,
                        help='Input FASTA file with multiple sequences')
    parser.add_argument('-w', '--weights', required=True,
                        help='Weights file (format: seqID weight)')
    parser.add_argument('-n', '--num', type=int, required=True,
                        help='Number of fragments to generate')
    parser.add_argument('-o', '--output', required=True,
                        help='Output FASTA file for fragments')
    parser.add_argument('-s', '--seed', type=int, default=None,
                        help='Random seed for reproducibility (optional)')
    parser.add_argument('--read-length-mean', type=int, default=175,
                        help='Mean read length (default: 175)')
    parser.add_argument('--read-length-min', type=int, default=150,
                        help='Minimum read length (default: 150)')
    parser.add_argument('--read-length-max', type=int, default=200,
                        help='Maximum read length (default: 200)')

    args = parser.parse_args()

    # Set random seed if provided
    if args.seed is not None:
        random.seed(args.seed)
        logger.info(f"Using random seed: {args.seed}")

    # Parse input files
    logger.info("Parsing FASTA file...")
    sequences = parse_fasta(args.fasta)

    logger.info("Parsing weights file...")
    weights = parse_weights(args.weights)

    # Generate fragments
    generate_fragments(
        sequences, weights, args.num, args.output,
        read_length_mean=args.read_length_mean,
        read_length_min=args.read_length_min,
        read_length_max=args.read_length_max
    )


if __name__ == '__main__':
    main()
