#!/usr/bin/env python3
"""
Metrics Calculator for BaitBench

Calculates True Positive, False Positive, True Negative, and False Negative rates
by comparing detected references against expected targets and distractors.
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path


def parse_targets(targets_file):
    """
    Parse targets file and return set of expected target IDs.
    Format: reference_id<tab>weight (or just reference_id per line)
    """
    targets = set()
    with open(targets_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if parts:
                    ref_id = parts[0]
                    targets.add(ref_id)
    return targets


def parse_distractors(distractors_file):
    """
    Parse distractors file and return set of distractor IDs.
    Format: reference_id<tab>weight (or just reference_id per line)
    """
    distractors = set()
    with open(distractors_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if parts:
                    ref_id = parts[0]
                    distractors.add(ref_id)
    return distractors


def parse_detected(detected_file):
    """
    Parse detection list file and return dict of reference_id -> read_count.
    Format: reference_id<tab>count
    """
    detected = {}
    with open(detected_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                parts = line.split('\t')
                if len(parts) >= 2:
                    ref_id = parts[0]
                    count = int(parts[1])
                    detected[ref_id] = count
    return detected


def count_sequences(fasta_file):
    """Count number of sequences in a FASTA file."""
    count = 0
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                count += 1
    return count


def calculate_metrics(targets, distractors, detected):
    """
    Calculate TP/FP/FN/TN metrics.

    Args:
        targets: set of target reference IDs (should be detected)
        distractors: set of distractor reference IDs (should NOT be detected)
        detected: dict of detected reference IDs -> read counts

    Returns:
        dict with metric values
    """
    detected_set = set(detected.keys())

    # Set operations for targets
    true_positives = targets & detected_set       # Targets that were detected
    false_negatives = targets - detected_set      # Targets that were NOT detected

    # Set operations for distractors
    false_positives = distractors & detected_set  # Distractors that were detected
    true_negatives = distractors - detected_set   # Distractors that were NOT detected

    # Also track any unknown detections (not in targets or distractors)
    known = targets | distractors
    unknown_detected = detected_set - known

    tp_count = len(true_positives)
    fp_count = len(false_positives)
    fn_count = len(false_negatives)
    tn_count = len(true_negatives)

    # Calculate rates
    # Sensitivity (True Positive Rate) = TP / (TP + FN)
    sensitivity = tp_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0

    # Specificity (True Negative Rate) = TN / (TN + FP)
    specificity = tn_count / (tn_count + fp_count) if (tn_count + fp_count) > 0 else 0.0

    # Precision = TP / (TP + FP)
    precision = tp_count / (tp_count + fp_count) if (tp_count + fp_count) > 0 else 0.0

    # False Negative Rate = FN / (TP + FN)
    fnr = fn_count / (tp_count + fn_count) if (tp_count + fn_count) > 0 else 0.0

    # False Positive Rate = FP / (FP + TN)
    fpr = fp_count / (fp_count + tn_count) if (fp_count + tn_count) > 0 else 0.0

    # F1 Score = 2 * (precision * sensitivity) / (precision + sensitivity)
    f1_score = 2 * (precision * sensitivity) / (precision + sensitivity) if (precision + sensitivity) > 0 else 0.0

    return {
        'true_positives': sorted(true_positives),
        'false_positives': sorted(false_positives),
        'false_negatives': sorted(false_negatives),
        'true_negatives': sorted(true_negatives),
        'unknown_detected': sorted(unknown_detected),
        'tp_count': tp_count,
        'fp_count': fp_count,
        'fn_count': fn_count,
        'tn_count': tn_count,
        'sensitivity': sensitivity,
        'specificity': specificity,
        'precision': precision,
        'fnr': fnr,
        'fpr': fpr,
        'f1_score': f1_score
    }


def write_summary_tsv(output_file, run_info, capture_stats, metrics):
    """Write summary results to TSV file."""
    headers = [
        'run_name', 'timestamp',
        'num_reads', 'seed',
        'reads_generated', 'reads_captured', 'capture_rate',
        'targets_total', 'distractors_total',
        'tp_count', 'fp_count', 'fn_count', 'tn_count',
        'sensitivity', 'specificity', 'precision', 'f1_score'
    ]

    values = [
        run_info['run_name'],
        run_info['timestamp'],
        run_info['num_reads'],
        run_info['seed'],
        capture_stats['reads_generated'],
        capture_stats['reads_captured'],
        f"{capture_stats['capture_rate']:.4f}",
        metrics['tp_count'] + metrics['fn_count'],  # total targets
        metrics['tn_count'] + metrics['fp_count'],  # total distractors
        metrics['tp_count'],
        metrics['fp_count'],
        metrics['fn_count'],
        metrics['tn_count'],
        f"{metrics['sensitivity']:.4f}",
        f"{metrics['specificity']:.4f}",
        f"{metrics['precision']:.4f}",
        f"{metrics['f1_score']:.4f}"
    ]

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        f.write('\t'.join(str(v) for v in values) + '\n')


def write_detail_tsv(output_file, targets, distractors, detected, metrics):
    """Write per-reference detail to TSV file."""
    headers = ['reference_id', 'category', 'expected', 'detected', 'read_count', 'classification']

    rows = []

    # Add detected references
    for ref_id in detected.keys():
        if ref_id in targets:
            category = 'target'
            classification = 'TP'
        elif ref_id in distractors:
            category = 'distractor'
            classification = 'FP'
        else:
            category = 'unknown'
            classification = 'UNKNOWN'

        rows.append({
            'reference_id': ref_id,
            'category': category,
            'expected': 'true' if ref_id in targets else 'false',
            'detected': 'true',
            'read_count': detected[ref_id],
            'classification': classification
        })

    # Add false negatives (targets not detected)
    for ref_id in metrics['false_negatives']:
        rows.append({
            'reference_id': ref_id,
            'category': 'target',
            'expected': 'true',
            'detected': 'false',
            'read_count': 0,
            'classification': 'FN'
        })

    # Sort by classification priority (TP, FP, FN) then by read_count descending
    classification_order = {'TP': 0, 'FP': 1, 'UNKNOWN': 2, 'FN': 3}
    rows.sort(key=lambda x: (classification_order.get(x['classification'], 99), -x['read_count']))

    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for row in rows:
            f.write('\t'.join([
                row['reference_id'],
                row['category'],
                row['expected'],
                row['detected'],
                str(row['read_count']),
                row['classification']
            ]) + '\n')


def write_json(output_file, run_info, capture_stats, metrics):
    """Write results to JSON file."""
    data = {
        'run_info': run_info,
        'capture_stats': capture_stats,
        'metrics': {
            'tp_count': metrics['tp_count'],
            'fp_count': metrics['fp_count'],
            'fn_count': metrics['fn_count'],
            'tn_count': metrics['tn_count'],
            'sensitivity': metrics['sensitivity'],
            'specificity': metrics['specificity'],
            'precision': metrics['precision'],
            'f1_score': metrics['f1_score'],
            'fnr': metrics['fnr'],
            'fpr': metrics['fpr']
        },
        'details': {
            'true_positives': metrics['true_positives'],
            'false_positives': metrics['false_positives'],
            'false_negatives': metrics['false_negatives'],
            'unknown_detected': metrics['unknown_detected']
        }
    }

    with open(output_file, 'w') as f:
        json.dump(data, f, indent=2)


def main():
    parser = argparse.ArgumentParser(
        description='Calculate TP/FP/FN/TN metrics for BaitBench capture pipeline'
    )

    parser.add_argument('--targets', required=True,
                        help='Targets file with expected reference IDs')
    parser.add_argument('--distractors', required=True,
                        help='Distractors file with background reference IDs')
    parser.add_argument('--detected', required=True,
                        help='Detection list file (reference_id<tab>count)')
    parser.add_argument('--reads', required=True,
                        help='Generated reads FASTA file')
    parser.add_argument('--captured', required=True,
                        help='Captured reads FASTA file')
    parser.add_argument('--run-name', required=True,
                        help='Name for this run')
    parser.add_argument('--num-reads', type=int, required=True,
                        help='Number of reads requested')
    parser.add_argument('--seed', default='NA',
                        help='Random seed used')
    parser.add_argument('--output-summary', required=True,
                        help='Output TSV file for summary metrics')
    parser.add_argument('--output-detail', required=True,
                        help='Output TSV file for per-reference detail')
    parser.add_argument('--output-json',
                        help='Optional: Output JSON file for programmatic access')

    args = parser.parse_args()

    # Parse input files
    print("Parsing targets file...", file=sys.stderr)
    targets = parse_targets(args.targets)
    print(f"  Found {len(targets)} target references", file=sys.stderr)

    print("Parsing distractors file...", file=sys.stderr)
    distractors = parse_distractors(args.distractors)
    print(f"  Found {len(distractors)} distractor references", file=sys.stderr)

    print("Parsing detection list...", file=sys.stderr)
    detected = parse_detected(args.detected)
    print(f"  Found {len(detected)} detected references", file=sys.stderr)

    # Count sequences in FASTA files
    print("Counting sequences...", file=sys.stderr)
    reads_generated = count_sequences(args.reads)
    reads_captured = count_sequences(args.captured)
    capture_rate = reads_captured / reads_generated if reads_generated > 0 else 0.0
    print(f"  Reads generated: {reads_generated}", file=sys.stderr)
    print(f"  Reads captured: {reads_captured}", file=sys.stderr)
    print(f"  Capture rate: {capture_rate:.4f}", file=sys.stderr)

    # Calculate metrics
    print("Calculating metrics...", file=sys.stderr)
    metrics = calculate_metrics(targets, distractors, detected)
    print(f"  True Positives: {metrics['tp_count']}", file=sys.stderr)
    print(f"  False Positives: {metrics['fp_count']}", file=sys.stderr)
    print(f"  False Negatives: {metrics['fn_count']}", file=sys.stderr)
    print(f"  True Negatives: {metrics['tn_count']}", file=sys.stderr)
    print(f"  Sensitivity: {metrics['sensitivity']:.4f}", file=sys.stderr)
    print(f"  Specificity: {metrics['specificity']:.4f}", file=sys.stderr)
    print(f"  Precision: {metrics['precision']:.4f}", file=sys.stderr)
    print(f"  F1 Score: {metrics['f1_score']:.4f}", file=sys.stderr)

    # Prepare run info
    run_info = {
        'run_name': args.run_name,
        'timestamp': datetime.now().isoformat(),
        'num_reads': args.num_reads,
        'seed': args.seed
    }

    capture_stats = {
        'reads_generated': reads_generated,
        'reads_captured': reads_captured,
        'capture_rate': capture_rate
    }

    # Write output files
    print(f"Writing summary to {args.output_summary}...", file=sys.stderr)
    write_summary_tsv(args.output_summary, run_info, capture_stats, metrics)

    print(f"Writing detail to {args.output_detail}...", file=sys.stderr)
    write_detail_tsv(args.output_detail, targets, distractors, detected, metrics)

    if args.output_json:
        print(f"Writing JSON to {args.output_json}...", file=sys.stderr)
        write_json(args.output_json, run_info, capture_stats, metrics)

    print("Done!", file=sys.stderr)


if __name__ == '__main__':
    main()
