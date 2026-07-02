#!/usr/bin/env python3
"""Generate tutorial example datasets for BaitBench documentation.

Creates two datasets:
  examples/tutorial/         -- standard mode (3 RNA viruses + human distractors)
  examples/tutorial-genome/  -- genome mode (2 bacteria: 16S targets + genome fragments)

Sequences are pseudorandom with per-organism GC content.
Probes are tiled directly from target sequences so they are guaranteed to match.
"""

import os
import random

SEED = 42
random.seed(SEED)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def rand_seq(length: int, gc: float) -> str:
    at = ['A', 'T']
    gc_ = ['G', 'C']
    return ''.join(random.choice(gc_) if random.random() < gc else random.choice(at)
                   for _ in range(length))


def mutate(seq: str, rate: float) -> str:
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice([b for b in bases if b != c]) if random.random() < rate else c
                   for c in seq)


def tile_probes(seq: str, length: int = 120, step: int = 60) -> list[str]:
    return [seq[i:i + length] for i in range(0, len(seq) - length + 1, step)]


def write_fasta(path: str, records: list[tuple[str, str, str]]) -> None:
    """records: list of (seq_id, description, sequence)"""
    with open(path, 'w') as fh:
        for seq_id, desc, seq in records:
            fh.write(f'>{seq_id} {desc}\n')
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + '\n')


def write_tsv(path: str, rows: list[tuple]) -> None:
    with open(path, 'w') as fh:
        for row in rows:
            fh.write('\t'.join(str(c) for c in row) + '\n')


# ---------------------------------------------------------------------------
# Dataset 1: standard mode — RNA viruses
# ---------------------------------------------------------------------------
# Three respiratory/arboviral pathogens. Dengue is in the target panel but
# *not* in the sample manifest, demonstrating discrimination.
#
# GC content references (approximate):
#   Influenza A H3N2 segments: ~42%
#   SARS-CoV-2 genome:         ~38%
#   Dengue virus 2:            ~48%
#   Human genomic DNA:         ~41%

SEQ_LEN = 2000
PROBE_LEN = 120
PROBE_STEP = 60

flu_seq    = rand_seq(SEQ_LEN, gc=0.42)
sars_seq   = rand_seq(SEQ_LEN, gc=0.38)
dengue_seq = rand_seq(SEQ_LEN, gc=0.48)

# Distractors: human genomic fragments, mutated away from targets
human1_seq = mutate(rand_seq(SEQ_LEN, gc=0.41), rate=0.30)
human2_seq = mutate(rand_seq(SEQ_LEN, gc=0.41), rate=0.30)

targets = [
    ('Influenza_A_H3N2',  'Influenza A virus (H3N2) hemagglutinin segment', flu_seq),
    ('SARS_CoV_2',        'Severe acute respiratory syndrome coronavirus 2 spike region', sars_seq),
    ('Dengue_virus_2',    'Dengue virus 2 envelope gene',                   dengue_seq),
]

distractors = [
    ('Human_chr1_frag',   'Homo sapiens chromosome 1 fragment', human1_seq),
    ('Human_chr22_frag',  'Homo sapiens chromosome 22 fragment', human2_seq),
]

probes = []
for seq_id, _, seq in targets:
    for j, p in enumerate(tile_probes(seq, PROBE_LEN, PROBE_STEP), 1):
        probes.append((f'{seq_id}_probe_{j:03d}', f'Probe {j} targeting {seq_id}', p))

os.makedirs('tutorial', exist_ok=True)
write_fasta('tutorial/targets.fa',    targets)
write_fasta('tutorial/distractors.fa', distractors)
write_fasta('tutorial/probes.fa',     probes)

# Sample manifest: Influenza and SARS-CoV-2 are the sample; Dengue is not.
# Weights reflect approximate relative abundance (equal here for simplicity).
write_tsv('tutorial/sample.tsv', [
    ('Influenza_A_H3N2', 1.0),
    ('SARS_CoV_2',       1.0),
])

print(f"tutorial/: {len(targets)} targets, {len(distractors)} distractors, {len(probes)} probes")
print(f"  sample.tsv: 2 of 3 targets in sample (Dengue_virus_2 is non-sample)")


# ---------------------------------------------------------------------------
# Dataset 2: genome mode — bacteria with 16S rRNA targets
# ---------------------------------------------------------------------------
# Two clinically relevant bacteria. 16S rRNA genes are the sequencing targets;
# full genome fragments are the simulation source.
# Each genome = upstream_flank + 16S_sequence + downstream_flank.
#
# GC content references (approximate):
#   M. tuberculosis genome:  ~65%   16S: ~54%
#   S. aureus genome:        ~33%   16S: ~37%
#   Human mtDNA (distractor):~44%

S16_LEN    = 1500   # typical 16S rRNA length
FLANK_LEN  = 3500   # genomic flanking sequence on each side
PROBE_STEP_16S = 40  # tighter tiling for 16S probes

mtb_16s    = rand_seq(S16_LEN, gc=0.54)
staph_16s  = rand_seq(S16_LEN, gc=0.37)

mtb_genome   = rand_seq(FLANK_LEN, gc=0.65) + mtb_16s   + rand_seq(FLANK_LEN, gc=0.65)
staph_genome = rand_seq(FLANK_LEN, gc=0.33) + staph_16s + rand_seq(FLANK_LEN, gc=0.33)

# Distractors: human mitochondrial-like fragment
human_mt_seq = mutate(rand_seq(SEQ_LEN, gc=0.44), rate=0.25)

genome_targets = [
    ('Mycobacterium_tuberculosis_H37Rv_16S', 'Mycobacterium tuberculosis H37Rv 16S rRNA gene', mtb_16s),
    ('Staphylococcus_aureus_MRSA252_16S',    'Staphylococcus aureus MRSA252 16S rRNA gene',    staph_16s),
]

genomes = [
    ('Mycobacterium_tuberculosis_H37Rv',
     'Mycobacterium tuberculosis H37Rv representative genome fragment',
     mtb_genome),
    ('Staphylococcus_aureus_MRSA252',
     'Staphylococcus aureus MRSA252 representative genome fragment',
     staph_genome),
]

genome_distractors = [
    ('Human_mtDNA_frag', 'Homo sapiens mitochondrial DNA fragment', human_mt_seq),
]

genome_probes = []
for seq_id, _, seq in genome_targets:
    for j, p in enumerate(tile_probes(seq, PROBE_LEN, PROBE_STEP_16S), 1):
        genome_probes.append((f'{seq_id}_probe_{j:03d}', f'Probe {j} targeting {seq_id}', p))

# sample-target-map: genome ID → 16S target ID (one-to-one here)
sample_target_map = [
    ('Mycobacterium_tuberculosis_H37Rv',  'Mycobacterium_tuberculosis_H37Rv_16S'),
    ('Staphylococcus_aureus_MRSA252',     'Staphylococcus_aureus_MRSA252_16S'),
]

# Sample manifest: both bacteria in the sample
genome_sample = [
    ('Mycobacterium_tuberculosis_H37Rv', 1.0),
    ('Staphylococcus_aureus_MRSA252',    1.0),
]

os.makedirs('tutorial-genome', exist_ok=True)
write_fasta('tutorial-genome/targets.fa',     genome_targets)
write_fasta('tutorial-genome/genomes.fa',     genomes)
write_fasta('tutorial-genome/distractors.fa', genome_distractors)
write_fasta('tutorial-genome/probes.fa',      genome_probes)
write_tsv('tutorial-genome/sample_target_map.tsv', sample_target_map)
write_tsv('tutorial-genome/sample.tsv',            genome_sample)

print(f"\ntutorial-genome/: {len(genomes)} genomes, {len(genome_targets)} 16S targets, "
      f"{len(genome_distractors)} distractors, {len(genome_probes)} probes")
print(f"  sample_target_map.tsv: {len(sample_target_map)} genome→target mappings")
