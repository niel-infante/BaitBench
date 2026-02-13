# BaitBench - Claude Code Guide

## Project Overview

BaitBench is a generic tool for testing probe capture efficiency via in-silico simulation. Users provide probesets, target genomes, and distractor genomes to evaluate how well probes capture intended targets while avoiding off-target sequences.

## Architecture

### Pipeline Flow
```
targets.fa + distractors.fa
         ↓
   PREPARE_REFERENCE (combine, generate weights)
         ↓
   GENERATE_READS (weighted random fragments)
         ↓
   CAPTURE (minimap2 or BLAST)
         ↓
   FILTER_HOST (optional)
         ↓
   MAP_READS (back to references)
         ↓
   GENERATE_LIST (count reads per reference)
         ↓
   CALCULATE_METRICS (TP/FP/FN/TN)
         ↓
   GENERATE_REPORT (HTML summary)
```

### Key Files

| File | Purpose |
|------|---------|
| `main.nf` | Nextflow DSL2 pipeline - orchestrates all processes |
| `nextflow.config` | Default parameters and execution profiles |
| `bin/fasta_sampler.py` | Generates weighted random fragments from FASTA |
| `bin/metrics.py` | Calculates TP/FP/FN/TN and derived metrics |
| `bin/prepare_reference.py` | Combines targets/distractors, generates weights |
| `bin/generate_report.py` | Creates HTML report from results |
| `environment.yml` | Conda environment with all dependencies |

### Metrics Definitions

- **Targets**: Genomes in `targets.fa` - probes SHOULD capture these
- **Distractors**: Genomes in `distractors.fa` - probes should NOT capture these
- **TP (True Positive)**: Target genome detected
- **FP (False Positive)**: Distractor genome detected
- **FN (False Negative)**: Target genome NOT detected
- **TN (True Negative)**: Distractor genome NOT detected

## Development Guidelines

### Nextflow Conventions
- Use DSL2 syntax with process definitions
- Processes should emit named outputs for clarity
- Use `publishDir` to copy outputs to results directory
- Scripts in `bin/` are automatically added to PATH by Nextflow

### Python Scripts
- All scripts should be executable (`chmod +x`)
- Use `#!/usr/bin/env python3` shebang
- Accept inputs via argparse CLI arguments
- Print progress/status to stderr, data to stdout or files
- No external dependencies beyond standard library + jinja2

### Testing Changes
```bash
# Activate environment
conda activate baitbench

# Run with minimal example
nextflow run main.nf \
  --targets examples/minimal/targets.fa \
  --distractors examples/minimal/distractors.fa \
  --probes examples/minimal/probes.fa \
  --num_reads 1000 \
  --outdir test_results

# Check outputs
ls test_results/
cat test_results/results.tsv
```

### Common Modifications

**Adding a new metric**: Edit `bin/metrics.py`, update `calculate_metrics()` function and TSV output headers.

**Changing capture parameters**: Edit `nextflow.config` default params or pass via CLI `--param_name value`.

**Adding a new process**: Add process definition in `main.nf`, wire into workflow at bottom of file.

**Modifying read generation**: Edit `bin/fasta_sampler.py` - fragment length, sampling strategy, etc.

## Origin

This tool was extracted from the MTEC_probes repository (arbovirus probe design project). The original in-silico pipeline was generalized to work with any organism/probe combination.

Key differences from MTEC_probes:
- Separate targets/distractors input model (vs single weights file)
- Generic reference_id naming (vs virus_id)
- Added TN/specificity metrics
- HTML report generation
- No hardcoded species or paths
