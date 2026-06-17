# Reference

## Metrics Definitions

### Genome-Level Metrics

These answer: "Was each genome detected?"

| Metric | Formula | Meaning |
|--------|---------|---------|
| **Sensitivity** | TP / (TP + FN) | Fraction of sample targets that were detected |
| **Specificity** | TN / (TN + FP) | Fraction of non-sample references correctly not detected |
| **Precision** | TP / (TP + FP) | Of detected references, fraction that are sample targets |
| **F1 Score** | 2 * (Precision * Sensitivity) / (Precision + Sensitivity) | Harmonic mean of precision and sensitivity |

Where:
- **TP** = sample targets detected (at least one read maps)
- **FN** = sample targets not detected
- **FP** = FP_target + FP_distractor (non-sample references incorrectly detected)
- **TN** = TN_target + TN_distractor (non-sample references correctly not detected)

### Read-Level Metrics

These track how fragments and reads flow through the pipeline:

| Metric | Description |
|--------|-------------|
| `sample_captured` | Fragments from sample targets that were captured by probes |
| `nonsample_target_captured` | Fragments from non-sample targets that were captured |
| `distractor_captured` | Fragments from distractors that were captured |
| `untargeted_captured` | Fragments from untargeted genomes that were captured (genome mode) |
| `reads_correctly_mapped` | Reads that map back to their source reference |
| `reads_incorrectly_mapped` | Reads that map to a different reference than their source |
| `reads_sequenced` | Number of reads after the sequencing step (may differ from captured if `--num-sequences` is used) |
| `reads_after_filter` | Number of reads after host filtering (0 if `--host-fasta` not provided) |
| `reads_mapped` | Total reads that mapped to any reference (= correctly + incorrectly mapped) |
| `reads_unmapped` | Reads that entered the mapping step but did not map to any reference |

In genome mode, a read from genome G mapping to target T is considered correctly mapped if T is linked to G in the sample-target-map.

Read source is determined from the fragment naming pattern `{seq_id}_fragment_{n}`, using the last `_fragment_` occurrence as the delimiter.

---

## Input File Formats

### FASTA Files

Standard FASTA format. Sequence IDs are the first whitespace-delimited word of the header:

```
>dengue_1 Dengue virus type 1
ATGCTAGCTAGCTAGC...
>zika_virus
GCTAGCTAGCTAGCTA...
```

**Requirements:**
- Sequence IDs must be unique within each file
- Sequence IDs must not contain spaces (use underscores)
- IDs must be consistent across input files (sample manifest IDs must match FASTA headers)

### Sample Manifest Format

The `--sample` flag accepts two formats:

**Inline IDs** (on the command line):

```bash
--sample id1 id2 id3
```

All IDs default to weight 1.0. A number following an ID sets that ID's weight:

```bash
--sample dengue_1 5 zika_virus chikungunya 0.5
# Result: dengue_1=5.0, zika_virus=1.0, chikungunya=0.5
```

**TSV file** (if a single argument that is an existing file):

```
# Optional comment lines starting with #
dengue_1	5.0
zika_virus
chikungunya	0.5
```

- First column: sequence ID (required)
- Second column: weight (optional, defaults to 1.0)
- Empty lines and lines starting with `#` are ignored

In standard mode, IDs must match target FASTA headers. In genome mode, IDs must match genome FASTA headers.

### Sample-Target Map Format

TSV file mapping genome IDs to target IDs (used with `--genomes` via `--sample-target-map`):

```
# genome_id	target_id
e_coli	e_coli_16S
e_coli	e_coli_gyrB
influenza_a	influenza_a
```

- One mapping per line: `genome_id<TAB>target_id`
- Multiple targets per genome supported (one line per mapping)
- Lines starting with `#` are ignored
- Empty lines are ignored

**Auto-linking:** When `--sample-target-map` is omitted (or for genomes not listed in the map), BaitBench auto-links genomes to targets by:

1. **Exact match**: genome ID equals a target ID (e.g., genome `influenza_a` → target `influenza_a`)
2. **Prefix match**: target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` → targets `Bartonella_grahamii|ompB`, `Bartonella_grahamii|16S`)

This means you can name targets using the `organism|gene` convention and genomes using just `organism`, and they will auto-link without needing an explicit map file:

```
# genomes.fa
>Bartonella_grahamii
ATGC...
>Rickettsia_montanensis
ATGC...

# targets.fa (organism|gene naming)
>Bartonella_grahamii|ompB
ATGC...
>Rickettsia_montanensis|ompA
ATGC...
>Rickettsia_montanensis|gltA
ATGC...
```

With this naming, `Bartonella_grahamii` auto-links to `Bartonella_grahamii|ompB`, and `Rickettsia_montanensis` auto-links to both `Rickettsia_montanensis|ompA` and `Rickettsia_montanensis|gltA`.

**Using `--sample-target-map` for non-standard naming:** If your genome and target IDs don't follow either naming convention (exact match or `organism|gene`), provide an explicit mapping file:

```
# mapping.tsv — needed when genome IDs don't match target IDs
NC_012846.1	bartonella_ompB
NC_012846.1	bartonella_16S
GCF_000022725.1	rickettsia_gltA
```

Explicit mappings take precedence over auto-linking for the same genome ID.

**Untargeted genomes:** Sample genomes with no target mapping (explicit or auto-linked) become "untargeted" -- they generate fragments but have no expected target to detect. This models unknown organisms.

**Validation:** BaitBench errors if the map references genome or target IDs not found in their respective FASTA files.

### Groups File Format

TSV file mapping sequence IDs to group names. Used by `--groups` (target grouping) and `--distractor-groups` (distractor grouping):

```
# Optional comment lines starting with #
# seq_id	group_name
West_Nile_virus_0001	West_Nile_virus
West_Nile_virus_0002	West_Nile_virus
West_Nile_virus_0003	West_Nile_virus
Dengue_virus_1_6275	Dengue_virus_1
Dengue_virus_1_2274	Dengue_virus_1
Dengue_virus_2_8773	Dengue_virus_2
```

- One mapping per line: `seq_id<TAB>group_name`
- Lines starting with `#` are ignored
- Empty lines are ignored
- Leading `>` characters on sequence IDs are stripped automatically

**Target groups (`--groups`):**
- Sequence IDs must exist in the targets FASTA (BaitBench errors on unknown IDs)
- Sequences not listed in the file form singleton groups (their own ID = their group name)
- Without `--groups`, each target sequence is its own singleton group (no behavioral change)

**Distractor groups (`--distractor-groups`):**
- Overrides the default automatic grouping by FASTA file stem
- Sequence IDs must exist in the distractor sequences (BaitBench errors on unknown IDs)
- Without `--distractor-groups`, all contigs from each `--distractors` FASTA file are automatically grouped together using the file stem as the group name (e.g., all contigs in `Aaegypti.fa` → group `"Aaegypti"`)

---

## Dependencies

### Rust (build-time, managed by Cargo)

| Crate | Purpose |
|-------|---------|
| clap | CLI argument parsing (derive macros) |
| anyhow | Error handling |
| serde, serde_json | Serialization (JSON output) |
| rand, rand_distr | Random sampling and normal distribution |
| chrono | Timestamps |
| log, env_logger | Logging |

### External (runtime, installed via conda)

| Tool | Purpose | Required? |
|------|---------|-----------|
| minimap2 | Alignment (simulate, mapping, filtering) | Yes |
| BLAST+ | Cross-reactivity analysis (xreact) | Only if `baitbench xreact` is used |
| cd-hit | Sequence clustering (build-probes, tool collapse) | Only if `baitbench build-probes` or `baitbench tool collapse` is used |
| R + packages | HTML report generation | Only if reports are enabled |

Install all via:

```bash
conda env create -f environment.yml
conda activate baitbench
```
