# Input Formats

Precise specifications for all files accepted by BaitBench.

---

## FASTA Files

Standard FASTA format. The sequence ID is the first whitespace-delimited word of the header:

```
>dengue_1 Dengue virus type 1
ATGCTAGCTAGCTAGC...
>zika_virus
GCTAGCTAGCTAGCTA...
```

**Requirements:**

- Sequence IDs must be unique within each file
- Sequence IDs must not contain spaces (use underscores instead)
- IDs must be consistent across related input files: sample manifest IDs must match FASTA headers, sample-target-map IDs must match genome and target FASTA headers

Accepted by `--targets`, `--genomes`, `--distractors`, `--probes`, `--host-fasta`, `--filter-genomes`.

---

## Sample Manifest

The `--sample` flag accepts two formats:

### Inline IDs (on the command line)

```bash
--sample id1 id2 id3
```

All IDs default to weight 1.0. A numeric value following an ID sets that ID's weight:

```bash
--sample dengue_1 5 zika_virus chikungunya 0.5
# Result: dengue_1=5.0, zika_virus=1.0, chikungunya=0.5
```

### TSV file

If a single argument is given that is the path to an existing file, it is parsed as a TSV manifest:

```
# Optional comment lines starting with #
dengue_1	5.0
zika_virus
chikungunya	0.5
```

- First column: sequence ID (required)
- Second column: weight (optional, defaults to 1.0)
- Empty lines and lines starting with `#` are ignored
- Tab-separated

**Standard mode:** IDs must match target FASTA headers.  
**Genome mode:** IDs must match genome FASTA headers.

---

## Sample-Target Map

TSV file mapping genome IDs to target IDs. Used with `--genomes` via `--sample-target-map`:

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

### Auto-linking

When `--sample-target-map` is omitted, BaitBench auto-links genomes to targets by:

1. **Exact match**: genome ID equals a target ID (e.g., genome `influenza_a` → target `influenza_a`)
2. **Prefix match**: target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` → targets `Bartonella_grahamii|ompB`, `Bartonella_grahamii|16S`)

The `organism|gene` naming convention enables auto-linking without an explicit map:

```
# genomes.fa
>Bartonella_grahamii
ATGC...

# targets.fa
>Bartonella_grahamii|ompB
ATGC...
>Bartonella_grahamii|16S
ATGC...
```

With this naming, `Bartonella_grahamii` auto-links to both `Bartonella_grahamii|ompB` and `Bartonella_grahamii|16S`.

### Untargeted genomes

Sample genomes with no matching target (via explicit map or auto-linking) become "untargeted" — they generate fragments but are not expected to produce reads mapping to any target. These are tracked separately and excluded from TP/FP/FN/TN counts.

**Validation:** BaitBench errors if the map references genome or target IDs not found in their respective FASTA files.

---

## Groups File

TSV file mapping sequence IDs to group names. Used by `--groups` (target grouping) and `--distractor-groups` (distractor grouping):

```
# Optional comment lines starting with #
# seq_id	group_name
West_Nile_virus_0001	West_Nile_virus
West_Nile_virus_0002	West_Nile_virus
West_Nile_virus_0003	West_Nile_virus
Dengue_virus_1_6275	Dengue_virus_1
Dengue_virus_2_8773	Dengue_virus_2
```

- One mapping per line: `seq_id<TAB>group_name`
- Lines starting with `#` are ignored
- Empty lines are ignored
- Leading `>` characters on sequence IDs are stripped automatically

### Target groups (`--groups`)

- Sequence IDs must exist in the targets FASTA (BaitBench errors on unknown IDs)
- Sequences not listed form singleton groups (their own ID is the group name)
- Without `--groups`, each target sequence is its own singleton group

### Distractor groups (`--distractor-groups`)

- Overrides the default automatic grouping by FASTA file stem
- Without `--distractor-groups`, all contigs from each `--distractors` FASTA file are grouped together using the file stem (e.g., all contigs in `Aaegypti.fa` → group `"Aaegypti"`)

---

## CT Calibration Points

The `--ct-calibration` flag accepts two quoted strings, each in `"CT,fraction"` format:

```bash
--ct-calibration "20,0.01" "25,0.0003"
```

- CT value and target fraction separated by a comma, no spaces
- Must provide exactly two points
- Both CT values must be distinct
- Fractions must be in (0, 1]

---

## Dependencies

### External (runtime, installed via conda)

| Tool | Purpose | Required? |
|------|---------|-----------|
| blastn (BLAST+) | Cross-reactivity analysis (`xreact`) | Only if `baitbench xreact` is used |
| cd-hit-est | Sequence clustering (`build-probes`, `tool collapse`) | Only if `baitbench build-probes` or `baitbench tool collapse` is used |
| R + packages | HTML report generation | Only if `--report full` is used |
| art_modern | Illumina read simulation | Only if `--read-simulator art` is used |
| badread | ONT / PacBio read simulation | Only if `--read-simulator badread` is used |

Alignment is handled by the rammap library compiled into the BaitBench binary — no external minimap2 installation required.

Install via:

```bash
conda env create -f environment.yml
conda activate baitbench
```

### Rust crates (build-time, managed by Cargo)

| Crate | Purpose |
|-------|---------|
| clap | CLI argument parsing |
| anyhow | Error handling |
| serde, serde_json | Serialization (JSON output) |
| rand, rand_distr | Random sampling and distributions |
| chrono | Timestamps |
| log, env_logger | Logging |
