# Prepare Input Files

This guide covers the format requirements for every input file BaitBench accepts. Get the formats right before running anything else.

---

## FASTA files (targets, distractors, probes, genomes)

All FASTA files follow standard format. The sequence ID is the **first whitespace-delimited word** of the header line — everything after `>` up to the first space:

```
>Dengue_virus_2 Dengue virus type 2 complete genome
ATGCTAGCTAGCTAGC...
```

Here the sequence ID is `Dengue_virus_2`. The rest of the header (`Dengue virus type 2 complete genome`) is ignored.

**Requirements for all FASTA files:**

- Sequence IDs must be **unique within each file**
- Sequence IDs must **not contain spaces** — use underscores: `>Zika_virus`, not `>Zika virus`
- IDs used in TSV files (sample manifest, groups files, sample-target-map) must exactly match the FASTA header IDs

**What goes in each file:**

| File | Flag | Contents |
|------|------|----------|
| Targets | `--targets` | Sequences your probes are designed to capture. Expected to be detected. |
| Distractors | `--distractors` | Background sequences that should NOT be captured. Can specify multiple FASTA files. |
| Probes | `--probes` | The probe sequences you want to evaluate. |
| Genomes | `--genomes` | Full genomes for fragment generation (genome mode only). See [Genome Mode tutorial](../tutorials/genome-mode.md). |

---

## Sample manifest (`--sample`)

The sample manifest specifies which targets (or genomes in genome mode) are present in the simulated specimen. Everything else becomes a non-sample target — present in the probe panel but not in the specimen, and therefore a true negative.

Without `--sample`, all targets are treated as sample with equal weight.

### TSV file format

```
# id	weight
Influenza_A_H3N2	1.0
SARS_CoV_2	1.0
Dengue_virus_2	0.5
```

- First column: sequence ID (must match a FASTA header)
- Second column: relative abundance weight (optional, defaults to 1.0)
- Lines starting with `#` and blank lines are ignored

Higher weights generate proportionally more fragments from that sequence. A weight of 2.0 generates twice as many fragments as weight 1.0.

### Inline IDs (on the command line)

```bash
--sample Influenza_A_H3N2 SARS_CoV_2 Dengue_virus_2
```

All IDs default to weight 1.0. Place a number after an ID to set its weight:

```bash
--sample Influenza_A_H3N2 2.0 SARS_CoV_2 Dengue_virus_2 0.5
# Result: Influenza_A_H3N2=2.0, SARS_CoV_2=1.0, Dengue_virus_2=0.5
```

---

## Sample-target map (`--sample-target-map`)

Required in genome mode when genome IDs and target IDs do not follow the auto-linking naming convention.

```
# genome_id	target_id
Mycobacterium_tuberculosis_H37Rv	Mycobacterium_tuberculosis_H37Rv_16S
Staphylococcus_aureus_MRSA252	Staphylococcus_aureus_MRSA252_16S
Bartonella_grahamii	Bartonella_grahamii|ompB
Bartonella_grahamii	Bartonella_grahamii|gltA
```

- One mapping per line: `genome_id<TAB>target_id`
- Multiple targets per genome are supported (one line each)
- Lines starting with `#` and blank lines are ignored

**Auto-linking (when `--sample-target-map` is omitted):** BaitBench links genomes to targets automatically if:

1. The genome ID exactly matches a target ID, OR
2. A target ID starts with `{genome_id}|` (e.g., genome `Bartonella_grahamii` → targets `Bartonella_grahamii|ompB`, `Bartonella_grahamii|gltA`)

If your IDs don't follow either convention, provide an explicit map.

---

## Groups file (`--groups`)

Maps individual target sequence IDs to logical group names. Use this when multiple sequence variants represent the same organism and should be counted as one entity.

```
# seq_id	group_name
West_Nile_virus_0001	West_Nile_virus
West_Nile_virus_0002	West_Nile_virus
West_Nile_virus_0003	West_Nile_virus
Dengue_virus_1_v1	Dengue_virus_1
Dengue_virus_1_v2	Dengue_virus_1
```

With this file, a detection of any `West_Nile_virus_*` variant counts as detecting the `West_Nile_virus` group. Sequences not listed in the file form singleton groups (their own ID = group name) and behave as before.

**When to use it:** Your target panel has many strain variants of the same species. Without `--groups`, detecting any one variant counts as TP, but detecting a different variant of the same species when only the first was in the sample would count as FP_target. With `--groups`, all variants collapse to one group.

---

## Distractor groups (`--distractor-groups`)

By default, all sequences from each `--distractors` FASTA file are automatically grouped under the file stem. All contigs in `Aedes_aegypti.fa` → group `Aedes_aegypti`.

If you have multiple organisms in one distractor FASTA, provide an explicit map:

```
# contig_id	group_name
scaffold_001	Aedes_aegypti
scaffold_002	Aedes_aegypti
chr1	Homo_sapiens
chr2	Homo_sapiens
```

Without `--distractor-groups`, every sequence in a multi-organism FASTA would be grouped as one entity, making it impossible to distinguish which organism produced the false positive.

---

## Sequence ID naming tips

- **No spaces:** `>Zika_virus` not `>Zika virus` — the space truncates the ID to `Zika`
- **No special characters that break shell globbing:** avoid `*`, `?`, `[`, `]`
- **Consistency across files:** if your FASTA says `>flu_H3N2`, your sample manifest must also say `flu_H3N2` — not `Flu_H3N2` or `flu-H3N2`
- **Fragment naming convention:** fragments are named `{seq_id}_fragment_{n}`. If your sequence ID contains `_fragment_`, BaitBench uses the *last* occurrence as the delimiter, so this is safe but worth knowing

---

## Checking your files before running

```bash
# Count sequences in a FASTA
grep -c ">" targets.fa

# List all sequence IDs
grep ">" targets.fa | awk '{print $1}' | sed 's/>//'

# Check for duplicate IDs
grep ">" targets.fa | awk '{print $1}' | sort | uniq -d

# Check that sample.tsv IDs exist in targets.fa
cut -f1 sample.tsv | grep -v "^#" | while read id; do
  grep -q "^>$id" targets.fa || echo "Missing: $id"
done
```

If BaitBench reports "Sample ID 'X' not found in targets FASTA", the ID in your manifest does not match any header in the targets file — check for spaces, case differences, or special characters.
