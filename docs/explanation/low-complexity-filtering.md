# Low-Complexity Filtering (sDUST)

BaitBench uses the sDUST algorithm ([Morgulis et al. 2006](https://pubmed.ncbi.nlm.nih.gov/16796549/)) to detect and remove probes whose sequences contain repetitive or low-complexity regions. This page explains why that matters, how the algorithm works, how the three parameters interact, and how to choose appropriate cutoffs for your panel.

---

## Why Filter Low-Complexity Probes?

Probe sequences like `AAAAAAAAAAAAAAAAAAAA` or `ATATATATATATATATAT` are chemically valid but problematic in capture experiments:

- **Non-specific binding**: Simple repeats are ubiquitous in eukaryotic genomes. A poly-A probe will bind anything with a poly-T stretch — host DNA, adapters, other viruses — producing false-positive capture and wasted sequencing capacity.
- **Poor thermodynamic discrimination**: Low-complexity probes tend to form stable but non-specific duplexes; thermodynamic scoring cannot reliably distinguish real targets from background.
- **Parasitic amplification**: In PCR-based library preparation, primer-like repeats can seed off-target amplification.

The goal of low-complexity filtering is to remove probes whose sequences are so repetitive that they cannot reliably discriminate between their intended target and everything else.

---

## How sDUST Works

sDUST scores the complexity of a sequence window by asking: **how non-uniformly distributed are the 3-mers in this window?**

There are 64 possible 3-mers (trinucleotides) over the four DNA bases. In a random sequence, all 64 appear roughly equally. In a poly-A stretch, only one appears — `AAA`. In a dinucleotide repeat like `ATATAT...`, only two appear — `ATA` and `TAT`. sDUST formalises this intuition into a score.

### The Score Formula

For a window containing *L* overlapping trinucleotides, let *c_i* be the count of the *i*-th trinucleotide type. The DUST score is:

```
score = Σ c_i × (c_i − 1) / 2
        ─────────────────────────
               L − 1
```

This is the number of **coincident trinucleotide pairs** (pairs of positions with the same 3-mer) divided by the maximum possible. Intuitively, a high score means many trinucleotide slots are occupied by the same few types — a hallmark of repetitiveness.

### Score Examples

The default window is 64 bases, which contains 62 overlapping trinucleotides (*L* = 62):

| Sequence type | Score | Interpretation |
|---------------|-------|----------------|
| `AAAA...` (poly-A) | ~31 | Only 1 distinct trinucleotide; maximally repetitive |
| `ATAT...` (dinucleotide repeat) | ~4.2 | 2 distinct trinucleotides; clearly low-complexity |
| `ACGACG...` (trinucleotide repeat) | ~2.1 | 3 distinct trinucleotides; borderline at default threshold |
| `ACGTACGT...` (4-mer repeat) | ~1.4 | Not masked at default threshold |
| Random DNA | ~0.5 | All 64 trinucleotides roughly equal; high-complexity |

A region is flagged as low-complexity when its score exceeds the threshold *T* (default 2.0). The algorithm slides this window across the probe and collects all flagged regions.

### The Sliding Window

sDUST does not score a fixed window and move on — it finds the highest-scoring **suffix** of each window position. This makes it sensitive to sub-window repetitive runs that are flanked by complex sequence. Concretely:

- A 20 bp poly-A tail on an otherwise complex 120 bp probe will be correctly identified as low-complexity even though the window that includes it also contains flanking complex sequence.
- The output is a set of masked intervals on the probe sequence, not a single per-probe score.

Non-ACGT bases (including N) reset the window: they are treated as sequence breaks, so an N-heavy stretch is scored independently on each side.

---

## From Masking to Filtering: Three Parameters

BaitBench uses three parameters to control how sDUST is applied during `build-probes`:

| Parameter | Default | What it controls |
|-----------|---------|-----------------|
| `--dust-threshold` | 2.0 | Score *T* above which a window is considered low-complexity |
| `--dust-window` | 64 | Window size *W* in bases |
| `--max-masked-frac` | 0.25 | Maximum fraction of a probe's bases that can be low-complexity before the probe is rejected |

### `--dust-threshold` (*T*)

This is the original DUST threshold parameter. The default of 2.0 comes from the Morgulis et al. (2006) paper and is widely used in bioinformatics tools (NCBI BLAST, RepeatMasker, seqtk).

| T value | Effect |
|---------|--------|
| 1.0 | Stricter — also catches weak repeats like 4-mer repeats |
| 2.0 | Default — catches mono-, di-, and trinucleotide repeats |
| 3.0 | Permissive — only strong repeats (mono- and dinucleotide) |
| 5.0 | Very permissive — only catches near-perfect mononucleotide runs |

**Recommendation**: Leave `--dust-threshold` at 2.0 for almost all use cases. The original authors calibrated this value empirically and it is the de facto standard across the field. Use `--max-masked-frac` to tune how strict the filtering is, not the threshold.

### `--dust-window` (*W*)

The window size controls how far back the algorithm looks when computing local complexity. Larger windows dilute short repetitive runs (a 10 bp poly-A inside a 200 bp window barely moves the score); smaller windows are more sensitive to short runs.

For typical probe lengths (80–150 bp), the default of 64 is well-matched. You might consider:

- **Longer probes (>200 bp)**: increase to 128 or 256 to score the full probe in one pass and avoid fragmented masking at window boundaries.
- **Very short probes (<50 bp)**: decrease to 32 so the window is not larger than the probe itself.

### `--max-masked-frac`

This is the primary knob for controlling stringency in BaitBench. It sets the maximum fraction of a probe's bases that can fall inside sDUST-flagged regions before the probe is discarded.

| `--max-masked-frac` | Effect |
|---------------------|--------|
| 0.0 | Reject any probe with even one masked base |
| 0.10 | Reject probes with more than 12 masked bases (of 120 bp) |
| 0.25 _(default)_ | Reject probes with more than 30 masked bases (of 120 bp) |
| 0.50 | Only reject probes that are mostly low-complexity |
| 1.0 | Disable filtering entirely |

---

## Practical Examples

### Example 1: Poly-A tail

A 120 bp probe with a 40 bp poly-A tail:

```
ACGTCGATCGGCATGCATCG...ACGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
└──────── 80 bp complex ──────────┘└────────── 40 bp poly-A ──────────┘
```

- sDUST flags the 40 bp poly-A region (score ≈ 31, far above T=2.0)
- Masked fraction = 40/120 = 0.33
- At `--max-masked-frac 0.25`: **rejected** (0.33 > 0.25)
- At `--max-masked-frac 0.50`: **kept** (0.33 < 0.50)

### Example 2: Short dinucleotide run in an otherwise complex probe

A 120 bp probe with a 20 bp `ATATAT...` run in the middle:

- Masked fraction = 20/120 = 0.17
- At default `--max-masked-frac 0.25`: **kept** (0.17 < 0.25)
- At `--max-masked-frac 0.10`: **rejected** (0.17 > 0.10)

### Example 3: Pure trinucleotide repeat probe

A 120 bp probe made entirely of `ACGACG...`:

- Score ≈ 2.1 at T=2.0 → entire probe flagged
- Masked fraction = 1.0
- **Rejected** at any `--max-masked-frac` below 1.0

### Example 4: `ACGTACGT...` repeat

A 120 bp probe made of repeating `ACGT`:

- Score ≈ 1.4 at T=2.0 → **not flagged** (below threshold)
- Masked fraction = 0
- **Kept** at any threshold

This shows the algorithm is not overly conservative: regular but non-trivial repeats pass through.

---

## Visualising Masking

The [`baitbench tool dustview`](../reference/commands.md#tool-dustview) command shows exactly which regions of a FASTA file would be masked by sDUST, using the same threshold and window parameters. This is useful for auditing probes that were unexpectedly removed, or for choosing `--max-masked-frac` based on the actual masking patterns in your probe set.

---

## Choosing Cutoffs

**Start with the defaults.** The defaults (T=2.0, W=64, max-masked-frac=0.25) are appropriate for most probe panels. They catch probes that would genuinely cause non-specific binding while keeping probes with minor repetitive content.

**Tighten `--max-masked-frac` if** you are targeting sequences that are intrinsically repetitive (e.g., viral genomes rich in dinucleotide repeats) and you want to maximise specificity, even at the cost of reduced coverage.

**Loosen `--max-masked-frac` (or set to 1.0 to disable) if** your targets require probes in regions that are unavoidably low-complexity — for example, tandem repeat regions that are nonetheless diagnostically relevant. Inspect the `assess-probes` coverage report to see whether those regions end up covered.

**Adjust `--dust-threshold` only if** you have specific reasons to deviate from the standard. Lower it (e.g., 1.5) to catch moderate 4-mer repeats; raise it (e.g., 3.0) if you are getting false positives from sequences that BLAST handles without issue.
