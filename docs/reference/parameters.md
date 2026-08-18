# Parameters

Complete parameter reference for `baitbench run` and related subcommands.

---

## Input Files

| Parameter | Flag | Default | Applies to | Description |
|-----------|------|---------|------------|-------------|
| Targets | `--targets`, `-t` | required | run, prepare, probe-coverage, coverage-curve | FASTA of target sequences the probes are designed to capture |
| Genomes | `--genomes`, `-g` | none | run, prepare, coverage-curve | FASTA of full genomes for fragment generation (genome mode) |
| Distractors | `--distractors`, `-d` | required | run, prepare, coverage-curve | FASTA of background sequences. Can be specified multiple times. |
| Probes | `--probes`, `-p` | required | run, coverage-curve | FASTA of probe sequences |
| Sample | `--sample` | all targets | run, coverage-curve | Sample targets or genomes: TSV file path OR inline IDs with optional weights. See [Sample Manifest](input-formats.md#sample-manifest) |
| Sample-target map | `--sample-target-map` | none | run, prepare, coverage-curve | TSV mapping genome IDs to target IDs (genome mode). See [Sample-Target Map](input-formats.md#sample-target-map) |
| Groups | `--groups` | none | run, prepare | TSV mapping target sequence IDs to group names. See [Groups File](input-formats.md#groups-file) |
| Distractor groups | `--distractor-groups` | none | run, prepare | TSV mapping distractor sequence IDs to group names (overrides default file-stem grouping). See [Groups File](input-formats.md#groups-file) |
| Host FASTA | `--host-fasta` | none | run, coverage-curve | Host genome for read filtering |

---

## Fragment Generation

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Num fragments | `--num-fragments` | 10000 | Fragments to generate. More = better statistical power, slower. |
| Fragment length mean | `--fragment-length-mean` | 175 | Mean fragment length in bp |
| Fragment length min | `--fragment-length-min` | 150 | Minimum fragment length in bp |
| Fragment length max | `--fragment-length-max` | 200 | Maximum fragment length in bp |

Fragment lengths follow a truncated normal distribution clamped to [min, max].

---

## Target Abundance

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Distractor fraction | `--distractor-fraction` (`-f` on `prepare` only) | 0.9 | Fraction of fragments from distractors (0–1). Higher = lower target abundance. **Mutually exclusive with `--ct`** |
| CT score | `--ct` | none | qPCR CT value. Converted to distractor fraction via calibration formula. Lower CT = more target. **Mutually exclusive with `--distractor-fraction`** |

If neither is specified, defaults to `--distractor-fraction 0.9`.

---

## CT Score Parameters

These calibrate the CT-to-fraction conversion. Only used when `--ct` is specified. See [CT Score Calculation](#ct-score-calculation) below.

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| CT baseline | `--ct-baseline` | 20.0 | CT value at which target fraction equals the baseline fraction |
| CT baseline fraction | `--ct-baseline-fraction` | 0.01 | Target fraction at the baseline CT |
| CT efficiency | `--ct-efficiency` | 1.0 | PCR amplification efficiency (0–1). Default 1.0 = 100% efficiency (doubling per cycle). Real assays typically 0.90–0.98. |
| CT calibration | `--ct-calibration` | none | Two-point calibration: `"CT1,FRAC1" "CT2,FRAC2"`. Derives efficiency automatically. Mutually exclusive with `--ct-baseline`, `--ct-baseline-fraction`, `--ct-efficiency`. |

---

## Simulation Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Simulate mode | `--simulate-mode` | thermodynamic | `thermodynamic` (SantaLucia TNN Boltzmann weighting) or `simple` (uniform probe-site weights) |
| Hybridisation temperature | `--hybridization-temperature` | 70.0 | Hybridisation temperature in °C. Only used in thermodynamic mode. |
| Capture fraction | `--capture-fraction` | 0.5 | Fraction of fragments drawn from probe binding sites (0.0–1.0); remainder are background. |

---

## Sequencing Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Read length | `--read-length` | 120 | Trim captured fragments to this length (bp). Used by `perfect` and `art`. Not applicable for `badread`. |
| Num sequences | `--num-sequences` | all | Number of reads to sample with replacement. If not set, all captured fragments become reads. |
| Read simulator | `--read-simulator` | `perfect` | `perfect` (trim, no errors), `art` (Illumina via ART-modern), `badread` (ONT/PacBio via Badread) |
| Sequencer profile | `--sequencer-profile` | `HiSeq2500_150bp` / `ont` | Chemistry/error model. Required when `--read-simulator` is `art` or `badread`. |
| Coverage depth | `--coverage-depth` | 1.0 | Reads generated per fragment (art/badread only). |
| Output format | `--output-format` | `fasta` | Output format: `fasta` or `fastq`. `fastq` preserves quality scores; perfect simulator writes dummy Q40 scores. |
| Paired-end | `--paired-end` | false | Paired-end output (art only). Produces reads.fa + reads_R2.fa. |
| PE fragment mean | `--pe-frag-len-mean` | 200 | Mean insert size for paired-end (art + `--paired-end` only) |
| PE fragment SD | `--pe-frag-len-sd` | 50 | Insert size std-dev for paired-end (art + `--paired-end` only) |

### Read Simulator Details

#### `--read-simulator perfect` (default)

Trims each captured fragment to `--read-length` bp from its start. No errors introduced. One read per fragment. Fragment names are preserved as-is.

#### `--read-simulator art` — Illumina short reads

Uses [ART-modern](https://github.com/YU-Zhejian/art_modern) to simulate Illumina-realistic base-call errors and quality scores. Requires `art_modern` on PATH (via conda: `conda install -c bioconda art_modern`).

`--sequencer-profile` selects the built-in quality profile. Common values:

| Profile | Platform | Read length |
|---------|----------|-------------|
| `HiSeq2500_150bp` (default) | Illumina HiSeq 2500 | 150 bp |
| `HiSeq2500_100bp` | Illumina HiSeq 2500 | 100 bp |
| `MiSeq_250bp` | Illumina MiSeq | 250 bp |

Run `art_modern --list-profiles` for all built-in profiles.

#### `--read-simulator badread` — ONT / PacBio CLR long reads

Uses [Badread](https://github.com/rrwick/Badread) to simulate long-read sequencing. Requires `badread` on PATH (via conda: `conda install -c conda-forge badread`).

`--sequencer-profile` selects the chemistry:

| Profile | Platform | Error model |
|---------|----------|-------------|
| `ont` (default) | ONT R10.4.1 / Kit14 | nanopore2023 |
| `ont-2020` | ONT R9.4.1 | nanopore2020 |
| `pacbio` | PacBio CLR | pacbio2016 |

`--read-length` is not used for badread — read length is bounded by fragment length and a per-profile lognormal distribution. Paired-end is not supported for long reads.

Additional flags for `badread` only:

| Flag | Default | Description |
|------|---------|-------------|
| `--long-read-length-mean` | profile default: 9000 (`ont`/`ont-2020`), 15000 (`pacbio`) | Mean read length for badread (bp) |
| `--long-read-length-sd` | profile default: 7000 (`ont`/`ont-2020`), 13000 (`pacbio`) | Standard deviation of read length for badread (bp) |
| `--badread-glitches` | badread's own default | Badread glitch parameters (`rate,size,skip`). |
| `--badread-junk-reads` | badread's own default (~1%) | Percentage of junk reads in badread output |
| `--badread-random-reads` | badread's own default (~1%) | Percentage of random reads in badread output |
| `--badread-chimeras` | badread's own default (~1%) | Percentage of chimeric reads in badread output |

---

## Execution Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Threads | `--threads` | 1 | Threads for external tools (BLAST, cd-hit-est) |
| Output dir | `--outdir`, `-o` | ./results | Output directory. A timestamped subdirectory is created for each run. |
| Output prefix | `--output-prefix` | (empty) | String prepended to every auto-generated output filename |
| Run name | `--run-name` | auto | Custom name for the run. Default: `run_YYYYMMDD_HHMMSS` |
| Report mode | `--report` | `both-r` | `full` (render HTML), `none` (skip), `rmd` (editable RMarkdown), `both-r` (HTML and RMarkdown) |
| Seed | `--seed`, `-s` | random | Random seed for reproducibility. Without this, results vary between runs. |
| Verbose | `--verbose` | false | Enable debug logging (global flag) |
| Minimap preset | `--minimap-preset` | auto | Alignment preset for read mapping. Auto-selected based on simulator: `sr` for perfect/art, `map-ont` for ONT badread, `map-pb` for PacBio badread. Override with explicit value if needed. |
| Host minimap preset | `--host-minimap-preset` | sr | Alignment preset for host read filtering |
| Cleanup | `--cleanup` | false | Delete intermediate files after completion |
| Identify | `--identify` | false | Enable species-level identification after metrics (genome mode, requires `--sample-target-map`) |
| Identity threshold | `--identity-threshold` | 90.0 | Minimum sequence identity % for target similarity in `run`, `panel-qc`, `identify` |
| Min unique targets | `--min-unique-targets` | 1 | Minimum unique target detections to call species PRESENT in `run`, `identify` |

---

## Probe Building Parameters

Used by `baitbench build-probes`. All method-specific flags are silently ignored when a different method is selected.

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Targets | `-t` / `--targets` | _(required)_ | Target sequences FASTA |
| Method | `--method` | `tile` | Probe design algorithm: `tile`, `catch-lite`, `catch`, `syotti-lite`, `probetools-lite` |
| Probe length | `--probe-length` | 120 | Probe length in bp |
| Step (tile) | `--step` | -60 | Step from end of each probe to start of next: negative = overlap, 0 = perfectly tiled, positive = gap (`tile` method only) |
| Max N fraction | `--max-n-frac` | 0.05 | Max fraction of ambiguous (N) bases in a target sequence; sequences above this are removed before design |
| No-N in probes | `--no-n-in-probes` | false | Replace N bases in designed probes: each N → T unless adjacent to T, then A/C/G |
| Min GC | `--min-gc` | 0.20 | Minimum GC fraction (0–1) to keep a probe |
| Max GC | `--max-gc` | 0.80 | Maximum GC fraction (0–1) to keep a probe |
| DUST threshold | `--dust-threshold` | 2.0 | DUST score threshold for low-complexity filtering (Morgulis et al. 2006) |
| DUST window | `--dust-window` | 64 | Window size in bases for DUST scoring |
| Max masked fraction | `--max-masked-frac` | 0.25 | Max fraction of DUST-masked bases to keep a probe; set to 1.0 to disable complexity filtering |
| Collapse threshold | `--collapse-threshold` | 0.95 | cd-hit-est identity threshold for initial target collapse |
| Dedup threshold | `--dedup-threshold` | 0.95 | cd-hit-est identity threshold for final probe deduplication |
| Threads | `--threads` | 5 | Number of threads for cd-hit-est, and for the blastn search when the assessment step uses `--aligner blast` |
| Output directory | `-o` / `--outdir` | `./build_probes_results` | Output directory |
| Output prefix | `--output-prefix` | _(none)_ | String prepended to every output filename |
| Report | `--report` | `both-r` | Report mode: `full` (HTML only), `rmd` (RMarkdown only), `both-r` (HTML + RMarkdown), `none` (skip) |
| Cleanup | `--cleanup` | false | Delete intermediate files after completion |
| Skip assess | `--skip-assess` | false | Skip the automatic probe assessment step (coverage + cross-reactivity) |

### Assessment Parameters

These apply to the `assess-probes` step that `build-probes` chains into automatically (unless `--skip-assess` is set).

| Flag | Default | Description |
|------|---------|-------------|
| `--genomes` | _(none)_ | Genome FASTA(s) for cross-reactivity checking; can be specified multiple times |
| `--threshold` | 80.0 | Minimum homology (%) for cross-reactivity detection |
| `--aligner` | `minimap2` | Alignment backend for the cross-reactivity step: `minimap2` (fast, embedded) or `blast` (blastn-short; more sensitive to weak/short homology, requires BLAST+ on PATH, uses `--threads`) |
| `--minimap-preset` | `sr` | Minimap2 alignment preset for coverage/assessment alignments (always used regardless of `--aligner`) |
| `--proximity` | 50 | Proximity distance (bp) for pull-down zone metric in coverage assessment |
| `--refine-threshold` | 80.0 | 1× coverage threshold (%) below which targets are re-analyzed in refinement iterations |
| `--refine-iterations` | _(none)_ | Number of refinement iterations to run on low-coverage targets |
| `--refine-until-stable` | false | Repeat refinement until no targets remain below `--refine-threshold` or the set stabilizes (conflicts with `--refine-iterations`) |

### `--method catch-lite` / `--method catch` Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--catch-probe-stride` | 60 | Step between candidate probe start positions during tiling |
| `--catch-mismatches` | 5 | Maximum Hamming-distance mismatches for a probe to cover a target window |
| `--catch-extension` | 0 | Extend covered interval by this many bp on each side of a probe match |
| `--catch-coverage` | 1.0 | Minimum fraction of each target sequence that must be covered |
| `--catch-minhash-threshold` | 0.6 | Jaccard similarity threshold for MinHash near-duplicate removal (0.0 = disabled) |

### `--method syotti-lite` Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--syotti-mismatches` | 40 | Maximum Hamming-distance mismatches for a bait to cover a reference window |
| `--syotti-seed-len` | 20 | Seed length (k-mer size) for the approximate search index |

### `--method probetools-lite` Parameters

| Flag | Default | Description |
|------|---------|-------------|
| `--pt-step` | 1 | Sliding window step between k-mer positions during enumeration (1 = every position) |
| `--pt-identity` | 0.9 | cd-hit-est identity threshold for k-mer clustering |
| `--pt-coverage` | 0.9 | Target coverage fraction (10th-percentile across all targets) to reach before stopping |
| `--pt-batch-size` | 100 | Number of probes to add per iteration |
| `--pt-max-panel-size` | none | Hard cap on total probes; no cap if omitted |
| `--pt-min-depth` | 1 | Minimum per-position depth to count a position as covered |
| `--pt-max-iterations` | 20 | Maximum iterations regardless of coverage progress |
| `--pt-min-coverage-gain` | 0.001 | Stop if 10th-percentile coverage improves by less than this between iterations |

---

## CT Score Calculation

CT (cycle threshold) from qPCR expresses target abundance. Lower CT = more target DNA. BaitBench converts CT to distractor fraction using:

```
target_fraction = ct_baseline_fraction × (1 + E)^(ct_baseline - ct)
distractor_fraction = 1 - target_fraction
```

Where:
- `ct_baseline` is a known reference CT (default: 20.0)
- `ct_baseline_fraction` is the target fraction at that CT (default: 0.01)
- `E` is the PCR efficiency (default: 1.0 = 100% = doubling per cycle)
- `ct` is the CT value you want to simulate

At default efficiency (E = 1.0), the formula simplifies to: `target_fraction = 0.01 × 2^(20 - ct)`.

### CT Reference Table

Default calibration: CT 20 = 1% target, 100% efficiency.

| CT | Target fraction | Distractor fraction | Interpretation |
|----|-----------------|---------------------|----------------|
| 10 | 100%* | 0% | Pure target (capped at 100%) |
| 15 | 32% | 68% | Very high abundance |
| 18 | 4% | 96% | High abundance |
| 20 | 1% | 99% | Moderate (baseline) |
| 22 | 0.25% | 99.75% | Low-moderate |
| 25 | 0.031% | 99.97% | Low abundance |
| 28 | 0.004% | 99.996% | Very low |
| 30 | 0.001% | 99.999% | Near limit of detection |
| 35 | 0.00003% | ~100% | Extremely low |
| 40 | 0.000001% | ~100% | At qPCR detection limit |

*Target fractions above 100% are capped at 100%.

### Custom One-Point Calibration

If your assay calibration differs, specify `--ct-baseline` and `--ct-baseline-fraction`:

```bash
# CT 25 = 0.1% target in your assay
baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --ct 30 \
  --ct-baseline 25 \
  --ct-baseline-fraction 0.001 \
  --num-fragments 10000 \
  --outdir results
```

### Two-Point Calibration

If you have two reference measurements, BaitBench derives efficiency automatically:

```bash
--ct-calibration "20,0.01" "25,0.0003"
```

Derived efficiency: `E = (f1/f2)^(1/(CT2 - CT1)) - 1`. The derived value is logged for verification. Mutually exclusive with `--ct-baseline`, `--ct-baseline-fraction`, and `--ct-efficiency`.

### Tips

- **Match your experimental system.** Use empirical CT-to-fraction data to set calibration parameters.
- **Sweep CT values.** Use `coverage-curve --ct-values 20 25 30` to find the limit of detection.
- **Remember the log scale.** Each 1-CT unit is a ~2× change; a 10-CT range spans ~1000×.
