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
| Fragment length sd | `--fragment-length-sd` | `(max - min) / 6` | Standard deviation of the length distribution in bp. The default reads the min/max span as a ±3σ interval, which is a guess about intent rather than a measurement — set it explicitly when the library's real spread is known. |
| Fragment length distribution | `--fragment-length-distribution` | `normal` | `normal` or `lognormal`. Sonication and enzymatic fragmentation give right-skewed lengths, which a Normal cannot represent: it is symmetric, so matching a long upper tail forces an equally long lower one. For `lognormal`, `--fragment-length-mean` and `--fragment-length-sd` are the mean and sd of the **lengths**, not of their logarithms. |

Both shapes are truncated to `[min, max]` by **redrawing**, not by clamping to the
bounds. Clamping puts a point mass on each bound, which at a wide `sd` becomes a large
spike at both ends of a distribution that should have none.

Note that truncation shifts the realised mean whenever `[min, max]` is asymmetric about
`--fragment-length-mean`, because one tail is cut more than the other. With mean 5800,
min 1000 and max 20000 the realised mean is about 6240: the lower bound sits 1.5σ below
the mean while the upper sits 4.5σ above, so the lower tail loses more. Centre the
window on the mean if you need the realised mean to match the requested one.

A contig shorter than `--fragment-length-mean` still contributes to the library in
proportion to its length, and yields fragments truncated at the contig end. `simulate`
warns when a non-trivial share of the weighted input sits on such contigs, since that
shifts the library's length distribution — draft assemblies with many short contigs are
the usual cause.

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
| Simulate mode | `--simulate-mode` | thermodynamic | `thermodynamic` (per-site occupancy from SantaLucia TNN ΔG) or `simple` (one shared occupancy for every placed probe, no temperature dependency) |
| Hybridisation temperature | `--hybridization-temperature` | 70.0 | Hybridisation temperature in °C. Only used in thermodynamic mode. |
| Capture fraction | `--capture-fraction` | 0.5 | Capture efficiency (0.0–1.0): of the fragments carrying enough of a probe's target site, the fraction the probes pull down. A preset over probe concentration — see [Capture efficiency is not output composition](#capture-efficiency-is-not-output-composition). |
| Background retention | `--background-retention` | 0.001 | Probability that a fragment overlapping no probe is retained anyway. Washing is imperfect; in a library that is mostly off-target this, not capture, supplies most of the background reads, and it sets how much enrichment the capture achieves. |
| Probe concentration | `--probe-concentration` | 0.13534 (e⁻²) | Concentration driving occupancy, in the arbitrary units absorbed into K. Fixed by default rather than back-solved — see [Why the concentration is fixed](#why-the-concentration-is-fixed). |
| Solve concentration | `--solve-concentration` | off | Back-solve the concentration to hit `--capture-fraction` instead of using a fixed one. Restores the pre-2026-11 behaviour. |
| Melting window | `--melting-window` | 10 | Pairs in the cooperative melting unit; occupancy is driven by the least stable window of this length rather than by the whole duplex. 0 uses the whole duplex. |
| Minimum site overlap | `--min-overlap-frac` | 0.5 | *Hidden/advanced.* Fraction of a probe's target site a fragment must carry before that probe can capture it. A threshold standing in for partial-duplex thermodynamics. |

### Why the concentration is fixed

`[probe]` used to be back-solved on every run so that the realized capture efficiency
matched `--capture-fraction`. That turns out to make the model insensitive to its own
thermodynamics.

With a single scalar free to move, the solver re-centres it around any change to the
affinity basis. Because occupancy saturates and the ΔG spread is enormous
(ΔΔG/RT = 34.6 between a perfect 120-mer and one with 9 mismatches at 58 °C), it hits
any target by choosing *which* probes saturate rather than by scaling them — a
rank-threshold device. Measured: changing the melting window produced bit-identical
output while the concentration was solved, with `[probe]` moving 22 orders of magnitude
to compensate.

Fixing the concentration alone does not help either, because the whole-duplex ΔG is too
sharp for any single value to place weak and strong loci both at intermediate
occupancy. The two changes are only useful together:

| | per-locus log-RMSE | worst locus |
|---|---|---|
| solved `[probe]`, whole-duplex ΔG | 2.71 | 6.98× off |
| fixed `[probe]`, 10 bp window | **1.02** | **1.74× off** |

The D6331 mock community improves at the same time (Spearman 0.904 → 0.950,
Lin's CCC 0.791 → 0.803), so the constant is not purely a Borrelia artifact.

**This is provisional.** It is one constant fitted to eight loci from one organism and
checked on one other dataset, not fitted jointly across all of them. Use
`--solve-concentration` for the old behaviour when a target on-target read fraction
matters more than relative capture between loci.

### Capture efficiency is not output composition

(Only consulted under `--solve-concentration`.) `--capture-fraction` describes **the panel and the protocol**: given a fragment that
carries enough of a probe's target site, how often the probes pull it down. It does
*not* set the on-target fraction of the output.

That distinction is deliberate. Output composition is
`abundance × capture efficiency`, so it emerges from the simulation rather than being
imposed:

* a dilution series **titrates**, instead of sitting on a plateau fixed by the
  parameter;
* `--distractor-fraction` stays meaningful across its whole range, rather than going
  inert above `1 - capture_fraction`.

Earlier versions filled a fixed quota of "capture" and "background" fragments, which
pinned the realized composition. A version during development instead solved for a
target *output* on-target fraction, which reintroduced the same problem by a different
route: adding distractor made the solver raise the probe concentration to compensate.
Both are gone.

The realized composition is reported rather than requested. `simulate` logs it and
writes it to `fragments.model_stats.tsv`, and `metrics` carries it into the summary as
`model_yield_fraction`, `model_on_target_fraction` and `model_enrichment`.
| Probe-hit cache | `--probe-hits` | none | Path to a probe-placement cache. Written on the first run, reused on every later run with the same panel, reference and aligner. See [Caching probe placements](#caching-probe-placements). |
| Probe aligner | `--probe-aligner` | seed-extend | How probes are placed on the reference for capture scoring: `seed-extend` (built-in local aligner), `minimap2` (embedded read mapper), or `blast` (blastn-short, requires BLAST+). See [Probe alignment sensitivity](#probe-alignment-sensitivity). |

### Caching probe placements

Placing a probe panel on a reference is the slowest step of a simulation — on a large
panel against a multi-gigabyte genome it can run for over an hour — and it depends only
on **the probes and the reference**. Abundance, temperature, salt, capture efficiency,
fragment length and seed do not move a single probe. A titration series, a replicate
set or a parameter sweep therefore recomputes an identical alignment on every run.

`--probe-hits PATH` stores that work at a location you choose:

```bash
# First run computes the alignment and writes the cache.
baitbench run ... --probe-hits ~/caches/mtec_vs_4path.tsv --outdir run_ct20

# Every later run with the same panel and reference reuses it, whatever else changes.
baitbench run ... --probe-hits ~/caches/mtec_vs_4path.tsv --outdir run_ct25 --ct 25
baitbench run ... --probe-hits ~/caches/mtec_vs_4path.tsv --outdir run_hot --hybridization-temperature 65
```

The cache holds the placements and their aligned base pairs, **not** the free energy.
ΔG is recomputed each run, which costs almost nothing, so one cache stays valid across
`--hybridization-temperature`, `--salt-concentration` and `--simulate-mode` sweeps.

**Validity.** The file is stamped with a content hash of the probes and of the
reference, plus the aligner used. If any differ it is recomputed and the reason is
logged. Content, not timestamps: `prepare` rebuilds `combined_reference.fa` in every run
directory, so a timestamp check would miss every time.

**One combination per file.** A cache holds one (probes, reference, aligner) triple.
Pointing one path at two different panels is safe but self-defeating — each run discards
the other's work — and BaitBench warns when it happens. Give each combination its own
path.

`coverage-curve` does this automatically, keeping one cache per prepared reference
inside its output directory, so a sweep pays for alignment once per CT value. Pass
`--probe-hits` to put it somewhere durable instead.

### Probe alignment sensitivity

Probes are placed on the reference before any thermodynamics are computed, so a probe
the aligner fails to place is invisible to the simulation. This matters because
per-locus capture efficiency is driven primarily by **probe density** — losing probes
at a locus directly distorts its predicted coverage.

Probe binding is a **local** phenomenon: 60 bp of perfect duplex out of a 120 bp probe
hybridizes and captures perfectly well. Read mappers solve the opposite problem —
placing a query end-to-end — and discard partial matches.

Measured on the MTEC bacterial panel against *Borrelia burgdorferi* B31, for the nifS
and clpA loci whose probes come only from related *Borrelia* species:

| `--probe-aligner` | Probes placed (of 228) | Time |
|---|---:|---:|
| `seed-extend` (default) | **77 (33.8%)** | 0.19 s |
| `blast` | 74 (32.5%) | 16.3 s |
| `minimap2` | 21 (9.2%) | 0.43 s |

- **`seed-extend`** — built-in ungapped local seed-and-extend. Reports partial
  placements, needs no external tools, and is the default. It is **not** a BLAST
  reimplementation: there is no gapped extension, no two-hit heuristic, no
  Karlin–Altschul statistics, no scoring matrices and no masking.
- **`minimap2`** — the embedded read mapper, retuned for divergent probes (the `sr`
  preset's `mismatch_penalty` of 8 rejects a 120-mer at 80% identity outright). Fastest
  on well-matched panels, but places probes end-to-end so partial matches are missed.
- **`blast`** — `blastn-short` with `-perc_identity 80 -qcov_hsp_perc 50`. Retained as
  an independent check; roughly 85× slower than `seed-extend` here. Without those two
  filters `blastn-short` returns a flood of tiny spurious HSPs (882k for this panel
  against one 0.9 Mb chromosome).

BaitBench reports the fraction of probes placed, and warns when it falls below 50% on
the read-mapper path.

#### How `seed-extend` works

1. Index every k-mer of both orientations of every probe (default seed length 10),
   2-bit encoded into a flat CSR table rather than a hash map.
2. Stream the reference, maintaining a rolling seed code, and collect candidate
   diagonals. Seeds occurring in more than `max_seed_occ` probe positions are skipped
   as uninformative.
3. For each candidate, take the best-scoring contiguous segment (+1 match, −2 mismatch)
   by Kadane's algorithm — this is the ungapped local extension.
4. Keep segments covering at least 50% of the probe at at least 80% identity.

The **probes** are indexed, not the reference. Indexing the reference costs one entry
per base: ~1.5 GB for a 93 Mb mock community and ~36 GB for a 2.25 GB host genome.
Probe panels are under a megabyte, so memory scales with the panel instead.

The reference scan is parallelised across reference sequences, bounded by `--threads`
(a scoped rayon pool, so it will not take every core on a shared machine). On the 93 Mb
D6331 mock community with the 16,007-probe RES panel it takes ~40 s on 10 threads,
against ~215 s on one.

`max_seed_occ` caps how many probe positions a seed may occur in before it is skipped
as uninformative, the analogue of minimap2's `mid_occ`. On typical panels it is
inactive — the average bucket holds only a few entries — but it bounds the cost on
pathological inputs such as low-complexity or homopolymer probes.


---

## Sequencing Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Read length | `--read-length` | 120 | Trim captured fragments to this length (bp). Used by `perfect` and `art`. Not applicable for `badread`. |
| Num sequences | `--num-sequences` | all | Number of reads to sample with replacement. If not set, all captured fragments become reads. |
| Read simulator | `--read-simulator` | `perfect` | `perfect` (trim, no errors), `art` (Illumina via ART-modern), `badread` (ONT/PacBio via Badread) |
| Sequencer profile | `--sequencer-profile` | `HiSeq2500_150bp` / `ont` | Chemistry/error model. Required when `--read-simulator` is `art` or `badread`. |
| Coverage depth | `--coverage-depth` | 1.0 | Reads generated per fragment (art/badread only). See the `art` section below if it produces 0 reads on short fragments. |
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

**Coverage-depth floor**: art_modern splits `--coverage-depth` across `--threads` internally
(BaitBench pins art_modern's `--parallel` to its own `--threads`, default 1). On short
(~150–350 bp) fragments, a high `--threads` combined with the default `--coverage-depth 1.0`
can push per-thread coverage below art_modern's internal minimum and produce 0 reads.
BaitBench detects this (0 reads from non-empty input) and errors out rather than silently
continuing with an empty read set. If you hit it, raise `--coverage-depth` (e.g. `5`) and use
`--num-sequences` to subsample back down to your intended read count — `--num-sequences` is
applied *after* sequencing, so it doesn't reintroduce the problem.

#### `--read-simulator badread` — ONT / PacBio CLR long reads

Uses [Badread](https://github.com/rrwick/Badread) to simulate long-read sequencing. Requires `badread` on PATH (via conda: `conda install -c conda-forge badread`).

`--sequencer-profile` selects the chemistry:

| Profile | Platform | Error model | Identity |
|---------|----------|-------------|----------|
| `ont` (default) | ONT R10.4.1 / Kit14 | nanopore2023 | badread default |
| `ont-2020` | ONT R9.4.1 | nanopore2020 | badread default |
| `pacbio` | PacBio CLR | pacbio2016 | badread default |
| `pacbio-hifi` | PacBio HiFi / CCS | pacbio2021 | `99.5,100,1.5` |

**Identity is a separate axis from the error model.** The error model shapes *what kind*
of errors appear; `--read-identity` sets *how many*. Badread's own default is
`95,99,2.5`, so a `pacbio` simulation runs at CLR-grade accuracy (~95% mean identity)
whichever error model is selected — which is why `--read-simulator perfect` was
previously the closer match to HiFi data. `pacbio-hifi` sets both a HiFi error model and
a HiFi identity.

To match a specific run rather than the generic default, measure it and pass
`--read-identity mean,max,stdev`. For reference, the D6331 mock community's real CCS
reads measure median 99.83%, mean 98.53%, sd 2.82, p5 91.53% — strongly left-skewed,
which badread's beta form reproduces.

`--read-length` is not used for badread — read length is bounded by fragment length and a per-profile lognormal distribution. Paired-end is not supported for long reads.

Additional flags for `badread` only:

| Flag | Default | Description |
|------|---------|-------------|
| `--long-read-length-mean` | profile default: 9000 (`ont`/`ont-2020`), 15000 (`pacbio`) | Mean read length for badread (bp) |
| `--read-identity` | profile default (see above) | Identity distribution: `mean,max,stdev` (beta) or `mean,stdev` (normal qscore) |
| `--long-read-length-sd` | profile default: 7000 (`ont`/`ont-2020`), 13000 (`pacbio`) | Standard deviation of read length for badread (bp) |
| `--badread-glitches` | badread's own default | Badread glitch parameters (`rate,size,skip`). |
| `--badread-junk-reads` | badread's own default (~1%) | Percentage of junk reads in badread output |
| `--badread-random-reads` | badread's own default (~1%) | Percentage of random reads in badread output |
| `--badread-chimeras` | badread's own default (~1%) | Percentage of chimeric reads in badread output |

---

## Execution Parameters

| Parameter | Flag | Default | Description |
|-----------|------|---------|-------------|
| Threads | `--threads` | 1 | Threads for probe placement, read mapping and host filtering, for external tools (BLAST, cd-hit-est), and for art_modern's `--parallel` (`--read-simulator art`). Mapping and probe placement are the two dominant costs of a run and both scale close to linearly, so raise this to your core count. The standalone `filter` and `map` subcommands take their own `--threads`. |
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
| `--catch-strand` | forward | Orientations counted as coverage: `forward` scores probes against the target as given (matches CATCH; panel sizes comparable to the original tool), `both` also credits reverse-complement matches. See [Strand handling](#strand-handling-catch-lite). |

#### Strand handling (catch-lite)

`--catch-strand` (and `--strand` on `baitbench tool catch`) controls which orientations
count when deciding whether a probe already covers a target region.

- **`forward`** (default) scores a probe against the target as given. This is what
  CATCH does, so panel sizes are directly comparable to `design.py` and to published
  panels built with it.
- **`both`** additionally credits a probe whose *reverse complement* matches the
  target. In hybrid capture of double-stranded DNA this is the more physical model — a
  bait equal to `revcomp(R)` hybridizes to the + strand of a fragment containing `R`
  and does capture it — so fewer probes are needed when a target set contains genes
  deposited on opposite strands.

Use `both` only for double-stranded targets. It is **not** valid for single-stranded
targets such as RNA-virus baits, where a bait captures one orientation only, and the
resulting panel sizes cannot be compared with CATCH's.

Because the greedy search also stops when no remaining candidate adds coverage,
`--catch-coverage` is a target rather than a guarantee. catch-lite logs the realized
coverage and warns when any target falls short of the requested fraction.

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
