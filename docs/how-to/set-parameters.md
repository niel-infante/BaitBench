# Set Parameters

This guide explains the practical effect of the main tunable parameters and when to change them. For a complete list of every flag see the [Parameters reference](../reference/parameters.md).

---

## Capture fraction (`--capture-fraction`, default 0.5)

Controls the fraction of simulated fragments that come from probe binding sites (the rest are background fragments drawn proportionally by sequence weight and length).

```bash
--capture-fraction 0.5   # 50% from probe sites, 50% background (default)
--capture-fraction 0.8   # 80% from probe sites — models highly efficient capture
--capture-fraction 0.2   # 20% from probe sites — models poor probe efficiency
```

**Effect:** Higher values concentrate reads on probe-covered regions and improve sensitivity. Lower values add more background noise, making discrimination harder.

**When to change it:**

- Increase toward 0.7–0.9 to model a well-optimised capture panel
- Decrease toward 0.1–0.3 to stress-test whether probes can work in a noisy background
- Sweep with `--capture-fraction-values 0.3 0.5 0.7` in `coverage-curve` to find the threshold where your panel breaks down

---

## Distractor fraction (`--distractor-fraction`, default 0.9)

Fraction of total fragments from distractor sequences. The default (0.9) represents 10% target DNA — a challenging but realistic clinical specimen.

```bash
--distractor-fraction 0.9    # 10% target (default)
--distractor-fraction 0.5    # 50% target — more target-rich
--distractor-fraction 0.999  # 0.1% target — near limit of detection
```

Alternatively use `--ct` to express abundance as a qPCR CT score. See [Choose a Simulation Mode](choose-simulation-mode.md) for the trade-offs.

**Effect on weights:** The weight assigned to each distractor sequence is derived from this fraction and the total sample weight. With more distractor sequences or lower sample weight, each distractor sequence gets a smaller weight per-sequence. See the [Simulation Modes explanation](../explanation/simulation-modes.md) for the formula.

---

## Hybridization temperature (`--hybridization-temperature`, default 70°C)

Temperature in °C used in the SantaLucia thermodynamic model. Affects how strongly Boltzmann weighting biases fragment sampling toward high-affinity probe sites.

```bash
--hybridization-temperature 65   # lower temp — stronger binding bias, higher enrichment
--hybridization-temperature 70   # default — standard hybridization conditions
--hybridization-temperature 75   # higher temp — weaker binding, more selective
```

**Effect:** Lower temperatures increase Boltzmann weights, biasing sampling more strongly toward high-affinity sites. Higher temperatures flatten the weights, making sampling more uniform.

**When to change it:** Match your actual hybridization conditions. A typical probe capture protocol uses 65–72°C. Only applies in `--simulate-mode thermodynamic` (the default).

---

## Number of fragments (`--num-fragments`, default 10000)

Total fragments to simulate before sequencing.

```bash
--num-fragments 2000    # fast, less statistical power
--num-fragments 10000   # default
--num-fragments 50000   # slower, more stable metrics
```

**Effect:** More fragments give smoother coverage profiles and more stable sensitivity/specificity estimates. For quick tests or parameter sweeps, 2000–5000 is usually sufficient. For final evaluations, 10000–50000 gives better confidence.

**Trade-off:** Runtime scales roughly linearly with fragment count.

---

## Read length (`--read-length`, default 120)

Trims simulated fragments to this length in base pairs (for the `perfect` and `art` simulators).

```bash
--read-length 120   # default, suitable for most Illumina short-read panels
--read-length 150   # longer Illumina reads
--read-length 75    # shorter reads, e.g. miRNA panels
```

**When to change it:** Match your sequencing platform. This parameter affects mapping uniqueness — shorter reads are harder to map unambiguously, especially in repetitive probe-target regions.

---

## Number of sequences (`--num-sequences`)

Subsamples the captured reads after sequencing. Models a fixed sequencing depth budget shared across all captured fragments.

```bash
# Default: all captured fragments become reads
# Override:
--num-sequences 5000   # subsample to 5000 reads regardless of how many were captured
```

**When to use it:** You want to compare panels at a fixed sequencing cost. Without this flag, panels with more fragments get proportionally more reads.

Sweep with `--num-sequences-values 500 1000 5000` in `coverage-curve` to find the minimum sequencing depth needed for adequate coverage.

---

## Fragment length (`--fragment-length-mean`, `--fragment-length-min`, `--fragment-length-max`)

Fragment lengths follow a truncated normal distribution:

```bash
--fragment-length-mean 175   # default
--fragment-length-min  150   # default
--fragment-length-max  200   # default
```

**When to change it:** Match your library preparation protocol. Sonication typically gives 150–300 bp; enzymatic fragmentation may give 100–200 bp. Fragment length affects how much of a probe binding site ends up in each read.

---

## Random seed (`--seed`)

Sets the random seed for reproducibility. Without `--seed`, each run uses a different seed and results will vary slightly.

```bash
--seed 42    # any integer works; results identical on re-run
```

Always use `--seed` when:

- Comparing two parameter settings (you want variation from the parameter, not from sampling noise)
- Generating results for a report or publication
- Debugging unexpected output

---

## Threads (`--threads`, default 1)

Number of threads for BLAST and cd-hit. Does not affect the main pipeline (alignment is single-threaded via the embedded rammap library).

```bash
--threads 8   # use 8 threads for xreact and build-probes
```

Only relevant when running `baitbench xreact` or `baitbench build-probes`.

---

## Parameter combinations to avoid

| Combination | Problem |
|-------------|---------|
| `--ct` and `--distractor-fraction` together | Mutually exclusive — BaitBench will error |
| `--ct-values` and `--distractor-fraction` in `coverage-curve` | Mutually exclusive |
| `--capture-fraction` very close to 0 | Almost no probe-biased fragments; sensitivity will be poor regardless of probe quality |
| `--num-fragments` very low (< 500) with many targets | Some targets may receive zero fragments by chance, producing spurious FN |
| `--read-length` > `--fragment-length-min` | Fragments shorter than `--read-length` are used as-is; reads will not all be the same length |
