# BaitBench

BaitBench is an in-silico probe capture simulation tool for evaluating how well a probeset performs. It answers questions like:

- Does the probeset capture all target sequences?
- Does it reject background (distractor) sequences?
- Can it discriminate between organisms within the target panel?
- How does performance change at different target abundances (CT values)?
- What sequencing depth is needed for adequate genome coverage?

---

## Documentation

BaitBench's documentation is organized around four kinds of content:

| Section | Purpose |
|---------|---------|
| [**Tutorials**](tutorials/index.md) | Step-by-step walkthroughs for newcomers. Start here. |
| [**How-To Guides**](how-to/index.md) | Goal-oriented guides for specific tasks (input prep, parameter tuning, interpreting results). |
| [**Reference**](reference/index.md) | Complete, precise information: every subcommand, flag, output column, and file format. |
| [**Explanation**](explanation/index.md) | Background and concepts: the thermodynamic model, classification system, pipeline design. |

---

## Quick Start

```bash
# Install dependencies and build
conda env create -f environment.yml
conda activate baitbench
cargo build --release

# Run a simulation
./target/release/baitbench run \
  --targets targets.fa \
  --distractors distractors.fa \
  --probes probes.fa \
  --outdir results/
```

See [Your First Simulation](tutorials/first-run.md) for a complete walkthrough, or the [README](https://github.com/niel-infante/BaitBench) for a quick-start overview.
