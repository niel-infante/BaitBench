# BaitBench

BaitBench is an in-silico probe capture simulation tool for evaluating how well a probe set performs. It answers questions like:

- Does the probe set capture all target sequences?
- Does it reject background (distractor) sequences?
- Can it discriminate between organisms within the target panel?
- How does performance change at different target abundances (CT values)?
- What sequencing depth is needed for adequate genome coverage?

The tool aligns probes to reference sequences, scores each binding site using thermodynamic nearest-neighbor free energy (SantaLucia 1998), and generates fragments biased toward high-affinity binding sites. Background fragments fill the remainder. Reads are then mapped back to references and detection metrics are computed.

## Installation

### Prerequisites

- [Rust toolchain](https://rustup.rs/) (for building)
- [Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/) (for runtime dependencies)

### Steps

```bash
# 1. Install runtime dependencies
conda env create -f environment.yml
conda activate baitbench

# 2. Build
cargo build --release

# 3. Verify
./target/release/baitbench --help
```

The binary is at `target/release/baitbench`. Copy it to a location on your PATH or use the full path.

### Runtime Dependencies

Installed via `environment.yml`:

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >= 2.24 | Sequence alignment (simulate, mapping, filtering) |
| BLAST+ | >= 2.12 | Cross-reactivity analysis (xreact) |
| R | >= 4.2 | Report generation (optional) |
| r-ggplot2 | >= 3.4 | Figures |
| r-rmarkdown | >= 2.20 | HTML report rendering |
| r-dplyr | >= 1.1 | Data manipulation |
| r-tidyr | >= 1.3 | Data reshaping |
| r-scales | >= 1.2 | Axis formatting |
| r-knitr | >= 1.40 | Report rendering |
| r-optparse | >= 1.7 | R script CLI parsing |
| r-DT | >= 0.27 | Interactive tables |
| pandoc | >= 2.19 | Document conversion |

R and its packages are only required for full HTML report generation (`--report full`). Use `--report none` to skip report generation entirely, or `--report rmd` to produce an editable RMarkdown file that can be rendered later without requiring R at pipeline run time.
