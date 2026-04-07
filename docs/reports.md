# Report Guide

## Report Modes

The `--report` flag controls report output for `run`, `probe-coverage`, `coverage-curve`, and `report` commands:

| Mode | Description |
|------|-------------|
| `full` (default) | Render the full HTML report using R/rmarkdown. Requires R and pandoc. |
| `none` | Skip report generation entirely. All other outputs (TSV, JSON) are still produced. |
| `rmd` | Write a parameterized RMarkdown (`.Rmd`) file with all file paths and parameters pre-filled. Does not require R at run time. |

**Using `--report rmd`:**

The `rmd` mode produces an `.Rmd` file in the output directory with all parameters baked in. You can then:

1. Open the `.Rmd` file in RStudio or any text editor
2. Customize the report -- add sections, change figures, adjust formatting
3. Render it when ready:

```bash
Rscript -e 'rmarkdown::render("results/run_20250101_120000/report.Rmd")'
```

This is useful when you want to:
- Customize the report before rendering
- Render on a different machine that has R installed
- Add project-specific analysis sections
- Iterate on the report without re-running the pipeline

## HTML Report Sections

The main pipeline report (`report.html`) includes:

1. **Run Parameters** -- Table of all configuration values. Also shows the reconstructed command line for reproducibility.

2. **Capture Summary** -- Bar chart comparing fragments generated vs captured, broken down by source (sample, non-sample target, distractor, untargeted).

3. **Detection Performance** -- Bar chart of sensitivity, specificity, precision, and F1 score.

4. **Read Mapping Accuracy** -- Correctly vs incorrectly mapped reads. Incorrect mapping indicates cross-reactivity (e.g., virus A reads mapping to virus B).

5. **Confusion Matrix** -- Heatmap showing TP, FN, FP, and TN counts.

6. **Detection Detail** -- Table of every reference sequence with detection status, fragment counts, read counts, and coverage statistics.

7. **Detection Lollipop** -- Reads per detected reference, colored by classification (TP, FP_target, FP_distractor).

8. **Coverage Plots** (if coverage data available) -- Per-position read depth plots for each detected reference, with faceted overview and expandable per-reference detail views.

## Coverage Curve Report

The coverage curve report (`coverage_curve_report.html`) shows:

- **Depth curves** -- % genome covered (Y-axis) vs depth of coverage threshold on log10 scale (X-axis), with one line per parameter combination
- **Faceting** -- With < 10 combinations, all lines on one plot. With >= 10 combinations, faceted by the parameter with the fewest levels
- **Per-target panels** -- If the sample contains multiple targets, each gets its own panel
- **Summary table** -- Key depth thresholds (1x, 5x, 10x, 20x, 50x, 100x) for each combination

## Probe Coverage Report

The probe coverage report (`probe_coverage_report.html`) shows:

- **Summary table** -- Per-target coverage statistics (adapts to dataset size: simple table for <= 20 targets, interactive DT table for > 20)
- **Coverage bar charts** -- Per-target 1x/2x/5x/10x coverage (switches to histograms/boxplots for > 100 targets)
- **Depth profiles** -- Per-position probe depth plots for each target (omitted for > 100 targets)
- **Gap analysis** -- Uncovered regions and gap statistics
- **Multi-mapping probes** -- Probes that align to multiple targets (specificity concerns)

## Panel QC Report

The panel QC report (`panel_qc_report.html`) shows:

- **Panel Summary** -- Total species, targets, similar pairs, and species with zero unique targets
- **Species Discriminability** -- Bar chart (≤50 species) or histogram (>50) of discriminability scores
- **Target Composition** -- Stacked bar chart of unique vs shared targets per species
- **Species Confusion Matrix** -- Heatmap (≤30 species) or distribution statistics (>30) of shared target counts
- **Discriminability Table** -- Full per-species discriminability data (simple table ≤20, interactive DT table >20)
- **Target Similarity Pairs** -- All pairwise target similarities above the threshold

## Species Identification in Main Report

When species calls are available (from `--identify` or standalone `baitbench identify`), the main HTML report includes a "Species Identification" section with:

- **Summary table** -- Species-level sensitivity and specificity (when ground truth is available via `--sample`)
- **Species call chart** -- Bar chart of PRESENT/ABSENT/AMBIGUOUS calls per species
- **Evidence detail table** -- Full breakdown with unique/shared detected counts, reads, and explanation

## Probe Assessment Report

The probe assessment report (`assess_probes_report.html`) combines coverage and cross-reactivity analysis into a single document:

- **Build Pipeline** (conditional, when chained from build-probes) -- Pipeline stats table, sequence/base count bar charts
- **Probe Coverage** -- Summary table, coverage breadth bar charts, tiered coverage, gap analysis, pangenome depth, per-target depth profiles, proximity coverage, multi-mapping probes
- **Self-Homology** -- Plotly heatmap (≤1000 probes), density plots, hits table
- **Cross-Reactivity vs Genomes** (conditional, when `--genomes` provided) -- Plotly heatmap, per-genome bar chart, density plots, hits table
- **Parameters** -- Run configuration under a collapsible fold
