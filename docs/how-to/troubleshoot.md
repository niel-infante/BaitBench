# Troubleshoot

Symptoms and solutions for common BaitBench problems.

---

## Input and ID problems

### "Sample ID 'X' not found in targets FASTA"

The ID in your sample manifest does not match any header in `--targets`.

**Causes and fixes:**

- **Space in FASTA header:** `>Zika virus` becomes ID `Zika`. Rename to `>Zika_virus`.
- **Case mismatch:** `Dengue_1` vs `dengue_1` — IDs are case-sensitive.
- **Trailing whitespace:** check the TSV file with `cat -A sample.tsv` (shows `$` at line ends and `^I` for tabs).
- **Wrong file:** you may be pointing `--sample` at a file from a different run whose IDs don't match the current targets.

```bash
# List all target IDs
grep ">" targets.fa | awk '{print $1}' | sed 's/>//'
# List all sample IDs
cut -f1 sample.tsv | grep -v "^#"
```

### "Sample-target-map references genome ID not found in genomes FASTA"

The genome ID in `--sample-target-map` does not match any header in `--genomes`. Apply the same ID checks as above to the genomes FASTA.

### Duplicate sequence IDs

If BaitBench reads or skips duplicates unexpectedly, check for repeated IDs:

```bash
grep ">" targets.fa | awk '{print $1}' | sort | uniq -d
```

---

## Pipeline failures

### "Cannot find R scripts directory"

Report generation requires R and the `R/` directory.

- **Fix 1:** Run from the BaitBench project root where `R/` lives.
- **Fix 2:** Set the environment variable: `export BAITBENCH_R_DIR=/path/to/BaitBench/R`
- **Fix 3:** Skip report generation: `--report none`

### "blastn: command not found" or "cd-hit-est: command not found"

`blastn` is only needed for `--aligner blast` on `xreact`, `assess-probes`, or `build-probes` (the default `--aligner minimap2` is embedded and needs no external install). `cd-hit-est` is needed for `build-probes` and `tool collapse`.

The conda environment is not active.

```bash
conda activate baitbench
```

If the error persists after activating, the tool may not be installed:

```bash
conda install -c bioconda blast   # for blastn
conda install -c bioconda cd-hit  # for cd-hit-est
```

### Build fails (Rust/cargo error)

```bash
# Check cargo is on PATH
which cargo
# If not: source ~/.cargo/env
```

If cargo is installed but the build fails with a linker error, ensure development tools are installed (on macOS: Xcode Command Line Tools via `xcode-select --install`).

---

## Simulation problems

### Very low capture rate

Symptoms: `capture_rate` near 0; almost no reads; `sample_captured` very low.

**Causes:**

- Probes don't actually match the targets — check that `--probes` and `--targets` are the correct files for each other
- `--capture-fraction` is very low (default is 0.5 — this is not the issue at default)
- Probes are too short to align reliably

**Diagnose:**

```bash
# Check that probes align to targets at all
baitbench assess-probes \
  --targets targets.fa \
  --probes probes.fa \
  --outdir check_coverage \
  --report none
cat check_coverage/cov_probe_coverage_summary.tsv
```

If `pct_covered_1x` is near 0%, the probes simply don't cover the targets.

### All targets showing as detected (no discrimination)

Symptoms: `tp_count` = total targets, `fp_target_count` = 0, `tn_target_count` = 0 — even when you provided `--sample`.

**Check:** Did you pass `--sample` on the command line? Without it, all targets are in the sample.

**Check:** Are non-sample targets very short or sharing sequence with sample targets? Probes for sample targets may cross-map to similar non-sample targets.

### High `reads_incorrectly_mapped`

Reads from one source are mapping to a different reference. This indicates cross-mapping, not missed detections.

- Run `baitbench xreact --self` on your probes to check probe self-homology
- Check `detected_detail.tsv` to see which references are collecting unexpected reads

### results.tsv shows `sensitivity: NaN` or division by zero

This happens when `sample_total = 0` — there are no sample targets. Check that `--sample` is correct and that the IDs match the targets FASTA.

---

## Unexpected metric values

### Sensitivity is 1.0 but specificity is 0.0

Every sample target was found (good), but every non-sample entity was also detected (bad). With small target counts this can happen easily:

- There may only be one non-sample entity (one distractor group) and it was detected
- Check `detected_detail.tsv` — look at the coverage pattern of the detected distractors (high 5× but low 20× breadth = background noise, not genuine capture)

### Sensitivity is 0.0

No sample targets detected at all.

- Verify that your sample manifest IDs match the FASTA — run the ID check commands above
- Check `detected_detail.tsv` — are any reads assigned? If not, there may be a problem at the mapping step
- Run `--report none` and check the log for errors during mapping

### F1 is NaN

Precision or sensitivity is 0 and the other is also 0 (or both undefined). This occurs when no sample targets exist or no detections were made. Check that the input files are not empty.

---

## Reproducibility problems

### Results differ between runs with the same settings

Without `--seed`, each run uses a different random seed. Add `--seed 42` (or any integer) to make runs reproducible:

```bash
baitbench run ... --seed 42
```

### Results differ even with `--seed`

The seed controls the fragment sampling random number generator, but not the order of alignment output from the embedded aligner, which can vary by thread count. Use `--threads 1` (the default) if strict reproducibility is required.

---

## Performance problems

### Runs are very slow

- Reduce `--num-fragments` for testing (2000–5000 is usually enough for parameter exploration)
- Use `--report none` to skip R-based report generation during iteration
- Use `--cleanup` to avoid accumulating large intermediate files across many runs

### Output directory fills up disk

Each run writes fragments, reads, SAM files, and other intermediates. Use `--cleanup` to automatically remove them after the run, keeping only the final results and report inputs:

```bash
baitbench run ... --cleanup
```

Or manually remove old run directories:

```bash
ls results/   # see what's there
rm -rf results/run_20260601_*   # remove specific runs
```
