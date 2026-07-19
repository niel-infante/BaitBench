# Release Checklist

Steps to complete before tagging a new BaitBench release.

---

## 1. Build and Test

- [ ] `cargo build --release` succeeds with no warnings
- [ ] Run the minimal example end-to-end and confirm results look correct:
  ```bash
  ./target/release/baitbench run \
    --targets examples/minimal/targets.fa \
    --distractors examples/minimal/distractors.fa \
    --probes examples/minimal/probes.fa \
    --num-fragments 1000 --seed 42 --outdir /tmp/bb_release_test
  ```
- [ ] GUI builds and launches (`cd gui && make copy-sidecar && make dev`)
- [ ] **GUI flags** — for each tool in `RunView.svelte`, verify every CLI flag is wired up and matches the current `src/cli.rs` definition. Check new flags added since last release are present in the GUI.
- [ ] **GUI tooltips** — for each tool, verify all non-obvious fields have a `data-tooltip` / `tooltip` prop. Any field added since last release should have a tooltip.

---

## 2. Reference Documentation Audit

Cross-check the reference docs against the actual source. These must be accurate at every release.

- [ ] **CLI flags** — compare `baitbench run --help` (and `--help` for each subcommand) against `docs/reference/parameters.md` and `docs/reference/commands.md`. Look for new flags, removed flags, changed defaults.
- [ ] **Output columns** — check `src/commands/metrics.rs` (`write_summary_tsv`, `write_detail_tsv`, `write_group_detail_tsv`) against `docs/reference/output-formats.md`. Confirm column names and data types match.
- [ ] **Report modes** — confirm `docs/reference/reports.md` lists all `ReportMode` variants from `src/cli.rs`.
- [ ] **Probe methods** — confirm `docs/reference/commands.md` and `docs/reference/parameters.md` list all `ProbeMethod` variants from `src/cli.rs`.
- [ ] **CLAUDE.md** — confirm the pipeline diagram, Metrics Definitions table, and CT Score Support section reflect the current code.

---

## 3. How-To Guides and Tutorials Review

These drift more slowly and don't need updating for every change, but should be reviewed before a release.

- [ ] **Tutorials** — work through `docs/tutorials/first-run.md` and `docs/tutorials/genome-mode.md` using the current binary. Confirm every command runs and produces the described output.
- [ ] **How-To Guides** — skim each guide in `docs/how-to/` for mentions of flags, defaults, or behaviors that changed since the last release. Update any that are stale.
- [ ] **Explanation pages** — verify thermodynamic formulas in `docs/explanation/thermodynamic-scoring.md` and CT formula in `docs/explanation/ct-scores.md` still match `src/thermodynamics.rs` and `src/main.rs`.

---

## 4. README

- [ ] Key Parameters table is current
- [ ] Quick Start examples still work with the current CLI
- [ ] No references to removed flags or subcommands

---

## 5. ARCHITECTURE.md

- [ ] All source files listed are still present
- [ ] New source files added since last release are documented
- [ ] Pipeline flow description matches current subcommand behavior

---

## 6. Version and Changelog

- [ ] Update version in `Cargo.toml`
- [ ] Summarize changes for the GitHub release notes

---

## 7. Tag and Release

- [ ] `git tag vX.Y.Z && git push origin vX.Y.Z`
- [ ] GUI release workflow triggers (`gui-release.yml`) — confirm artifacts build for macOS ARM64, macOS x64, Windows
- [ ] Attach release notes to the GitHub release
