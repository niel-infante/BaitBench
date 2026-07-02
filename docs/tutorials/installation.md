# Installation

This tutorial walks you through installing BaitBench from source. By the end you will have a working `baitbench` binary and all runtime dependencies ready to use.

## Prerequisites

You will need two tools before starting:

**[Conda](https://docs.conda.io/) or [Mamba](https://mamba.readthedocs.io/)** — manages the runtime dependencies (BLAST, R, read simulators). If you do not have either installed, Mamba is recommended for faster environment creation.

**[Rust toolchain](https://rustup.rs/)** — compiles the BaitBench binary. Install it with:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

After the installer completes, open a new terminal (or run `source ~/.cargo/env`) so that `cargo` is on your PATH.

---

## 1. Get the code

```bash
git clone https://github.com/niel-infante/BaitBench.git
cd BaitBench
```

---

## 2. Install runtime dependencies

```bash
conda env create -f environment.yml
```

This creates a conda environment named `baitbench` containing:

| Package | Purpose | Required? |
|---------|---------|-----------|
| BLAST+ | Cross-reactivity analysis (`xreact`) | Yes |
| cd-hit | Sequence clustering (`build-probes`) | Yes |
| R + packages | HTML report generation | Optional |
| pandoc | HTML report rendering | Optional |
| art_modern | Read error simulation (`--read-simulator art`) | Optional |
| badread | Nanopore read simulation (`--read-simulator badread`) | Optional |
| catch | CATCH probe design (`build-probes --method catch`) | Optional |

Alignment is handled by the [rammap](https://github.com/lh3/rammap) library compiled into the BaitBench binary — no external alignment tool required.

---

## 3. Build

Activate the conda environment first so the build can find any needed system libraries, then compile:

```bash
conda activate baitbench
cargo build --release
```

The first build downloads and compiles Rust dependencies, which takes a few minutes. Subsequent builds are much faster.

The binary is produced at `target/release/baitbench`.

---

## 4. Verify

```bash
./target/release/baitbench --version
./target/release/baitbench --help
```

You should see the version string and the list of subcommands. If you see a "command not found" error, check that you are in the BaitBench directory.

---

## Add to PATH (optional)

To run `baitbench` from anywhere without the `./target/release/` prefix, choose one of these approaches:

**Symlink into a directory already on your PATH:**

```bash
ln -s "$(pwd)/target/release/baitbench" ~/.local/bin/baitbench
```

**Copy the binary:**

```bash
cp target/release/baitbench ~/.local/bin/
```

**Add the build directory to your PATH** (in `~/.bashrc` or `~/.zshrc`):

```bash
export PATH="$HOME/path/to/BaitBench/target/release:$PATH"
```

After adding to PATH, restart your shell or run `source ~/.bashrc`.

!!! note "Always activate the conda environment"
    Whatever method you choose for PATH, remember to run `conda activate baitbench`
    before using BaitBench. The conda environment provides BLAST, cd-hit, R, and other
    runtime dependencies that the binary calls out to.

---

## Reports (optional)

HTML reports require R and pandoc, which are included in `environment.yml`. They are generated automatically when you run `baitbench run` (or any report-producing command) with the conda environment active.

To skip report generation entirely, add `--report none` to any command:

```bash
baitbench run ... --report none
```

---

## Desktop GUI (optional)

A desktop GUI is available for macOS and Windows. See the [GUI installer page](https://github.com/niel-infante/BaitBench/releases/latest) for pre-built downloads, or the developer guide in `gui/README.md` if you want to build it from source.

---

## Next step

Continue to [Your First Simulation](first-run.md) to run BaitBench on the included tutorial dataset.
