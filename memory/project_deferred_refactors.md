---
name: project-deferred-refactors
description: Deferred code improvements to revisit at specific trigger points
metadata:
  type: project
---

## sequence.rs — split when adding a new read simulator

`src/commands/sequence.rs` (~660 lines) mixes the `ReadSimulator` enum dispatch with per-simulator setup and read-renaming logic. Per-simulator code for ART and Badread already lives in `src/external/`, but the dispatch and renaming stays in `sequence.rs`.

**Trigger**: next time a new read simulator is added. At that point, convert to a `src/commands/sequence/` sub-module with one file per simulator (mirroring the `src/external/` split).

**Why**: avoids growing the file further with each new simulator; makes per-simulator logic independently testable.

## GUI — CLI flag validation when GUI is overhauled

The GUI (`gui/`) passes all CLI flags as opaque strings in a `HashMap<String, String>` built in `RunView.svelte`. There is no compile-time or runtime check that GUI flag names match the actual CLI. Renames, removals, or new required arguments silently break the GUI with no visible error.

**Trigger**: the GUI major redo. At that point, add either:
- A `tests/cli_flags.rs` integration test that runs `baitbench run --help` and asserts all GUI-visible flag names appear in the output (lightweight, catches renames/removals)
- A build script that parses `cli.rs` and generates the GUI flag list at compile time (stronger, eliminates the class of bug entirely)

**Why**: the failure mode is invisible — the GUI appears to run but the flag is silently ignored or causes a wrong-value error.
