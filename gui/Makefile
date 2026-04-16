.PHONY: all dev build install copy-sidecar release clean-release

TRIPLE := $(shell rustc -vV 2>/dev/null | grep '^host:' | cut -d' ' -f2)
RELEASE_DIR := dist-release

# ── Development ────────────────────────────────────────────────────────────

# Build the CLI and start Tauri in dev mode
dev: copy-sidecar
	npm run tauri:dev

# Install frontend dependencies
install:
	npm install

# Build the CLI binary and copy it as a Tauri sidecar
copy-sidecar:
	cd .. && cargo build --release
	node scripts/copy-sidecar.mjs

all: install copy-sidecar

# ── Distribution ───────────────────────────────────────────────────────────

# Full release build: CLI + GUI app bundle, signed, packaged into dist-release/
release: copy-sidecar
	npm run tauri:build
	node scripts/package-release.mjs

clean-release:
	rm -rf $(RELEASE_DIR)
