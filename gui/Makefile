.PHONY: all dev build install copy-sidecar

# Build the CLI and start Tauri in dev mode
dev: copy-sidecar
	npm run tauri:dev

# Build everything for distribution
build: copy-sidecar
	npm run tauri:build

# Install frontend dependencies
install:
	npm install

# Build the CLI binary and copy it as a Tauri sidecar
copy-sidecar:
	cd .. && cargo build --release
	node scripts/copy-sidecar.mjs

all: install copy-sidecar
