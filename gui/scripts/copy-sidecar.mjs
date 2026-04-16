#!/usr/bin/env node
/**
 * Copies the compiled baitbench binary to the Tauri sidecar location.
 * Run this before `npm run tauri:dev` or `npm run tauri:build`.
 *
 * Usage:
 *   node scripts/copy-sidecar.mjs
 */

import { execSync } from 'child_process';
import { cpSync, mkdirSync, existsSync } from 'fs';
import { resolve, dirname } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));
const repoRoot = resolve(__dirname, '../..');

// Get the current host triple from rustc
const triple = execSync('rustc -vV', { encoding: 'utf8' })
  .split('\n')
  .find(l => l.startsWith('host:'))
  ?.split(': ')[1]
  ?.trim();

if (!triple) {
  console.error('Could not determine Rust target triple from `rustc -vV`');
  process.exit(1);
}

const src = resolve(repoRoot, 'target/release/baitbench');
const destDir = resolve(__dirname, '../src-tauri/binaries');
const dest = resolve(destDir, `baitbench-${triple}`);

if (!existsSync(src)) {
  console.error(`Source binary not found: ${src}`);
  console.error('Run `cargo build --release` in the repo root first.');
  process.exit(1);
}

mkdirSync(destDir, { recursive: true });
cpSync(src, dest);
// Make executable
execSync(`chmod +x "${dest}"`);

console.log(`Copied  ${src}`);
console.log(`    →   ${dest}`);
